"""Pluggable external-scorer hook, VS enrichment metrics, and a leak-free baseline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Scoring import Scoring
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

from hitpilot.diversity import sphere_exclusion
from hitpilot.fingerprints import fingerprint, morgan_generator
from hitpilot.hygiene import read_compound_table, standardize_smiles
from hitpilot.io_utils import build_audit, load_compound_frame, sha256_file, write_audit
from hitpilot.models import HygieneConfig

MODEL_NAME = "ecfp4_rf"


class ScorerHookConfig(BaseModel):
    """Join rules, ranking direction, shortlist size, and optional hard gates."""

    id_column: str = "compound_id"
    smiles_column: str | None = None
    score_column: str = "score"
    higher_is_better: bool = True
    top_n: int = 25
    diversity_cutoff: float = 0.65
    require_lipinski: bool = False
    require_veber: bool = False
    exclude_known_matches: bool = False
    label_column: str = "label"
    radius: int = 2
    fp_size: int = 2048


class EnrichmentMetrics(BaseModel):
    """Standard virtual-screening enrichment metrics on a labeled, scored set."""

    n_compounds: int
    n_actives: int
    n_inactives: int
    roc_auc: float
    pr_auc: float
    ef_1pct: float
    ef_5pct: float
    ef_10pct: float
    bedroc_20: float
    bedroc_80_5: float
    hits_top_25: int
    hits_top_50: int
    hits_top_100: int
    k_1pct: int
    k_5pct: int
    k_10pct: int
    hits_at_1pct: int
    hits_at_5pct: int
    hits_at_10pct: int


class ScoredShortlistRecord(BaseModel):
    """One compound in the score-ordered diverse med-chem shortlist."""

    compound_id: str
    smiles: str = ""
    score: float
    score_rank: int
    shortlist_rank: int
    label: str | None = None
    mw: float | None = None
    clogp: float | None = None
    qed: float | None = None
    sa_score: float | None = None
    novelty_tier: str | None = None
    max_tanimoto: float | None = None
    murcko_scaffold_smiles: str | None = None


class ScoringSummary(BaseModel):
    """Match counts, optional enrichment, and shortlist statistics."""

    n_compounds: int
    n_score_rows: int
    n_matched: int
    n_matched_by_id: int
    n_matched_by_smiles: int
    n_unmatched_compounds: int
    n_unmatched_scores: int
    n_duplicate_score_ids: int = 0
    n_actives_scored: int = 0
    n_inactives_scored: int = 0
    n_gated_out: int = 0
    gates_skipped: list[str] = Field(default_factory=list)
    n_shortlist: int = 0
    n_shortlist_actives: int | None = None
    n_shortlist_scaffolds: int = 0
    diverse_hits_top_n: int | None = None
    enrichment: EnrichmentMetrics | None = None
    higher_is_better: bool = True
    top_n: int = 25
    diversity_cutoff: float = 0.65


def _as_str_series(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].astype(str).str.strip()


def _is_active(value: Any) -> bool:
    text = str(value).strip().lower()
    return text in {"active", "1", "true", "yes"}


def _is_inactive(value: Any) -> bool:
    text = str(value).strip().lower()
    return text in {"inactive", "0", "false", "no"}


def _truthy(value: Any) -> bool:
    text = str(value).strip().lower()
    return text in {"true", "1", "yes"}


def _canonical_smiles(value: Any) -> str | None:
    raw = "" if value is None else str(value).strip()
    if not raw:
        return None
    result = standardize_smiles(raw, HygieneConfig())
    return result.smiles


def _pick_smiles_column(frame: pd.DataFrame, explicit: str | None) -> str | None:
    if explicit and explicit in frame.columns:
        return explicit
    for name in (
        "standardized_smiles",
        "canonical_smiles",
        "smiles",
        "SMILES",
        "input_smiles",
    ):
        if name in frame.columns:
            return name
    return None


def _topk_size(n_compounds: int, fraction: float) -> int:
    if n_compounds <= 0:
        return 1
    return max(1, int(round(fraction * n_compounds)))


def _hits_in_top(labels: list[bool], k: int) -> int:
    if k <= 0:
        return 0
    return int(sum(labels[: min(k, len(labels))]))


def compute_enrichment_metrics(
    df: pd.DataFrame,
    score_col: str,
    label_col: str = "label",
    higher_is_better: bool = True,
) -> EnrichmentMetrics | None:
    """Compute ROC-AUC, PR-AUC, EF at 1/5/10%, BEDROC, and raw hit counts.

    Returns None unless the table has at least one active and one inactive
    with a numeric score. Ranking puts better scores first (high scores when
    ``higher_is_better`` is True; low docking energies when False).
    """
    if score_col not in df.columns or label_col not in df.columns:
        return None
    work = df[[score_col, label_col]].copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work.dropna(subset=[score_col])
    if work.empty:
        return None
    # Build boolean masks in Python so empty/string-dtype frames stay boolean.
    act_mask = pd.Series(
        [_is_active(v) for v in work[label_col]], index=work.index, dtype=bool
    )
    inact_mask = pd.Series(
        [_is_inactive(v) for v in work[label_col]], index=work.index, dtype=bool
    )
    labeled = work[act_mask | inact_mask]
    n_actives = int(labeled[label_col].map(_is_active).sum())
    n_inactives = int(len(labeled) - n_actives)
    if n_actives < 1 or n_inactives < 1:
        return None

    ranked = labeled.sort_values(
        score_col, ascending=not higher_is_better, kind="mergesort"
    )
    y_true = ranked[label_col].map(_is_active).astype(int).to_numpy()
    y_score = ranked[score_col].to_numpy(dtype=float)
    if not higher_is_better:
        y_score = -y_score
    labels = [bool(v) for v in y_true]
    n_compounds = len(labels)
    rate = n_actives / n_compounds

    def enrichment_factor(fraction: float) -> tuple[int, int, float]:
        k = _topk_size(n_compounds, fraction)
        hits = _hits_in_top(labels, k)
        ef = (hits / k) / rate if rate > 0 else 0.0
        return k, hits, float(ef)

    k1, hits1, ef1 = enrichment_factor(0.01)
    k5, hits5, ef5 = enrichment_factor(0.05)
    k10, hits10, ef10 = enrichment_factor(0.10)
    bedroc_rows = [[int(flag)] for flag in labels]
    return EnrichmentMetrics(
        n_compounds=n_compounds,
        n_actives=n_actives,
        n_inactives=n_inactives,
        roc_auc=float(roc_auc_score(y_true, y_score)),
        pr_auc=float(average_precision_score(y_true, y_score)),
        ef_1pct=ef1,
        ef_5pct=ef5,
        ef_10pct=ef10,
        bedroc_20=float(Scoring.CalcBEDROC(bedroc_rows, 0, 20.0)),
        bedroc_80_5=float(Scoring.CalcBEDROC(bedroc_rows, 0, 80.5)),
        hits_top_25=_hits_in_top(labels, 25),
        hits_top_50=_hits_in_top(labels, 50),
        hits_top_100=_hits_in_top(labels, 100),
        k_1pct=k1,
        k_5pct=k5,
        k_10pct=k10,
        hits_at_1pct=hits1,
        hits_at_5pct=hits5,
        hits_at_10pct=hits10,
    )


def merge_external_scores(
    compounds_df: pd.DataFrame,
    scores_df: pd.DataFrame,
    config: ScorerHookConfig | None = None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Join an external scorer CSV onto a Phase 1/2 compound table and rank it."""
    config = config or ScorerHookConfig()
    compounds = compounds_df.copy()
    scores = scores_df.copy()
    id_col = config.id_column
    score_col = config.score_column
    if id_col not in compounds.columns:
        raise ValueError(f"Compounds table is missing id column {id_col!r}")
    if id_col not in scores.columns:
        raise ValueError(f"Scores table is missing id column {id_col!r}")
    if score_col not in scores.columns:
        raise ValueError(f"Scores table is missing score column {score_col!r}")

    compounds[id_col] = _as_str_series(compounds, id_col)
    scores[id_col] = _as_str_series(scores, id_col)
    scores[score_col] = pd.to_numeric(scores[score_col], errors="coerce")
    n_duplicate_score_ids = int(scores[id_col].duplicated().sum())
    scores = scores.drop_duplicates(subset=[id_col], keep="first")

    score_cols = [id_col, score_col]
    for extra in ("model_name", "split"):
        if extra in scores.columns and extra not in score_cols:
            score_cols.append(extra)
    score_slim = scores[score_cols].rename(columns={score_col: "score"})
    # Phase 1/2 tables may already have an empty or legacy "score" column.
    # Keep it under a different name so the external scorer values win.
    if "score" in compounds.columns:
        dest = "input_score" if "input_score" not in compounds.columns else "legacy_score"
        compounds = compounds.rename(columns={"score": dest})

    merged = compounds.merge(score_slim, on=id_col, how="left", suffixes=("", "_score"))
    matched_by_id = merged["score"].notna()
    n_matched_by_id = int(matched_by_id.sum())

    n_matched_by_smiles = 0
    used_score_ids: set[str] = set()
    smiles_col_c = _pick_smiles_column(compounds, config.smiles_column)
    smiles_col_s = _pick_smiles_column(scores, config.smiles_column)
    if smiles_col_c and smiles_col_s:
        unmatched_idx = merged.index[~matched_by_id]
        if len(unmatched_idx) > 0:
            score_by_smi: dict[str, dict[str, Any]] = {}
            for _, row in scores.iterrows():
                canon = _canonical_smiles(row[smiles_col_s])
                if canon and canon not in score_by_smi and pd.notna(row[score_col]):
                    payload = {"score": float(row[score_col]), "_score_id": str(row[id_col])}
                    for extra in ("model_name", "split"):
                        if extra in scores.columns:
                            payload[extra] = row[extra]
                    score_by_smi[canon] = payload
            for idx in unmatched_idx:
                canon = _canonical_smiles(merged.at[idx, smiles_col_c])
                if canon and canon in score_by_smi:
                    payload = score_by_smi[canon]
                    merged.at[idx, "score"] = payload["score"]
                    for extra in ("model_name", "split"):
                        if extra in payload and extra in merged.columns:
                            merged.at[idx, extra] = payload[extra]
                        elif extra in payload:
                            merged.at[idx, extra] = payload[extra]
                    n_matched_by_smiles += 1
                    used_score_ids.add(payload["_score_id"])

    has_score = merged["score"].notna()
    n_matched = int(has_score.sum())
    used_score_ids |= set(merged.loc[matched_by_id, id_col].astype(str))
    score_ids = set(scores[id_col].astype(str))
    n_unmatched_scores = len(score_ids - used_score_ids)

    ranked = merged.loc[has_score].copy()
    ranked = ranked.sort_values(
        "score", ascending=not config.higher_is_better, kind="mergesort"
    )
    ranked["score_rank"] = range(1, len(ranked) + 1)
    stats = {
        "n_compounds": int(len(compounds)),
        "n_score_rows": int(len(scores_df)),
        "n_matched": n_matched,
        "n_matched_by_id": n_matched_by_id,
        "n_matched_by_smiles": n_matched_by_smiles,
        "n_unmatched_compounds": int(len(compounds) - n_matched),
        "n_unmatched_scores": int(n_unmatched_scores),
        "n_duplicate_score_ids": n_duplicate_score_ids,
    }
    return ranked, stats


def _apply_gates(
    frame: pd.DataFrame, config: ScorerHookConfig
) -> tuple[pd.DataFrame, int, list[str]]:
    kept = frame.copy()
    skipped: list[str] = []
    n_before = len(kept)
    if config.require_lipinski:
        if "lipinski_pass" in kept.columns:
            kept = kept[kept["lipinski_pass"].map(_truthy)]
        else:
            skipped.append("require_lipinski")
    if config.require_veber:
        if "veber_pass" in kept.columns:
            kept = kept[kept["veber_pass"].map(_truthy)]
        else:
            skipped.append("require_veber")
    if config.exclude_known_matches:
        if "novelty_tier" in kept.columns:
            kept = kept[
                kept["novelty_tier"].astype(str).str.strip() != "known_or_close_match"
            ]
        else:
            skipped.append("exclude_known_matches")
    return kept, int(n_before - len(kept)), skipped


def select_diverse_shortlist(
    ranked_df: pd.DataFrame,
    top_n: int = 25,
    exclusion_similarity: float = 0.65,
    *,
    smiles_column: str | None = None,
    id_column: str = "compound_id",
    radius: int = 2,
    fp_size: int = 2048,
) -> pd.DataFrame:
    """Score-ordered sphere exclusion so the best hits seed diverse spheres first."""
    if ranked_df.empty or top_n <= 0:
        out = ranked_df.iloc[0:0].copy()
        out["shortlist_rank"] = []
        return out

    smiles_col = _pick_smiles_column(ranked_df, smiles_column)
    if smiles_col is None:
        raise ValueError("No SMILES column available for diverse shortlist selection")

    # Preserve the incoming ranked order (already best-first).
    work = ranked_df.reset_index(drop=True)
    generator = morgan_generator(radius=radius, fp_size=fp_size)
    fps = []
    valid_idx: list[int] = []
    ids: list[str] = []
    for i, row in work.iterrows():
        mol = Chem.MolFromSmiles(str(row[smiles_col]).strip())
        if mol is None:
            continue
        fps.append(fingerprint(mol, generator))
        valid_idx.append(int(i))
        ids.append(str(row[id_column]))

    if not valid_idx:
        out = work.iloc[0:0].copy()
        out["shortlist_rank"] = []
        return out

    order = list(range(len(valid_idx)))  # already score-ranked
    selected, ranks, _excluded = sphere_exclusion(
        fps, ids, order, exclusion_similarity
    )
    keep_local = [i for i, flag in enumerate(selected) if flag and ranks[i] is not None]
    keep_local = sorted(keep_local, key=lambda i: int(ranks[i] or 0))[:top_n]
    rows = work.iloc[[valid_idx[i] for i in keep_local]].copy()
    rows["shortlist_rank"] = [int(ranks[i] or 0) for i in keep_local]
    return rows.reset_index(drop=True)


def _murcko_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return ""
    return Chem.MolToSmiles(scaffold, isomericSmiles=True, canonical=True)


def _shortlist_records(
    shortlist: pd.DataFrame, smiles_col: str | None, id_col: str
) -> list[ScoredShortlistRecord]:
    records: list[ScoredShortlistRecord] = []
    for _, row in shortlist.iterrows():
        smiles = str(row[smiles_col]).strip() if smiles_col else ""
        label = row["label"] if "label" in shortlist.columns else None
        if label is not None and str(label).strip() == "":
            label = None
        records.append(
            ScoredShortlistRecord(
                compound_id=str(row[id_col]),
                smiles=smiles,
                score=float(row["score"]),
                score_rank=int(row["score_rank"]),
                shortlist_rank=int(row["shortlist_rank"]),
                label=None if label is None else str(label),
                mw=_optional_float(row, "mw"),
                clogp=_optional_float(row, "clogp"),
                qed=_optional_float(row, "qed"),
                sa_score=_optional_float(row, "sa_score"),
                novelty_tier=_optional_str(row, "novelty_tier"),
                max_tanimoto=_optional_float(row, "max_tanimoto"),
                murcko_scaffold_smiles=_murcko_smiles(smiles) if smiles else "",
            )
        )
    return records


def _optional_float(row: Any, name: str) -> float | None:
    if name not in row.index:
        return None
    try:
        if pd.isna(row[name]):
            return None
        return float(row[name])
    except (TypeError, ValueError):
        return None


def _optional_str(row: Any, name: str) -> str | None:
    if name not in row.index:
        return None
    value = row[name]
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    return text or None


def run_scorer_hook(
    compounds_df: pd.DataFrame,
    scores_df: pd.DataFrame,
    config: ScorerHookConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, ScoringSummary, list[ScoredShortlistRecord]]:
    """Merge scores, compute enrichment, and build a diverse top-N shortlist."""
    config = config or ScorerHookConfig()
    ranked, stats = merge_external_scores(compounds_df, scores_df, config)
    gated, n_gated, skipped = _apply_gates(ranked, config)
    shortlist = select_diverse_shortlist(
        gated,
        top_n=config.top_n,
        exclusion_similarity=config.diversity_cutoff,
        smiles_column=config.smiles_column,
        id_column=config.id_column,
        radius=config.radius,
        fp_size=config.fp_size,
    )
    enrichment = compute_enrichment_metrics(
        ranked,
        "score",
        label_col=config.label_column,
        higher_is_better=config.higher_is_better,
    )
    smiles_col = _pick_smiles_column(shortlist, config.smiles_column)
    records = _shortlist_records(shortlist, smiles_col, config.id_column)
    n_actives = 0
    n_inactives = 0
    if config.label_column in ranked.columns:
        n_actives = int(ranked[config.label_column].map(_is_active).sum())
        n_inactives = int(ranked[config.label_column].map(_is_inactive).sum())
    n_short_actives: int | None = None
    if config.label_column in shortlist.columns:
        n_short_actives = int(shortlist[config.label_column].map(_is_active).sum())
    scaffolds = {
        rec.murcko_scaffold_smiles
        for rec in records
        if rec.murcko_scaffold_smiles
    }
    summary = ScoringSummary(
        n_compounds=stats["n_compounds"],
        n_score_rows=stats["n_score_rows"],
        n_matched=stats["n_matched"],
        n_matched_by_id=stats["n_matched_by_id"],
        n_matched_by_smiles=stats["n_matched_by_smiles"],
        n_unmatched_compounds=stats["n_unmatched_compounds"],
        n_unmatched_scores=stats["n_unmatched_scores"],
        n_duplicate_score_ids=stats["n_duplicate_score_ids"],
        n_actives_scored=n_actives,
        n_inactives_scored=n_inactives,
        n_gated_out=n_gated,
        gates_skipped=skipped,
        n_shortlist=int(len(shortlist)),
        n_shortlist_actives=n_short_actives,
        n_shortlist_scaffolds=len(scaffolds),
        diverse_hits_top_n=n_short_actives,
        enrichment=enrichment,
        higher_is_better=config.higher_is_better,
        top_n=config.top_n,
        diversity_cutoff=config.diversity_cutoff,
    )
    return ranked, shortlist, summary, records


def run_scorer_hook_file(
    compounds_path: str | Path,
    scores_path: str | Path,
    output_dir: str | Path,
    config: ScorerHookConfig | None = None,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, ScoringSummary, Any]:
    """Run the pluggable hook on two CSVs and write ranked/shortlist/audit outputs."""
    config = config or ScorerHookConfig()
    updates: dict[str, Any] = {}
    if smiles_column:
        updates["smiles_column"] = smiles_column
    if id_column:
        updates["id_column"] = id_column
    if label_column:
        updates["label_column"] = label_column
    if updates:
        config = config.model_copy(update=updates)

    compounds_path = Path(compounds_path)
    scores_path = Path(scores_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    compounds = load_compound_frame(
        compounds_path,
        smiles_column=config.smiles_column,
        id_column=config.id_column,
        label_column=config.label_column,
    )
    suffix = scores_path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        scores = pd.read_csv(scores_path, sep="\t")
    else:
        scores = pd.read_csv(scores_path)

    ranked, shortlist, summary, records = run_scorer_hook(compounds, scores, config)
    ranked_path = output_dir / "ranked_compounds.csv"
    short_path = output_dir / "shortlist_topN_diverse.csv"
    audit_path = output_dir / "audit_log.json"
    ranked.to_csv(ranked_path, index=False)
    shortlist.to_csv(short_path, index=False)
    manifest = build_audit(
        command="score",
        config={
            **config.model_dump(mode="json"),
            "compounds_path": str(compounds_path.resolve()),
            "compounds_sha256": sha256_file(compounds_path),
            "scores_path": str(scores_path.resolve()),
            "scores_sha256": sha256_file(scores_path),
        },
        input_path=compounds_path,
        output_dir=output_dir,
        summary=summary,
        outputs=[
            str(ranked_path.resolve()),
            str(short_path.resolve()),
            str(audit_path.resolve()),
        ],
    )
    write_audit(audit_path, manifest)
    return ranked, shortlist, summary, manifest


def _discover_split_files(root: Path) -> dict[tuple[str, str], Path]:
    found: dict[tuple[str, str], Path] = {}
    for path in root.rglob("*.smi"):
        name = path.name.lower()
        if name.endswith("_t.smi") or "_t." in name:
            split = "T"
        elif name.endswith("_v.smi") or "_v." in name:
            split = "V"
        else:
            continue
        if "inactive" in name:
            kind = "inactive"
        elif "active" in name:
            kind = "active"
        else:
            continue
        found[(kind, split)] = path
    return found


def _fp_matrix(smiles_list: list[str], radius: int = 2, fp_size: int = 2048) -> np.ndarray:
    generator = morgan_generator(radius=radius, fp_size=fp_size)
    matrix = np.zeros((len(smiles_list), fp_size), dtype=np.uint8)
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        DataStructs.ConvertToNumpyArray(fingerprint(mol, generator), matrix[i])
    return matrix


def generate_litpcba_baseline_scores(
    target: str = "ESR1_ant",
    data_dir: str | Path = "data/litpcba",
    output_path: str | Path | None = None,
    *,
    random_state: int = 42,
    n_estimators: int = 200,
    radius: int = 2,
    fp_size: int = 2048,
) -> dict[str, Any]:
    """Train an ECFP4 random forest on LIT-PCBA ``_T`` and score held-out ``_V``.

    Any standardized SMILES shared by train and validation is removed from
    the validation set before scoring so the evaluation is leak-free.
    """
    data_dir = Path(data_dir)
    target_dir = data_dir / target if (data_dir / target).is_dir() else data_dir
    if not target_dir.is_dir():
        raise FileNotFoundError(f"LIT-PCBA target directory not found: {target_dir}")
    files = _discover_split_files(target_dir)
    required = (("active", "T"), ("inactive", "T"), ("active", "V"), ("inactive", "V"))
    missing = [f"{k}_{s}" for k, s in required if (k, s) not in files]
    if missing:
        raise FileNotFoundError(
            f"Missing LIT-PCBA split files under {target_dir}: {missing}"
        )

    def load_split(kind: str, split: str) -> pd.DataFrame:
        table = read_compound_table(files[(kind, split)])
        table = table.copy()
        table["label"] = kind
        table["split"] = split
        return table

    train = pd.concat([load_split("active", "T"), load_split("inactive", "T")], ignore_index=True)
    valid = pd.concat([load_split("active", "V"), load_split("inactive", "V")], ignore_index=True)

    def standardize_table(frame: pd.DataFrame) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        n_invalid = 0
        for row in frame.to_dict(orient="records"):
            result = standardize_smiles(str(row.get("smiles", "")), HygieneConfig())
            if result.smiles is None:
                n_invalid += 1
                continue
            item = dict(row)
            item["canonical_smiles"] = result.smiles
            rows.append(item)
        out = pd.DataFrame(rows)
        out.attrs["n_invalid"] = n_invalid
        return out

    train_std = standardize_table(train)
    valid_std = standardize_table(valid)
    train_smiles = set(train_std["canonical_smiles"])
    overlap_mask = valid_std["canonical_smiles"].isin(train_smiles)
    n_overlap = int(overlap_mask.sum())
    overlap_ids = valid_std.loc[overlap_mask, "compound_id"].astype(str).tolist()
    valid_clean = valid_std.loc[~overlap_mask].copy()

    if train_std.empty or valid_clean.empty:
        raise ValueError("Train or leak-free validation set is empty after standardization")

    x_train = _fp_matrix(train_std["canonical_smiles"].tolist(), radius, fp_size)
    y_train = train_std["label"].map(_is_active).astype(int).to_numpy()
    model = RandomForestClassifier(
        n_estimators=n_estimators,
        class_weight="balanced",
        random_state=random_state,
        n_jobs=1,
    )
    model.fit(x_train, y_train)
    x_valid = _fp_matrix(valid_clean["canonical_smiles"].tolist(), radius, fp_size)
    proba = model.predict_proba(x_valid)
    # Class 1 is active if both classes were seen; otherwise fall back to column order.
    classes = list(model.classes_)
    active_col = classes.index(1) if 1 in classes else len(classes) - 1
    valid_clean = valid_clean.copy()
    valid_clean["score"] = proba[:, active_col]
    valid_clean["model_name"] = MODEL_NAME
    valid_clean["split"] = "V"

    out_frame = valid_clean[
        ["compound_id", "score", "model_name", "split"]
    ].copy()
    out_frame["compound_id"] = out_frame["compound_id"].astype(str)
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_frame.to_csv(output_path, index=False)

    def _counts(frame: pd.DataFrame) -> dict[str, int]:
        return {
            "n": int(len(frame)),
            "n_active": int(frame["label"].map(_is_active).sum()) if len(frame) else 0,
            "n_inactive": int(frame["label"].map(_is_inactive).sum()) if len(frame) else 0,
        }

    return {
        "target": target,
        "model_name": MODEL_NAME,
        "random_state": random_state,
        "n_estimators": n_estimators,
        "output_path": str(output_path.resolve()) if output_path is not None else None,
        "train_raw": _counts(train),
        "valid_raw": _counts(valid),
        "train_standardized": _counts(train_std),
        "valid_standardized": _counts(valid_std),
        "n_invalid_train": int(train_std.attrs.get("n_invalid", 0)),
        "n_invalid_valid": int(valid_std.attrs.get("n_invalid", 0)),
        "n_overlap_removed": n_overlap,
        "overlap_ids": overlap_ids,
        "valid_leakfree": _counts(valid_clean),
        "n_scores": int(len(out_frame)),
        "scores": out_frame,
        "train_files": {f"{k}_{s}": str(files[(k, s)]) for k, s in required},
    }
