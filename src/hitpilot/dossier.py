"""Publication-grade med-chem dossier: figures, Markdown, and self-contained HTML."""

from __future__ import annotations

import base64
import html
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import Draw
from sklearn.metrics import roc_curve

from hitpilot.io_utils import rdkit_version, sha256_file, utc_timestamp, write_audit
from hitpilot.scoring import ScoringSummary, _is_active, _is_inactive
from hitpilot.version import __version__

FIGURE_NAMES = (
    "top_shortlist_grid.png",
    "flagged_alerts_grid.png",
    "physchem_distributions.png",
    "enrichment_and_triage_funnel.png",
    "scaffold_and_clusters.png",
)

COLOR_INACTIVE = "#4C78A8"
COLOR_ACTIVE = "#E45756"
COLOR_SHORTLIST = "#E9C46A"
COLOR_SHORTLIST_EDGE = "#1D3557"
COLOR_FUNNEL = ("#1B4F72", "#1F6F8B", "#2A9D8F", "#E9C46A", "#E76F51")
NOVELTY_COLORS = {
    "known_or_close_match": "#9B2226",
    "close_analog": "#CA6702",
    "moderate_novelty": "#0A9396",
    "novel_chemotype": "#005F73",
}
NOVELTY_LABELS = {
    "known_or_close_match": "Known / close match (>=0.85)",
    "close_analog": "Close analog (0.65-0.85)",
    "moderate_novelty": "Moderate novelty (0.45-0.65)",
    "novel_chemotype": "Novel chemotype (<0.45)",
}


class DossierConfig(BaseModel):
    """Options for the publication-grade med-chem dossier."""

    target: str = "ESR1_ant"
    target_chembl_id: str = "CHEMBL206"
    target_description: str = ""
    top_n: int = 25
    n_grid: int = 12
    n_alert_examples: int = 8
    mols_per_row: int = 4
    sub_img_width: int = 380
    sub_img_height: int = 300
    dpi: int = 180
    higher_is_better: bool = True


class DossierManifest(BaseModel):
    """Auditable manifest of dossier inputs, artifacts, and key metrics."""

    timestamp_utc: str
    hitpilot_version: str
    rdkit_version: str
    command: str = "dossier"
    target: str
    target_chembl_id: str
    config: dict[str, Any] = Field(default_factory=dict)
    inputs: dict[str, dict[str, str]] = Field(default_factory=dict)
    artifacts: dict[str, dict[str, Any]] = Field(default_factory=dict)
    summary: dict[str, Any] = Field(default_factory=dict)


class DossierResult(BaseModel):
    """Paths and manifest returned by a dossier run."""

    output_dir: str
    markdown_path: str
    html_path: str
    manifest_path: str
    figures: dict[str, str] = Field(default_factory=dict)
    manifest: DossierManifest


def _setup_style(dpi: int = 180) -> None:
    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "font.size": 11,
            "axes.titlesize": 12.5,
            "axes.labelsize": 11,
            "axes.titleweight": "bold",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.28,
            "grid.linestyle": "--",
            "grid.linewidth": 0.7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "DejaVu Sans",
            "legend.frameon": False,
            "legend.fontsize": 9.5,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.28,
        }
    )


def _as_frame(source: str | Path | pd.DataFrame | None) -> pd.DataFrame:
    if source is None:
        return pd.DataFrame()
    if isinstance(source, pd.DataFrame):
        return source.copy()
    path = Path(source)
    if not path.is_file():
        return pd.DataFrame()
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _as_json(source: str | Path | dict | BaseModel | None) -> dict[str, Any]:
    if source is None:
        return {}
    if isinstance(source, BaseModel):
        return source.model_dump(mode="json")
    if isinstance(source, dict):
        return source
    path = Path(source)
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _num(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _fmt(value: Any, digits: int = 3, default: str = "\u2014") -> str:
    if value is None or (isinstance(value, float) and (math.isnan(value) or math.isinf(value))):
        return default
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none"}:
        return default
    try:
        number = float(text)
    except (TypeError, ValueError):
        return text
    if float(number).is_integer() and digits == 0:
        return str(int(number))
    return f"{number:.{digits}f}"


def _fmt_int(value: Any) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return "\u2014"


def _fmt_pct(part: int | float, whole: int | float, digits: int = 2) -> str:
    if not whole:
        return "\u2014"
    return f"{100.0 * float(part) / float(whole):.{digits}f}%"


def _is_true(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _smiles_of(row: pd.Series) -> str:
    for column in ("standardized_smiles", "smiles", "canonical_smiles", "input_smiles"):
        if column in row.index and str(row[column]).strip():
            return str(row[column]).strip()
    return ""


def _mol(smiles: str) -> Chem.Mol | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.rdDepictor.Compute2DCoords(mol)
    except Exception:  # noqa: BLE001
        pass
    return mol


def _save_fig(fig: plt.Figure, path: Path, dpi: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)
    return path


def _empty_figure(path: Path, message: str, dpi: int, size: tuple[float, float] = (10, 4)) -> Path:
    _setup_style(dpi)
    fig, ax = plt.subplots(figsize=size, constrained_layout=True)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=13, color="#4A4A4A")
    return _save_fig(fig, path, dpi)


def _legend_shortlist(row: pd.Series, rank: int) -> str:
    label = str(row["label"]).strip() if "label" in row.index else ""
    tag = f" ({label})" if label else ""
    cid = str(row.get("compound_id", "")).strip()
    score = _fmt(row["score"] if "score" in row.index else None, 3)
    qed = _fmt(row["qed"] if "qed" in row.index else None, 2)
    sa = _fmt(row["sa_score"] if "sa_score" in row.index else None, 2)
    tan = _fmt(row["max_tanimoto"] if "max_tanimoto" in row.index else None, 2)
    return f"#{rank} {cid}{tag}\nScore: {score} | QED: {qed} | SA: {sa} | Tan: {tan}"


def render_shortlist_grid(
    shortlist: pd.DataFrame,
    path: Path,
    config: DossierConfig,
) -> Path:
    """2D structure grid of the top score-ordered diverse shortlist members."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = shortlist.head(config.n_grid)
    mols: list[Chem.Mol] = []
    legends: list[str] = []
    for _i, row in rows.iterrows():
        mol = _mol(_smiles_of(row))
        if mol is None:
            continue
        rank = int(row["shortlist_rank"]) if "shortlist_rank" in row.index else len(mols) + 1
        mols.append(mol)
        legends.append(_legend_shortlist(row, rank))
    if not mols:
        return _empty_figure(path, "No shortlist structures available.", config.dpi, (12, 4))
    image = Draw.MolsToGridImage(
        mols,
        molsPerRow=config.mols_per_row,
        subImgSize=(config.sub_img_width, config.sub_img_height),
        legends=legends,
        returnPNG=False,
        legendFontSize=15,
    )
    image.save(path, format="PNG", dpi=(config.dpi, config.dpi))
    return path


def _alert_family(catalog: str) -> str:
    name = catalog.strip().upper()
    if name.startswith("PAINS"):
        return "PAINS"
    if name.startswith("BRENK"):
        return "BRENK"
    if name.startswith("NIH"):
        return "NIH"
    return name or "OTHER"


def _parse_alerts(text: str) -> list[tuple[str, str]]:
    hits: list[tuple[str, str]] = []
    for part in str(text).split("|"):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            catalog, rule = part.split(":", 1)
            hits.append((catalog.strip(), rule.strip()))
        else:
            hits.append((part, part))
    return hits


def select_alert_examples(all_compounds: pd.DataFrame, limit: int = 8) -> pd.DataFrame:
    """Pick diverse PAINS / Brenk / NIH examples, preferring actives."""
    if all_compounds.empty or "alerts" not in all_compounds.columns:
        return all_compounds.iloc[0:0].copy()
    flagged = all_compounds[all_compounds["alerts"].astype(str).str.strip() != ""].copy()
    if flagged.empty:
        return flagged
    flagged["_has_mol"] = flagged.apply(lambda row: _mol(_smiles_of(row)) is not None, axis=1)
    flagged = flagged[flagged["_has_mol"]].copy()
    if flagged.empty:
        return flagged

    buckets: dict[str, list[pd.Series]] = {"PAINS": [], "BRENK": [], "NIH": [], "OTHER": []}
    for _, row in flagged.iterrows():
        hits = _parse_alerts(row["alerts"])
        if not hits:
            continue
        families = {_alert_family(catalog) for catalog, _rule in hits}
        family = next((name for name in ("PAINS", "BRENK", "NIH") if name in families), "OTHER")
        buckets[family].append(row)

    def _sort_key(row: pd.Series) -> tuple[int, str, str]:
        active_rank = 0 if _is_active(row.get("label")) else 1
        hits = _parse_alerts(row.get("alerts", ""))
        rule = hits[0][1] if hits else ""
        return (active_rank, rule, str(row.get("compound_id", "")))

    for family in buckets:
        buckets[family] = sorted(buckets[family], key=_sort_key)

    chosen: list[pd.Series] = []
    seen_ids: set[str] = set()
    seen_rules: set[str] = set()
    families = ["PAINS", "BRENK", "NIH", "OTHER"]
    indices = {name: 0 for name in families}
    while len(chosen) < limit:
        progressed = False
        for family in families:
            pool = buckets[family]
            idx = indices[family]
            while idx < len(pool):
                row = pool[idx]
                idx += 1
                cid = str(row.get("compound_id", ""))
                if cid in seen_ids:
                    continue
                hits = _parse_alerts(row.get("alerts", ""))
                rule_key = "|".join(f"{c}:{r}" for c, r in hits)
                unused_remains = any(indices[name] < len(buckets[name]) for name in families)
                if rule_key in seen_rules and unused_remains and idx < len(pool):
                    continue
                chosen.append(row)
                seen_ids.add(cid)
                seen_rules.add(rule_key)
                progressed = True
                break
            indices[family] = idx
            if len(chosen) >= limit:
                break
        if not progressed:
            break
    if not chosen:
        return flagged.iloc[0:0].copy()
    return pd.DataFrame(chosen).reset_index(drop=True)


def render_alerts_grid(
    all_compounds: pd.DataFrame,
    path: Path,
    config: DossierConfig,
) -> tuple[Path, pd.DataFrame]:
    """2D structure grid of real Phase 1 structural-alert-flagged molecules."""
    examples = select_alert_examples(all_compounds, limit=config.n_alert_examples)
    path.parent.mkdir(parents=True, exist_ok=True)
    mols: list[Chem.Mol] = []
    legends: list[str] = []
    for _, row in examples.iterrows():
        mol = _mol(_smiles_of(row))
        if mol is None:
            continue
        cid = str(row.get("compound_id", "")).strip()
        label = str(row.get("label", "")).strip()
        hits = _parse_alerts(row.get("alerts", ""))
        parts = [f"{catalog}: {rule}" for catalog, rule in hits[:2]]
        extra = f" (+{len(hits) - 2} more)" if len(hits) > 2 else ""
        head = f"{cid}" + (f" ({label})" if label else "")
        legends.append(head + "\n" + "; ".join(parts) + extra)
        mols.append(mol)
    if not mols:
        return (
            _empty_figure(path, "No flagged alert examples available.", config.dpi, (12, 4)),
            examples,
        )
    image = Draw.MolsToGridImage(
        mols,
        molsPerRow=config.mols_per_row,
        subImgSize=(config.sub_img_width, config.sub_img_height),
        legends=legends,
        returnPNG=False,
        legendFontSize=14,
    )
    image.save(path, format="PNG", dpi=(config.dpi, config.dpi))
    return path, examples


def _scatter_groups(
    phase2: pd.DataFrame, shortlist_ids: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labels = phase2["label"] if "label" in phase2.columns else pd.Series("", index=phase2.index)
    actives = phase2[labels.map(_is_active)].copy()
    inactives = phase2[~labels.map(_is_active)].copy()
    if "compound_id" in phase2.columns:
        short = phase2[phase2["compound_id"].astype(str).isin(shortlist_ids)].copy()
    else:
        short = phase2.iloc[0:0].copy()
    return inactives, actives, short


def render_physchem_figure(
    phase2: pd.DataFrame,
    shortlist: pd.DataFrame,
    path: Path,
    config: DossierConfig,
) -> Path:
    """2x2 physchem figure: Lipinski space, QED/SA, TPSA, and novelty Tanimoto."""
    if phase2.empty:
        return _empty_figure(path, "No Phase 2 compounds available.", config.dpi, (12, 9))
    _setup_style(config.dpi)
    short_ids = set()
    if "compound_id" in shortlist.columns:
        short_ids = set(shortlist["compound_id"].astype(str))
    inactives, actives, short = _scatter_groups(phase2, short_ids)
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 10.6), constrained_layout=True)

    def _points(ax, xcol, ycol, xlabel, ylabel, title):
        xi, yi = _num(inactives, xcol), _num(inactives, ycol)
        xa, ya = _num(actives, xcol), _num(actives, ycol)
        xs, ys = _num(short, xcol), _num(short, ycol)
        ax.scatter(
            xi,
            yi,
            s=18,
            c=COLOR_INACTIVE,
            alpha=0.38,
            linewidths=0,
            label="Passed inactives",
            zorder=2,
        )
        ax.scatter(
            xa,
            ya,
            s=42,
            c=COLOR_ACTIVE,
            alpha=0.92,
            linewidths=0.3,
            edgecolors="white",
            label="Passed actives",
            zorder=3,
        )
        if len(short):
            ax.scatter(
                xs,
                ys,
                s=86,
                facecolors="none",
                edgecolors=COLOR_SHORTLIST_EDGE,
                linewidths=1.6,
                label=f"Diverse shortlist (n={len(short)})",
                zorder=4,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", pad=8)

    ax = axes[0, 0]
    _points(
        ax, "mw", "clogp", "Molecular weight (Da)", "cLogP", "A.  MW vs cLogP (Lipinski region)"
    )
    mw = _num(phase2, "mw")
    clogp = _num(phase2, "clogp")
    x_min = min(0.0, float(mw.min()) - 10) if mw.notna().any() else 0.0
    y_min = (float(clogp.min()) - 0.4) if clogp.notna().any() else -2.0
    y_min = min(y_min, 0.0)
    width = 500.0 - x_min
    height = 5.0 - y_min
    if width > 0 and height > 0:
        ax.add_patch(
            Rectangle(
                (x_min, y_min),
                width,
                height,
                fill=True,
                facecolor="#2A9D8F",
                alpha=0.08,
                edgecolor="#2A9D8F",
                linewidth=1.3,
                linestyle="--",
                zorder=1,
                label="Lipinski box (MW ≤ 500, cLogP ≤ 5)",
            )
        )
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=2, borderaxespad=0.0, fontsize=8.8)
    ax.set_xlim(left=max(-5, x_min - 15))
    if clogp.notna().any():
        ax.set_ylim(y_min - 0.15, max(float(clogp.max()) + 0.5, 5.6))

    ax = axes[0, 1]
    _points(ax, "qed", "sa_score", "QED", "Ertl SA score", "B.  QED vs Ertl SA score")
    ax.set_xlim(-0.02, 1.05)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=2, borderaxespad=0.0, fontsize=8.8)

    ax = axes[1, 0]
    tpsa_i = _num(inactives, "tpsa").dropna()
    tpsa_a = _num(actives, "tpsa").dropna()
    tpsa_s = _num(short, "tpsa").dropna()
    upper = 180.0
    if len(tpsa_i) + len(tpsa_a):
        upper = max(180.0, float(pd.concat([tpsa_i, tpsa_a]).max()) + 10)
    bins = np.linspace(0, upper, 28)
    if len(tpsa_i):
        ax.hist(
            tpsa_i, bins=bins, color=COLOR_INACTIVE, alpha=0.55, label="Passed inactives", zorder=2
        )
    if len(tpsa_a):
        ax.hist(tpsa_a, bins=bins, color=COLOR_ACTIVE, alpha=0.85, label="Passed actives", zorder=3)
    if len(tpsa_s):
        ymin = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1
        ax.scatter(
            tpsa_s,
            np.full(len(tpsa_s), max(ymin * 0.04, 0.2)),
            s=55,
            facecolors=COLOR_SHORTLIST,
            edgecolors=COLOR_SHORTLIST_EDGE,
            linewidths=1.1,
            zorder=4,
            label="Shortlist",
            clip_on=False,
        )
    ax.axvline(
        140, color="#9B2226", linestyle="--", linewidth=1.4, label="Veber TPSA ≤ 140 Å²", zorder=5
    )
    ax.set_xlabel("TPSA (Å²)")
    ax.set_ylabel("Count")
    ax.set_title("C.  TPSA distribution", loc="left", pad=8)
    ax.legend(loc="upper right", fontsize=8.8)
    ax.set_xlim(left=0)

    ax = axes[1, 1]
    tan_i = _num(inactives, "max_tanimoto").dropna()
    tan_a = _num(actives, "max_tanimoto").dropna()
    bins = np.linspace(0, 1.0, 26)
    if len(tan_i):
        ax.hist(
            tan_i, bins=bins, color=COLOR_INACTIVE, alpha=0.55, label="Passed inactives", zorder=2
        )
    if len(tan_a):
        ax.hist(tan_a, bins=bins, color=COLOR_ACTIVE, alpha=0.85, label="Passed actives", zorder=3)
    for threshold, name, color in (
        (0.45, "0.45  novel / moderate", "#005F73"),
        (0.65, "0.65  moderate / analog", "#CA6702"),
        (0.85, "0.85  analog / known", "#9B2226"),
    ):
        ax.axvline(threshold, color=color, linestyle="--", linewidth=1.3, label=name, zorder=4)
    ax.set_xlabel("Max ChEMBL Morgan Tanimoto")
    ax.set_ylabel("Count")
    ax.set_title("D.  Novelty (max Tanimoto vs ChEMBL actives)", loc="left", pad=8)
    ax.set_xlim(0, 1.02)
    ax.legend(loc="upper right", fontsize=8.5)
    return _save_fig(fig, path, config.dpi)


def _funnel_stages(
    all_compounds: pd.DataFrame,
    hygiene_audit: dict[str, Any],
    scoring_summary: dict[str, Any],
    shortlist: pd.DataFrame,
) -> list[tuple[str, int]]:
    hygiene = hygiene_audit.get("summary", hygiene_audit)
    n_raw = int(hygiene.get("n_total") or len(all_compounds) or 0)
    if not all_compounds.empty and "standardized_smiles" in all_compounds.columns:
        std = all_compounds["standardized_smiles"].astype(str).str.strip()
        n_unique = int(std[std != ""].nunique())
    else:
        n_unique = n_raw
    n_passed = int(hygiene.get("n_passed") or 0)
    n_scored = int(
        (scoring_summary.get("enrichment") or {}).get("n_compounds")
        or scoring_summary.get("n_matched")
        or 0
    )
    n_short = int(
        scoring_summary.get("n_shortlist") or (len(shortlist) if shortlist is not None else 0) or 0
    )
    return [
        ("Raw input", n_raw),
        ("Unique standardized", n_unique),
        ("Hygiene passed", n_passed),
        ("Leak-free V scored", n_scored),
        ("Diverse top-N shortlist", n_short),
    ]


def _roc_and_gain(
    ranked: pd.DataFrame, higher_is_better: bool
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    if ranked.empty or "score" not in ranked.columns or "label" not in ranked.columns:
        return None
    work = ranked.copy()
    work["score"] = pd.to_numeric(work["score"], errors="coerce")
    work = work.dropna(subset=["score"])
    y_true = work["label"].map(_is_active).astype(int).to_numpy()
    if y_true.sum() < 1 or (1 - y_true).sum() < 1:
        return None
    y_score = work["score"].to_numpy(dtype=float)
    if not higher_is_better:
        y_score = -y_score
    fpr, tpr, _thresholds = roc_curve(y_true, y_score)
    order = np.argsort(-y_score, kind="mergesort")
    recovered = np.cumsum(y_true[order]) / y_true.sum()
    screened = (np.arange(len(order)) + 1) / len(order)
    return fpr, tpr, screened, recovered


def render_enrichment_funnel_figure(
    all_compounds: pd.DataFrame,
    ranked: pd.DataFrame,
    shortlist: pd.DataFrame,
    hygiene_audit: dict[str, Any],
    scoring_summary: dict[str, Any],
    path: Path,
    config: DossierConfig,
) -> Path:
    """1x3 figure: triage funnel, ROC curve, and cumulative enrichment."""
    _setup_style(config.dpi)
    fig, axes = plt.subplots(1, 3, figsize=(16.8, 5.7), constrained_layout=True)
    stages = _funnel_stages(all_compounds, hygiene_audit, scoring_summary, shortlist)
    raw = stages[0][1] or 1

    ax = axes[0]
    labels = [name for name, _count in stages]
    counts = [count for _name, count in stages]
    y_pos = np.arange(len(stages))[::-1]
    bars = ax.barh(y_pos, counts, color=list(COLOR_FUNNEL), height=0.68, zorder=3)
    ax.set_yticks(y_pos, labels)
    ax.set_xlabel("Compounds")
    ax.set_title("A.  HitPilot triage funnel", loc="left", pad=8)
    xmax = max(counts) if counts else 1
    ax.set_xlim(0, xmax * 1.28)
    for bar, count in zip(bars, counts, strict=True):
        ax.text(
            bar.get_width() + xmax * 0.02,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}  ({_fmt_pct(count, raw, 1)})",
            va="center",
            ha="left",
            fontsize=9.5,
            color="#1D3557",
        )
    ax.grid(axis="y", visible=False)

    curves = _roc_and_gain(ranked, config.higher_is_better)
    enrichment = scoring_summary.get("enrichment") or {}
    ax = axes[1]
    if curves is None:
        ax.text(0.5, 0.5, "ROC requires labeled actives and inactives.", ha="center", va="center")
        ax.set_axis_off()
    else:
        fpr, tpr, _screened, _recovered = curves
        ax.plot(fpr, tpr, color="#1B4F72", linewidth=2.2, label="ECFP4 RF rescorer", zorder=3)
        ax.plot(
            [0, 1], [0, 1], color="#8D99AE", linestyle="--", linewidth=1.2, label="Random", zorder=2
        )
        auc = enrichment.get("roc_auc")
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        title = "B.  ROC curve"
        if auc is not None:
            title = f"B.  ROC curve  (AUC = {_fmt(auc, 3)})"
        ax.set_title(title, loc="left", pad=8)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-0.01, 1.03)
        ax.legend(loc="lower right")
        ax.set_aspect("equal", adjustable="box")

    ax = axes[2]
    if curves is None:
        ax.text(0.5, 0.5, "Enrichment curve requires labeled scores.", ha="center", va="center")
        ax.set_axis_off()
    else:
        _fpr, _tpr, screened, recovered = curves
        ax.plot(
            100 * screened,
            100 * recovered,
            color="#E76F51",
            linewidth=2.2,
            label="Observed recovery",
            zorder=3,
        )
        ax.plot(
            [0, 100],
            [0, 100],
            color="#8D99AE",
            linestyle="--",
            linewidth=1.2,
            label="Random",
            zorder=2,
        )
        notes = []
        for frac, ef, name in (
            (1.0, enrichment.get("ef_1pct"), "EF 1%"),
            (5.0, enrichment.get("ef_5pct"), "EF 5%"),
        ):
            if ef is None:
                continue
            idx = int(np.searchsorted(screened, frac / 100.0, side="left"))
            idx = min(max(idx, 0), len(recovered) - 1)
            y = 100 * recovered[idx]
            ax.scatter([frac], [y], color="#1D3557", s=28, zorder=4)
            notes.append(f"{name} = {_fmt(ef, 2)}")
        bedroc = enrichment.get("bedroc_20")
        if bedroc is not None:
            notes.append(f"BEDROC \u03b1=20 = {_fmt(bedroc, 3)}")
        if notes:
            ax.text(
                0.98,
                0.22,
                "\n".join(notes),
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=9.5,
                color="#1D3557",
                linespacing=1.45,
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "#D6DEE8",
                    "linewidth": 0.8,
                    "alpha": 0.92,
                },
                zorder=5,
            )
        ax.set_xlabel("% of library screened")
        ax.set_ylabel("% of actives recovered")
        ax.set_title("C.  Cumulative enrichment", loc="left", pad=8)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 105)
        ax.legend(loc="upper left")
    return _save_fig(fig, path, config.dpi)


def render_scaffold_cluster_figure(
    phase2: pd.DataFrame,
    scaffold_summary: pd.DataFrame,
    path: Path,
    config: DossierConfig,
) -> Path:
    """Cluster sizes plus novelty tiers and top Murcko frequencies."""
    _setup_style(config.dpi)
    fig = plt.figure(figsize=(14.6, 6.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.05, 1.15])
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b1 = fig.add_subplot(grid[0, 1])
    ax_b2 = fig.add_subplot(grid[1, 1])

    sizes = pd.Series(dtype=int)
    if "cluster_id" in phase2.columns and "cluster_size" in phase2.columns:
        sizes = (
            phase2.assign(_sz=_num(phase2, "cluster_size"))
            .dropna(subset=["_sz"])
            .drop_duplicates(subset=["cluster_id"])["_sz"]
            .astype(int)
        )
    elif "cluster_size" in phase2.columns:
        sizes = _num(phase2, "cluster_size").dropna().astype(int)
    buckets = {"Singletons (1)": 0, "Size 2–3": 0, "Size 4–5": 0, "Size 6+": 0}
    for size in sizes:
        if size <= 1:
            buckets["Singletons (1)"] += 1
        elif size <= 3:
            buckets["Size 2–3"] += 1
        elif size <= 5:
            buckets["Size 4–5"] += 1
        else:
            buckets["Size 6+"] += 1
    colors = ["#1B4F72", "#2A9D8F", "#E9C46A", "#E76F51"]
    bars = ax_a.bar(list(buckets.keys()), list(buckets.values()), color=colors, width=0.7, zorder=3)
    for bar, count in zip(bars, buckets.values(), strict=True):
        ax_a.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{count:,}",
            ha="center",
            va="bottom",
            fontsize=10,
            color="#1D3557",
        )
    ax_a.set_ylabel("Number of clusters")
    ax_a.set_title("A.  Butina cluster size distribution", loc="left", pad=8)
    ax_a.set_ylim(0, max(list(buckets.values()) or [1]) * 1.15)
    ax_a.tick_params(axis="x", labelsize=9.5)
    ax_a.grid(axis="x", visible=False)

    if "novelty_tier" in phase2.columns:
        counts = phase2["novelty_tier"].replace("", np.nan).dropna().value_counts()
    else:
        counts = pd.Series(dtype=int)
    order = ["novel_chemotype", "moderate_novelty", "close_analog", "known_or_close_match"]
    labels = [NOVELTY_LABELS.get(key, key) for key in order if key in counts.index]
    values = [int(counts.get(key, 0)) for key in order if key in counts.index]
    bar_colors = [NOVELTY_COLORS.get(key, "#4A4A4A") for key in order if key in counts.index]
    if labels:
        y_pos = np.arange(len(labels))[::-1]
        bars = ax_b1.barh(y_pos, values, color=bar_colors, height=0.64, zorder=3)
        ax_b1.set_yticks(y_pos, labels, fontsize=9)
        for bar, value in zip(bars, values, strict=True):
            ax_b1.text(
                bar.get_width() + max(values) * 0.02,
                bar.get_y() + bar.get_height() / 2,
                f"{value:,}",
                va="center",
                fontsize=9.5,
            )
        ax_b1.set_xlim(0, max(values) * 1.2)
    ax_b1.set_xlabel("Compounds")
    ax_b1.set_title("B.  Novelty tier breakdown", loc="left", pad=8)
    ax_b1.grid(axis="y", visible=False)

    scaf = scaffold_summary.copy() if scaffold_summary is not None else pd.DataFrame()
    if scaf.empty and "murcko_scaffold_smiles" in phase2.columns:
        tmp = phase2.copy()
        tmp["murcko_scaffold_smiles"] = tmp["murcko_scaffold_smiles"].replace("", "(acyclic)")
        scaf = (
            tmp.groupby("murcko_scaffold_smiles", dropna=False)
            .size()
            .reset_index(name="count")
            .sort_values("count", ascending=False)
        )
    if not scaf.empty:
        if "count" not in scaf.columns:
            scaf = scaf.copy()
            scaf["count"] = 1
        scaf["count"] = pd.to_numeric(scaf["count"], errors="coerce").fillna(0).astype(int)
        if "murcko_scaffold_smiles" in scaf.columns:
            smiles = scaf["murcko_scaffold_smiles"].fillna("").astype(str)
            if "acyclic" in scaf.columns:
                acyclic = scaf["acyclic"].map(_is_true)
            else:
                acyclic = smiles.isin(["", "nan"])
            labels_s = []
            for smi, acy in zip(smiles, acyclic, strict=True):
                if acy or not smi or smi.lower() == "nan":
                    labels_s.append("(acyclic)")
                elif len(smi) <= 32:
                    labels_s.append(smi)
                else:
                    labels_s.append(smi[:30] + "...")
            scaf = scaf.assign(_label=labels_s)
            top = scaf.sort_values("count", ascending=False).head(8)
            y_pos = np.arange(len(top))[::-1]
            ax_b2.barh(y_pos, top["count"], color="#1F6F8B", height=0.66, zorder=3)
            ax_b2.set_yticks(y_pos, list(top["_label"]), fontsize=8.2)
            for y, value in zip(y_pos, top["count"], strict=True):
                ax_b2.text(
                    value + max(int(top["count"].max()), 1) * 0.02,
                    y,
                    f"{int(value):,}",
                    va="center",
                    fontsize=9,
                )
            ax_b2.set_xlim(0, max(int(top["count"].max()), 1) * 1.22)
        else:
            ax_b2.text(0.5, 0.5, "No scaffold summary available.", ha="center", va="center")
            ax_b2.set_axis_off()
    else:
        ax_b2.text(0.5, 0.5, "No scaffold summary available.", ha="center", va="center")
        ax_b2.set_axis_off()
    ax_b2.set_xlabel("Compounds")
    ax_b2.set_title("B.  Top Bemis-Murcko scaffolds", loc="left", pad=8)
    ax_b2.grid(axis="y", visible=False)
    return _save_fig(fig, path, config.dpi)


def _md_escape(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _html_esc(value: Any) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def _badge_label(label: str) -> str:
    if _is_active(label):
        return '<span class="badge badge-active">active</span>'
    if _is_inactive(label):
        return '<span class="badge badge-inactive">inactive</span>'
    return _html_esc(label) if label else "\u2014"


def _badge_tier(tier: str) -> str:
    key = (tier or "").strip()
    if not key:
        return "\u2014"
    css = {
        "known_or_close_match": "badge-known",
        "close_analog": "badge-analog",
        "moderate_novelty": "badge-moderate",
        "novel_chemotype": "badge-novel",
    }.get(key, "badge-muted")
    return f'<span class="badge {css}">{_html_esc(key)}</span>'


def _img_md(name: str, title: str) -> str:
    return f"![{title}](figures/{name})"


def _b64_uri(path: Path) -> str:
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def _input_hashes(paths: dict[str, Path | None]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for name, path in paths.items():
        if path is None or not Path(path).is_file():
            continue
        resolved = Path(path).resolve()
        out[name] = {"path": str(resolved), "sha256": sha256_file(resolved)}
    return out


def _artifact_record(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved),
        "bytes": int(resolved.stat().st_size),
    }


def compute_funnel_and_summary(
    all_compounds: pd.DataFrame,
    phase2: pd.DataFrame,
    ranked: pd.DataFrame,
    shortlist: pd.DataFrame,
    hygiene_audit: dict[str, Any],
    phase2_audit: dict[str, Any],
    scoring_summary: dict[str, Any],
    baseline_stats: dict[str, Any],
    config: DossierConfig,
    alert_examples: pd.DataFrame,
) -> dict[str, Any]:
    """Assemble the exact counts and metrics used in the dossier text and figures."""
    stages = _funnel_stages(all_compounds, hygiene_audit, scoring_summary, shortlist)
    hygiene = hygiene_audit.get("summary", hygiene_audit)
    phase2_sum = phase2_audit.get("summary", phase2_audit)
    enrichment = scoring_summary.get("enrichment")
    n_unique_inchikey = 0
    if "inchikey" in all_compounds.columns:
        keys = all_compounds["inchikey"].astype(str).str.strip()
        n_unique_inchikey = int(keys[keys != ""].nunique())
    n_short_scaf = scoring_summary.get("n_shortlist_scaffolds")
    if n_short_scaf is None and "murcko_scaffold_smiles" in shortlist.columns and len(shortlist):
        n_short_scaf = int(shortlist["murcko_scaffold_smiles"].replace("", np.nan).nunique())
    return {
        "funnel": [{"stage": name, "count": count} for name, count in stages],
        "hygiene": hygiene,
        "phase2": phase2_sum,
        "scoring": scoring_summary,
        "enrichment": enrichment,
        "baseline": baseline_stats,
        "n_phase2": int(len(phase2)),
        "n_ranked": int(len(ranked)),
        "n_shortlist": int(len(shortlist)),
        "n_shortlist_scaffolds": int(n_short_scaf or 0),
        "diverse_hits_top_n": scoring_summary.get("diverse_hits_top_n"),
        "n_unique_inchikey": n_unique_inchikey,
        "n_alert_examples": int(len(alert_examples)),
        "target": config.target,
        "target_chembl_id": config.target_chembl_id,
        "target_description": config.target_description,
        "top_n": config.top_n,
    }


def _scaffold_markdown_rows(scaffold_summary: pd.DataFrame) -> list[str]:
    rows = [
        "| Rank | Murcko scaffold | Generic framework | Acyclic | Count | Actives | Inactives | Examples |",  # noqa: E501
        "| ---: | --- | --- | --- | ---: | ---: | ---: | --- |",
    ]
    if scaffold_summary is None or scaffold_summary.empty:
        rows.append("| \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 |")
        return rows
    scaf = scaffold_summary.copy()
    if "count" in scaf.columns:
        scaf["count"] = pd.to_numeric(scaf["count"], errors="coerce").fillna(0).astype(int)
        scaf = scaf.sort_values("count", ascending=False)
    for rank, (_, row) in enumerate(scaf.head(15).iterrows(), start=1):
        smi = str(row.get("murcko_scaffold_smiles") or "").strip()
        acy = _is_true(row.get("acyclic")) if "acyclic" in row.index else not bool(smi)
        if acy or not smi or smi.lower() == "nan":
            smi = "(acyclic)"
            acy = True
        gen = str(row.get("generic_scaffold_smiles") or "").strip() or "\u2014"
        if gen.lower() == "nan":
            gen = "\u2014"
        examples = str(row.get("example_compound_ids") or "").replace("|", ", ")
        n_act = row.get("n_active", "")
        n_inact = row.get("n_inactive", "")
        rows.append(
            f"| {rank} | `{_md_escape(smi)}` | `{_md_escape(gen)}` | {acy} | "
            f"{_fmt_int(row.get('count', 0))} | "
            f"{_fmt_int(n_act) if str(n_act) != '' else chr(0x2014)} | "
            f"{_fmt_int(n_inact) if str(n_inact) != '' else chr(0x2014)} | "
            f"{_md_escape(examples)} |"
        )
    return rows


def _write_markdown(
    path: Path,
    config: DossierConfig,
    summary: dict[str, Any],
    shortlist: pd.DataFrame,
    scaffold_summary: pd.DataFrame,
    alert_examples: pd.DataFrame,
    inputs: dict[str, dict[str, str]],
    timestamp: str,
) -> Path:
    hygiene = summary.get("hygiene") or {}
    enrichment = summary.get("enrichment") or {}
    baseline = summary.get("baseline") or {}
    phase2_sum = summary.get("phase2") or {}
    title_bits = f"{config.target}"
    if config.target_description:
        title_bits += f" \u2014 {config.target_description}"
    lines: list[str] = [
        f"# HitPilot Med-Chem Dossier: {title_bits}",
        "",
        f"**Target:** `{config.target}` (LIT-PCBA)  \u00b7  **ChEMBL reference:** `{config.target_chembl_id}`  ",  # noqa: E501
        f"**Generated:** {timestamp} (UTC)  \u00b7  **HitPilot** {__version__}  \u00b7  **RDKit** {rdkit_version()}",  # noqa: E501
        "",
        "All chemistry and metrics in this dossier are computed from public LIT-PCBA and ChEMBL",
        "data with verified RDKit / scikit-learn code. Scores, SMILES, and enrichment values are",
        "never invented.",
        "",
        "## Input audit (SHA-256)",
        "",
        "| Input | SHA-256 | Path |",
        "| --- | --- | --- |",
    ]
    for name, meta in inputs.items():
        lines.append(
            f"| `{_md_escape(name)}` | `{meta['sha256']}` | `{_md_escape(meta['path'])}` |"
        )
    if not inputs:
        lines.append("| \u2014 | \u2014 | \u2014 |")

    raw = (summary.get("funnel") or [{"count": 0}])[0].get("count", 0)
    lines += [
        "",
        "## Executive summary and triage funnel",
        "",
        "HitPilot applies a strict, auditable funnel: standardize and deduplicate the public",
        "LIT-PCBA library, drop PAINS / Brenk / NIH liabilities, attach physicochemical",
        "descriptors and ChEMBL nearest-neighbor novelty, then rescore only the leak-free",
        "held-out `_V` split (training is confined to `_T`) and pick a score-ordered diverse",
        "shortlist by Morgan sphere exclusion.",
        "",
        "| Stage | Count | % of raw input |",
        "| --- | ---: | ---: |",
    ]
    for item in summary.get("funnel") or []:
        count = int(item["count"])
        lines.append(f"| {item['stage']} | {_fmt_int(count)} | {_fmt_pct(count, raw, 2)} |")

    if baseline:
        train = baseline.get("train_standardized") or {}
        valid = baseline.get("valid_standardized") or {}
        leak = baseline.get("valid_leakfree") or {}
        n_overlap = baseline.get("n_overlap_removed", 0)
        lines += [
            "",
            "### Leak-free `_T` -> `_V` audit",
            "",
            "The official LIT-PCBA train/validation files are SID-level, so the same standardized",
            "structure can appear in both `_T` and `_V`. Before the baseline scorer is trained,",
            "every `_V` structure whose canonical SMILES is present in `_T` is removed.",
            "",
            "| Split | Standardized | Actives | Inactives |",
            "| --- | ---: | ---: | ---: |",
            f"| Train `_T` | {_fmt_int(train.get('n', ''))} | {_fmt_int(train.get('n_active', ''))} | {_fmt_int(train.get('n_inactive', ''))} |",  # noqa: E501
            f"| Valid `_V` (before leak removal) | {_fmt_int(valid.get('n', ''))} | {_fmt_int(valid.get('n_active', ''))} | {_fmt_int(valid.get('n_inactive', ''))} |",  # noqa: E501
            f"| Held-out `_V` scored (leak-free) | {_fmt_int(leak.get('n', ''))} | {_fmt_int(leak.get('n_active', ''))} | {_fmt_int(leak.get('n_inactive', ''))} |",  # noqa: E501
            "",
            f"**Standardized SMILES overlap removed from `_V`:** {_fmt_int(n_overlap)} "
            f"(structures shared with `_T`; evaluation is leak-free).",
        ]

    by_label = hygiene.get("by_label") or {}
    if by_label:
        lines += [
            "",
            "Hygiene pass rates by label:",
            "",
            "| Label | Input | Passed | Dropped | Pass rate |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
        for label, breakdown in by_label.items():
            if isinstance(breakdown, dict):
                lines.append(
                    f"| {label} | {_fmt_int(breakdown.get('n_total', 0))} | "
                    f"{_fmt_int(breakdown.get('n_passed', 0))} | "
                    f"{_fmt_int(breakdown.get('n_dropped', 0))} | "
                    f"{_fmt_pct(breakdown.get('n_passed', 0), breakdown.get('n_total', 0) or 1, 2)} |"  # noqa: E501
                )

    props = (phase2_sum.get("properties") or {}) if isinstance(phase2_sum, dict) else {}
    nov = (phase2_sum.get("novelty") or {}) if isinstance(phase2_sum, dict) else {}
    div = (phase2_sum.get("diversity") or {}) if isinstance(phase2_sum, dict) else {}
    if props or nov or div:
        lines += ["", "Phase 2 context (hygiene-passed compounds):", ""]
        if props:
            lines.append(
                f"- Median MW {_fmt(props.get('median_mw'), 1)}, cLogP {_fmt(props.get('median_clogp'), 2)}, "  # noqa: E501
                f"TPSA {_fmt(props.get('median_tpsa'), 1)}, QED {_fmt(props.get('median_qed'), 3)}, "  # noqa: E501
                f"SA {_fmt(props.get('median_sa_score'), 3)}; "
                f"Lipinski pass {_fmt_int(props.get('lipinski_pass_count'))} "
                f"({_fmt_pct(props.get('lipinski_pass_count') or 0, props.get('n_total') or 0, 2)}), "  # noqa: E501
                f"Veber pass {_fmt_int(props.get('veber_pass_count'))} "
                f"({_fmt_pct(props.get('veber_pass_count') or 0, props.get('n_total') or 0, 2)})."
            )
        if nov:
            tiers = nov.get("counts_by_tier") or {}
            lines.append(
                f"- Novelty vs `{nov.get('target_chembl_id', config.target_chembl_id)}` "
                f"({_fmt_int(nov.get('n_reference_actives'))} unique actives): "
                + ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in tiers.items())
                + f"; median max Tanimoto {_fmt(nov.get('median_max_tanimoto'), 3)}."
            )
        if div:
            lines.append(
                f"- {_fmt_int(div.get('n_unique_murcko_scaffolds'))} unique Murcko scaffolds, "
                f"{_fmt_int(div.get('n_clusters'))} Butina clusters, "
                f"{_fmt_int(div.get('n_sphere_selected'))} sphere-exclusion representatives "
                f"({_fmt_int(div.get('n_actives_selected'))}/"
                f"{_fmt_int(div.get('n_actives_total'))} actives captured)."
            )

    lines += [
        "",
        "## Figures",
        "",
        "### Top diverse shortlist",
        "",
        _img_md("top_shortlist_grid.png", "Top diverse shortlist structures"),
        "",
        "### Phase 1 structural-alert examples",
        "",
        _img_md("flagged_alerts_grid.png", "Flagged PAINS, Brenk, and NIH examples"),
        "",
        "### Physicochemical and novelty distributions",
        "",
        _img_md("physchem_distributions.png", "Physicochemical distributions"),
        "",
        "### Enrichment and triage funnel",
        "",
        _img_md("enrichment_and_triage_funnel.png", "Triage funnel, ROC, and enrichment"),
        "",
        "### Scaffolds and clusters",
        "",
        _img_md("scaffold_and_clusters.png", "Butina clusters and Murcko scaffolds"),
        "",
        "## Virtual-screening enrichment metrics",
        "",
        "Metrics are computed only on hygiene-passed compounds that also belong to the leak-free",
        "held-out `_V` set (matched to baseline or user-supplied scores). Ranking respects",
        "`higher_is_better`. ROC-AUC and PR-AUC use scikit-learn; EF uses the standard VS",
        "definition with k = max(1, round(fN)); BEDROC uses RDKit CalcBEDROC at alpha = 20",
        "and alpha = 80.5.",
        "",
    ]
    if enrichment:
        rows = [
            ("Scored compounds (N)", _fmt_int(enrichment.get("n_compounds"))),
            (
                "Actives / inactives",
                f"{_fmt_int(enrichment.get('n_actives'))} / {_fmt_int(enrichment.get('n_inactives'))}",  # noqa: E501
            ),
            ("ROC-AUC", _fmt(enrichment.get("roc_auc"), 4)),
            ("PR-AUC", _fmt(enrichment.get("pr_auc"), 4)),
            (
                "EF 1%",
                f"{_fmt(enrichment.get('ef_1pct'), 3)} (k={enrichment.get('k_1pct')}, hits={enrichment.get('hits_at_1pct')})",  # noqa: E501
            ),
            (
                "EF 5%",
                f"{_fmt(enrichment.get('ef_5pct'), 3)} (k={enrichment.get('k_5pct')}, hits={enrichment.get('hits_at_5pct')})",  # noqa: E501
            ),
            (
                "EF 10%",
                f"{_fmt(enrichment.get('ef_10pct'), 3)} (k={enrichment.get('k_10pct')}, hits={enrichment.get('hits_at_10pct')})",  # noqa: E501
            ),
            ("BEDROC alpha=20", _fmt(enrichment.get("bedroc_20"), 4)),
            ("BEDROC alpha=80.5", _fmt(enrichment.get("bedroc_80_5"), 4)),
            (
                "Hits in top 25 / 50 / 100",
                f"{enrichment.get('hits_top_25')} / {enrichment.get('hits_top_50')} / {enrichment.get('hits_top_100')}",  # noqa: E501
            ),
            (
                "Diverse shortlist hits / unique Murcko scaffolds",
                f"{summary.get('diverse_hits_top_n')} / {summary.get('n_shortlist_scaffolds')}",
            ),
        ]
        lines += ["| Metric | Value |", "| --- | --- |"]
        for key, value in rows:
            lines.append(f"| {key} | {value} |")
    else:
        lines.append(
            "_No enrichment metrics were available (need both actives and inactives with scores)._"
        )

    lines += [
        "",
        f"## Top-{config.top_n} diverse hit shortlist",
        "",
        "Score-ordered Morgan sphere exclusion (pairwise Tanimoto < 0.65 by default) so the",
        "best-scoring hits seed diverse spheres first.",
        "",
        "| Rank | ID | Label | Score | MW | cLogP | TPSA | QED | SA | Lipinski | Veber | Cluster | Max Tan | Nearest ChEMBL (pChEMBL) | Novelty | Canonical SMILES |",  # noqa: E501
        "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: | --- | --- | --- |",  # noqa: E501
    ]
    table = shortlist.head(config.top_n)
    if table.empty:
        lines.append(
            "| \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 | \u2014 |"  # noqa: E501
        )
    for _, row in table.iterrows():
        rank = row["shortlist_rank"] if "shortlist_rank" in row.index else ""
        lip = row["lipinski_pass"] if "lipinski_pass" in row.index else ""
        veber = row["veber_pass"] if "veber_pass" in row.index else ""
        cluster = row["cluster_id"] if "cluster_id" in row.index else ""
        if "cluster_size" in row.index and str(row.get("cluster_size", "")).strip():
            cluster = f"{cluster} (n={row['cluster_size']})"
        nearest = str(row["nearest_ref_id"]).strip() if "nearest_ref_id" in row.index else ""
        pchembl = (
            _fmt(row["nearest_ref_pchembl"], 2) if "nearest_ref_pchembl" in row.index else "\u2014"
        )
        nearest_cell = f"{nearest} ({pchembl})" if nearest else "\u2014"
        smiles = _smiles_of(row)
        rank_s = str(int(rank)) if str(rank).strip() and str(rank) != "nan" else "\u2014"
        lines.append(
            "| "
            + " | ".join(
                [
                    rank_s,
                    _md_escape(row.get("compound_id", "")),
                    _md_escape(row.get("label", "")),
                    _fmt(row.get("score"), 3),
                    _fmt(row.get("mw"), 1),
                    _fmt(row.get("clogp"), 2),
                    _fmt(row.get("tpsa"), 1),
                    _fmt(row.get("qed"), 3),
                    _fmt(row.get("sa_score"), 3),
                    _md_escape(lip) if str(lip).strip() else "\u2014",
                    _md_escape(veber) if str(veber).strip() else "\u2014",
                    _md_escape(cluster) if str(cluster).strip() else "\u2014",
                    _fmt(row.get("max_tanimoto"), 3),
                    _md_escape(nearest_cell),
                    _md_escape(row.get("novelty_tier", "")),
                    f"`{_md_escape(smiles)}`" if smiles else "\u2014",
                ]
            )
            + " |"
        )

    lines += ["", "## Top Butina clusters and Bemis-Murcko scaffolds", ""]
    lines += _scaffold_markdown_rows(scaffold_summary)

    drop = hygiene.get("counts_by_drop_reason") or {}
    catalogs = hygiene.get("counts_by_catalog") or {}
    lines += [
        "",
        "## Phase 1 structural-alert breakdown",
        "",
        "Drop-reason counts are not mutually exclusive: a compound may match more than one",
        "catalog or also be marked as a duplicate.",
        "",
        "| Drop reason or catalog | Compounds |",
        "| --- | ---: |",
    ]
    for key, value in drop.items():
        lines.append(f"| `{_md_escape(key)}` | {_fmt_int(value)} |")
    if catalogs:
        lines.append("| | |")
        for key, value in catalogs.items():
            lines.append(f"| Catalog `{_md_escape(key)}` | {_fmt_int(value)} |")
    lines += [
        "",
        "### Example flagged liabilities",
        "",
        "Examples below are the same structures drawn in `flagged_alerts_grid.png`, selected to",
        "cover PAINS, Brenk, and NIH (actives first when available).",
        "",
        "| ID | Label | Catalog + rule | Canonical SMILES |",
        "| --- | --- | --- | --- |",
    ]
    if alert_examples.empty:
        lines.append("| \u2014 | \u2014 | \u2014 | \u2014 |")
    for _, row in alert_examples.iterrows():
        hits = _parse_alerts(row.get("alerts", ""))
        rule = "; ".join(f"{c}: {r}" for c, r in hits)
        lines.append(
            f"| {_md_escape(row.get('compound_id', ''))} | {_md_escape(row.get('label', ''))} | "
            f"{_md_escape(rule)} | `{_md_escape(_smiles_of(row))}` |"
        )
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


HTML_CSS = """
:root { --ink:#1D3557; --muted:#5C6B7A; --line:#D6DEE8; --bg:#F7F9FC; --card:#FFFFFF; }
*{box-sizing:border-box;}
html,body{margin:0;padding:0;background:var(--bg);color:var(--ink);
  font:16px/1.55 "Source Sans 3","Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
.wrap{max-width:1180px;margin:0 auto;padding:32px 24px 64px;}
header{margin-bottom:28px;}
h1{font-size:1.85rem;line-height:1.25;margin:0 0 8px;}
h2{font-size:1.35rem;margin:2.2rem 0 0.7rem;border-bottom:1px solid var(--line);padding-bottom:6px;}
h3{font-size:1.08rem;margin:1.4rem 0 0.5rem;}
.meta{color:var(--muted);font-size:0.95rem;}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;
  padding:18px 20px;margin:14px 0 22px;box-shadow:0 1px 2px rgba(29,53,87,0.04);}
figure{margin:0 0 8px;}
figure img{width:100%;height:auto;border:1px solid var(--line);border-radius:10px;background:#fff;}
figcaption{color:var(--muted);font-size:0.9rem;margin-top:6px;}
table{width:100%;border-collapse:collapse;font-size:0.86rem;}
th,td{border-bottom:1px solid var(--line);padding:7px 8px;text-align:left;vertical-align:top;}
th{background:#EEF3F8;font-weight:600;}
tr:nth-child(even) td{background:#FBFCFE;}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:0.82em;word-break:break-all;}
.badge{display:inline-block;padding:1px 8px;border-radius:999px;font-size:0.75rem;font-weight:650;}
.badge-active,.badge-known{background:#FDECEC;color:#9B2226;}
.badge-inactive{background:#E7F0F8;color:#1B4F72;}
.badge-analog{background:#FFF1E0;color:#9A4A00;}
.badge-moderate{background:#E4F5F5;color:#0A6B6E;}
.badge-novel{background:#E4F0F3;color:#005F73;}
.badge-muted{background:#EEE;color:#444;}
p{margin:0.6rem 0;}
.note{color:var(--muted);font-size:0.92rem;}
"""


def _html_table(headers: list[str], rows: list[list[str]], numeric: set[int] | None = None) -> str:
    numeric = numeric or set()
    head = "".join(
        f'<th class="num">{_html_esc(h)}</th>' if i in numeric else f"<th>{_html_esc(h)}</th>"
        for i, h in enumerate(headers)
    )
    body_rows = []
    for row in rows:
        cells = []
        for i, cell in enumerate(row):
            klass = ' class="num"' if i in numeric else ""
            cells.append(f"<td{klass}>{cell}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    body = "".join(body_rows)
    return (
        "<div class='card'><table><thead><tr>"
        + head
        + "</tr></thead><tbody>"
        + body
        + "</tbody></table></div>"
    )


def _write_html(
    path: Path,
    config: DossierConfig,
    summary: dict[str, Any],
    shortlist: pd.DataFrame,
    scaffold_summary: pd.DataFrame,
    alert_examples: pd.DataFrame,
    figures: dict[str, Path],
    inputs: dict[str, dict[str, str]],
    timestamp: str,
) -> Path:
    enrichment = summary.get("enrichment") or {}
    baseline = summary.get("baseline") or {}
    hygiene = summary.get("hygiene") or {}
    title = f"HitPilot Med-Chem Dossier: {config.target}"
    if config.target_description:
        title += f" \u2014 {config.target_description}"
    parts: list[str] = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        f"<title>{_html_esc(title)}</title>",
        f"<style>{HTML_CSS}</style></head><body><div class='wrap'>",
        f"<header><h1>{_html_esc(title)}</h1>",
        "<p class='meta'>",
        f"Target <code>{_html_esc(config.target)}</code> (LIT-PCBA) \u00b7 ChEMBL ",
        f"<code>{_html_esc(config.target_chembl_id)}</code> \u00b7 Generated {_html_esc(timestamp)} UTC \u00b7 ",  # noqa: E501
        f"HitPilot {__version__} \u00b7 RDKit {_html_esc(rdkit_version())}",
        "</p>",
        "<p class='note'>All chemistry and metrics are computed from public LIT-PCBA and ChEMBL ",
        "data with verified RDKit / scikit-learn code. Scores, SMILES, and enrichment values are never invented.</p>",  # noqa: E501
        "</header>",
        "<h2>Input audit (SHA-256)</h2>",
    ]
    audit_rows = [
        [
            _html_esc(name),
            f"<code>{_html_esc(meta['sha256'])}</code>",
            f"<code>{_html_esc(meta['path'])}</code>",
        ]
        for name, meta in inputs.items()
    ]
    parts.append(
        _html_table(["Input", "SHA-256", "Path"], audit_rows or [["\u2014", "\u2014", "\u2014"]])
    )
    parts.append("<h2>Executive summary and triage funnel</h2>")
    parts.append(
        "<p>HitPilot applies a strict, auditable funnel: standardize and deduplicate the public "
        "LIT-PCBA library, drop PAINS / Brenk / NIH liabilities, attach physicochemical descriptors "  # noqa: E501
        "and ChEMBL nearest-neighbor novelty, then rescore only the leak-free held-out <code>_V</code> "  # noqa: E501
        "split (training is confined to <code>_T</code>) and pick a score-ordered diverse shortlist "  # noqa: E501
        "by Morgan sphere exclusion.</p>"
    )
    raw = (summary.get("funnel") or [{"count": 0}])[0]["count"]
    funnel_rows = [
        [_html_esc(item["stage"]), _fmt_int(item["count"]), _fmt_pct(item["count"], raw, 2)]
        for item in summary.get("funnel") or []
    ]
    parts.append(_html_table(["Stage", "Count", "% of raw input"], funnel_rows, numeric={1, 2}))

    if baseline:
        train = baseline.get("train_standardized") or {}
        valid = baseline.get("valid_standardized") or {}
        leak = baseline.get("valid_leakfree") or {}
        parts += [
            "<h3>Leak-free <code>_T</code> \u2192 <code>_V</code> audit</h3>",
            "<p>The official LIT-PCBA train/validation files are SID-level, so the same standardized "  # noqa: E501
            "structure can appear in both splits. Every <code>_V</code> structure whose canonical SMILES "  # noqa: E501
            "is present in <code>_T</code> is removed before scoring.</p>",
            _html_table(
                ["Split", "Standardized", "Actives", "Inactives"],
                [
                    [
                        "Train <code>_T</code>",
                        _fmt_int(train.get("n")),
                        _fmt_int(train.get("n_active")),
                        _fmt_int(train.get("n_inactive")),
                    ],
                    [
                        "Valid <code>_V</code> (before leak removal)",
                        _fmt_int(valid.get("n")),
                        _fmt_int(valid.get("n_active")),
                        _fmt_int(valid.get("n_inactive")),
                    ],
                    [
                        "Held-out <code>_V</code> scored (leak-free)",
                        _fmt_int(leak.get("n")),
                        _fmt_int(leak.get("n_active")),
                        _fmt_int(leak.get("n_inactive")),
                    ],
                ],
                numeric={1, 2, 3},
            ),
            "<p><strong>Standardized SMILES overlap removed from <code>_V</code>:</strong> "
            f"{_fmt_int(baseline.get('n_overlap_removed', 0))} "
            "(structures shared with <code>_T</code>; evaluation is leak-free).</p>",
        ]

    captions = [
        (
            "top_shortlist_grid.png",
            "Top diverse shortlist",
            "Score-ordered diverse shortlist structures with QED, SA, and nearest-neighbor Tanimoto.",  # noqa: E501
        ),
        (
            "flagged_alerts_grid.png",
            "Phase 1 structural-alert examples",
            "Diverse PAINS, Brenk, and NIH FilterCatalog matches (actives first when available).",
        ),
        (
            "physchem_distributions.png",
            "Physicochemical and novelty distributions",
            "Lipinski space, QED vs Ertl SA, Veber TPSA, and ChEMBL Morgan Tanimoto novelty.",
        ),
        (
            "enrichment_and_triage_funnel.png",
            "Enrichment and triage funnel",
            "HitPilot funnel, ROC curve, and cumulative enrichment on the leak-free held-out set.",
        ),
        (
            "scaffold_and_clusters.png",
            "Scaffolds and clusters",
            "Butina cluster sizes, novelty tiers, and the most frequent Bemis-Murcko scaffolds.",
        ),
    ]
    parts.append("<h2>Figures</h2>")
    for name, heading, caption in captions:
        fig_path = figures.get(name)
        if fig_path is None or not Path(fig_path).is_file():
            continue
        parts += [
            f"<h3>{_html_esc(heading)}</h3>",
            "<div class='card'><figure>",
            f"<img src='{_b64_uri(Path(fig_path))}' alt='{_html_esc(heading)}'>",
            f"<figcaption>{_html_esc(caption)}</figcaption>",
            "</figure></div>",
        ]

    parts.append("<h2>Virtual-screening enrichment metrics</h2>")
    if enrichment:
        metric_rows = [
            ["Scored compounds (N)", _fmt_int(enrichment.get("n_compounds"))],
            [
                "Actives / inactives",
                f"{_fmt_int(enrichment.get('n_actives'))} / {_fmt_int(enrichment.get('n_inactives'))}",  # noqa: E501
            ],
            ["ROC-AUC", _fmt(enrichment.get("roc_auc"), 4)],
            ["PR-AUC", _fmt(enrichment.get("pr_auc"), 4)],
            [
                "EF 1%",
                f"{_fmt(enrichment.get('ef_1pct'), 3)} (k={enrichment.get('k_1pct')}, hits={enrichment.get('hits_at_1pct')})",  # noqa: E501
            ],
            [
                "EF 5%",
                f"{_fmt(enrichment.get('ef_5pct'), 3)} (k={enrichment.get('k_5pct')}, hits={enrichment.get('hits_at_5pct')})",  # noqa: E501
            ],
            [
                "EF 10%",
                f"{_fmt(enrichment.get('ef_10pct'), 3)} (k={enrichment.get('k_10pct')}, hits={enrichment.get('hits_at_10pct')})",  # noqa: E501
            ],
            ["BEDROC alpha=20", _fmt(enrichment.get("bedroc_20"), 4)],
            ["BEDROC alpha=80.5", _fmt(enrichment.get("bedroc_80_5"), 4)],
            [
                "Hits in top 25 / 50 / 100",
                f"{enrichment.get('hits_top_25')} / {enrichment.get('hits_top_50')} / {enrichment.get('hits_top_100')}",  # noqa: E501
            ],
            [
                "Diverse shortlist hits / unique Murcko scaffolds",
                f"{summary.get('diverse_hits_top_n')} / {summary.get('n_shortlist_scaffolds')}",
            ],
        ]
        parts.append(_html_table(["Metric", "Value"], metric_rows, numeric={1}))
    else:
        parts.append("<p class='note'>No enrichment metrics were available.</p>")

    parts.append(f"<h2>Top-{config.top_n} diverse hit shortlist</h2>")
    parts.append(
        "<p class='note'>Score-ordered Morgan sphere exclusion (pairwise Tanimoto &lt; 0.65 by default).</p>"  # noqa: E501
    )
    sl_rows = []
    for _, row in shortlist.head(config.top_n).iterrows():
        rank = row["shortlist_rank"] if "shortlist_rank" in row.index else ""
        cluster = str(row["cluster_id"]).strip() if "cluster_id" in row.index else ""
        if "cluster_size" in row.index and str(row.get("cluster_size", "")).strip():
            cluster = f"{cluster} (n={row['cluster_size']})"
        nearest = str(row["nearest_ref_id"]).strip() if "nearest_ref_id" in row.index else ""
        pchembl = (
            _fmt(row["nearest_ref_pchembl"], 2) if "nearest_ref_pchembl" in row.index else "\u2014"
        )
        nearest_cell = f"{_html_esc(nearest)} ({pchembl})" if nearest else "\u2014"
        rank_s = str(int(rank)) if str(rank).strip() and str(rank) != "nan" else "\u2014"
        sl_rows.append(
            [
                rank_s,
                _html_esc(row.get("compound_id", "")),
                _badge_label(str(row.get("label", ""))),
                _fmt(row.get("score"), 3),
                _fmt(row.get("mw"), 1),
                _fmt(row.get("clogp"), 2),
                _fmt(row.get("tpsa"), 1),
                _fmt(row.get("qed"), 3),
                _fmt(row.get("sa_score"), 3),
                _html_esc(row.get("lipinski_pass", "")) or "\u2014",
                _html_esc(row.get("veber_pass", "")) or "\u2014",
                _html_esc(cluster) if cluster else "\u2014",
                _fmt(row.get("max_tanimoto"), 3),
                nearest_cell,
                _badge_tier(str(row.get("novelty_tier", ""))),
                f"<code>{_html_esc(_smiles_of(row))}</code>" if _smiles_of(row) else "\u2014",
            ]
        )
    parts.append(
        _html_table(
            [
                "Rank",
                "ID",
                "Label",
                "Score",
                "MW",
                "cLogP",
                "TPSA",
                "QED",
                "SA",
                "Lipinski",
                "Veber",
                "Cluster",
                "Max Tan",
                "Nearest ChEMBL (pChEMBL)",
                "Novelty",
                "Canonical SMILES",
            ],
            sl_rows or [["\u2014"] * 16],
            numeric={0, 3, 4, 5, 6, 7, 8, 12},
        )
    )

    parts.append("<h2>Top Butina clusters and Bemis-Murcko scaffolds</h2>")
    scaf = scaffold_summary.copy() if scaffold_summary is not None else pd.DataFrame()
    scaf_rows: list[list[str]] = []
    if not scaf.empty:
        if "count" in scaf.columns:
            scaf["count"] = pd.to_numeric(scaf["count"], errors="coerce").fillna(0).astype(int)
            scaf = scaf.sort_values("count", ascending=False)
        for rank, (_, row) in enumerate(scaf.head(15).iterrows(), start=1):
            smi = str(row.get("murcko_scaffold_smiles") or "").strip()
            acy = _is_true(row.get("acyclic")) if "acyclic" in row.index else not bool(smi)
            if acy or not smi or smi.lower() == "nan":
                smi = "(acyclic)"
                acy = True
            gen = str(row.get("generic_scaffold_smiles") or "").strip() or "\u2014"
            if gen.lower() == "nan":
                gen = "\u2014"
            examples = str(row.get("example_compound_ids") or "").replace("|", ", ")
            n_act = row.get("n_active", "")
            n_inact = row.get("n_inactive", "")
            scaf_rows.append(
                [
                    str(rank),
                    f"<code>{_html_esc(smi)}</code>",
                    f"<code>{_html_esc(gen)}</code>",
                    str(acy),
                    _fmt_int(row.get("count", 0)),
                    _fmt_int(n_act) if str(n_act) != "" else "\u2014",
                    _fmt_int(n_inact) if str(n_inact) != "" else "\u2014",
                    _html_esc(examples),
                ]
            )
    parts.append(
        _html_table(
            [
                "Rank",
                "Murcko scaffold",
                "Generic framework",
                "Acyclic",
                "Count",
                "Actives",
                "Inactives",
                "Examples",
            ],
            scaf_rows or [["\u2014"] * 8],
            numeric={0, 4, 5, 6},
        )
    )

    parts.append("<h2>Phase 1 structural-alert breakdown</h2>")
    parts.append(
        "<p class='note'>Drop-reason counts are not mutually exclusive: a compound may match more "
        "than one catalog or also be marked as a duplicate.</p>"
    )
    drop_rows = [
        [f"<code>{_html_esc(k)}</code>", _fmt_int(v)]
        for k, v in (hygiene.get("counts_by_drop_reason") or {}).items()
    ]
    drop_rows += [
        [f"Catalog <code>{_html_esc(k)}</code>", _fmt_int(v)]
        for k, v in (hygiene.get("counts_by_catalog") or {}).items()
    ]
    parts.append(
        _html_table(
            ["Drop reason or catalog", "Compounds"],
            drop_rows or [["\u2014", "\u2014"]],
            numeric={1},
        )
    )
    parts.append("<h3>Example flagged liabilities</h3>")
    ex_rows = []
    for _, row in alert_examples.iterrows():
        hits = _parse_alerts(row.get("alerts", ""))
        rule = "; ".join(f"{c}: {r}" for c, r in hits)
        ex_rows.append(
            [
                _html_esc(row.get("compound_id", "")),
                _badge_label(str(row.get("label", ""))),
                _html_esc(rule),
                f"<code>{_html_esc(_smiles_of(row))}</code>",
            ]
        )
    parts.append(
        _html_table(
            ["ID", "Label", "Catalog + rule", "Canonical SMILES"],
            ex_rows or [["\u2014", "\u2014", "\u2014", "\u2014"]],
        )
    )
    parts.append("</div></body></html>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(parts), encoding="utf-8")
    return path


def generate_medchem_dossier(
    output_dir: str | Path,
    *,
    all_compounds: str | Path | pd.DataFrame,
    phase2_enriched: str | Path | pd.DataFrame,
    ranked: str | Path | pd.DataFrame | None = None,
    shortlist: str | Path | pd.DataFrame | None = None,
    hygiene_audit: str | Path | dict | None = None,
    phase2_audit: str | Path | dict | None = None,
    scoring_audit: str | Path | dict | ScoringSummary | None = None,
    baseline_stats: str | Path | dict | None = None,
    scaffold_summary: str | Path | pd.DataFrame | None = None,
    config: DossierConfig | None = None,
    input_paths: dict[str, str | Path | None] | None = None,
) -> DossierResult:
    """Write figures, Markdown, self-contained HTML, and an auditable manifest."""
    config = config or DossierConfig()
    output_dir = Path(output_dir)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    all_df = _as_frame(all_compounds)
    phase2_df = _as_frame(phase2_enriched)
    ranked_df = _as_frame(ranked)
    short_df = _as_frame(shortlist)
    scaf_df = _as_frame(scaffold_summary)
    hygiene = _as_json(hygiene_audit)
    phase2_js = _as_json(phase2_audit)
    scoring = _as_json(scoring_audit)
    if scoring and "summary" in scoring and "enrichment" not in scoring:
        scoring = {**scoring["summary"], "config": scoring.get("config", {})}
    baseline = _as_json(baseline_stats)
    if "higher_is_better" in scoring:
        config.higher_is_better = bool(scoring["higher_is_better"])

    fig_paths: dict[str, Path] = {name: figures_dir / name for name in FIGURE_NAMES}
    render_shortlist_grid(short_df, fig_paths["top_shortlist_grid.png"], config)
    _alerts_path, alert_examples = render_alerts_grid(
        all_df, fig_paths["flagged_alerts_grid.png"], config
    )
    render_physchem_figure(phase2_df, short_df, fig_paths["physchem_distributions.png"], config)
    render_enrichment_funnel_figure(
        all_df,
        ranked_df,
        short_df,
        hygiene,
        scoring,
        fig_paths["enrichment_and_triage_funnel.png"],
        config,
    )
    render_scaffold_cluster_figure(
        phase2_df, scaf_df, fig_paths["scaffold_and_clusters.png"], config
    )

    summary = compute_funnel_and_summary(
        all_df,
        phase2_df,
        ranked_df,
        short_df,
        hygiene,
        phase2_js,
        scoring,
        baseline,
        config,
        alert_examples,
    )
    timestamp = utc_timestamp()
    resolved_inputs = _input_hashes(
        {
            key: Path(value) if value is not None else None
            for key, value in (input_paths or {}).items()
        }
    )
    if not resolved_inputs:
        auto: dict[str, Path | None] = {}
        for name, source in {
            "all_compounds": all_compounds,
            "phase2_enriched": phase2_enriched,
            "ranked_compounds": ranked,
            "shortlist": shortlist,
            "hygiene_audit": hygiene_audit,
            "phase2_audit": phase2_audit,
            "scoring_audit": scoring_audit,
            "baseline_stats": baseline_stats,
            "scaffold_summary": scaffold_summary,
        }.items():
            if isinstance(source, (str, Path)) and Path(source).is_file():
                auto[name] = Path(source)
        resolved_inputs = _input_hashes(auto)

    md_path = output_dir / "medchem_dossier.md"
    html_path = output_dir / "medchem_dossier.html"
    manifest_path = output_dir / "dossier_manifest.json"
    _write_markdown(
        md_path,
        config,
        summary,
        short_df,
        scaf_df,
        alert_examples,
        resolved_inputs,
        timestamp,
    )
    _write_html(
        html_path,
        config,
        summary,
        short_df,
        scaf_df,
        alert_examples,
        fig_paths,
        resolved_inputs,
        timestamp,
    )

    artifacts = {name: _artifact_record(path) for name, path in fig_paths.items() if path.is_file()}
    artifacts["medchem_dossier.md"] = _artifact_record(md_path)
    artifacts["medchem_dossier.html"] = _artifact_record(html_path)
    manifest = DossierManifest(
        timestamp_utc=timestamp,
        hitpilot_version=__version__,
        rdkit_version=rdkit_version(),
        command="dossier",
        target=config.target,
        target_chembl_id=config.target_chembl_id,
        config=config.model_dump(mode="json"),
        inputs=resolved_inputs,
        artifacts=artifacts,
        summary=summary,
    )
    write_audit(manifest_path, manifest)
    artifacts["dossier_manifest.json"] = _artifact_record(manifest_path)
    manifest.artifacts = artifacts
    write_audit(manifest_path, manifest)

    return DossierResult(
        output_dir=str(output_dir.resolve()),
        markdown_path=str(md_path.resolve()),
        html_path=str(html_path.resolve()),
        manifest_path=str(manifest_path.resolve()),
        figures={name: str(path.resolve()) for name, path in fig_paths.items()},
        manifest=manifest,
    )


def generate_medchem_dossier_from_dirs(
    output_dir: str | Path,
    *,
    phase1_dir: str | Path,
    phase2_dir: str | Path,
    phase3_dir: str | Path | None = None,
    config: DossierConfig | None = None,
) -> DossierResult:
    """Build a dossier from Phase 1-3 output directories."""
    phase1_dir = Path(phase1_dir)
    phase2_dir = Path(phase2_dir)
    phase3_dir = Path(phase3_dir) if phase3_dir is not None else None
    baseline_stats = None
    if phase3_dir is not None and (phase3_dir / "baseline_stats.json").is_file():
        baseline_stats = phase3_dir / "baseline_stats.json"
    return generate_medchem_dossier(
        output_dir,
        all_compounds=phase1_dir / "all_compounds.csv",
        phase2_enriched=phase2_dir / "phase2_enriched.csv",
        ranked=None if phase3_dir is None else phase3_dir / "ranked_compounds.csv",
        shortlist=None if phase3_dir is None else phase3_dir / "shortlist_topN_diverse.csv",
        hygiene_audit=phase1_dir / "audit_log.json",
        phase2_audit=phase2_dir / "audit_log.json",
        scoring_audit=None if phase3_dir is None else phase3_dir / "audit_log.json",
        baseline_stats=baseline_stats,
        scaffold_summary=phase2_dir / "scaffold_summary.csv",
        config=config,
    )
