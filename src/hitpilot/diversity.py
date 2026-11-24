"""Bemis-Murcko scaffolds, Butina clustering, and sphere-exclusion selection."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

from hitpilot.fingerprints import (
    fingerprint,
    morgan_generator,
    tanimoto_distance_matrix,
)
from hitpilot.io_utils import build_audit, load_compound_frame, write_audit, write_csv
from hitpilot.models import DiversityConfig, DiversityRecord, DiversitySummary


def murcko_scaffolds(mol: Chem.Mol) -> tuple[str, str, bool]:
    """Return (murcko_smiles, generic_smiles, acyclic) for one molecule.

    Acyclic molecules (no rings) have empty scaffold SMILES and ``acyclic=True``.
    """
    if mol is None or mol.GetNumAtoms() == 0:
        return "", "", True
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return "", "", True
    murcko = Chem.MolToSmiles(scaffold, isomericSmiles=True, canonical=True)
    try:
        generic_mol = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    except Exception:  # noqa: BLE001
        generic_mol = None
    generic = ""
    if generic_mol is not None and generic_mol.GetNumAtoms() > 0:
        generic = Chem.MolToSmiles(generic_mol, isomericSmiles=True, canonical=True)
    return murcko, generic, False


def _parse_score(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "na"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _selection_order(
    rows: list[dict[str, Any]],
    config: DiversityConfig,
) -> list[int]:
    indices = list(range(len(rows)))
    if not config.rank_by_score:
        return indices
    scored: list[tuple[float | None, int]] = []
    for idx, row in enumerate(rows):
        raw = row.get(config.score_column) if config.score_column else row.get("score")
        score = _parse_score(raw)
        scored.append((score, idx))

    def sort_key(item: tuple[float | None, int]) -> tuple[int, float, int]:
        score, idx = item
        if score is None:
            return (1, 0.0, idx)
        # Ranked compounds first; direction depends on score_ascending.
        directed = score if config.score_ascending else -score
        return (0, directed, idx)

    return [idx for _, idx in sorted(scored, key=sort_key)]


def butina_clusters(fps: list, tanimoto_cutoff: float) -> list[tuple[int, ...]]:
    """Cluster fingerprints with Butina; *tanimoto_cutoff* is a similarity cutoff."""
    n_pts = len(fps)
    if n_pts == 0:
        return []
    if n_pts == 1:
        return [(0,)]
    dist_cutoff = 1.0 - float(tanimoto_cutoff)
    dists = tanimoto_distance_matrix(fps)
    clusters = Butina.ClusterData(dists, n_pts, dist_cutoff, isDistData=True, reordering=True)
    return [tuple(cluster) for cluster in clusters]


def sphere_exclusion(
    fps: list,
    compound_ids: list[str],
    order: list[int],
    exclusion_similarity: float,
) -> tuple[list[bool], list[int | None], list[str | None]]:
    """Greedy leader-follower selection. Selected pairs all have Tanimoto < threshold."""
    n_pts = len(fps)
    selected = [False] * n_pts
    ranks: list[int | None] = [None] * n_pts
    excluded_by: list[str | None] = [None] * n_pts
    taken = [False] * n_pts
    next_rank = 1
    for idx in order:
        if taken[idx] or fps[idx] is None:
            continue
        selected[idx] = True
        ranks[idx] = next_rank
        next_rank += 1
        taken[idx] = True
        # Exclude remaining neighbors of this leader.
        neighbors = [
            (j, fp) for j, fp in enumerate(fps) if not taken[j] and fp is not None
        ]
        if not neighbors:
            continue
        neighbor_idx = [j for j, _ in neighbors]
        neighbor_fps = [fp for _, fp in neighbors]
        sims = DataStructs.BulkTanimotoSimilarity(fps[idx], neighbor_fps)
        for j, sim in zip(neighbor_idx, sims, strict=True):
            if float(sim) >= exclusion_similarity:
                taken[j] = True
                excluded_by[j] = compound_ids[idx]
    return selected, ranks, excluded_by


def run_diversity(
    rows: Iterable[dict[str, Any]],
    config: DiversityConfig | None = None,
    *,
    smiles_field: str = "smiles",
    id_field: str = "compound_id",
    label_field: str = "label",
) -> tuple[list[DiversityRecord], DiversitySummary]:
    config = config or DiversityConfig()
    row_list = list(rows)
    generator = morgan_generator(radius=config.radius, fp_size=config.fp_size)

    records: list[DiversityRecord] = []
    valid_indices: list[int] = []
    valid_fps = []
    valid_ids: list[str] = []
    murcko_list: list[str] = []

    for index, row in enumerate(row_list):
        compound_id = str(row.get(id_field) or row.get("compound_id") or f"cmpd_{index + 1:06d}")
        smiles = str(row.get(smiles_field, row.get("smiles", "")) or "")
        raw_label = row.get(label_field, row.get("label"))
        label = None if raw_label is None or str(raw_label).strip() == "" else str(raw_label)
        mol = Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
        if mol is None:
            records.append(
                DiversityRecord(
                    compound_id=compound_id,
                    smiles=smiles,
                    label=label,
                    error="invalid SMILES" if smiles.strip() else "empty SMILES",
                )
            )
            continue
        murcko, generic, acyclic = murcko_scaffolds(mol)
        rec = DiversityRecord(
            compound_id=compound_id,
            smiles=smiles,
            murcko_scaffold_smiles=murcko,
            generic_scaffold_smiles=generic,
            acyclic=acyclic,
            label=label,
        )
        records.append(rec)
        valid_indices.append(index)
        valid_fps.append(fingerprint(mol, generator))
        valid_ids.append(compound_id)
        murcko_list.append(murcko)

    frequencies = Counter(murcko_list)
    for rec in records:
        if rec.error is None:
            rec.scaffold_frequency = int(frequencies.get(rec.murcko_scaffold_smiles, 0))

    # Butina clustering on valid fingerprints only.
    clusters = butina_clusters(valid_fps, config.tanimoto_cutoff)
    cluster_of: dict[int, tuple[int, int, bool]] = {}
    for cluster_id, members in enumerate(clusters, start=1):
        size = len(members)
        for position, member in enumerate(members):
            cluster_of[member] = (cluster_id, size, position == 0)
    for local_idx, rec_idx in enumerate(valid_indices):
        assignment = cluster_of.get(local_idx)
        if assignment is None:
            continue
        cluster_id, size, is_centroid = assignment
        records[rec_idx].cluster_id = cluster_id
        records[rec_idx].cluster_size = size
        records[rec_idx].is_cluster_centroid = is_centroid

    # Sphere exclusion uses the same valid subset, mapped back to full rows.
    order_full = _selection_order(row_list, config)
    # Restrict order to valid local indices, preserving the requested ranking.
    rec_index_by_id_row = {row_i: local for local, row_i in enumerate(valid_indices)}
    order_local = [rec_index_by_id_row[i] for i in order_full if i in rec_index_by_id_row]
    selected, ranks, excluded_by = sphere_exclusion(
        valid_fps,
        valid_ids,
        order_local,
        config.exclusion_similarity,
    )
    for local_idx, rec_idx in enumerate(valid_indices):
        records[rec_idx].sphere_selected = selected[local_idx]
        records[rec_idx].sphere_rank = ranks[local_idx]
        records[rec_idx].excluded_by = excluded_by[local_idx]

    return records, summarize_diversity(records, config)


def summarize_diversity(
    records: list[DiversityRecord],
    config: DiversityConfig,
) -> DiversitySummary:
    valid = [rec for rec in records if rec.error is None]
    murcko = {rec.murcko_scaffold_smiles for rec in valid if rec.murcko_scaffold_smiles}
    generic = {rec.generic_scaffold_smiles for rec in valid if rec.generic_scaffold_smiles}
    sizes = sorted(
        {rec.cluster_id: rec.cluster_size for rec in valid if rec.cluster_id is not None}.values(),
        reverse=True,
    )
    n_selected = sum(1 for rec in valid if rec.sphere_selected)
    actives = [rec for rec in valid if rec.label == "active"]
    n_actives_total = len(actives) if any(rec.label for rec in records) else None
    n_actives_selected = (
        sum(1 for rec in actives if rec.sphere_selected) if n_actives_total is not None else None
    )
    n_clusters = len({rec.cluster_id for rec in valid if rec.cluster_id is not None})
    return DiversitySummary(
        n_total=len(records),
        n_valid=len(valid),
        n_invalid=len(records) - len(valid),
        n_unique_murcko_scaffolds=len(murcko),
        n_unique_generic_scaffolds=len(generic),
        n_acyclic=sum(1 for rec in valid if rec.acyclic),
        n_clusters=n_clusters,
        top_cluster_sizes=sizes[:10],
        n_sphere_selected=n_selected,
        n_actives_selected=n_actives_selected,
        n_actives_total=n_actives_total,
        tanimoto_cutoff=config.tanimoto_cutoff,
        exclusion_similarity=config.exclusion_similarity,
    )


def scaffold_summary_rows(records: list[DiversityRecord]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    for rec in records:
        if rec.error is not None:
            continue
        key = (rec.murcko_scaffold_smiles, rec.generic_scaffold_smiles)
        bucket = grouped.setdefault(
            key,
            {
                "murcko_scaffold_smiles": rec.murcko_scaffold_smiles,
                "generic_scaffold_smiles": rec.generic_scaffold_smiles,
                "acyclic": rec.acyclic,
                "count": 0,
                "n_active": 0,
                "n_inactive": 0,
                "example_compound_ids": [],
            },
        )
        bucket["count"] += 1
        if rec.label == "active":
            bucket["n_active"] += 1
        elif rec.label == "inactive":
            bucket["n_inactive"] += 1
        if len(bucket["example_compound_ids"]) < 5:
            bucket["example_compound_ids"].append(rec.compound_id)
    rows = []
    ordered = sorted(
        grouped.values(),
        key=lambda item: (-item["count"], item["murcko_scaffold_smiles"]),
    )
    for bucket in ordered:
        row = dict(bucket)
        row["example_compound_ids"] = "|".join(row["example_compound_ids"])
        rows.append(row)
    return rows


def run_diversity_file(
    input_path: str | Path,
    output_dir: str | Path,
    config: DiversityConfig | None = None,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
    score_column: str | None = None,
):
    """Cluster and select a diverse subset; write CSVs plus an audit log."""
    config = config or DiversityConfig()
    if score_column and not config.score_column:
        config = config.model_copy(update={"score_column": score_column, "rank_by_score": True})
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = load_compound_frame(
        input_path,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
        score_column=score_column or config.score_column,
    )
    records, summary = run_diversity(frame.to_dict(orient="records"), config)
    all_path = output_dir / "diversity.csv"
    subset_path = output_dir / "diverse_subset.csv"
    scaffold_path = output_dir / "scaffold_summary.csv"
    write_csv(all_path, [rec.to_csv_row() for rec in records])
    write_csv(subset_path, [rec.to_csv_row() for rec in records if rec.sphere_selected])
    write_csv(scaffold_path, scaffold_summary_rows(records))
    manifest = build_audit(
        command="diversity",
        config=config,
        input_path=input_path,
        output_dir=output_dir,
        summary=summary,
        outputs=[
            str(all_path.resolve()),
            str(subset_path.resolve()),
            str(scaffold_path.resolve()),
            str((output_dir / "audit_log.json").resolve()),
        ],
    )
    write_audit(output_dir / "audit_log.json", manifest)
    return records, summary, manifest
