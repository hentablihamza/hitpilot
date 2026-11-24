"""Shared Morgan fingerprint helpers (ECFP4-like, 2048-bit)."""

from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def morgan_generator(radius: int = 2, fp_size: int = 2048):
    return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)


def fingerprint(
    mol: Chem.Mol,
    generator: rdFingerprintGenerator.FingerprintGenerator64 | None = None,
    *,
    radius: int = 2,
    fp_size: int = 2048,
):
    gen = generator or morgan_generator(radius=radius, fp_size=fp_size)
    return gen.GetFingerprint(mol)


def fingerprints_from_mols(
    mols: list[Chem.Mol | None],
    generator: rdFingerprintGenerator.FingerprintGenerator64 | None = None,
    *,
    radius: int = 2,
    fp_size: int = 2048,
) -> list:
    gen = generator or morgan_generator(radius=radius, fp_size=fp_size)
    return [None if mol is None else gen.GetFingerprint(mol) for mol in mols]


def bulk_tanimoto(query_fp, ref_fps: list) -> list[float]:
    if not ref_fps:
        return []
    return list(DataStructs.BulkTanimotoSimilarity(query_fp, ref_fps))


def max_tanimoto(query_fp, ref_fps: list) -> tuple[float, int]:
    """Return (max_similarity, index). Index is -1 when *ref_fps* is empty."""
    if not ref_fps:
        return 0.0, -1
    sims = DataStructs.BulkTanimotoSimilarity(query_fp, ref_fps)
    best_idx = max(range(len(sims)), key=lambda i: sims[i])
    return float(sims[best_idx]), int(best_idx)


def tanimoto_distance_matrix(fps: list) -> list[float]:
    """Condensed lower-triangle distance matrix (1 - Tanimoto) for Butina."""
    dists: list[float] = []
    for i, fp in enumerate(fps):
        if i == 0:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps[:i])
        dists.extend(1.0 - float(sim) for sim in sims)
    return dists
