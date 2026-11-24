"""Tests for the public LIT-PCBA downloader and loader."""

from __future__ import annotations

import io
import tarfile
from pathlib import Path

import pandas as pd
import pytest

from hitpilot.data.litpcba import (
    KNOWN_TARGETS,
    LitPcbaError,
    candidate_urls,
    fetch_litpcba,
    load_litpcba,
    official_url,
    write_unified_csv,
)

IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CATECHOL = "Oc1ccccc1O"


def _make_target_tarball(tmp_path: Path) -> bytes:
    actives = tmp_path / "actives.smi"
    inactives = tmp_path / "inactives.smi"
    actives.write_text(
        f"{IBUPROFEN} CHEMBL521\nCC(=O)Oc1ccccc1C(=O)O CHEMBL25\n",
        encoding="utf-8",
    )
    inactives.write_text(
        f"{CATECHOL} DECOY-1\nCCO DECOY-2\n",
        encoding="utf-8",
    )
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        tar.add(actives, arcname="PPARG/SMILES/actives.smi")
        tar.add(inactives, arcname="PPARG/SMILES/inactives.smi")
    return buffer.getvalue()


def test_known_targets_and_official_urls():
    assert "ESR1_ant" in KNOWN_TARGETS
    assert "PPARG" in KNOWN_TARGETS
    assert "MAPK1" in KNOWN_TARGETS
    assert official_url("esr1_ant") == "https://drugdesign.unistra.fr/LIT-PCBA/Files/ESR1_ant.tgz"
    urls = candidate_urls("PPARG")
    assert urls[0].endswith("/PPARG.tgz")


def test_load_litpcba_from_extracted_smi(tmp_path: Path):
    target_dir = tmp_path / "ESR1_ant" / "raw"
    smiles_dir = target_dir / "ESR1_ant"
    smiles_dir.mkdir(parents=True)
    (smiles_dir / "actives.smi").write_text(f"{IBUPROFEN} PUBCHEM-A\n", encoding="utf-8")
    (smiles_dir / "inactives.smi").write_text(
        f"{CATECHOL} PUBCHEM-I1\nCCO PUBCHEM-I2\n",
        encoding="utf-8",
    )
    frame = load_litpcba(tmp_path / "ESR1_ant", target="ESR1_ant")
    assert list(frame.columns[:3]) == ["compound_id", "smiles", "label"]
    assert set(frame["label"]) == {"active", "inactive"}
    assert int((frame["label"] == "active").sum()) == 1
    assert int((frame["label"] == "inactive").sum()) == 2
    assert "PUBCHEM-A" in set(frame["compound_id"])
    assert IBUPROFEN in set(frame["smiles"])


def test_load_litpcba_from_csv_and_unified_write(tmp_path: Path):
    target_dir = tmp_path / "MAPK1"
    target_dir.mkdir()
    (target_dir / "actives.csv").write_text(
        f"compound_id,smiles\nA1,{IBUPROFEN}\n",
        encoding="utf-8",
    )
    (target_dir / "inactives.tsv").write_text(
        f"id\tsmiles\nI1\t{CATECHOL}\n",
        encoding="utf-8",
    )
    frame = load_litpcba(target_dir, target="MAPK1")
    path = write_unified_csv(frame, target_dir / "mapk1_litpcba.csv")
    loaded = pd.read_csv(path)
    assert len(loaded) == 2
    assert set(loaded["label"]) == {"active", "inactive"}


def test_fetch_litpcba_uses_archive_and_writes_csv(tmp_path: Path, monkeypatch):
    from hitpilot.data import litpcba as litpcba_mod

    payload = _make_target_tarball(tmp_path)

    def fake_download(url: str, timeout: int = 15) -> bytes:
        assert url.endswith("PPARG.tgz")
        return payload

    monkeypatch.setattr(litpcba_mod, "download_bytes", fake_download)
    dest = tmp_path / "data" / "litpcba"
    result = fetch_litpcba("PPARG", dest=dest)
    assert result["target"] == "PPARG"
    assert result["n_active"] == 2
    assert result["n_inactive"] == 2
    assert result["n_compounds"] == 4
    csv_path = Path(result["csv_path"])
    assert csv_path.name == "pparg_litpcba.csv"
    frame = pd.read_csv(csv_path)
    assert set(frame["compound_id"]) >= {"CHEMBL521", "CHEMBL25", "DECOY-1", "DECOY-2"}


def test_fetch_litpcba_raises_on_bad_payload(tmp_path: Path, monkeypatch):
    from hitpilot.data import litpcba as litpcba_mod

    monkeypatch.setattr(litpcba_mod, "download_bytes", lambda url, timeout=15: b"<html>nope</html>")
    with pytest.raises(LitPcbaError):
        fetch_litpcba("ESR1_ant", dest=tmp_path / "dest")


def test_normalize_aliases_and_prefer_complete_over_splits(tmp_path: Path):
    from hitpilot.data.litpcba import load_litpcba, normalize_target

    assert normalize_target("ESR_antago") == "ESR1_ant"
    target_dir = tmp_path / "ESR1_ant"
    smiles_dir = target_dir / "ESR_antago"
    smiles_dir.mkdir(parents=True)
    (smiles_dir / "active.smi").write_text(f"{IBUPROFEN} FULL-A\n", encoding="utf-8")
    (smiles_dir / "inactive.smi").write_text(f"{CATECHOL} FULL-I\n", encoding="utf-8")
    (smiles_dir / "active_T.smi").write_text(f"{IBUPROFEN} SPLIT-A\n", encoding="utf-8")
    (smiles_dir / "inactive_T.smi").write_text(f"{CATECHOL} SPLIT-I\n", encoding="utf-8")
    frame = load_litpcba(target_dir, target="ESR1_ant")
    assert len(frame) == 2
    assert set(frame["compound_id"]) == {"FULL-A", "FULL-I"}
