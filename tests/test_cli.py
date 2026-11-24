"""CLI tests for fetch-litpcba and hygiene commands."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from hitpilot.cli import app
from hitpilot.data.litpcba import write_unified_csv

runner = CliRunner()

IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CATECHOL = "Oc1ccccc1O"


def test_hygiene_cli_writes_outputs_and_summary(tmp_path: Path):
    input_path = tmp_path / "compounds.csv"
    input_path.write_text(
        "compound_id,smiles,label\n"
        f"ibu,{IBUPROFEN},active\n"
        f"cat,{CATECHOL},inactive\n"
        "bad,not_a_smiles,inactive\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "phase1"
    result = runner.invoke(
        app,
        ["hygiene", str(input_path), "--output-dir", str(output_dir)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "Hygiene summary" in result.output
    assert "Drop reasons" in result.output
    assert (output_dir / "passed_compounds.csv").is_file()
    assert (output_dir / "all_compounds.csv").is_file()
    audit = json.loads((output_dir / "audit_log.json").read_text(encoding="utf-8"))
    assert audit["summary"]["n_total"] == 3
    assert audit["summary"]["n_passed"] == 1
    assert audit["summary"]["counts_by_drop_reason"]["invalid_smiles"] == 1
    assert audit["summary"]["counts_by_drop_reason"]["pains_alert"] == 1


def test_hygiene_cli_flag_only_nih_and_brenk(tmp_path: Path):
    input_path = tmp_path / "aldehyde.smi"
    input_path.write_text("O=Cc1ccccc1 benzaldehyde\n", encoding="utf-8")
    output_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "hygiene",
            str(input_path),
            "--output-dir",
            str(output_dir),
            "--no-drop-brenk",
            "--no-drop-nih",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    audit = json.loads((output_dir / "audit_log.json").read_text(encoding="utf-8"))
    assert audit["config"]["drop_brenk"] is False
    assert audit["config"]["drop_nih"] is False
    assert audit["summary"]["n_passed"] == 1


def test_fetch_litpcba_cli_from_local_archive(tmp_path: Path, monkeypatch):
    from hitpilot.data import litpcba as litpcba_mod

    raw = tmp_path / "raw_source"
    actives = raw / "actives.smi"
    inactives = raw / "inactives.smi"
    raw.mkdir()
    actives.write_text(f"{IBUPROFEN} PUBCHEM-1\n", encoding="utf-8")
    inactives.write_text(f"{CATECHOL} PUBCHEM-2\nCCO PUBCHEM-3\n", encoding="utf-8")

    import io
    import tarfile

    archive_bytes = io.BytesIO()
    with tarfile.open(fileobj=archive_bytes, mode="w:gz") as tar:
        tar.add(actives, arcname="ESR1_ant/actives.smi")
        tar.add(inactives, arcname="ESR1_ant/inactives.smi")
    payload = archive_bytes.getvalue()

    def fake_download(url: str, timeout: int = 15) -> bytes:
        assert "ESR1_ant" in url
        return payload

    monkeypatch.setattr(litpcba_mod, "download_bytes", fake_download)
    dest = tmp_path / "data" / "litpcba"
    result = runner.invoke(
        app,
        ["fetch-litpcba", "ESR1_ant", "--dest", str(dest)],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "ESR1_ant" in result.output
    csv_path = dest / "ESR1_ant" / "esr1_ant_litpcba.csv"
    assert csv_path.is_file()
    frame = pd.read_csv(csv_path)
    assert set(frame.columns) >= {"compound_id", "smiles", "label"}
    assert set(frame["label"]) == {"active", "inactive"}
    assert len(frame) == 3
    assert int((frame["label"] == "active").sum()) == 1


def test_write_unified_csv_helper(tmp_path: Path):
    frame = pd.DataFrame(
        {
            "compound_id": ["a", "b"],
            "smiles": [IBUPROFEN, CATECHOL],
            "label": ["active", "inactive"],
            "extra": [1, 2],
        }
    )
    path = write_unified_csv(frame, tmp_path / "out.csv")
    loaded = pd.read_csv(path)
    assert list(loaded.columns) == ["compound_id", "smiles", "label"]
