"""Download and load public LIT-PCBA target datasets."""

from __future__ import annotations

import io
import re
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import pandas as pd
import requests

DEFAULT_BASE_URL = "https://drugdesign.unistra.fr/LIT-PCBA/Files"
DEFAULT_TIMEOUT = 15

# Official LIT-PCBA target names (Tran-Nguyen et al., J. Chem. Inf. Model. 2020)
# mapped to the human protein they represent.
KNOWN_TARGETS: dict[str, dict[str, str]] = {
    "ADRB2": {
        "protein": "ADRB2",
        "query": "ADRB2",
        "description": "Beta-2 adrenergic receptor",
    },
    "ALDH1": {
        "protein": "ALDH1A1",
        "query": "ALDH1A1",
        "description": "Aldehyde dehydrogenase 1A1",
    },
    "ESR1_ago": {
        "protein": "ESR1",
        "query": "ESR1",
        "description": "Estrogen receptor alpha (agonist)",
    },
    "ESR1_ant": {
        "protein": "ESR1",
        "query": "ESR1",
        "description": "Estrogen receptor alpha (antagonist)",
    },
    "FEN1": {"protein": "FEN1", "query": "FEN1", "description": "Flap endonuclease 1"},
    "GBA": {"protein": "GBA", "query": "GBA", "description": "Glucosylceramidase"},
    "IDH1": {
        "protein": "IDH1",
        "query": "IDH1",
        "description": "Isocitrate dehydrogenase 1",
    },
    "KAT2A": {
        "protein": "KAT2A",
        "query": "KAT2A",
        "description": "Histone acetyltransferase KAT2A",
    },
    "MAPK1": {"protein": "MAPK1", "query": "MAPK1", "description": "MAP kinase ERK2"},
    "MTORC1": {
        "protein": "MTOR",
        "query": "MTOR",
        "description": "Serine/threonine-protein kinase mTOR",
    },
    "OPRK1": {
        "protein": "OPRK1",
        "query": "OPRK1",
        "description": "Kappa-type opioid receptor",
    },
    "PKM2": {"protein": "PKM", "query": "PKM", "description": "Pyruvate kinase PKM"},
    "PPARG": {
        "protein": "PPARG",
        "query": "PPARG",
        "description": "Peroxisome proliferator-activated receptor gamma",
    },
    "TP53": {
        "protein": "TP53",
        "query": "TP53",
        "description": "Cellular tumor antigen p53",
    },
    "VDR": {"protein": "VDR", "query": "VDR", "description": "Vitamin D receptor"},
}

# Alternative folder / archive names used by official Unistra releases.
TARGET_ALIASES: dict[str, list[str]] = {
    "ESR1_ago": ["ESR1_ago", "ESR_ago", "ESR-ago"],
    "ESR1_ant": ["ESR1_ant", "ESR_antago", "ESR_ant", "ESR-antago", "ESRantago"],
}

# Legacy per-target tarballs (often 404 now) plus official combined archives
# and a public per-target SMILES directory mirror.
MIRROR_URL_TEMPLATES: list[str] = [
    "https://drugdesign.unistra.fr/LIT-PCBA/Files/{target}.tgz",
    "https://drugdesign.unistra.fr/LIT-PCBA/Files/{alias}.tgz",
    "https://storage.valencelabs.com/datasets-public-research/LIT_PCBA/{target}/",
]

OFFICIAL_COMBINED_ARCHIVES: list[str] = [
    "https://drugdesign.unistra.fr/LIT-PCBA/Files/AVE_unbiased.tgz",
    "https://drugdesign.unistra.fr/LIT-PCBA/Files/full_data.tgz",
    "https://drugdesign.unistra.fr/downloads/datasets/LIT-PCBA_AVE_unbiased.tar.gz",
    "https://drugdesign.unistra.fr/downloads/datasets/LIT-PCBA_full.tar.gz",
]

_COMPLETE_STEMS = {
    "active": {"active", "actives"},
    "inactive": {"inactive", "inactives"},
}


class LitPcbaError(RuntimeError):
    """Raised when a LIT-PCBA target cannot be downloaded or parsed."""


def official_url(target: str) -> str:
    return f"{DEFAULT_BASE_URL}/{normalize_target(target)}.tgz"


def _alias_lookup() -> dict[str, str]:
    lookup = {key.lower(): key for key in KNOWN_TARGETS}
    for canonical, aliases in TARGET_ALIASES.items():
        for alias in aliases:
            lookup[alias.lower()] = canonical
    return lookup


def normalize_target(target: str) -> str:
    name = target.strip()
    return _alias_lookup().get(name.lower(), name)


def target_aliases(target: str) -> list[str]:
    name = normalize_target(target)
    aliases = TARGET_ALIASES.get(name, [name])
    if name not in aliases:
        aliases = [name, *aliases]
    return aliases


def candidate_urls(target: str, override: str | None = None) -> list[str]:
    if override:
        return [override]
    name = normalize_target(target)
    urls: list[str] = []
    for template in MIRROR_URL_TEMPLATES:
        if "{alias}" in template:
            for alias in target_aliases(name):
                url = template.format(target=name, alias=alias)
                if url not in urls:
                    urls.append(url)
        else:
            url = template.format(target=name, alias=name)
            if url not in urls:
                urls.append(url)
    return urls


def _is_probably_archive(content: bytes) -> bool:
    if len(content) < 4:
        return False
    return content.startswith(b"\x1f\x8b") or content.startswith(b"PK") or b"ustar" in content[:512]


def _is_probably_html(content: bytes) -> bool:
    head = content[:256].lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html")


def download_bytes(url: str, timeout: int = DEFAULT_TIMEOUT) -> bytes:
    """Download a URL with a short connect timeout and streamed reads."""
    response = requests.get(url, timeout=(timeout, max(timeout, 60)), stream=True)
    response.raise_for_status()
    chunks: list[bytes] = []
    for chunk in response.iter_content(chunk_size=1024 * 256):
        if chunk:
            chunks.append(chunk)
    return b"".join(chunks)


def _member_matches_target(name: str, aliases: set[str]) -> bool:
    parts = Path(name).parts
    return any(part.lower() in aliases for part in parts)


def _extract_archive(data: bytes, dest: Path, target: str | None = None) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    aliases = {alias.lower() for alias in target_aliases(target)} if target else set()
    if data.startswith(b"PK"):
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            names = archive.namelist()
            if aliases:
                selected = [name for name in names if _member_matches_target(name, aliases)]
            else:
                selected = names
            if selected:
                for name in selected:
                    archive.extract(name, dest)
                return
            archive.extractall(dest)
        return
    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:*") as archive:
            members = archive.getmembers()
            if aliases:
                selected = [m for m in members if _member_matches_target(m.name, aliases)]
                if selected:
                    archive.extractall(dest, members=selected, filter="data")
                    return
            archive.extractall(dest, filter="data")
            return
    except tarfile.TarError as exc:
        raise LitPcbaError(f"Could not extract LIT-PCBA archive: {exc}") from exc


def _iter_smiles_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() in {".smi", ".smiles", ".csv", ".tsv", ".txt", ".sdf"}:
            files.append(path)
    return files


def _classify_file(path: Path) -> str | None:
    text = str(path).lower()
    name = path.name.lower()
    if "inactive" in text or name.startswith("inact"):
        return "inactive"
    if "active" in text or name.startswith("act"):
        return "active"
    return None


def _is_split_file(path: Path) -> bool:
    stem = path.stem.lower()
    return stem.endswith(("_t", "_v", "_train", "_val", "_validation"))


def select_source_files(files: list[Path], kind: str) -> list[Path]:
    """Prefer complete active/inactive files over train/validation splits."""
    classified = [path for path in files if _classify_file(path) == kind]
    complete = [path for path in classified if path.stem.lower() in _COMPLETE_STEMS[kind]]
    if complete:
        return complete
    splits = [path for path in classified if _is_split_file(path)]
    if splits:
        return splits
    return classified


def _scope_root(target_dir: Path, target: str | None) -> Path:
    if not target:
        return target_dir
    aliases = {alias.lower() for alias in target_aliases(target)}
    matches = [
        path
        for path in target_dir.rglob("*")
        if path.is_dir() and path.name.lower() in aliases
    ]
    if not matches:
        return target_dir
    return sorted(matches, key=lambda path: (len(path.parts), str(path)))[0]


def _load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".sdf":
        from rdkit import Chem

        rows: list[dict[str, str]] = []
        supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        for index, mol in enumerate(supplier):
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol) if mol.GetNumAtoms() else ""
            props = {key: str(mol.GetProp(key)) for key in mol.GetPropNames()}
            compound_id = (
                props.get("compound_id")
                or props.get("ID")
                or props.get("id")
                or props.get("_Name")
            )
            if not compound_id and mol.HasProp("_Name"):
                compound_id = mol.GetProp("_Name")
            if not compound_id:
                compound_id = f"{path.stem}_{index + 1}"
            if not smiles and mol.HasProp("SMILES"):
                smiles = mol.GetProp("SMILES")
            rows.append({"compound_id": str(compound_id), "smiles": smiles})
        return pd.DataFrame(rows)

    from hitpilot.hygiene import read_compound_table

    return read_compound_table(path)


def load_litpcba(target_dir: str | Path, target: str | None = None) -> pd.DataFrame:
    """Load actives and inactives from an extracted LIT-PCBA target directory."""
    target_dir = Path(target_dir)
    if not target_dir.is_dir():
        raise FileNotFoundError(target_dir)

    root = _scope_root(target_dir, target)
    files = _iter_smiles_files(root)
    grouped: dict[str, list[pd.DataFrame]] = {"active": [], "inactive": []}
    for kind in ("active", "inactive"):
        for path in select_source_files(files, kind):
            table = _load_table(path)
            if table.empty or "smiles" not in table.columns:
                continue
            table = table.copy()
            table["label"] = kind
            table["source_file"] = str(path.relative_to(target_dir))
            grouped[kind].append(table)

    frames = grouped["active"] + grouped["inactive"]
    if not frames:
        raise LitPcbaError(
            f"No LIT-PCBA actives/inactives SMILES files found under {target_dir}"
        )
    combined = pd.concat(frames, ignore_index=True)
    if target:
        combined["target"] = normalize_target(target)
    if "compound_id" not in combined.columns:
        combined["compound_id"] = [f"cmpd_{i + 1:06d}" for i in range(len(combined))]
    combined["compound_id"] = combined["compound_id"].astype(str)
    duplicated = combined["compound_id"].duplicated(keep=False)
    if duplicated.any():
        combined.loc[duplicated, "compound_id"] = [
            f"{cid}_{label}"
            for cid, label in zip(
                combined.loc[duplicated, "compound_id"],
                combined.loc[duplicated, "label"],
                strict=True,
            )
        ]
    columns = ["compound_id", "smiles", "label"]
    extras = [col for col in combined.columns if col not in columns]
    return combined[columns + extras]


def write_unified_csv(frame: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = frame[["compound_id", "smiles", "label"]].copy()
    out.to_csv(path, index=False)
    return path


def unified_csv_path(target_dir: Path, target: str) -> Path:
    slug = normalize_target(target).lower()
    return target_dir / f"{slug}_litpcba.csv"


def _parse_directory_links(page_url: str, html: str) -> list[str]:
    hrefs = re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.I)
    urls: list[str] = []
    for href in hrefs:
        if href.endswith("/"):
            continue
        lower = href.lower()
        exts = (".smi", ".smiles", ".csv", ".tsv", ".txt", ".sdf")
        if not any(lower.endswith(ext) for ext in exts):
            continue
        urls.append(urljoin(page_url if page_url.endswith("/") else page_url + "/", href))
    return urls


def _fetch_smiles_directory(url: str, dest: Path, timeout: int) -> list[Path]:
    """Download active/inactive SMILES files listed on a public directory page."""
    page = requests.get(url, timeout=(timeout, max(timeout, 30)))
    page.raise_for_status()
    links = _parse_directory_links(url, page.text)
    saved: list[Path] = []
    dest.mkdir(parents=True, exist_ok=True)
    for link in links:
        kind = _classify_file(Path(link))
        if kind is None:
            continue
        data = download_bytes(link, timeout=timeout)
        if _is_probably_html(data):
            continue
        path = dest / Path(link).name
        path.write_bytes(data)
        saved.append(path)
    if not saved:
        raise LitPcbaError(f"No SMILES files found at directory {url}")
    return saved


def fetch_litpcba(
    target: str,
    dest: str | Path = "data/litpcba",
    *,
    url: str | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    session_urls: list[str] | None = None,
) -> dict[str, Any]:
    """Download one LIT-PCBA target archive, extract it, and write a unified CSV.

    Returns a metadata dictionary with paths, URL used, and active/inactive counts.
    """
    target_name = normalize_target(target)
    dest = Path(dest)
    target_dir = dest / target_name
    target_dir.mkdir(parents=True, exist_ok=True)
    extract_dir = target_dir / "raw"

    urls = session_urls or candidate_urls(target_name, override=url)
    last_error: Exception | None = None
    used_url: str | None = None
    archive_path = target_dir / f"{target_name}.tgz"

    for candidate in urls:
        try:
            if candidate.endswith("/"):
                if extract_dir.exists():
                    shutil.rmtree(extract_dir)
                extract_dir.mkdir(parents=True, exist_ok=True)
                _fetch_smiles_directory(candidate, extract_dir, timeout=timeout)
                used_url = candidate
                break
            data = download_bytes(candidate, timeout=timeout)
            if _is_probably_html(data) or not _is_probably_archive(data):
                preview = data[:200].decode("utf-8", errors="replace")
                raise LitPcbaError(
                    f"{candidate} did not return an archive (got {preview!r})"
                )
            archive_path.write_bytes(data)
            if extract_dir.exists():
                shutil.rmtree(extract_dir)
            extract_dir.mkdir(parents=True, exist_ok=True)
            _extract_archive(data, extract_dir, target=target_name)
            used_url = candidate
            break
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            continue

    if used_url is None:
        raise LitPcbaError(
            f"Failed to download LIT-PCBA target {target_name} from {urls}: {last_error}"
        )

    frame = load_litpcba(target_dir, target=target_name)
    csv_path = write_unified_csv(frame, unified_csv_path(target_dir, target_name))
    n_active = int((frame["label"] == "active").sum())
    n_inactive = int((frame["label"] == "inactive").sum())
    return {
        "target": target_name,
        "url": used_url,
        "target_dir": str(target_dir.resolve()),
        "archive_path": str(archive_path.resolve()),
        "csv_path": str(csv_path.resolve()),
        "n_compounds": int(len(frame)),
        "n_active": n_active,
        "n_inactive": n_inactive,
        "protein": KNOWN_TARGETS.get(target_name, {}).get("protein"),
        "description": KNOWN_TARGETS.get(target_name, {}).get("description"),
    }
