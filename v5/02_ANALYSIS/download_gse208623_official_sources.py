"""Resumable acquisition of official GSE208623 beta and GPL21145 annotation files."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "01_DATA_AUDIT" / "gse208623"
SOURCE_MANIFEST_PATH = INPUT_DIR / "gse208623_official_source_manifest.json"
BETA_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208623/suppl/GSE208623_Norm_beta.csv.gz"
ANNOTATION_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
SOURCES = {
    "beta": (BETA_URL, INPUT_DIR / "GSE208623_Norm_beta.csv.gz", 363_136_110),
    "annotation": (ANNOTATION_URL, INPUT_DIR / "GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz", 169_013_044),
}


def parse_content_range(value: str) -> tuple[int, int, int]:
    """Parse one RFC 7233 byte Content-Range response header."""
    prefix, separator, span = value.strip().partition(" ")
    range_text, slash, total_text = span.partition("/")
    start_text, dash, end_text = range_text.partition("-")
    if prefix.casefold() != "bytes" or not separator or not slash or not dash:
        raise ValueError(f"invalid Content-Range header: {value!r}")
    try:
        start, end, total = int(start_text), int(end_text), int(total_text)
    except ValueError as error:
        raise ValueError(f"invalid Content-Range header: {value!r}") from error
    if start < 0 or end < start or total <= end:
        raise ValueError(f"invalid Content-Range header: {value!r}")
    return start, end, total


def next_range(offset: int, expected_bytes: int, chunk_bytes: int) -> tuple[int, int] | None:
    """Return the next inclusive byte span, or None once the file is complete."""
    if offset < 0 or expected_bytes <= 0 or chunk_bytes <= 0:
        raise ValueError("offset, expected_bytes and chunk_bytes must be positive where applicable")
    if offset > expected_bytes:
        raise ValueError("partial file exceeds official expected byte count")
    if offset == expected_bytes:
        return None
    return offset, min(offset + chunk_bytes - 1, expected_bytes - 1)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_completed_gzip(path: Path, expected_bytes: int) -> dict[str, object]:
    if not path.exists() or path.stat().st_size != expected_bytes:
        raise ValueError(
            f"official source byte count mismatch for {path.name}: "
            f"expected {expected_bytes}, observed {path.stat().st_size if path.exists() else 0}"
        )
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
            header = handle.readline().rstrip("\n\r")
    except (OSError, UnicodeError) as error:
        raise ValueError(f"official source is not a valid readable gzip file: {path.name}") from error
    if not header:
        raise ValueError(f"official source gzip header is empty: {path.name}")
    return {"path": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path), "header": header}


def write_completed_source_manifest() -> dict[str, object]:
    """Lock the two fully downloaded official sources without storing sample-level data."""
    if SOURCE_MANIFEST_PATH.exists():
        raise ValueError(f"official GSE208623 source manifest is immutable: {SOURCE_MANIFEST_PATH.name}")
    records: dict[str, dict[str, object]] = {}
    for logical_name, (url, path, expected_bytes) in SOURCES.items():
        completed = validate_completed_gzip(path, expected_bytes)
        records[logical_name] = {
            "logical_name": logical_name,
            "url": url,
            "expected_bytes": expected_bytes,
            "sha256": completed["sha256"],
            "header": completed["header"],
        }
    manifest = {"status": "COMPLETE", "sources": records}
    SOURCE_MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = SOURCE_MANIFEST_PATH.with_suffix(SOURCE_MANIFEST_PATH.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(SOURCE_MANIFEST_PATH)
    return manifest


def download_resumable(
    url: str,
    destination: Path,
    expected_bytes: int,
    chunk_bytes: int = 8 * 1024 * 1024,
    max_chunks: int | None = None,
) -> dict[str, object]:
    """Resume a verified byte-range download, atomically promoting a complete file."""
    if chunk_bytes <= 0:
        raise ValueError("chunk_bytes must be positive")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size == expected_bytes:
        return validate_completed_gzip(destination, expected_bytes)
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists() and partial.stat().st_size > expected_bytes:
        raise ValueError(f"partial file exceeds official expected byte count: {partial.name}")
    chunks_downloaded = 0
    while not partial.exists() or partial.stat().st_size < expected_bytes:
        if max_chunks is not None and chunks_downloaded >= max_chunks:
            break
        offset = partial.stat().st_size if partial.exists() else 0
        next_span = next_range(offset, expected_bytes, chunk_bytes)
        if next_span is None:
            break
        _, end = next_span
        request = Request(url, headers={"Range": f"bytes={offset}-{end}", "User-Agent": "MCI-public-reliability-audit/1.0"})
        try:
            with urlopen(request, timeout=120) as response:
                content_range = response.headers.get("Content-Range", "")
                response_start, response_end, response_total = parse_content_range(content_range)
                if response.status != 206 or response_start != offset or response_end != end or response_total != expected_bytes:
                    raise ValueError(f"unexpected official range response for {destination.name}: {content_range!r}")
                payload = response.read()
        except HTTPError as error:
            raise ValueError(f"official download request failed for {destination.name}: HTTP {error.code}") from error
        if len(payload) != end - offset + 1:
            raise ValueError(f"official range payload length mismatch for {destination.name}")
        with partial.open("ab") as handle:
            handle.write(payload)
        chunks_downloaded += 1
    observed_bytes = partial.stat().st_size if partial.exists() else 0
    if observed_bytes == expected_bytes:
        partial.replace(destination)
        return {"status": "COMPLETE", "chunks_downloaded": chunks_downloaded, **validate_completed_gzip(destination, expected_bytes)}
    return {
        "status": "PARTIAL",
        "chunks_downloaded": chunks_downloaded,
        "path": partial.name,
        "bytes": observed_bytes,
        "expected_bytes": expected_bytes,
    }


def main(selected: list[str] | None = None, chunk_bytes: int = 8 * 1024 * 1024, max_chunks: int | None = None) -> int:
    names = selected or list(SOURCES)
    unknown = set(names) - set(SOURCES)
    if unknown:
        raise ValueError(f"unknown GSE208623 source selection: {sorted(unknown)}")
    results = {}
    for name in names:
        url, path, expected_bytes = SOURCES[name]
        results[name] = {"url": url, **download_resumable(url, path, expected_bytes, chunk_bytes, max_chunks)}
    if not SOURCE_MANIFEST_PATH.exists() and all(path.exists() and path.stat().st_size == expected_bytes for _, path, expected_bytes in SOURCES.values()):
        results["source_manifest"] = write_completed_source_manifest()
    print(json.dumps(results, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", action="append", choices=sorted(SOURCES))
    parser.add_argument("--chunk-bytes", type=int, default=8 * 1024 * 1024)
    parser.add_argument("--max-chunks", type=int)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.source, arguments.chunk_bytes, arguments.max_chunks))
