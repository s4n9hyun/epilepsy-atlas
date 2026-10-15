#!/usr/bin/env python3
"""Helpers for timestamped CSV output names."""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path


TIMESTAMP_RE = re.compile(r"_(\d{8}_\d{6})$")


def make_timestamp(value: str | None = None) -> str:
    if value:
        if not re.fullmatch(r"\d{8}_\d{6}", value):
            raise ValueError("timestamp must use YYYYMMDD_HHMMSS format")
        return value
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def timestamp_path(path: Path, timestamp: str | None) -> Path:
    if not timestamp:
        return path
    return path.with_name(f"{path.stem}_{timestamp}{path.suffix}")


def csv_timestamp(path: Path) -> str:
    match = TIMESTAMP_RE.search(path.stem)
    return match.group(1) if match else ""


def resolve_csv_path(directory: Path, stem: str) -> Path | None:
    """Return the newest timestamped CSV for stem, falling back to the legacy name."""
    exact = directory / f"{stem}.csv"
    timestamped = [
        path for path in directory.glob(f"{stem}_*.csv")
        if path.is_file() and csv_timestamp(path)
    ]
    if timestamped:
        return max(timestamped, key=lambda path: (csv_timestamp(path), path.stat().st_mtime))
    if exact.exists():
        return exact
    return None


def resolve_latest_path(path: Path) -> Path | None:
    if path.exists():
        return path
    if path.suffix.lower() == ".csv":
        return resolve_csv_path(path.parent, path.stem)
    return None
