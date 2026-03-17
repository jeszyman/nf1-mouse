"""One-off metadata ingest script.

Reads existing data/metadata.xlsx and 3 WashU source files from inputs/metadata/,
normalizes animal IDs, builds unified dataset, writes updated workbook.

Usage: conda run -n basecamp python scripts/ingest_metadata.py
"""
import re

import pandas as pd

# Canonical model names
SHEET_TO_MODEL = {
    "MN2": "MN2",
    "WU225": "WU-225",
    "WU487": "WU-487",
    "JH-2-055": "JH-2-055",
}

TREATMENT_MAP = {
    "Vehicle": "vehicle",
    "A 0.3mg/kg qdx5days": "mirdametinib",
    "VS 50mg/kg bidx5days": "volasertib",
    "A+VS": "mirdametinib_volasertib",
    "untreated": "untreated",
}


def parse_model_from_plasma_bank(sample_id: str) -> str:
    """Extract canonical model name from a Plasma Bank sample ID."""
    s = sample_id.strip()
    if s.startswith("No PDX"):
        return "No PDX"
    if s.startswith("WU-487"):
        return "WU-487"
    if s.startswith("WU-225"):
        return "WU-225"
    if s.startswith("MN2"):
        return "MN2"
    if s.startswith("2-055"):
        return "JH-2-055"
    raise ValueError(f"Unknown model in Plasma Bank ID: {sample_id!r}")


def normalize_animal_id(raw_id: str, source: str, model_sheet: str = None) -> str:
    """Normalize an animal ID to canonical form: '<MODEL> <group>-<animal>'.

    Args:
        raw_id: Raw animal ID string from source.
        source: One of 'plasma_bank', 'tumor_vol', 'existing'.
        model_sheet: Sheet name (required for tumor_vol source).

    Returns:
        Normalized ID string, e.g. 'WU-487 2-1'.
    """
    s = raw_id.strip()

    if source == "plasma_bank":
        model = parse_model_from_plasma_bank(s)
        if model == "No PDX":
            return s
        if model == "JH-2-055":
            animal_part = s[len("2-055 "):]
        else:
            animal_part = s[len(model) + 1:]
        return f"{model} {animal_part}"

    elif source == "tumor_vol":
        if model_sheet is None:
            raise ValueError("model_sheet required for tumor_vol source")
        model = SHEET_TO_MODEL[model_sheet]
        animal_part = s.lstrip("#")
        return f"{model} {animal_part}"

    elif source == "existing":
        s = re.sub(r"^WU487\b", "WU-487", s)
        s = re.sub(r"^WU225\b", "WU-225", s)
        return s

    else:
        raise ValueError(f"Unknown source: {source!r}")
