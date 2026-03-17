"""One-off metadata ingest script.

Reads existing data/metadata.xlsx and 3 WashU source files from inputs/metadata/,
normalizes animal IDs, builds unified dataset, writes updated workbook.

Usage: conda run -n basecamp python scripts/ingest_metadata.py
"""
import re
from pathlib import Path

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


# ---------------------------------------------------------------------------
# Source file readers
# ---------------------------------------------------------------------------

def read_plasma_bank(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read Murine Plasma Sample Bank xlsx.

    Returns:
        (animals_df, samples_df):
          - animals_df: unique animals with normalized_id, model, treatment_raw, tumor_injection_date
          - samples_df: one row per bleed with normalized_id, collection_date, bleed_type, notes, cfdna_conc_ng_ul
    """
    df = pd.read_excel(path, sheet_name="Sheet1")
    df.columns = [
        "plasma_sample", "tumor_injection_date", "treatment",
        "collection_date", "bleed_type", "storage_note", "cfdna_conc_ng_ul",
    ]

    df["normalized_id"] = df["plasma_sample"].apply(
        lambda x: normalize_animal_id(str(x).strip(), source="plasma_bank") if pd.notna(x) else None
    )
    df["model"] = df["plasma_sample"].apply(
        lambda x: parse_model_from_plasma_bank(str(x).strip()) if pd.notna(x) else None
    )
    df = df.dropna(subset=["normalized_id"])

    # Unique animals — take first treatment and injection date per animal
    animals = (
        df.groupby("normalized_id")
        .agg({"model": "first", "treatment": "first", "tumor_injection_date": "first"})
        .reset_index()
        .rename(columns={"treatment": "treatment_raw"})
    )

    # Sample records
    samples = df[["normalized_id", "collection_date", "bleed_type", "storage_note", "cfdna_conc_ng_ul"]].copy()
    samples["bleed_type"] = samples["bleed_type"].str.strip().str.lower()
    samples = samples.rename(columns={"storage_note": "notes"})

    return animals, samples


def read_tumor_volume(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read tumor volume and body weight xlsx.

    Returns:
        (animals_df, measurements_df):
          - animals_df: unique animals with normalized_id, model
          - measurements_df: one row per (mouse, date) with tumor_volume_mm3, days_post_implant, measurement_date
    """
    wb = pd.ExcelFile(path)
    all_animals = []
    all_measurements = []

    for sheet_name in ["MN2", "WU225", "WU487", "JH-2-055"]:
        ws = pd.read_excel(wb, sheet_name=sheet_name, header=None)
        model = SHEET_TO_MODEL[sheet_name]

        # Row 1 (0-indexed) has dates starting around col 5
        date_row = ws.iloc[1]
        dates = {}
        for col_idx in range(5, len(date_row)):
            val = date_row.iloc[col_idx]
            if pd.notna(val):
                try:
                    dates[col_idx] = pd.Timestamp(val)
                except (ValueError, TypeError):
                    pass

        # Row 2 has days_post_implant
        days_row = ws.iloc[2]

        # Animals start around row 4, column B (index 1)
        for row_idx in range(4, len(ws)):
            animal_raw = ws.iloc[row_idx, 1]
            if pd.isna(animal_raw) or not str(animal_raw).strip().startswith("#"):
                continue

            normalized = normalize_animal_id(
                str(animal_raw).strip(), source="tumor_vol", model_sheet=sheet_name
            )
            all_animals.append({"normalized_id": normalized, "model": model})

            for col_idx, mdate in dates.items():
                vol = ws.iloc[row_idx, col_idx]
                if pd.notna(vol):
                    try:
                        vol_float = float(vol)
                    except (ValueError, TypeError):
                        continue
                    day_val = days_row.iloc[col_idx] if col_idx < len(days_row) else None
                    all_measurements.append({
                        "normalized_id": normalized,
                        "measurement_date": mdate,
                        "days_post_implant": int(day_val) if pd.notna(day_val) else None,
                        "tumor_volume_mm3": vol_float,
                    })

    animals_df = pd.DataFrame(all_animals).drop_duplicates(subset=["normalized_id"])
    measurements_df = pd.DataFrame(all_measurements)
    return animals_df, measurements_df


def read_qubit_data(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read JH 2-055 Serial Samples WashU Qubit Data xlsx.

    Returns:
        (qubit_df, bioanalyzer_df)
    """
    df = pd.read_excel(path, sheet_name="Sheet1")

    qubit = pd.DataFrame({
        "source_sample": df.iloc[:, 0],
        "treatment_group": df.iloc[:, 1],
        "collection_date": df.iloc[:, 2],
        "day_post_treatment": df.iloc[:, 3],
        "tumor_volume": df.iloc[:, 4],
        "qubit_conc_ng_ul": df.iloc[:, 5],
        "total_cfdna_ng": df.iloc[:, 6],
        "plasma_volume_ul": df.iloc[:, 7],
        "cfdna_per_ul_plasma": df.iloc[:, 8],
    })
    qubit = qubit.dropna(subset=["source_sample"])

    # Normalize: sample names like '#5-2 ' or '5-2' -> 'JH-2-055 5-2'
    def _norm_qubit_id(x):
        s = str(x).strip().lstrip("#")
        # Remove trailing space/text after the group-animal part
        s = s.strip()
        return normalize_animal_id(f"2-055 {s}", source="plasma_bank")

    qubit["normalized_id"] = qubit["source_sample"].apply(_norm_qubit_id)

    # Sheet 2: Bioanalyzer — complex layout, placeholder for now
    bioanalyzer = pd.DataFrame()

    return qubit, bioanalyzer


# ---------------------------------------------------------------------------
# Roster and sample builders
# ---------------------------------------------------------------------------

def build_mouse_roster(
    existing_mice: pd.DataFrame,
    plasma_animals: pd.DataFrame,
    tumor_animals: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build unified mouse roster from all sources.

    Returns:
        (roster_df, audit_df)
    """
    audit_rows = []

    # Normalize existing alt_subj for matching
    existing_norm = {}
    if not existing_mice.empty:
        for _, row in existing_mice.iterrows():
            if pd.notna(row.get("alt_subj")):
                norm = normalize_animal_id(str(row["alt_subj"]), source="existing")
                existing_norm[norm] = row["mouse_id"]

    # Collect all unique animals from plasma bank
    plasma_lookup = {}
    for _, row in plasma_animals.iterrows():
        plasma_lookup[row["normalized_id"]] = row

    # Animals only in tumor volume
    all_animals = set(plasma_lookup.keys())
    tumor_only = set()
    for _, row in tumor_animals.iterrows():
        nid = row["normalized_id"]
        if nid not in all_animals:
            all_animals.add(nid)
            tumor_only.add(nid)

    # Determine next ID
    next_id = 19
    if not existing_mice.empty:
        max_existing = existing_mice["mouse_id"].str.extract(r"(\d+)")[0].astype(int).max()
        next_id = max(next_id, max_existing + 1)

    roster_rows = []

    # Preserve existing mice, backfill from plasma bank
    for _, row in existing_mice.iterrows():
        new_row = row.to_dict()
        norm = None
        if pd.notna(row.get("alt_subj")):
            norm = normalize_animal_id(str(row["alt_subj"]), source="existing")
            new_row["alt_subj"] = norm

        if norm and norm in plasma_lookup:
            pb = plasma_lookup[norm]
            raw_tx = pb.get("treatment_raw")
            if pd.notna(raw_tx):
                new_row["treatment"] = TREATMENT_MAP.get(str(raw_tx).strip())
            new_row["tumor_injection_date"] = pb.get("tumor_injection_date")
            audit_rows.append({
                "source_file": "plasma_bank",
                "source_id": norm,
                "normalized_id": norm,
                "mou_id": row["mouse_id"],
                "match_type": "existing",
                "anomaly_flag": "",
            })
        else:
            # Existing mouse not in any source — preserve as-is
            if "treatment" not in new_row:
                new_row["treatment"] = None
            if "tumor_injection_date" not in new_row:
                new_row["tumor_injection_date"] = None

        roster_rows.append(new_row)

    # New animals
    assigned = set(existing_norm.keys())
    new_animals = sorted(all_animals - assigned)
    for nid in new_animals:
        mou_id = f"mou{next_id:04d}"
        next_id += 1

        if nid.startswith("No PDX"):
            model = None
        elif nid in plasma_lookup:
            model = plasma_lookup[nid]["model"]
        else:
            model = nid.split(" ")[0] if " " in nid else None

        treatment = None
        tumor_injection_date = None
        anomaly = ""

        if nid in plasma_lookup:
            pb = plasma_lookup[nid]
            raw_tx = pb.get("treatment_raw")
            if pd.notna(raw_tx):
                treatment = TREATMENT_MAP.get(str(raw_tx).strip())
            tumor_injection_date = pb.get("tumor_injection_date")
            source = "plasma_bank"
        else:
            source = "tumor_vol"
            anomaly = "tumor_vol_only"

        # Flag three-part IDs
        parts = nid.split(" ", 1)
        if len(parts) > 1 and parts[1].count("-") > 1:
            anomaly = "three-part ID" if not anomaly else f"{anomaly}; three-part ID"

        roster_rows.append({
            "mouse_id": mou_id,
            "strain": "nod",
            "stock": "nod_cg_rag1_null_il2rg_null",
            "sex": None,
            "birth_date": None,
            "institution_housed": "washu",
            "study_protocol_id": None,
            "pdx_line_id": model if model != "No PDX" else None,
            "alt_subj": nid,
            "notes": None,
            "treatment": treatment,
            "tumor_injection_date": tumor_injection_date,
        })
        audit_rows.append({
            "source_file": source,
            "source_id": nid,
            "normalized_id": nid,
            "mou_id": mou_id,
            "match_type": "new",
            "anomaly_flag": anomaly,
        })

    roster_df = pd.DataFrame(roster_rows)
    audit_df = pd.DataFrame(audit_rows)
    return roster_df, audit_df


def build_samples(
    existing_samples: pd.DataFrame,
    plasma_samples: pd.DataFrame,
    mouse_roster: pd.DataFrame,
) -> pd.DataFrame:
    """Build unified samples table. Renames sam->smp prefix and adds new samples."""
    existing = existing_samples.copy()
    existing["sample_id"] = existing["sample_id"].str.replace("^sam", "smp", regex=True)

    norm_to_mou = dict(zip(mouse_roster["alt_subj"], mouse_roster["mouse_id"]))

    if not existing.empty:
        max_num = existing["sample_id"].str.extract(r"(\d+)")[0].astype(int).max()
    else:
        max_num = 624
    next_num = max_num + 1

    new_rows = []
    for _, row in plasma_samples.iterrows():
        nid = row["normalized_id"]
        mouse_id = norm_to_mou.get(nid)
        if mouse_id is None:
            continue

        # Skip if already exists (same mouse + date)
        date_val = row["collection_date"]
        if not existing.empty and mouse_id in existing["mouse_id"].values:
            existing_dates = existing.loc[
                existing["mouse_id"] == mouse_id, "collection_date"
            ]
            if pd.notna(date_val):
                date_str = str(date_val)[:10]
                match = existing_dates.apply(lambda d: str(d)[:10] == date_str if pd.notna(d) else False)
                if match.any():
                    continue

        smp_id = f"smp{next_num:04d}"
        next_num += 1
        new_rows.append({
            "sample_id": smp_id,
            "mouse_id": mouse_id,
            "bleed_type": row.get("bleed_type"),
            "timepoint_day": None,
            "collection_date": date_val,
            "volume_ul": None,
            "visual_hemolysis": None,
            "alt_spec_id": None,
            "notes": row.get("notes"),
        })

    new_df = pd.DataFrame(new_rows)
    return pd.concat([existing, new_df], ignore_index=True)


def build_tumor_measurements(
    tv_measurements: pd.DataFrame,
    mouse_roster: pd.DataFrame,
) -> pd.DataFrame:
    """Build tumor_measurements sheet from parsed tumor volume data."""
    norm_to_mou = dict(zip(mouse_roster["alt_subj"], mouse_roster["mouse_id"]))

    rows = []
    for _, row in tv_measurements.iterrows():
        mouse_id = norm_to_mou.get(row["normalized_id"])
        if mouse_id is None:
            continue
        rows.append({
            "mouse_id": mouse_id,
            "measurement_date": row["measurement_date"],
            "days_post_implant": row.get("days_post_implant"),
            "tumor_volume_mm3": row.get("tumor_volume_mm3"),
            "body_weight_g": row.get("body_weight_g"),
        })

    return pd.DataFrame(rows)


def fix_fk_prefix(
    cfdna_isolations: pd.DataFrame,
    seq_libraries: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fix sam -> smp prefix in FK references."""
    iso = cfdna_isolations.copy()
    iso["sample_id"] = iso["sample_id"].str.replace("^sam", "smp", regex=True)
    lib = seq_libraries.copy()
    return iso, lib


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Run the one-off metadata ingest."""
    repo = Path(__file__).parent.parent
    inputs = Path("/mnt/gcs/jeszyman/projects/nf1-mouse/inputs/metadata")
    existing_path = repo / "data" / "metadata.xlsx"
    output_path = repo / "data" / "metadata.xlsx"
    audit_path = repo / "data" / "ingest_audit.csv"

    print("Reading existing metadata...")
    existing_mice = pd.read_excel(existing_path, sheet_name="mice")
    existing_samples = pd.read_excel(existing_path, sheet_name="samples")
    existing_iso = pd.read_excel(existing_path, sheet_name="cfdna_isolations")
    existing_lib = pd.read_excel(existing_path, sheet_name="sequencing_libraries")
    existing_pdx = pd.read_excel(existing_path, sheet_name="pdx_lines")

    print("Reading source files...")
    plasma_animals, plasma_samples = read_plasma_bank(
        inputs / "Murine Plasma Sample Bank.xlsx"
    )
    tumor_animals, tumor_measurements_raw = read_tumor_volume(
        inputs / "Tumor volume and body weight-A+VS-MN2, Shim8, PP4T, JH-2-055-110524.xlsx"
    )
    qubit_data, bioanalyzer_data = read_qubit_data(
        inputs / "JH 2-055 Serial Samples WashU Qubit Data.xlsx"
    )

    print("Building mouse roster...")
    mouse_roster, audit_df = build_mouse_roster(existing_mice, plasma_animals, tumor_animals)

    print("Building samples...")
    samples = build_samples(existing_samples, plasma_samples, mouse_roster)

    print("Building tumor measurements...")
    tumor_meas = build_tumor_measurements(tumor_measurements_raw, mouse_roster)

    print("Fixing FK prefixes...")
    iso_fixed, lib_fixed = fix_fk_prefix(existing_iso, existing_lib)

    # Add new pdx_lines
    new_pdx_rows = []
    for line_id in ["WU-225", "MN2"]:
        if line_id not in existing_pdx["pdx_line_id"].values:
            new_pdx_rows.append({
                "pdx_line_id": line_id,
                "chr8q_status": None,
                "institution": "washu",
                "notes": "chr8q status TBD",
            })
    pdx_lines = pd.concat([existing_pdx, pd.DataFrame(new_pdx_rows)], ignore_index=True)

    # Write audit
    print(f"Writing audit to {audit_path}...")
    audit_df.to_csv(audit_path, index=False)

    # Write metadata
    print(f"Writing metadata to {output_path}...")
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        mouse_roster.to_excel(writer, sheet_name="mice", index=False)
        pdx_lines.to_excel(writer, sheet_name="pdx_lines", index=False)
        samples.to_excel(writer, sheet_name="samples", index=False)
        iso_fixed.to_excel(writer, sheet_name="cfdna_isolations", index=False)
        lib_fixed.to_excel(writer, sheet_name="sequencing_libraries", index=False)
        tumor_meas.to_excel(writer, sheet_name="tumor_measurements", index=False)

    # Summary
    print(f"\nDone!")
    print(f"  Mice: {len(mouse_roster)}")
    print(f"  PDX lines: {len(pdx_lines)}")
    print(f"  Samples: {len(samples)}")
    print(f"  Tumor measurements: {len(tumor_meas)}")
    print(f"  Audit entries: {len(audit_df)}")
    anomalies = audit_df[audit_df["anomaly_flag"].str.len() > 0]
    if len(anomalies) > 0:
        print(f"  Anomalies flagged: {len(anomalies)}")
        for _, row in anomalies.iterrows():
            print(f"    {row['mou_id']} ({row['source_id']}): {row['anomaly_flag']}")


if __name__ == "__main__":
    main()
