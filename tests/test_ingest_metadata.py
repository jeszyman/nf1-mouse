"""Tests for metadata ingest script."""
import sys

sys.path.insert(0, "scripts")
from ingest_metadata import normalize_animal_id, parse_model_from_plasma_bank


class TestNormalizeAnimalId:
    """Test animal ID normalization across source formats."""

    def test_plasma_bank_wu487(self):
        assert normalize_animal_id("WU-487 2-1", source="plasma_bank") == "WU-487 2-1"

    def test_plasma_bank_2055_prepends_jh(self):
        assert normalize_animal_id("2-055 1-4", source="plasma_bank") == "JH-2-055 1-4"

    def test_plasma_bank_mn2(self):
        assert normalize_animal_id("MN2 5-2", source="plasma_bank") == "MN2 5-2"

    def test_plasma_bank_wu225(self):
        assert normalize_animal_id("WU-225 1-1", source="plasma_bank") == "WU-225 1-1"

    def test_plasma_bank_no_pdx(self):
        assert normalize_animal_id("No PDX 1", source="plasma_bank") == "No PDX 1"

    def test_tumor_vol_strips_hash(self):
        assert normalize_animal_id("#2-1", source="tumor_vol", model_sheet="WU487") == "WU-487 2-1"

    def test_tumor_vol_wu225_adds_hyphen(self):
        assert normalize_animal_id("#7-1", source="tumor_vol", model_sheet="WU225") == "WU-225 7-1"

    def test_tumor_vol_jh2055(self):
        assert normalize_animal_id("#4-5", source="tumor_vol", model_sheet="JH-2-055") == "JH-2-055 4-5"

    def test_tumor_vol_mn2(self):
        assert normalize_animal_id("#3-1", source="tumor_vol", model_sheet="MN2") == "MN2 3-1"

    def test_existing_metadata_wu487_no_hyphen(self):
        assert normalize_animal_id("WU487 2-1", source="existing") == "WU-487 2-1"

    def test_existing_metadata_jh2055(self):
        assert normalize_animal_id("JH-2-055 1-4", source="existing") == "JH-2-055 1-4"

    def test_existing_metadata_mn2(self):
        assert normalize_animal_id("MN2 5-2", source="existing") == "MN2 5-2"

    def test_three_part_id(self):
        assert normalize_animal_id("2-055 4-1-2", source="plasma_bank") == "JH-2-055 4-1-2"

    def test_whitespace_stripped(self):
        assert normalize_animal_id("  WU-487 2-1  ", source="plasma_bank") == "WU-487 2-1"


class TestParseModelFromPlasmaBank:
    """Test model name extraction from Plasma Bank sample IDs."""

    def test_wu487(self):
        assert parse_model_from_plasma_bank("WU-487 2-1") == "WU-487"

    def test_wu225(self):
        assert parse_model_from_plasma_bank("WU-225 1-1") == "WU-225"

    def test_mn2(self):
        assert parse_model_from_plasma_bank("MN2 5-2") == "MN2"

    def test_2055(self):
        assert parse_model_from_plasma_bank("2-055 1-4") == "JH-2-055"

    def test_no_pdx(self):
        assert parse_model_from_plasma_bank("No PDX 1") == "No PDX"
