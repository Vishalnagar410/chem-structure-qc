import json
from pathlib import Path

import pandas as pd
import pytest

from app.ui_helpers import analyze_dataframe, build_download_payloads, render_smiles


PROJECT_ROOT = Path(__file__).parents[1]


def test_example_dataset_analysis_matches_documented_summary():
    compounds = pd.read_csv(PROJECT_ROOT / "data" / "example_compounds.csv")

    analysis = analyze_dataframe(compounds)

    assert analysis.report.input_rows == 6
    assert analysis.report.valid_structures == 5
    assert analysis.report.invalid_structures == 1
    assert analysis.report.unique_structures == 4
    assert analysis.report.duplicate_rows == 2
    assert analysis.report.duplicate_groups == 1
    assert analysis.invalid.loc[0, "compound_id"] == "CMP006"
    assert analysis.duplicate_groups.loc[0, "compound_ids"] == '["CMP001", "CMP002"]'


def test_analysis_accepts_a_custom_smiles_column():
    compounds = pd.DataFrame({"compound_id": ["A"], "structure": ["CCO"]})

    analysis = analyze_dataframe(compounds, "structure")

    assert analysis.results.loc[0, "canonical_smiles"] == "CCO"


def test_analysis_uses_core_required_column_validation():
    with pytest.raises(ValueError, match="Missing required columns: compound_id, structure"):
        analyze_dataframe(pd.DataFrame({"smiles": ["CCO"]}), "structure")


def test_download_payloads_contain_all_expected_outputs():
    compounds = pd.read_csv(PROJECT_ROOT / "data" / "example_compounds.csv")
    payloads = build_download_payloads(analyze_dataframe(compounds))

    assert set(payloads) == {
        "qc_results.csv",
        "invalid_structures.csv",
        "duplicate_groups.csv",
        "qc_report.json",
    }
    assert b"CMP006" in payloads["invalid_structures.csv"]
    assert b"CMP001" in payloads["duplicate_groups.csv"]
    assert json.loads(payloads["qc_report.json"]) ["input_rows"] == 6


def test_render_smiles_returns_image_for_valid_input_and_none_for_invalid():
    image = render_smiles("CCO")

    assert image is not None
    assert image.size == (680, 480)
    assert render_smiles("not-a-smiles") is None