import json

import pandas as pd
import pytest

from chem_structure_qc.qc import descriptor_row, mol_from_smiles, process_dataframe


def test_invalid_smiles_returns_none():
    assert mol_from_smiles("not-a-smiles") is None


def test_process_dataframe_counts_invalid_and_duplicate_groups():
    df = pd.DataFrame(
        {
            "compound_id": ["A;1", "B", "C", "D"],
            "smiles": ["CCO", "OCC", "c1ccccc1", "invalid"],
        }
    )
    results, invalid, duplicate_groups, report = process_dataframe(df)

    assert len(results) == 4
    assert results["is_valid"].tolist() == [True, True, True, False]
    assert len(invalid) == 1
    assert invalid.loc[0, "error_reason"] == "Could not parse SMILES"
    assert len(duplicate_groups) == 1
    assert json.loads(duplicate_groups.loc[0, "compound_ids"]) == ["A;1", "B"]
    assert report.unique_structures == 2
    assert report.duplicate_rows == 2
    assert report.duplicate_groups == 1
    assert results.loc[0, "canonical_smiles"] == results.loc[1, "canonical_smiles"]
    assert results["is_duplicate"].tolist() == [True, True, False, False]


def test_missing_smiles_is_reported_without_crashing():
    df = pd.DataFrame({"compound_id": ["A", "B"], "smiles": [None, "   "]})

    results, invalid, duplicate_groups, report = process_dataframe(df)

    assert results["is_valid"].tolist() == [False, False]
    assert invalid["error_reason"].tolist() == ["Missing or empty SMILES"] * 2
    assert duplicate_groups.empty
    assert report.invalid_structures == 2
    assert report.unique_structures == 0


def test_descriptor_row_calculates_expected_ethanol_properties():
    descriptors = descriptor_row(mol_from_smiles("CCO"))

    assert descriptors["molecular_weight"] == pytest.approx(46.069, abs=0.001)
    assert descriptors["logp"] == pytest.approx(-0.0014, abs=0.001)
    assert descriptors["hbd"] == 1
    assert descriptors["hba"] == 1
    assert descriptors["heavy_atoms"] == 3
    assert descriptors["ring_count"] == 0


def test_descriptor_failure_returns_null_without_interrupting_calculation(monkeypatch):
    molecule = mol_from_smiles("CCO")
    assert molecule is not None

    def fail_molecular_weight(_molecule):
        raise RuntimeError("descriptor unavailable")

    monkeypatch.setattr("chem_structure_qc.qc.Descriptors.MolWt", fail_molecular_weight)

    descriptors = descriptor_row(molecule)

    assert descriptors["molecular_weight"] is None
    assert descriptors["heavy_atoms"] == 3


def test_required_columns_are_validated():
    with pytest.raises(ValueError, match="compound_id, smiles"):
        process_dataframe(pd.DataFrame({"structure": ["CCO"]}))


def test_empty_dataset_returns_empty_tables_and_zero_counts():
    empty = pd.DataFrame(columns=["compound_id", "smiles"])

    results, invalid, duplicate_groups, report = process_dataframe(empty)

    assert results.empty
    assert invalid.empty
    assert duplicate_groups.empty
    assert report.input_rows == 0
    assert report.valid_structures == 0
    assert report.invalid_rate == 0.0
    assert report.duplicate_rate == 0.0
