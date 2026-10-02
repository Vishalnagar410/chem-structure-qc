import pandas as pd

from chem_structure_qc.qc import mol_from_smiles, process_dataframe


def test_invalid_smiles_returns_none():
    assert mol_from_smiles("not-a-smiles") is None


def test_process_dataframe_counts_invalid_and_duplicates():
    df = pd.DataFrame(
        {
            "compound_id": ["A", "B", "C", "D"],
            "smiles": ["CCO", "OCC", "c1ccccc1", "invalid"],
        }
    )
    valid, invalid, duplicates, report = process_dataframe(df)

    assert len(valid) == 3
    assert len(invalid) == 1
    assert len(duplicates) == 2
    assert report.unique_structures == 2
    assert report.duplicate_rows == 2
    assert "canonical_smiles" in valid.columns
    assert "molecular_weight" in valid.columns
