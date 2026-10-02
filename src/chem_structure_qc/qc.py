from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit import rdBase


DESCRIPTOR_COLUMNS = (
    "molecular_weight",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atoms",
    "ring_count",
)


@dataclass(frozen=True)
class QCResult:
    input_rows: int
    valid_structures: int
    invalid_structures: int
    unique_structures: int
    duplicate_rows: int
    duplicate_groups: int

    @property
    def invalid_rate(self) -> float:
        return self.invalid_structures / self.input_rows if self.input_rows else 0.0

    @property
    def duplicate_rate(self) -> float:
        return self.duplicate_rows / self.valid_structures if self.valid_structures else 0.0


def _parse_smiles(smiles: object) -> tuple[Chem.Mol | None, str | None]:
    if smiles is None or pd.isna(smiles):
        return None, "Missing or empty SMILES"
    text = str(smiles).strip()
    if not text:
        return None, "Missing or empty SMILES"
    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(text)
    except Exception:
        return None, "Could not parse SMILES"
    if molecule is None:
        return None, "Could not parse SMILES"
    return molecule, None


def mol_from_smiles(smiles: object) -> Chem.Mol | None:
    """Parse a SMILES value, returning None when it is missing or invalid."""
    molecule, _ = _parse_smiles(smiles)
    return molecule


def _calculate(value: Any, transform: Any) -> float | int | None:
    try:
        result = transform(value)
        return round(result, 4) if isinstance(result, float) else int(result)
    except Exception:
        return None


def descriptor_row(mol: Chem.Mol) -> dict[str, float | int | None]:
    """Calculate the supported molecular descriptors for one valid molecule."""
    return {
        "molecular_weight": _calculate(mol, Descriptors.MolWt),
        "logp": _calculate(mol, Crippen.MolLogP),
        "tpsa": _calculate(mol, rdMolDescriptors.CalcTPSA),
        "hbd": _calculate(mol, Lipinski.NumHDonors),
        "hba": _calculate(mol, Lipinski.NumHAcceptors),
        "rotatable_bonds": _calculate(mol, Lipinski.NumRotatableBonds),
        "heavy_atoms": _calculate(mol, lambda molecule: molecule.GetNumHeavyAtoms()),
        "ring_count": _calculate(mol, rdMolDescriptors.CalcNumRings),
    }


def process_dataframe(
    df: pd.DataFrame,
    smiles_column: str = "smiles",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, QCResult]:
    """Validate input rows and return results, invalid rows, duplicate groups, and counts."""
    required_columns = ["compound_id", smiles_column]
    missing_columns = [column for column in dict.fromkeys(required_columns) if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    result_rows: list[dict[str, object]] = []
    invalid_rows: list[dict[str, object]] = []
    for row in df.to_dict(orient="records"):
        compound_id = row["compound_id"]
        input_smiles = row[smiles_column]
        molecule, error_reason = _parse_smiles(input_smiles)
        canonical_smiles = None
        descriptors: dict[str, float | int | None] = dict.fromkeys(DESCRIPTOR_COLUMNS)
        if molecule is not None:
            try:
                canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
            except Exception:
                error_reason = "Could not generate canonical SMILES"
            if error_reason is None:
                descriptors = descriptor_row(molecule)
        is_valid = error_reason is None
        result_rows.append(
            {
                "compound_id": compound_id,
                "input_smiles": input_smiles,
                "canonical_smiles": canonical_smiles if is_valid else None,
                "is_valid": is_valid,
                **descriptors,
                "is_duplicate": False,
            }
        )
        if not is_valid:
            invalid_rows.append(
                {
                    "compound_id": compound_id,
                    "input_smiles": input_smiles,
                    "error_reason": error_reason,
                }
            )

    result_columns = [
        "compound_id", "input_smiles", "canonical_smiles", "is_valid",
        *DESCRIPTOR_COLUMNS, "is_duplicate",
    ]
    results = pd.DataFrame(result_rows, columns=result_columns)
    invalid = pd.DataFrame(invalid_rows, columns=["compound_id", "input_smiles", "error_reason"])
    duplicate_group_rows: list[dict[str, object]] = []
    if not results.empty:
        valid_mask = results["is_valid"]
        duplicate_counts = results.loc[valid_mask, "canonical_smiles"].value_counts()
        duplicate_smiles = duplicate_counts[duplicate_counts > 1].index
        results.loc[results["canonical_smiles"].isin(duplicate_smiles), "is_duplicate"] = True
        for canonical_smiles in duplicate_smiles:
            members = results.loc[results["canonical_smiles"] == canonical_smiles, "compound_id"]
            duplicate_group_rows.append(
                {
                    "canonical_smiles": canonical_smiles,
                    "duplicate_count": len(members),
                    "compound_ids": json.dumps([str(compound_id) for compound_id in members]),
                }
            )
    duplicate_groups = pd.DataFrame(
        duplicate_group_rows,
        columns=["canonical_smiles", "duplicate_count", "compound_ids"],
    )
    valid_count = int(results["is_valid"].sum())
    duplicate_row_count = int(results["is_duplicate"].sum())
    report = QCResult(
        input_rows=len(df),
        valid_structures=valid_count,
        invalid_structures=len(invalid),
        unique_structures=valid_count - duplicate_row_count + len(duplicate_groups),
        duplicate_rows=duplicate_row_count,
        duplicate_groups=len(duplicate_groups),
    )
    return results, invalid, duplicate_groups, report


def write_outputs(
    results: pd.DataFrame,
    invalid: pd.DataFrame,
    duplicate_groups: pd.DataFrame,
    report: QCResult,
    output_dir: Path,
) -> None:
    """Write the machine-readable QC tables and summary report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "qc_results.csv", index=False)
    invalid.to_csv(output_dir / "invalid_structures.csv", index=False)
    duplicate_groups.to_csv(output_dir / "duplicate_groups.csv", index=False)

    payload = asdict(report)
    payload["invalid_rate"] = report.invalid_rate
    payload["duplicate_rate"] = report.duplicate_rate
    (output_dir / "qc_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
