from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


@dataclass(frozen=True)
class QCResult:
    input_rows: int
    valid_structures: int
    invalid_structures: int
    unique_structures: int
    duplicate_rows: int

    @property
    def invalid_rate(self) -> float:
        return self.invalid_structures / self.input_rows if self.input_rows else 0.0

    @property
    def duplicate_rate(self) -> float:
        return self.duplicate_rows / self.valid_structures if self.valid_structures else 0.0


def mol_from_smiles(smiles: object) -> Chem.Mol | None:
    """Parse a SMILES value, returning None for missing/invalid input."""
    if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)):
        return None
    text = str(smiles).strip()
    if not text:
        return None
    try:
        return Chem.MolFromSmiles(text)
    except Exception:
        return None


def descriptor_row(mol: Chem.Mol) -> dict[str, float | int | str]:
    """Return a compact, recruiter-readable descriptor set."""
    return {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "molecular_weight": round(Descriptors.MolWt(mol), 4),
        "logp": round(Crippen.MolLogP(mol), 4),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 4),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
    }


def process_dataframe(df: pd.DataFrame, smiles_column: str = "smiles") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, QCResult]:
    """Validate and enrich a dataframe containing a SMILES column."""
    if smiles_column not in df.columns:
        raise ValueError(f"Missing required column: {smiles_column}")

    working = df.copy()
    working["_mol"] = working[smiles_column].map(mol_from_smiles)
    invalid_mask = working["_mol"].isna()

    invalid = working.loc[invalid_mask].drop(columns=["_mol"]).copy()
    valid = working.loc[~invalid_mask].copy()

    if not valid.empty:
        descriptors = valid["_mol"].map(descriptor_row).apply(pd.Series)
        valid = pd.concat([valid.drop(columns=["_mol"]), descriptors], axis=1)
        valid["is_duplicate"] = valid.duplicated(subset=["canonical_smiles"], keep=False)
    else:
        valid = valid.drop(columns=["_mol"])
        valid["is_duplicate"] = pd.Series(dtype=bool)

    duplicate_rows = valid.loc[valid["is_duplicate"]].copy()
    unique_structures = int(valid["canonical_smiles"].nunique()) if not valid.empty else 0

    report = QCResult(
        input_rows=len(df),
        valid_structures=len(valid),
        invalid_structures=len(invalid),
        unique_structures=unique_structures,
        duplicate_rows=len(duplicate_rows),
    )
    return valid, invalid, duplicate_rows, report


def write_outputs(valid: pd.DataFrame, invalid: pd.DataFrame, duplicates: pd.DataFrame, report: QCResult, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    valid.to_csv(output_dir / "clean_compounds.csv", index=False)
    invalid.to_csv(output_dir / "invalid_structures.csv", index=False)
    duplicates.to_csv(output_dir / "duplicate_structures.csv", index=False)

    payload = asdict(report)
    payload["invalid_rate"] = report.invalid_rate
    payload["duplicate_rate"] = report.duplicate_rate
    (output_dir / "qc_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
