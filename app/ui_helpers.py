from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Mapping

import pandas as pd
from PIL import Image
from rdkit.Chem import Draw

from chem_structure_qc.qc import QCResult, mol_from_smiles, process_dataframe


@dataclass(frozen=True)
class DashboardAnalysis:
    results: pd.DataFrame
    invalid: pd.DataFrame
    duplicate_groups: pd.DataFrame
    report: QCResult


def analyze_dataframe(df: pd.DataFrame, smiles_column: str = "smiles") -> DashboardAnalysis:
    results, invalid, duplicate_groups, report = process_dataframe(df, smiles_column)
    return DashboardAnalysis(results, invalid, duplicate_groups, report)


def build_download_payloads(analysis: DashboardAnalysis) -> Mapping[str, bytes]:
    report = asdict(analysis.report)
    report["invalid_rate"] = analysis.report.invalid_rate
    report["duplicate_rate"] = analysis.report.duplicate_rate
    return {
        "qc_results.csv": analysis.results.to_csv(index=False).encode("utf-8"),
        "invalid_structures.csv": analysis.invalid.to_csv(index=False).encode("utf-8"),
        "duplicate_groups.csv": analysis.duplicate_groups.to_csv(index=False).encode("utf-8"),
        "qc_report.json": json.dumps(report, indent=2).encode("utf-8"),
    }


def render_smiles(smiles: object, size: tuple[int, int] = (680, 480)) -> Image.Image | None:
    molecule = mol_from_smiles(smiles)
    if molecule is None:
        return None
    return Draw.MolToImage(molecule, size=size)