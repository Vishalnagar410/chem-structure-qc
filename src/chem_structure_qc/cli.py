from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .qc import process_dataframe, write_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate, standardize and profile chemical structures in CSV format."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input CSV path")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--smiles-column", default="smiles", help="Name of the SMILES column")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    df = pd.read_csv(args.input)
    valid, invalid, duplicates, report = process_dataframe(df, args.smiles_column)
    write_outputs(valid, invalid, duplicates, report, args.output)

    print(f"Input rows:       {report.input_rows}")
    print(f"Valid structures: {report.valid_structures}")
    print(f"Invalid:          {report.invalid_structures}")
    print(f"Unique structures:{report.unique_structures}")
    print(f"Duplicate rows:   {report.duplicate_rows}")
    print(f"Outputs written:  {args.output.resolve()}")


if __name__ == "__main__":
    main()
