from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence

import pandas as pd

from .qc import process_dataframe, write_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate SMILES, generate RDKit canonical SMILES, identify duplicate "
            "structures, and calculate molecular descriptors from a CSV file. "
            "Canonicalization does not perform complete chemical standardization."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input CSV with compound_id and the selected SMILES column",
    )
    parser.add_argument("--output", required=True, type=Path, help="Directory for QC CSV and JSON outputs")
    parser.add_argument("--smiles-column", default="smiles", help="Name of the SMILES column")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run dataset QC and return a process-compatible status code."""
    args = build_parser().parse_args(argv)
    try:
        df = pd.read_csv(args.input)
        results, invalid, duplicate_groups, report = process_dataframe(df, args.smiles_column)
        write_outputs(results, invalid, duplicate_groups, report, args.output)
    except (OSError, UnicodeError, ValueError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    print("Chemical Structure QC")
    print("---------------------")
    print(f"Input rows          : {report.input_rows}")
    print(f"Valid structures    : {report.valid_structures}")
    print(f"Invalid structures  : {report.invalid_structures}")
    print(f"Unique structures   : {report.unique_structures}")
    print(f"Duplicate rows      : {report.duplicate_rows}")
    print(f"Duplicate groups    : {report.duplicate_groups}")
    print("\nOutputs written to:")
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
