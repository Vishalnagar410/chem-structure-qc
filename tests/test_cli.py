import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from chem_structure_qc.cli import main


PROJECT_ROOT = Path(__file__).parents[1]


def test_cli_writes_expected_outputs_and_summary(tmp_path, capsys):
    input_path = tmp_path / "compounds.csv"
    output_dir = tmp_path / "nested" / "qc-output"
    pd.DataFrame(
        {
            "compound_id": ["A", "B", "C"],
            "smiles": ["CCO", "OCC", "invalid"],
        }
    ).to_csv(input_path, index=False)

    exit_code = main(["--input", str(input_path), "--output", str(output_dir)])

    assert exit_code == 0
    assert "Duplicate rows      : 2" in capsys.readouterr().out
    assert {path.name for path in output_dir.iterdir()} == {
        "qc_results.csv",
        "invalid_structures.csv",
        "duplicate_groups.csv",
        "qc_report.json",
    }
    report = json.loads((output_dir / "qc_report.json").read_text(encoding="utf-8"))
    assert report["input_rows"] == 3
    assert report["valid_structures"] == 2
    assert report["invalid_structures"] == 1
    assert report["unique_structures"] == 1
    assert report["duplicate_rows"] == 2
    assert report["duplicate_groups"] == 1
    assert pd.read_csv(output_dir / "qc_results.csv").shape[0] == 3


def test_cli_missing_input_returns_error_without_traceback(tmp_path, capsys):
    exit_code = main(
        ["--input", str(tmp_path / "missing.csv"), "--output", str(tmp_path / "out")]
    )

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "Error:" in captured.err
    assert "Traceback" not in captured.err


def test_cli_missing_required_columns_returns_error(tmp_path, capsys):
    input_path = tmp_path / "wrong-columns.csv"
    pd.DataFrame({"structure": ["CCO"]}).to_csv(input_path, index=False)

    exit_code = main(["--input", str(input_path), "--output", str(tmp_path / "out")])

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "Missing required columns: compound_id, smiles" in captured.err
    assert "Traceback" not in captured.err


def test_cli_module_help_succeeds():
    completed = subprocess.run(
        [sys.executable, "-m", "chem_structure_qc.cli", "--help"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0
    assert "canonicalization does not perform complete chemical standardization" in completed.stdout.lower()
    assert "--input" in completed.stdout