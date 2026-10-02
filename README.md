# Chem Structure QC

A small, reproducible cheminformatics pipeline for validating and standardizing molecular structures from CSV files.

## Why this project

This project demonstrates practical chemical-data handling rather than a large application. It takes a table containing compound identifiers and SMILES, validates the structures, generates canonical SMILES and common molecular descriptors, detects duplicate structures, and writes separate quality-control outputs.

The implementation is intentionally independent and uses public/open-source tooling. It is not derived from any proprietary laboratory or ReySci application.

## Workflow

```text
CSV (compound_id, smiles)
        |
        v
  SMILES validation
        |
        +------ invalid_structures.csv
        |
        v
  Canonical SMILES
        |
        +------ duplicate structures
        |
        v
  Molecular descriptors
        |
        v
 clean_compounds.csv + qc_report.json
```

## Calculated descriptors

- Molecular weight
- LogP (Crippen)
- Topological polar surface area (TPSA)
- Hydrogen-bond donors (HBD)
- Hydrogen-bond acceptors (HBA)
- Rotatable bonds
- Heavy atoms
- Ring count

## Input

CSV with at least these columns:

```csv
compound_id,smiles
CMP001,CCO
CMP002,CC(=O)Oc1ccccc1C(=O)O
```

## Run

Create an environment and install the package:

```bash
python -m venv .venv
# Windows PowerShell
.venv\\Scripts\\Activate.ps1

pip install -e .
```

Run the CLI:

```bash
chem-structure-qc --input data/example_compounds.csv --output outputs
```

The output directory contains:

- `clean_compounds.csv` — valid, standardized structures and descriptors
- `invalid_structures.csv` — rows that could not be parsed as valid molecules
- `duplicate_structures.csv` — duplicate canonical structures
- `qc_report.json` — summary statistics

## Test

```bash
pytest
```

## Scope and limitations

This is a structure-QC demonstration, not a medicinal-chemistry registration standard. Standardization choices such as tautomer handling, charge normalization, stereochemistry policy, salt stripping and isotope treatment can be domain-specific and should be made explicit for a production workflow.

## License

MIT
