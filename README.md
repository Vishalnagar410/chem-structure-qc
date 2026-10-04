# Chemical Structure Validation, Canonicalization & QC

A reproducible **Python/RDKit cheminformatics toolkit** for validating chemical structures, canonicalizing SMILES, detecting duplicate compounds, calculating molecular descriptors, and assessing chemical dataset quality.

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-green)
![Tests](https://img.shields.io/badge/tests-16%20passed-success)
![License](https://img.shields.io/badge/license-MIT-blue)

---

## Overview

`chem-structure-qc` is a small Python package and command-line tool for performing reproducible quality control on chemical datasets stored as CSV files.

The package:

- validates chemical structures represented as SMILES
- generates RDKit canonical SMILES
- identifies duplicate structures
- calculates selected molecular descriptors
- reports invalid or missing structures
- produces machine-readable QC outputs in CSV and JSON formats

An optional Streamlit dashboard provides an interactive interface over the same core QC package.

---

## Why this matters in chemical data workflows

Chemical datasets frequently contain:

- invalid or incomplete structures
- different SMILES representations of the same molecule
- duplicate compound records
- inconsistent chemical records that can affect downstream analysis

A repeatable QC pass makes these issues visible before data is used for analysis, visualization, or model preparation.

This project performs **data-quality assessment**; it does not make scientific decisions about which compound record should be retained.

---

## Key Features

### Chemical Structure QC

- Validate required input columns
- Parse SMILES using RDKit
- Report invalid and missing structures with reasons
- Generate canonical SMILES for valid molecules
- Detect duplicate structures using exact canonical-SMILES matching
- Retain one result row per input record

### Molecular Descriptors

For valid molecules, the package calculates:

- Molecular weight
- Crippen LogP
- TPSA
- H-bond donors (HBD)
- H-bond acceptors (HBA)
- Rotatable bonds
- Heavy atoms
- Ring count

### Outputs

The workflow generates:

- QC results
- invalid-structure report
- duplicate-group report
- JSON QC summary

---

## Workflow

```text
CSV Input
   ↓
SMILES Validation
   ↓
RDKit Canonicalization
   ↓
Duplicate Detection
   ↓
Molecular Descriptor Calculation
   ↓
QC Results + JSON Summary
