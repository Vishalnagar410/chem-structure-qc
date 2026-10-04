from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from ui_helpers import analyze_dataframe, build_download_payloads, render_smiles


st.set_page_config(
    page_title="Chemical Structure Validation & QC",
    layout="wide",
)
st.markdown(
    """
    <style>
    .block-container {max-width: 1920px; padding: 1.5rem clamp(1rem, 2vw, 2.25rem) 3rem;}
    h1 {color: #173d35; font-size: 2.4rem; font-weight: 680; line-height: 1.2; letter-spacing: 0; margin-bottom: 0.25rem;}
    h2, h3 {letter-spacing: 0; margin-top: 0.85rem;}
    .stApp {background: #f8faf9;}
    [data-testid="stMetric"] {border-left: 3px solid #76a995; padding: 0.55rem 0.7rem; background: rgba(255,255,255,0.65); min-height: 92px;}
    [data-testid="stMetricLabel"] p {font-size: 0.82rem; margin-bottom: 0.3rem;}
    [data-testid="stMetricValue"] {font-size: 1.5rem; line-height: 1.2;}
    .scope-note {border-left: 3px solid #d19a47; padding: 0.8rem 1rem; background: #f7f5ef; border-radius: 0.4rem;}
    .stDataFrame {border: 1px solid rgba(49, 65, 61, 0.12); border-radius: 0.5rem;}
    div[data-testid="stVerticalBlock"] > div:has(> [data-testid="stMetric"]) {margin-bottom: 0.5rem;}
    @media (max-width: 900px) {
        h1 {font-size: 2rem;}
        [data-testid="stMetric"] {padding: 0.45rem 0.5rem;}
        [data-testid="stMetricValue"] {font-size: 1.3rem;}
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("Chemical Structure Validation & QC")
st.caption("RDKit-powered quality control for chemical datasets containing compound identifiers and SMILES.")

st.header("Data input")
source = st.radio("Dataset source", ["Example dataset", "Upload CSV"], horizontal=True)
uploaded_file = None
if source == "Upload CSV":
    uploaded_file = st.file_uploader("Choose a CSV file", type=["csv"])
    if uploaded_file is None:
        st.info("Upload a CSV with a `compound_id` column and a SMILES column to begin.")
        st.stop()

try:
    if source == "Example dataset":
        example_path = Path(__file__).resolve().parents[1] / "data" / "example_compounds.csv"
        input_df = pd.read_csv(example_path)
        st.caption("Using the bundled educational example dataset.")
    else:
        input_df = pd.read_csv(uploaded_file)
except (OSError, UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:
    st.error(f"Could not read the CSV file: {error}")
    st.stop()

st.write("Detected columns:", ", ".join(str(column) for column in input_df.columns) or "none")
if "compound_id" not in input_df.columns:
    st.error("Required column `compound_id` was not found.")
    st.stop()
smiles_options = [column for column in input_df.columns if column != "compound_id"]
if not smiles_options:
    st.error("No SMILES column is available. Include a column containing SMILES strings.")
    st.stop()
smiles_column = st.selectbox(
    "SMILES column",
    smiles_options,
    index=smiles_options.index("smiles") if "smiles" in smiles_options else 0,
)

try:
    analysis = analyze_dataframe(input_df, str(smiles_column))
except ValueError as error:
    st.error(str(error))
    st.stop()

st.header("QC summary")
summary = analysis.report
metric_columns = st.columns(6, gap="medium")
for column, label, value in zip(
    metric_columns,
    ("Input Structures", "Valid Structures", "Invalid Structures", "Unique Structures", "Duplicate Rows", "Duplicate Groups"),
    (
        summary.input_rows,
        summary.valid_structures,
        summary.invalid_structures,
        summary.unique_structures,
        summary.duplicate_rows,
        summary.duplicate_groups,
    ),
):
    column.metric(label, value)

st.header("Visual summary")
chart_columns = st.columns(2, gap="large")
with chart_columns[0]:
    st.subheader("Structure Validity")
    st.bar_chart(
        pd.DataFrame({"Structures": [summary.valid_structures, summary.invalid_structures]}, index=["Valid", "Invalid"]),
        width="stretch",
        height=400,
    )
with chart_columns[1]:
    st.subheader("Duplicate Status / unique vs duplicate")
    st.bar_chart(
        pd.DataFrame({"Structures": [summary.unique_structures, summary.duplicate_rows]}, index=["Unique", "Duplicate rows"]),
        width="stretch",
        height=400,
    )
    st.caption("Duplicate rows are valid records beyond one representative per exact canonical-SMILES group.")

st.header("Structure explorer")
if analysis.results.empty:
    st.info("There are no structures to explore.")
else:
    result_rows = analysis.results.reset_index(drop=True)
    choices = result_rows.index.tolist()

    def compound_label(index: int) -> str:
        row = result_rows.loc[index]
        return f"{row['compound_id']}  |  {row['input_smiles']}"

    selected_index = st.selectbox("Select compound", choices, format_func=compound_label)
    selected = result_rows.loc[selected_index]
    molecule_column, detail_column = st.columns([1.55, 1.0], gap="large", vertical_alignment="top")
    with molecule_column:
        image = render_smiles(selected["input_smiles"], size=(900, 620))
        if image is None:
            st.warning("No structure image is available for this invalid SMILES.")
        else:
            st.image(image, caption="RDKit 2D depiction", width="stretch")
    with detail_column:
        st.markdown("### Compound details")
        st.write(f"**Compound ID:** {selected['compound_id']}")
        st.write(f"**Original SMILES:** `{selected['input_smiles']}`")
        canonical = selected["canonical_smiles"]
        st.write(f"**Canonical SMILES:** `{canonical if pd.notna(canonical) else 'Unavailable'}`")
        if bool(selected["is_valid"]):
            descriptors = {
                "MW": selected["molecular_weight"],
                "LogP": selected["logp"],
                "TPSA": selected["tpsa"],
                "HBD": selected["hbd"],
                "HBA": selected["hba"],
                "Rotatable Bonds": selected["rotatable_bonds"],
                "Heavy Atoms": selected["heavy_atoms"],
                "Ring Count": selected["ring_count"],
            }
            descriptor_columns = st.columns(2, gap="small")
            for column, (label, value) in zip(descriptor_columns * 4, descriptors.items()):
                column.metric(label, value)
        else:
            issue = analysis.invalid.loc[analysis.invalid["compound_id"] == selected["compound_id"], "error_reason"]
            st.warning(f"Invalid structure: {issue.iloc[0] if not issue.empty else 'Unable to parse SMILES'}")

st.header("Results table")
filters, toggles = st.columns([3, 1], gap="small")
with filters:
    validity_filter = st.multiselect("Validity", ["Valid", "Invalid"], default=["Valid", "Invalid"])
with toggles:
    duplicates_only = st.checkbox("Duplicates only")
filtered_results = analysis.results.copy()
validity_values = {"Valid": True, "Invalid": False}
if validity_filter:
    filtered_results = filtered_results[filtered_results["is_valid"].isin([validity_values[value] for value in validity_filter])]
else:
    filtered_results = filtered_results.iloc[0:0]
if duplicates_only:
    filtered_results = filtered_results[filtered_results["is_duplicate"]]
st.dataframe(filtered_results, width="stretch", hide_index=True, height=520)

st.header("Invalid structures")
if analysis.invalid.empty:
    st.caption("No invalid structures.")
else:
    st.dataframe(analysis.invalid, width="stretch", hide_index=True, height=240)

st.header("Duplicate groups")
if analysis.duplicate_groups.empty:
    st.caption("No duplicate groups.")
else:
    display_groups = analysis.duplicate_groups.rename(
        columns={"canonical_smiles": "Canonical SMILES", "duplicate_count": "Duplicate count", "compound_ids": "Compound IDs"}
    )
    st.dataframe(display_groups, width="stretch", hide_index=True, height=240)

st.header("Downloads")
payloads = build_download_payloads(analysis)
download_columns = st.columns(4, gap="small")
download_labels = {
    "qc_results.csv": "QC results",
    "invalid_structures.csv": "Invalid structures",
    "duplicate_groups.csv": "Duplicate groups",
    "qc_report.json": "QC report",
}
for column, filename in zip(download_columns, payloads):
    mime_type = "application/json" if filename.endswith(".json") else "text/csv"
    column.download_button(
        label=download_labels[filename],
        data=payloads[filename],
        file_name=filename,
        mime=mime_type,
        width="stretch",
    )

st.markdown(
    "<div class='scope-note'><strong>Scientific scope:</strong> canonicalization is not complete chemical standardization. This tool does not apply comprehensive salt/fragment handling, charge normalization, tautomer normalization, isotope normalization, or a stereochemistry policy. Duplicate groups use exact canonical-SMILES matches under the installed RDKit behavior.</div>",
    unsafe_allow_html=True,
)