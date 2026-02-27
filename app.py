# app.py
# Streamlit app to build MS1 MassQL compendia grouped by Class / Subclass and export as ZIP

import io
import re
import zipfile
from typing import Dict, List
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="GPL MS1 MassQL Compendium Builder", layout="wide")

MONO_MASS = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
}

ADDUCT_SHIFT = {
    "H": 1.007276466812,
    "Na": 22.989218,
    "K": 38.963158,
    "NH4": 18.033823,
}

FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for elem, num in FORMULA_TOKEN_RE.findall(str(formula).strip()):
        if elem not in MONO_MASS:
            raise ValueError(f"Unsupported element '{elem}' in formula '{formula}'")
        counts[elem] = counts.get(elem, 0) + (int(num) if num else 1)
    return counts


def exact_mass_from_formula(formula: str) -> float:
    counts = parse_formula(formula)
    return sum(MONO_MASS[e] * n for e, n in counts.items())


def safe_str(x) -> str:
    return "" if x is None else str(x).strip()


def build_title(row: pd.Series) -> str:
    cls = safe_str(row.get("Class", ""))
    fa1 = safe_str(row.get("FA1", ""))
    fa2 = safe_str(row.get("FA2", ""))
    fa3 = safe_str(row.get("FA3", "")) if "FA3" in row.index else ""
    sumdb = safe_str(row.get("SumDB", ""))

    parts = [cls]
    if fa1:
        parts.append(fa1)
    if fa2:
        parts.append(fa2)
    if fa3:
        parts.append(fa3)
    parts.append(f"DB{sumdb}")

    return "_".join(parts).replace(" ", "")


def build_ms1_query(mz_list: List[float], tol: float, intensity_percent: int) -> str:
    mz_str = " OR ".join([f"{mz:.6f}" for mz in mz_list])
    return (
        f"QUERY scaninfo(MS1DATA) WHERE MS1MZ=({mz_str})"
        f":TOLERANCEMZ={tol}:INTENSITYPERCENT={intensity_percent}"
    )


def infer_subclass(row: pd.Series) -> str:
    # If FA2 is empty -> lyso (monoacyl), else diacyl
    fa2 = safe_str(row.get("FA2", ""))
    return "Lyso" if fa2 == "" or fa2.lower() == "nan" else "Diacyl"


def generate_queries(df: pd.DataFrame, adducts: List[str], tol: float, intensity_percent: int) -> pd.DataFrame:
    required = {"Class", "FA1", "FA2", "SumDB", "Formula"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = df.copy()

    # Normalize types
    for c in ["Class", "FA1", "FA2", "SumDB", "Formula"]:
        df[c] = df[c].astype(str).str.strip()
    if "FA3" in df.columns:
        df["FA3"] = df["FA3"].astype(str).str.strip()

    df["NeutralMass"] = df["Formula"].apply(exact_mass_from_formula)
    df["Subclass"] = df.apply(infer_subclass, axis=1)
    df["Title"] = df.apply(build_title, axis=1)

    # Ensure unique titles globally
    if df["Title"].duplicated().any():
        seen = {}
        new_titles = []
        for t in df["Title"].tolist():
            seen[t] = seen.get(t, 0) + 1
            new_titles.append(t if seen[t] == 1 else f"{t}_{seen[t]}")
        df["Title"] = new_titles

    # Build query
    def mzs(neutral_mass: float) -> List[float]:
        return [neutral_mass + ADDUCT_SHIFT[a] for a in adducts]

    df["Adducts"] = ",".join([f"[M+{a}]+" for a in adducts])
    df["AdductMZs"] = df["NeutralMass"].apply(lambda m: ",".join([f"{x:.6f}" for x in mzs(m)]))
    df["Query"] = df["NeutralMass"].apply(lambda m: build_ms1_query(mzs(m), tol, intensity_percent))

    return df


def to_compendium_tsv(df: pd.DataFrame) -> bytes:
    # Minimal columns for a compendium file:
    # Title \t Query
    out = df[["Title", "Query"]].copy()
    return out.to_csv(sep="\t", index=False).encode("utf-8")


def build_zip(compendia: Dict[str, pd.DataFrame]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for filename, subdf in compendia.items():
            zf.writestr(filename, to_compendium_tsv(subdf))
    return buf.getvalue()

# -----------------------------
# LOGOs
# -----------------------------

# Optional logos (won't crash if missing)
STATIC_DIR = Path(__file__).parent / "static"
for logo_name in ["logo_massQL.png", "LAABio.png"]:
    p = STATIC_DIR / logo_name
    try:
        from PIL import Image
        st.sidebar.image(Image.open(p), use_container_width=True)
    except Exception:
        pass

# Helpful links
st.markdown(
    """
---
"""
)


# -----------------------------
# UI
# -----------------------------
st.title("GPL MS1 MassQL Compendium Builder (Class/Subclass organized)")

# Helpful links
st.markdown(
    """
---
"""
)

with st.sidebar:
    st.header("Settings")

    tol = st.number_input("TOLERANCEMZ", min_value=0.0001, max_value=1.0, value=0.01, step=0.001, format="%.4f")
    intensity_percent = st.number_input("INTENSITYPERCENT", min_value=1, max_value=100, value=40, step=10)#, format="%.4f")

    adducts = st.multiselect("Positive-mode adducts", ["H", "Na", "K", "NH4"], default=["H", "Na", "K"])
    group_mode = st.radio(
        "Compendium organization",
        ["Class", 
        #"Class + Subclass (Diacyl vs Lyso)"
        ],
        index=0
    )

uploaded = st.file_uploader("Upload the GPL formula CSV", type=["csv"])

if uploaded and adducts:
    try:
        df_in = pd.read_csv(uploaded)
        st.success(f"Loaded {len(df_in):,} rows.")

        if st.button("Generate compendia", type="primary"):
            with st.spinner("Generating queries and compendia..."):
                df_out = generate_queries(df_in, adducts, tol, intensity_percent)

            st.success(f"Generated {len(df_out):,} queries.")

            st.markdown("### Preview")
            st.dataframe(df_out[["Class", "Subclass", "Title", "Query"]].head(30), use_container_width=True)

            # Build compendia dict: filename -> df
            compendia = {}

            if group_mode == "Class only":
                for cls, subdf in df_out.groupby(["Class"]):
                    cls_name = cls if isinstance(cls, str) else cls[0]
                    fname = f"GPL_{cls_name}_MS1.tsv"
                    compendia[fname] = subdf
            else:
                for (cls, subcls), subdf in df_out.groupby(["Class", "Subclass"]):
                    fname = f"GPL_{cls}_{subcls}_MS1.tsv"
                    compendia[fname] = subdf

            zip_bytes = build_zip(compendia)

            st.download_button(
                "Download ZIP (one compendium TSV per group)",
                data=zip_bytes,
                file_name="GPL_MS1_compendia.zip",
                mime="application/zip"
            )

            # Optional: also download a full metadata table
            meta_csv = df_out.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download full metadata CSV (all queries)",
                data=meta_csv,
                file_name="GPL_MS1_queries_full.csv",
                mime="text/csv"
            )

    except Exception as e:
        st.error(f"Error: {e}")
else:
    st.info("Upload a CSV and select at least one adduct.")
