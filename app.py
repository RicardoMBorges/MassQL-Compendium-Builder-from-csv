import io
import re
import zipfile
from typing import Dict, List
from pathlib import Path

import pandas as pd
import streamlit as st


st.set_page_config(page_title="MS1/MS2 MassQL Compendium Builder", layout="wide")

MONO_MASS = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "P": 30.97376199842,
    "S": 31.9720711744,
    "F": 18.99840316273,
    "Cl": 34.968852682,
    "Br": 78.9183376,
    "I": 126.9044719,
    "Na": 22.9897692820,
    "K": 38.9637064864,
    "Mg": 23.985041697,
    "Ca": 39.962590863,
    "Fe": 55.9349375,
    "B": 11.00930536
}

ADDUCT_SHIFT = {
    "M+H": lambda M: M + 1.007276466812,
    "M+Na": lambda M: M + 22.989218,
    "M+K": lambda M: M + 38.963158,
    "M+NH4": lambda M: M + 18.033823,
    "2M+H": lambda M: 2.0 * M + 1.007276466812,
}

FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
FRAG_SPLIT_RE = re.compile(r"[,\s;|]+")


def read_table_autodetect(uploaded_file) -> pd.DataFrame:
    """
    Accepts comma-, semicolon-, or tab-separated text files.
    Uses pandas sep=None + python engine to sniff delimiter,
    and falls back to common separators if needed.
    """
    data = uploaded_file.getvalue()
    text = data.decode("utf-8", errors="replace")

    try:
        return pd.read_csv(io.StringIO(text), sep=None, engine="python")
    except Exception:
        pass

    for sep in ["\t", ";", ","]:
        try:
            return pd.read_csv(io.StringIO(text), sep=sep)
        except Exception:
            continue

    raise ValueError("Could not read the file. Please upload a CSV/TSV with separators: tab, ';', or ','.")


def safe_str(x) -> str:
    return "" if x is None else str(x).strip()


def _normalize_name(s: str) -> str:
    return re.sub(r"[\s_\-]+", "", str(s).strip().lower())


def _find_column_case_insensitive(df: pd.DataFrame, target: str) -> str | None:
    """
    target: canonical name like 'Fragments' or 'neutral loss'
    Returns the actual column name in df if found (case-insensitive,
    ignores spaces/underscores/hyphens), else None.
    """
    t = _normalize_name(target)
    for c in df.columns:
        if _normalize_name(c) == t:
            return c
    return None



def find_first_existing_column(df: pd.DataFrame, candidates: List[str]) -> str | None:
    """
    Return the first matching column name in df using case-insensitive matching
    and ignoring spaces, underscores, and hyphens.
    """
    normalized_df_cols = {_normalize_name(c): c for c in df.columns}
    for cand in candidates:
        c = normalized_df_cols.get(_normalize_name(cand))
        if c is not None:
            return c
    return None



def parse_formula(formula: str) -> Dict[str, int]:
    # Remove trailing charge symbols
    formula = str(formula).strip()
    formula = re.sub(r"[+-]+$", "", formula)

    counts: Dict[str, int] = {}

    for elem, num in FORMULA_TOKEN_RE.findall(formula):
        if elem not in MONO_MASS:
            raise ValueError(f"Unsupported element '{elem}' in formula '{formula}'")

        counts[elem] = counts.get(elem, 0) + (int(num) if num else 1)

    return counts



def exact_mass_from_formula(formula: str) -> float:
    counts = parse_formula(formula)
    return sum(MONO_MASS[e] * n for e, n in counts.items())



def build_generic_title(row: pd.Series, id_col: str) -> str:
    val = safe_str(row.get(id_col, ""))
    if not val or val.lower() == "nan":
        return "Unknown"
    return val.replace(" ", "_")



def build_ms1_clause(mz_list: List[float], tol_ms1: float, ip_ms1: int) -> str:
    mz_str = " OR ".join([f"{mz:.6f}" for mz in mz_list])
    return f"MS2PREC=({mz_str}):TOLERANCEMZ={tol_ms1}:INTENSITYPERCENT={ip_ms1}"



def parse_Fragments_cell(cell: str) -> List[float]:
    """
    Accepts Fragments like:
      "184.0733, 104.1069; 86.0964"
      "184.0733 104.1069 86.0964"
      "[184.0733,104.1069]"
    Returns sorted unique floats.
    """
    s = safe_str(cell)
    if not s or s.lower() == "nan":
        return []
    s = s.strip().strip("[](){}")
    parts = [p for p in FRAG_SPLIT_RE.split(s) if p]
    frags = []
    for p in parts:
        try:
            frags.append(float(p))
        except ValueError:
            pass
    return sorted(set(frags))



def build_ms2_clause(frag_mzs: List[float], tol_ms2: float, ip_ms2: int) -> str:
    """
    AND MS2PROD=(m1 OR m2 OR ...):TOLERANCEMZ=xx:INTENSITYPERCENT=xxx
    """
    if not frag_mzs:
        return ""
    mz_str = " OR ".join([f"{mz:.6f}" for mz in frag_mzs])
    return f" AND MS2PROD=({mz_str}):TOLERANCEMZ={tol_ms2}:INTENSITYPERCENT={ip_ms2}"



def _parse_float_cell(x) -> float | None:
    s = safe_str(x)
    if not s or s.lower() == "nan":
        return None
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None



def build_ms2nl_clause(nl_value: float | None, nl_tol: float, nl_ip: int) -> str:
    """
    AND MS2NL=(xxxxx):TOLERANCEMZ=0.03:INTENSITYPERCENT=5
    """
    if nl_value is None:
        return ""
    return f" AND MS2NL=({nl_value:.6f}):TOLERANCEMZ={nl_tol}:INTENSITYPERCENT={nl_ip}"



def build_full_query(ms1_clause: str, ms2_clause: str, ms1data_key: str = "MS1DATA") -> str:
    return f"QUERY scaninfo({ms1data_key}) WHERE {ms1_clause}{ms2_clause}"



def generate_queries(
    df: pd.DataFrame,
    adducts: List[str],
    tol_ms1: float,
    ip_ms1: int,
    use_ms2: bool,
    tol_ms2: float,
    ip_ms2: int,
    use_nl: bool,
    nl_tol: float,
    nl_ip: int,
) -> pd.DataFrame:
    df = df.copy()

    id_col = find_first_existing_column(df, [
        "Title", "Compound_name", "Compound Name", "Name", "ID", "Identifier"
    ])
    if id_col is None:
        raise ValueError(
            "No identifier column found. Please include one of: "
            "'Title', 'Compound_name', 'Compound Name', 'Name', 'ID', or 'Identifier'."
        )

    formula_col = find_first_existing_column(df, [
        "Formula", "Molecular Formula", "Molecular_Formula", "MolFormula"
    ])
    if formula_col is None:
        raise ValueError(
            "No formula column found. Please include one of: "
            "'Formula', 'Molecular Formula', 'Molecular_Formula', or 'MolFormula'."
        )

    class_col = find_first_existing_column(df, ["Class", "Category", "Group", "Family"])

    df[id_col] = df[id_col].astype(str).str.strip()
    df[formula_col] = df[formula_col].astype(str).str.strip()
    if class_col is not None:
        df[class_col] = df[class_col].astype(str).str.strip()

    df["NeutralMass"] = df[formula_col].apply(exact_mass_from_formula)
    df["Title"] = df.apply(lambda r: build_generic_title(r, id_col), axis=1)

    if df["Title"].duplicated().any():
        seen = {}
        new_titles = []
        for t in df["Title"].tolist():
            seen[t] = seen.get(t, 0) + 1
            new_titles.append(t if seen[t] == 1 else f"{t}_{seen[t]}")
        df["Title"] = new_titles

    def mzs(neutral_mass: float) -> List[float]:
        bad = [a for a in adducts if a not in ADDUCT_SHIFT]
        if bad:
            raise ValueError(f"Unsupported adduct(s): {bad}. Allowed: {sorted(ADDUCT_SHIFT.keys())}")
        return [float(ADDUCT_SHIFT[a](neutral_mass)) for a in adducts]

    df["Adducts"] = ",".join([f"[{a}]+" for a in adducts])
    df["AdductMZs"] = df["NeutralMass"].apply(lambda m: ",".join([f"{x:.6f}" for x in mzs(m)]))

    df["_ms1_mzs"] = df["NeutralMass"].apply(mzs)
    df["_ms1_clause"] = df["_ms1_mzs"].apply(lambda L: build_ms1_clause(L, tol_ms1, ip_ms1))

    frag_col = None
    if use_ms2:
        frag_col = find_first_existing_column(df, ["Fragments", "Fragment", "MS2 Fragments"])
        if frag_col is None:
            raise ValueError("MS2PROD enabled, but file has no 'Fragments' column.")
        df["_frags"] = df[frag_col].apply(parse_Fragments_cell)
        df["_ms2prod_clause"] = df["_frags"].apply(lambda fr: build_ms2_clause(fr, tol_ms2, ip_ms2))
    else:
        df["_ms2prod_clause"] = ""

    nl_col = None
    if use_nl:
        nl_col = find_first_existing_column(df, [
            "Neutral Loss", "Neutral_Loss", "NeutralLoss", "NL"
        ])
        if nl_col is None:
            raise ValueError("MS2NL enabled, but file has no 'Neutral Loss' column.")
        df["_nl_val"] = df[nl_col].apply(_parse_float_cell)
        df["_ms2nl_clause"] = df["_nl_val"].apply(lambda v: build_ms2nl_clause(v, nl_tol, nl_ip))
    else:
        df["_ms2nl_clause"] = ""

    df["_ms2_clause"] = df["_ms2prod_clause"] + df["_ms2nl_clause"]
    df["Query"] = df.apply(lambda r: build_full_query(r["_ms1_clause"], r["_ms2_clause"]), axis=1)

    keep_cols = [
        "Title", "Query", id_col, formula_col, "NeutralMass", "Adducts", "AdductMZs"
    ]
    if class_col is not None and class_col not in keep_cols:
        keep_cols.append(class_col)
    if use_ms2 and frag_col is not None and frag_col not in keep_cols:
        keep_cols.append(frag_col)
    if use_nl and nl_col is not None and nl_col not in keep_cols:
        keep_cols.append(nl_col)

    return df[keep_cols]



def to_compendium_tsv(df: pd.DataFrame) -> bytes:
    out = df[["Title", "Query"]].copy()
    return out.to_csv(sep="\t", index=False).encode("utf-8")



def build_zip(compendia: Dict[str, pd.DataFrame], readme_text: str | None = None) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        if readme_text:
            zf.writestr("README.txt", readme_text)
        for filename, subdf in compendia.items():
            zf.writestr(filename, to_compendium_tsv(subdf))
    return buf.getvalue()


# -----------------------------
# LOGOs (optional)
# -----------------------------
STATIC_DIR = Path(__file__).parent / "static"
for logo_name in ["logo_massQL.png", "LAABio.png"]:
    p = STATIC_DIR / logo_name
    try:
        from PIL import Image
        st.sidebar.image(Image.open(p), use_container_width=True)
    except Exception:
        pass

# -----------------------------
# UI
# -----------------------------
st.title("MS1/MS2 MassQL Compendium Builder")
st.markdown("---")

with st.sidebar:
    st.header("MS1 Settings")
    tol_ms1 = st.number_input("MS1 TOLERANCEMZ", min_value=0.0001, max_value=1.0, value=0.01, step=0.001, format="%.4f")
    ip_ms1 = st.number_input("MS1 INTENSITYPERCENT", min_value=1, max_value=100, value=40, step=10)

    adducts = st.multiselect(
        "Positive-mode adducts",
        ["M+H", "M+Na", "M+K", "M+NH4", "2M+H"],
        default=["M+H", "M+Na", "M+K"]
    )

    st.markdown("---")
    st.header("MS2 (optional)")
    use_ms2 = st.toggle("Add MS2PROD constraints from CSV column 'Fragments'", value=False)
    tol_ms2 = st.number_input("MS2 TOLERANCEMZ", min_value=0.0001, max_value=1.0, value=0.02, step=0.001, format="%.4f", disabled=not use_ms2)
    ip_ms2 = st.number_input("MS2 INTENSITYPERCENT", min_value=1, max_value=100, value=10, step=5, disabled=not use_ms2)

    st.markdown("---")
    st.header("MS2 Neutral Loss (optional)")
    use_nl = st.toggle("Add MS2NL constraint from column 'Neutral Loss'", value=False)
    nl_tol = st.number_input(
        "MS2NL TOLERANCEMZ",
        min_value=0.0001, max_value=1.0, value=0.03, step=0.001, format="%.4f",
        disabled=not use_nl
    )
    nl_ip = st.number_input(
        "MS2NL INTENSITYPERCENT",
        min_value=1, max_value=100, value=5, step=1,
        disabled=not use_nl
    )

    st.markdown("---")
    group_mode = st.radio("Compendium organization", ["Single file", "By Class (if present)"], index=0)

uploaded = st.file_uploader(
    "Upload the formula file (csv/tsv)",
    type=["csv", "tsv", "txt"],
    help=(
        "Upload a CSV/TSV containing at least one identifier column "
        "('Title', 'Compound_name', 'Compound Name', 'Name', 'ID', or 'Identifier') "
        "and one formula column ('Formula' or 'Molecular Formula'). "
        "Optional columns: 'Fragments', 'Neutral Loss', and 'Class'."
    ),
)

if uploaded and adducts:
    try:
        df_in = read_table_autodetect(uploaded)
        st.success(f"Loaded {len(df_in):,} rows.")

        id_in_col = find_first_existing_column(df_in, [
            "Title", "Compound_name", "Compound Name", "Name", "ID", "Identifier"
        ])
        formula_in_col = find_first_existing_column(df_in, [
            "Formula", "Molecular Formula", "Molecular_Formula", "MolFormula"
        ])
        class_in_col = find_first_existing_column(df_in, ["Class", "Category", "Group", "Family"])
        frag_in_col = find_first_existing_column(df_in, ["Fragments", "Fragment", "MS2 Fragments"])
        nl_in_col = find_first_existing_column(df_in, ["Neutral Loss", "Neutral_Loss", "NeutralLoss", "NL"])

        if id_in_col:
            st.info(f"Detected identifier column: '{id_in_col}'.")
        else:
            st.warning("No identifier column found yet.")

        if formula_in_col:
            st.info(f"Detected formula column: '{formula_in_col}'.")
        else:
            st.warning("No formula column found yet.")

        if class_in_col:
            st.info(f"Detected class/group column: '{class_in_col}'.")

        if frag_in_col:
            st.info(f"Detected fragments column: '{frag_in_col}' (MS2PROD constraints available).")
        else:
            st.warning("No fragments column found (MS2PROD toggle will error if enabled).")

        if nl_in_col:
            st.info(f"Detected neutral-loss column: '{nl_in_col}' (MS2NL constraints available).")
        else:
            st.warning("No neutral-loss column found (MS2NL toggle will error if enabled).")

        if st.button("Generate compendia", type="primary"):
            with st.spinner("Generating queries and compendia..."):
                df_out = generate_queries(
                    df_in,
                    adducts=adducts,
                    tol_ms1=tol_ms1,
                    ip_ms1=ip_ms1,
                    use_ms2=use_ms2,
                    tol_ms2=tol_ms2,
                    ip_ms2=ip_ms2,
                    use_nl=use_nl,
                    nl_tol=nl_tol,
                    nl_ip=nl_ip,
                )

            st.success(f"Generated {len(df_out):,} queries.")

            st.markdown("### Preview")

            frag_out_col = find_first_existing_column(df_out, ["Fragments", "Fragment", "MS2 Fragments"])
            nl_out_col = find_first_existing_column(df_out, ["Neutral Loss", "Neutral_Loss", "NeutralLoss", "NL"])

            preview_cols = ["Title", "NeutralMass"]
            if class_in_col and class_in_col in df_out.columns:
                preview_cols.insert(1, class_in_col)
            if use_ms2 and frag_out_col:
                preview_cols.append(frag_out_col)
            if use_nl and nl_out_col:
                preview_cols.append(nl_out_col)
            preview_cols.append("Query")

            st.dataframe(df_out[preview_cols].head(30), use_container_width=True)

            compendia: Dict[str, pd.DataFrame] = {}

            suffix = "MS1"
            if use_ms2:
                suffix += "_MS2PROD"
            if use_nl:
                suffix += "_MS2NL"
            suffix += ".tsv"

            out_class_col = find_first_existing_column(df_out, ["Class", "Category", "Group", "Family"])
            if group_mode == "By Class (if present)" and out_class_col is not None:
                for cls, subdf in df_out.groupby(out_class_col):
                    safe_cls = str(cls).replace(" ", "_")
                    fname = f"Compendium_{safe_cls}_{suffix}"
                    compendia[fname] = subdf
            else:
                compendia[f"MassQL_compendium_all_{suffix}"] = df_out

            adducts_readme = [f"[{a}]+" for a in adducts]
            readme = (
                "MS1/MS2 MassQL Compendia\n"
                f"- MS1: TOLERANCEMZ={tol_ms1}, INTENSITYPERCENT={ip_ms1}\n"
                f"- Adducts={','.join(adducts_readme)}\n"
                f"- Identifier source column: '{id_in_col or 'NOT FOUND'}'\n"
                f"- Formula source column: '{formula_in_col or 'NOT FOUND'}'\n"
                f"- Class/group source column: '{class_in_col or 'NOT FOUND'}'\n"
                f"- MS2PROD enabled={use_ms2}\n"
                f"- MS2NL enabled={use_nl}\n"
            )
            if use_ms2:
                readme += f"- MS2PROD: TOLERANCEMZ={tol_ms2}, INTENSITYPERCENT={ip_ms2}\n"
                readme += f"- MS2PROD fragments source column: '{frag_in_col or 'NOT FOUND'}'\n"
            if use_nl:
                readme += f"- MS2NL: TOLERANCEMZ={nl_tol}, INTENSITYPERCENT={nl_ip}\n"
                readme += f"- MS2NL neutral-loss source column: '{nl_in_col or 'NOT FOUND'}'\n"

            zip_bytes = build_zip(compendia, readme_text=readme)

            st.download_button(
                "Download ZIP (one compendium TSV per group)",
                data=zip_bytes,
                file_name="MassQL_compendia.zip",
                mime="application/zip",
            )

            meta_csv = df_out.to_csv(sep="\t", index=False).encode("utf-8")
            st.download_button(
                "Download full metadata CSV (all queries)",
                data=meta_csv,
                file_name="MassQL_queries_full.csv",
                mime="text/csv",
            )

    except Exception as e:
        st.error(f"Error: {e}")
else:
    st.info("Upload a CSV and select at least one adduct.")
