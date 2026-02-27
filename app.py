# --- ADD THIS FEATURE TO YOUR EXISTING app.py ---
# MS2 fragment constraints (optional) read from CSV column "fragments"
# Builds an additional clause:
#   AND MS2PROD=(m1 OR m2 OR ...):TOLERANCEMZ=xx:INTENSITYPERCENT=xxx
# MS2 tolerances are independent from MS1 tolerances.

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
}

ADDUCT_SHIFT = {
    "M+H": lambda M: M + 1.007276466812,
    "M+Na": lambda M: M + 22.989218,
    "M+K": lambda M: M + 38.963158,
    "M+NH4": lambda M: M + 18.033823,
    "2M+H": lambda M: 2.0 * M + 1.007276466812,
}

FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)")

def read_table_autodetect(uploaded_file) -> pd.DataFrame:
    """
    Accepts comma-, semicolon-, or tab-separated text files.
    Uses pandas sep=None + python engine to sniff delimiter,
    and falls back to common separators if needed.
    """
    data = uploaded_file.getvalue()  # bytes
    text = data.decode("utf-8", errors="replace")

    # 1) Try delimiter sniffing
    try:
        return pd.read_csv(io.StringIO(text), sep=None, engine="python")
    except Exception:
        pass

    # 2) Fallbacks (ordered)
    for sep in ["\t", ";", ","]:
        try:
            return pd.read_csv(io.StringIO(text), sep=sep)
        except Exception:
            continue

    raise ValueError("Could not read the file. Please upload a CSV/TSV with separators: tab, ';', or ','.")

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
    if fa2 and fa2.lower() != "nan":
        parts.append(fa2)
    if fa3 and fa3.lower() != "nan":
        parts.append(fa3)
    parts.append(f"DB{sumdb}")

    return "_".join(parts).replace(" ", "")


def build_ms1_clause(mz_list: List[float], tol_ms1: float, ip_ms1: int) -> str:
    mz_str = " OR ".join([f"{mz:.6f}" for mz in mz_list])
    return f"MS1MZ=({mz_str}):TOLERANCEMZ={tol_ms1}:INTENSITYPERCENT={ip_ms1}"


# ---- MS2 parsing + clause ----
FRAG_SPLIT_RE = re.compile(r"[,\s;|]+")

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
            # ignore tokens that are not numbers
            pass
    # unique, stable sort
    frags = sorted(set(frags))
    return frags


def build_ms2_clause(frag_mzs: List[float], tol_ms2: float, ip_ms2: int) -> str:
    """
    AND MS2PROD=(m1 OR m2 OR ...):TOLERANCEMZ=xx:INTENSITYPERCENT=xxx
    """
    if not frag_mzs:
        return ""
    mz_str = " OR ".join([f"{mz:.6f}" for mz in frag_mzs])
    return f" AND MS2PROD=({mz_str}):TOLERANCEMZ={tol_ms2}:INTENSITYPERCENT={ip_ms2}"

def build_full_query(ms1_clause: str, ms2_clause: str, ms1data_key: str = "MS1DATA") -> str:
    # Keep your existing structure:
    return f"QUERY scaninfo({ms1data_key}) WHERE {ms1_clause}{ms2_clause}"


def infer_subclass(row: pd.Series) -> str:
    fa2 = safe_str(row.get("FA2", ""))
    return "Lyso" if fa2 == "" or fa2.lower() == "nan" else "Diacyl"


def generate_queries(
    df: pd.DataFrame,
    adducts: List[str],      # e.g. ["M+H","M+Na","M+K","M+NH4","2M+H"]
    tol_ms1: float,
    ip_ms1: int,
    use_ms2: bool,
    tol_ms2: float,
    ip_ms2: int,
    use_nl: bool,
    nl_tol: float,
    nl_ip: int,
) -> pd.DataFrame:
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

    # --- MS1 adduct handling (2M+H is just another selected adduct) ---
    # IMPORTANT: ADDUCT_SHIFT must be a dict of callables:
    #   ADDUCT_SHIFT["M+H"]  = lambda M: M + 1.007276466812
    #   ADDUCT_SHIFT["2M+H"] = lambda M: 2*M + 1.007276466812
    def mzs(neutral_mass: float) -> List[float]:
        bad = [a for a in adducts if a not in ADDUCT_SHIFT]
        if bad:
            raise ValueError(f"Unsupported adduct(s): {bad}. Allowed: {sorted(ADDUCT_SHIFT.keys())}")
        return [float(ADDUCT_SHIFT[a](neutral_mass)) for a in adducts]

    df["Adducts"] = ",".join([f"[{a}]+" for a in adducts])
    df["AdductMZs"] = df["NeutralMass"].apply(lambda m: ",".join([f"{x:.6f}" for x in mzs(m)]))

    # MS1 clause
    df["_ms1_mzs"] = df["NeutralMass"].apply(mzs)
    df["_ms1_clause"] = df["_ms1_mzs"].apply(lambda L: build_ms1_clause(L, tol_ms1, ip_ms1))

    # -------------------------
    # MS2 clauses (independent)
    # -------------------------

    # (A) MS2PROD from fragments
    frag_col = None
    if use_ms2:
        frag_col = _find_column_case_insensitive(df, "fragments")
        if frag_col is None:
            raise ValueError("MS2PROD enabled, but file has no 'Fragments/fragments' column (case-insensitive).")
        df["_frags"] = df[frag_col].apply(parse_Fragments_cell)
        df["_ms2prod_clause"] = df["_frags"].apply(lambda fr: build_ms2_clause(fr, tol_ms2, ip_ms2))
    else:
        df["_ms2prod_clause"] = ""

    # (B) MS2NL from neutral loss
    nl_col = None
    if use_nl:
        nl_col = _find_column_case_insensitive(df, "neutral loss")
        if nl_col is None:
            nl_col = (
                _find_column_case_insensitive(df, "neutral_loss")
                or _find_column_case_insensitive(df, "neutralloss")
            )
        if nl_col is None:
            raise ValueError("MS2NL enabled, but file has no 'Neutral Loss' column (case-insensitive).")

        df["_nl_val"] = df[nl_col].apply(_parse_float_cell)
        df["_ms2nl_clause"] = df["_nl_val"].apply(lambda v: build_ms2nl_clause(v, nl_tol, nl_ip))
    else:
        df["_ms2nl_clause"] = ""

    # Combine
    df["_ms2_clause"] = df["_ms2prod_clause"] + df["_ms2nl_clause"]

    # Final query
    df["Query"] = df.apply(lambda r: build_full_query(r["_ms1_clause"], r["_ms2_clause"]), axis=1)

    # Output columns
    keep_cols = [
        "Title", "Query", "Formula", "NeutralMass", "Adducts", "AdductMZs",
        "Class", "Subclass", "FA1", "FA2", "SumDB"
    ]
    if "FA3" in df.columns:
        keep_cols.insert(10, "FA3")
    if use_ms2 and frag_col is not None:
        keep_cols.append(frag_col)
    if use_nl and nl_col is not None:
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

# ---- MS2 / NL2 parsing + clause ----
def _find_column_case_insensitive(df: pd.DataFrame, target: str) -> str | None:
    """
    target: canonical name like 'Fragments' or 'neutral loss'
    Returns the actual column name in df if found (case-insensitive, ignores spaces/underscores/hyphens), else None.
    """
    def norm(s: str) -> str:
        return re.sub(r"[\s_\-]+", "", str(s).strip().lower())

    t = norm(target)
    for c in df.columns:
        if norm(c) == t:
            return c
    return None


def _parse_float_cell(x) -> float | None:
    s = safe_str(x)
    if not s or s.lower() == "nan":
        return None
    # allow "141.02" or "NL=141.02" (extract first number)
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


# --- In your sidebar UI, add these settings (independent) ---
# Put near your MS2 block:
# use_nl = st.toggle("Add MS2NL constraint from 'Neutral Loss' column", value=False)
# nl_tol = st.number_input("MS2NL TOLERANCEMZ", min_value=0.0001, max_value=1.0, value=0.03, step=0.001, format="%.4f", disabled=not use_nl)
# nl_ip  = st.number_input("MS2NL INTENSITYPERCENT", min_value=1, max_value=100, value=5, step=1, disabled=not use_nl)


# --- Then update your generate_queries signature to include: ---
#   use_nl: bool, nl_tol: float, nl_ip: int
#
# and replace ONLY the MS2 clause assembly part with this block:

# And in keep_cols, add these if enabled:
# if use_ms2 and frag_col is not None: keep_cols.append(frag_col)
# if use_nl and nl_col is not None: keep_cols.append(nl_col)


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
    group_mode = st.radio("Compendium organization", ["Class", "Class + Subclass (Diacyl vs Lyso)"], index=0)

uploaded = st.file_uploader("Upload the formula file (csv/tsv)", type=["csv", "tsv", "txt"])

if uploaded and adducts:
    try:
        # -----------------------------
        # Load input table
        # -----------------------------
        df_in = read_table_autodetect(uploaded)
        st.success(f"Loaded {len(df_in):,} rows.")

        # -----------------------------
        # Detect optional columns (case-insensitive)
        # -----------------------------
        frag_in_col = _find_column_case_insensitive(df_in, "fragments")

        nl_in_col = (
            _find_column_case_insensitive(df_in, "neutral loss")
            or _find_column_case_insensitive(df_in, "neutral_loss")
            or _find_column_case_insensitive(df_in, "neutralloss")
        )

        if frag_in_col:
            st.info(f"Detected fragments column: '{frag_in_col}' (MS2PROD constraints available).")
        else:
            st.warning("No fragments column found (MS2PROD toggle will error if enabled).")

        if nl_in_col:
            st.info(f"Detected neutral-loss column: '{nl_in_col}' (MS2NL constraints available).")
        else:
            st.warning("No neutral-loss column found (MS2NL toggle will error if enabled).")

        # -----------------------------
        # Generate
        # -----------------------------
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

            # -----------------------------
            # Preview
            # -----------------------------
            st.markdown("### Preview")

            frag_out_col = _find_column_case_insensitive(df_out, "fragments")
            nl_out_col = (
                _find_column_case_insensitive(df_out, "neutral loss")
                or _find_column_case_insensitive(df_out, "neutral_loss")
                or _find_column_case_insensitive(df_out, "neutralloss")
            )

            preview_cols = ["Class", "Subclass"]
            if use_ms2 and frag_out_col:
                preview_cols.append(frag_out_col)
            if use_nl and nl_out_col:
                preview_cols.append(nl_out_col)
            preview_cols += ["Title", "Query"]

            st.dataframe(df_out[preview_cols].head(30), use_container_width=True)

            # -----------------------------
            # Build compendia per group
            # -----------------------------
            compendia: Dict[str, pd.DataFrame] = {}

            suffix = "MS1"
            if use_ms2:
                suffix += "_MS2PROD"
            if use_nl:
                suffix += "_MS2NL"
            suffix += ".tsv"

            if group_mode == "Class":
                for cls, subdf in df_out.groupby("Class"):
                    fname = f"Compendium_{cls}_{suffix}"
                    compendia[fname] = subdf
            else:
                for (cls, subcls), subdf in df_out.groupby(["Class", "Subclass"]):
                    fname = f"MassQL_compendia_{cls}_{subcls}_{suffix}"
                    compendia[fname] = subdf

            # -----------------------------
            # README (FIX: define adducts_readme)
            # -----------------------------
            adducts_readme = [f"[{a}]+" for a in adducts]

            readme = (
                "MS1/MS2 MassQL Compendia\n"
                f"- MS1: TOLERANCEMZ={tol_ms1}, INTENSITYPERCENT={ip_ms1}\n"
                f"- Adducts={','.join(adducts_readme)}\n"
                f"- MS2PROD enabled={use_ms2}\n"
                f"- MS2NL enabled={use_nl}\n"
            )
            if use_ms2:
                readme += f"- MS2PROD: TOLERANCEMZ={tol_ms2}, INTENSITYPERCENT={ip_ms2}\n"
                readme += f"- MS2PROD fragments source column: '{frag_in_col or 'NOT FOUND'}'\n"
            if use_nl:
                readme += f"- MS2NL: TOLERANCEMZ={nl_tol}, INTENSITYPERCENT={nl_ip}\n"
                readme += f"- MS2NL neutral-loss source column: '{nl_in_col or 'NOT FOUND'}'\n"

            # -----------------------------
            # Downloads
            # -----------------------------
            zip_bytes = build_zip(compendia, readme_text=readme)

            st.download_button(
                "Download ZIP (one compendium TSV per group)",
                data=zip_bytes,
                file_name="MassQL_compendia.zip",
                mime="application/zip",
            )

            meta_csv = df_out.to_csv(index=False).encode("utf-8")
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