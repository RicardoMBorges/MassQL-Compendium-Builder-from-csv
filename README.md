# MassQL-Compendium-Builder-from-csv
Streamlit app to generate MS1-only MassQL compendia from a CSV containing glycerophospholipid molecular formulas.

## Input CSV requirements
Required columns:
- Class
- FA1
- FA2
- SumDB
- Formula

Optional:
- FA3

## What it generates
- One MS1 MassQL query per row/formula
- Positive-mode adduct support: [M+H]+, [M+Na]+, [M+K]+ (optional [M+NH4]+)
- TOLERANCEMZ and INTENSITYPERCENT are configurable
- Export as a ZIP with compendia organized by:
  - Class only, or
  - Class + Subclass (Diacyl vs Lyso inferred from FA2)

Each compendium TSV contains two columns:
- Title
- Query

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
