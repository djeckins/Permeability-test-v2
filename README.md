# Epidermal Barrier Screen

A Streamlit + CLI screening tool for estimating whether compounds fit a physicochemical window associated with epidermal barrier passage.

## Supported inputs

- single SMILES
- SMILES list
- single SDF
- ZIP archive with SDF files

SMILES are always canonicalized before screening.

## What changed in this build

- manual **pH** is respected everywhere (UI and CLI)
- ionization is handled by a **site-aware fallback engine** instead of a single molecule-wide pKa
- the app performs a **live DrugBank lookup** on public DrugBank search/detail pages when a compound name is available
- `logD_7_4` from input files is only used as an override when the selected pH is approximately **7.4**
- overall result logic is:
  - **PASS** = at most 1 non-optimal criterion
  - **BORDERLINE** = more than 1 non-optimal criterion and at most 2 poor criteria
  - **FAIL** = 3 or more poor criteria

## pKa / charge priority order

1. explicit pKa values from the input file
   - `input_pka_acidic`
   - `input_pka_basic`
   - `input_pka`
2. **live DrugBank lookup** by compound name
3. internal **site-aware heuristic fallback**

The live lookup parses public DrugBank drug pages for:
- `pKa (Strongest Acidic)`
- `pKa (Strongest Basic)`
- `Physiological Charge`

If live lookup fails or no confident match is found, the app falls back to the internal site-aware heuristic rules.

### Notes on live DrugBank lookup

- the lookup is designed for a **test version** and depends on the current public DrugBank page structure
- it works best when the input record has a good **compound name**
- it is intentionally conservative: if the match quality is weak, the code falls back instead of forcing a bad pKa
- you can disable the live lookup with:

```bash
export EPIDERMAL_DISABLE_DRUGBANK_LOOKUP=1
```

You can also tune the timeout / candidate depth:

```bash
export EPIDERMAL_DRUGBANK_TIMEOUT=8
export EPIDERMAL_DRUGBANK_MAX_CANDIDATES=5
```

## Optional SDF properties

These are read automatically when present:

- `pKa`, `PKA`, `predicted_pKa`, `input_pka`
- `pKa_strongest_acidic`, `acidic_pKa`, `input_pka_acidic`
- `pKa_strongest_basic`, `basic_pKa`, `input_pka_basic`
- `LogD`, `LOGD`, `logD_7_4`, `input_logd_7_4`
- `_Name`, `Name`, `TITLE`, `ID`, `Compound`

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## CLI examples

```bash
python -m epidermal_barrier_screen.cli --mode smiles --input "CC(=O)Oc1ccccc1C(=O)O" --ph 7.0 --output-prefix results/aspirin
python -m epidermal_barrier_screen.cli --mode smiles_list --input examples/example_smiles_list.txt --ph 5.5 --output-prefix results/list_run
python -m epidermal_barrier_screen.cli --mode sdf --input path/to/file.sdf --ph 7.4 --output-prefix results/sdf_run
python -m epidermal_barrier_screen.cli --mode sdf_zip --input path/to/archive.zip --ph 6.8 --output-prefix results/zip_run
```

The CLI writes both CSV and XLSX.
