# Epidermal Barrier Screen

A Streamlit + CLI screening tool for estimating whether compounds fit a physicochemical window associated with epidermal barrier passage.

## Supported inputs

- single SMILES
- SMILES list
- single SDF
- ZIP archive with SDF files
- compound names and synonyms (auto-resolved to canonical SMILES)
- InChIKey (resolved to canonical SMILES)

## pKa workflow (production)

1. Resolve input identity to canonical SMILES (for names/synonyms/InChIKey this is automatic).
2. Query databases in order: **PubChem → DrugBank → ChEMBL**.
3. Keep **all database pKa values** with provenance.
4. If and only if no database pKa exists, run **QupKake**.
5. Compute ionization + logD at the user-entered pH.

### Strict policy

- No heuristic/rule-based pKa prediction is used.
- QupKake is the only prediction backend.
- If QupKake is unavailable and no database pKa exists, output remains explicit with warnings (no silent fallback model).

## Output compatibility

The main final table preserves legacy columns and appends extended provenance fields, including:

- `input_name`, `input_type`, `matched_name`
- `canonical_smiles`, `inchikey`, `molecular_formula`
- `source_name`, `source_identifier`, `resolution_confidence`, `resolution_notes`
- `logp_value`, `logp_source`, `logp_evidence_type`
- `pka_values`, `pka_details`, `used_qupkake_fallback`, `warnings`

An exploded per-pKa table is also produced with one row per pKa entry.

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
