from __future__ import annotations

import argparse
from pathlib import Path

from epidermal_barrier_screen.io import parse_input
from epidermal_barrier_screen.screen import screen_records


def main() -> None:
    parser = argparse.ArgumentParser(description="Epidermal barrier screening CLI")
    parser.add_argument("--mode", required=True, choices=["smiles", "smiles_list", "sdf", "sdf_zip"])
    parser.add_argument("--input", required=True, help="SMILES string for smiles mode, file path otherwise")
    parser.add_argument("--output-prefix", required=True, help="Output path prefix without extension")
    parser.add_argument("--ph", type=float, default=5.5, help="Manual pH for ionization/logD calculations")
    args = parser.parse_args()

    if args.mode == "smiles":
        payload: str | bytes = args.input
        filename = None
    else:
        path = Path(args.input)
        payload = path.read_bytes()
        filename = path.name

    records = parse_input(args.mode, payload, filename=filename)
    df = screen_records(records, ph=float(args.ph))

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_prefix.with_suffix(".csv")
    xlsx_path = out_prefix.with_suffix(".xlsx")
    df.to_csv(csv_path, index=False)
    df.to_excel(xlsx_path, index=False)

    print(f"Saved: {csv_path}")
    print(f"Saved: {xlsx_path}")


if __name__ == "__main__":
    main()
