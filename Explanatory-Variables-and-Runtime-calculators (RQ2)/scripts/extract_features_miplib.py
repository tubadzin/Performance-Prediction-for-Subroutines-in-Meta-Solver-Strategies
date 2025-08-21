#!/usr/bin/env python3
"""
Run `feature-extractor` on every *.mps.gz file in MIP_DIR
and write the features to one CSV (header included).
"""

from __future__ import annotations
from pathlib import Path
import subprocess, csv, concurrent.futures, tqdm, sys, os, datetime
from typing import Optional, Tuple, List

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
MIP_DIR       = Path("../miplib_collection_set")
BIN           = "../binaries/MIPLIB-feature-extractor"
SETTINGS_FILE = Path("../MIPLIB-feature_extractor/settings/presolv-trivial.set")
TIMEOUT_SEC   = 60
N_WORKERS     = max(1, (os.cpu_count() or 2) - 1)

RESULT_DIR    = Path("./results")
dir_tag       = MIP_DIR.name                           # e.g. miplib_benchmark_set
timestamp     = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
OUT_CSV       = RESULT_DIR / f"{dir_tag}_{timestamp}.csv"

# ─────────────────────────────────────────────────────────────────────────────
# Worker helpers
# ─────────────────────────────────────────────────────────────────────────────
def run_extractor(inst_path: Path) -> Tuple[List[str], List[str]]:
    """
    Run `feature-extractor` and return (header, data_row).

    The script's stdout is expected to be two CSV lines: header, values.
    """
    cmd = [BIN, "-p", str(inst_path), "-s", str(SETTINGS_FILE)]
    out = subprocess.check_output(cmd, text=True, timeout=TIMEOUT_SEC)
    reader = list(csv.reader(out.splitlines()))
    if len(reader) < 2:
        raise ValueError(f"Unexpected extractor output for {inst_path.name}")
    header, row = reader[0], reader[1]
    # first column should be the instance name, but add it if extractor didn't
    if header[0].lower() != "instance":
        header = ["instance", *header]
        row    = [inst_path.name, *row]
    return header, row

def process_instance(inst_path: Path) -> Tuple[List[str], List[str]]:
    return run_extractor(inst_path)

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    instances = sorted(MIP_DIR.glob("*.mps.gz"))
    if not instances:
        sys.exit(f"No .mps.gz files found in {MIP_DIR}")

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    header_written = False

    with OUT_CSV.open("w", newline="") as csv_file, \
         concurrent.futures.ProcessPoolExecutor(max_workers=N_WORKERS) as pool:

        fut2inst = {pool.submit(process_instance, p): p for p in instances}
        writer: csv.writer = csv.writer(csv_file)

        for fut in tqdm.tqdm(concurrent.futures.as_completed(fut2inst),
                             total=len(fut2inst),
                             desc="Extracting", unit="inst"):
            inst = fut2inst[fut]
            try:
                header, row = fut.result()

                if not header_written:
                    writer.writerow(header)
                    header_written = True

                writer.writerow(row)
                csv_file.flush()
                rows_written += 1

            except subprocess.TimeoutExpired:
                print(f"[TIMEOUT]     {inst.name} (> {TIMEOUT_SEC}s)", file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"[ERROR]       {inst.name}: {e}", file=sys.stderr)
            except Exception as e:
                print(f"[UNEXPECTED]  {inst.name}: {e}", file=sys.stderr)

    if rows_written == 0:
        OUT_CSV.unlink(missing_ok=True)
        sys.exit("No features extracted – aborting.")

    print(f"✓ {rows_written} rows written to {OUT_CSV}")
