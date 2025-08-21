#!/usr/bin/env python3
"""
Extract MIPfeature output from *.mps.gz instances and store one row per
instance in a CSV, continuing from an earlier run.

* Uses ThreadPoolExecutor (I/O bound).
* Retries MIPfeature on the .mps if the .lp run segfaults.
* Failed instances are skipped (not written to the CSV).
"""

from __future__ import annotations
from pathlib import Path
import subprocess, csv, concurrent.futures, tqdm, sys, os, datetime, gzip, shutil, tempfile, resource
from typing import Tuple, List, Optional

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
MIP_DIR     = Path("./MIPLIB/miplib_collection_set")
BIN         = Path("./binaries/MIPfeature")
TIMEOUT_SEC = 120           # How long should each explanatory-variable be calculated
N_WORKERS   = 8             # threads
MAX_MPS_MB  = 500           # maximal size of the calculated instances
MEM_LIMIT_MB = 4096         # soft limit inside each worker (macOS/Linux)


CHECKPOINT = Path("PATH TO PREVIOUSLY CALCULATED .csv WITH FEATURES")

RESULT_DIR  = Path("./results")
dir_tag     = MIP_DIR.name
timestamp   = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
OUT_CSV     = RESULT_DIR / f"{dir_tag}_{timestamp}.csv"

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def parse_output(txt: str, inst_name: str) -> Tuple[List[str], List[str]]:
    """turn MIPfeature stdout into CSV header + row (no status column)"""
    lines  = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    h1_idx = next(i for i, ln in enumerate(lines) if ln.startswith("probtype"))
    h2_idx = next(i for i, ln in enumerate(lines) if ln.startswith("time_relax"))

    header1 = next(csv.reader([lines[h1_idx]]))
    data1   = next(csv.reader([lines[h1_idx + 1]]))
    header2 = next(csv.reader([lines[h2_idx]]))
    data2   = next(csv.reader([lines[h2_idx + 1]]))

    header  = ["instance"] + header1 + header2
    row     = [inst_name]  + data1   + data2
    return header, row

def unzip_mps(gz_path: Path, tmpdir: Path) -> Path:
    """*.mps.gz → *.mps"""
    mps_path = tmpdir / gz_path.stem
    with gzip.open(gz_path, "rb") as fin, mps_path.open("wb") as fout:
        shutil.copyfileobj(fin, fout)
    return mps_path

def mps_to_lp(mps_path: Path) -> Path:
    """*.mps → *.lp using SCIP (generic names to silence warnings)"""
    from pyscipopt import Model
    model = Model()
    model.readProblem(str(mps_path))
    lp_path = mps_path.with_suffix(".lp")
    model.writeProblem(str(lp_path), genericnames=True)
    return lp_path

def set_soft_memory_limit(mb: int) -> None:
    """Protect a worker from runaway allocations (best-effort)"""
    limit = mb * 1024 ** 2
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    except ValueError:
        pass  # not supported on some platforms

def run_mipfeature(input_path: Path) -> subprocess.CompletedProcess[str]:
    """call MIPfeature, raising on error/timeout"""
    return subprocess.run(
        [BIN, str(input_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=TIMEOUT_SEC,
        check=True,
    )

# ─────────────────────────────────────────────────────────────────────────────
# Worker: one instance
# ─────────────────────────────────────────────────────────────────────────────
def process_instance(inst_path: Path) -> Optional[Tuple[List[str], List[str]]]:
    """
    Return (header,row) on success.
    Return None if the instance is skipped or fails.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir   = Path(tmp)
        mps_path = unzip_mps(inst_path, tmpdir)

        # Size guard
        if mps_path.stat().st_size > MAX_MPS_MB * 1024 ** 2:
            print(f"[SKIP] {inst_path.name}: >{MAX_MPS_MB} MB", file=sys.stderr)
            return None

        # 1) try LP
        try:
            lp_path = mps_to_lp(mps_path)
            cp      = run_mipfeature(lp_path)
            return parse_output(cp.stdout, inst_path.name)

        # 2) LP segfault → retry original MPS
        except subprocess.CalledProcessError as e:
            if e.returncode < 0:  # signal
                try:
                    cp = run_mipfeature(mps_path)
                    return parse_output(cp.stdout, inst_path.name)
                except Exception:
                    pass
            print(f"[ERROR] {inst_path.name}: RC={e.returncode}", file=sys.stderr)
            return None

        except subprocess.TimeoutExpired:
            print(f"[TIMEOUT] {inst_path.name} (> {TIMEOUT_SEC}s)", file=sys.stderr)
            return None

        except Exception as err:
            print(f"[FAIL] {inst_path.name}: {type(err).__name__}", file=sys.stderr)
            return None

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. checkpoint (optional) ────────────────────────────────────────────
    processed: set[str] = set()
    cp_header: List[str] | None = None
    cp_rows:   List[List[str]]  = []
    if CHECKPOINT.is_file():
        with CHECKPOINT.open(newline="") as f_cp:
            reader = csv.reader(f_cp)
            try:
                cp_header = next(reader)
            except StopIteration:
                cp_header = None
            if cp_header and cp_header[-1].lower() == "status":
                # drop the old status column for compatibility
                cp_header = cp_header[:-1]
            for r in reader:
                if r:
                    if len(r) > len(cp_header):
                        r = r[:-1]            # strip status value
                    processed.add(r[0])
                    cp_rows.append(r)
        print(f"Loaded {len(processed)} rows from checkpoint.")

    # ── 2. decide which instances still need work ───────────────────────────
    all_instances  = sorted(MIP_DIR.glob("*.mps.gz"))
    todo_instances = [p for p in all_instances if p.name not in processed]
    print(f"{len(todo_instances)} of {len(all_instances)} instances still to do.")

    # ── 3. run the work pool and write CSV ──────────────────────────────────
    with OUT_CSV.open("w", newline="") as f_csv:
        writer = csv.writer(f_csv)

        # write checkpoint rows (if any)
        if cp_header:
            writer.writerow(cp_header)
            writer.writerows(cp_rows)
        header_written = bool(cp_header)

        with concurrent.futures.ThreadPoolExecutor(
                max_workers=N_WORKERS,
                initializer=set_soft_memory_limit,
                initargs=(MEM_LIMIT_MB,)) as pool:

            fut2inst = {pool.submit(process_instance, p): p for p in todo_instances}

            for fut in tqdm.tqdm(concurrent.futures.as_completed(fut2inst),
                                 total=len(fut2inst), desc="Extracting", unit="inst"):
                result = fut.result()
                if result is None:
                    continue                      # failed / skipped
                header, row = result

                if not header_written:
                    writer.writerow(header)
                    header_written = True

                writer.writerow(row)
                f_csv.flush()

    print(f"✓ CSV written to {OUT_CSV}")

# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
