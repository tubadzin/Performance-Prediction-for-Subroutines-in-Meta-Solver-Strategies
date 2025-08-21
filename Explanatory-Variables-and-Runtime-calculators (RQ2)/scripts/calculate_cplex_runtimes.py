"""
Run CPLEX sequentially on every *.mps.gz in MIP_DIR, except those listed
in OPEN_LIST.  Results stream into a CSV while the script runs; each run’s
logs live in results/logs/<timestamp>/.
CSV columns: instance , runtime_sec
"""

from pathlib import Path
import subprocess, csv, tqdm, sys, time, datetime

MIP_DIR   = Path("../miplib_collection_set")
OPEN_LIST = Path("../instances_list/open-v28.test")
CPLEX_BIN = "cplex"          # on PATH add manualy
TIME_LIMIT = 3600            # seconds

RESULT_DIR = Path("./results/calculate_cplex_runtimes_results")

stamp    = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
LOG_DIR  = RESULT_DIR / "logs" / stamp
CSV_PATH = RESULT_DIR / f"{MIP_DIR.name}_{stamp}.csv"

LOG_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)


def load_open_set(txt: Path) -> set[str]:
    return {ln.strip() for ln in txt.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")}


def run_cplex(inst: Path) -> float:
    cmd = [
        CPLEX_BIN, "-c",
        f"read {inst}",
        f"set timelimit {TIME_LIMIT}",
        "optimize",
        "quit"
    ]
    log_file = LOG_DIR / f"{inst.name}.log"
    t0 = time.perf_counter()
    try:
        with log_file.open("w") as log:
            subprocess.run(cmd,
                           stdout=log, stderr=subprocess.STDOUT,
                           text=True, timeout=TIME_LIMIT + 10, check=False)
    except subprocess.TimeoutExpired:
        return TIME_LIMIT
    return time.perf_counter() - t0


if __name__ == "__main__":
    instances = sorted(MIP_DIR.glob("*.mps.gz"))
    if not instances:
        sys.exit(f"No .mps.gz files found in {MIP_DIR}")

    open_set   = load_open_set(OPEN_LIST)

    with CSV_PATH.open("w", newline="") as csv_f:
        writer = csv.writer(csv_f)
        writer.writerow(["instance", "runtime_sec"])      # header

        for p in instances:
            if p.name in open_set:
                writer.writerow([p.name, TIME_LIMIT])
        csv_f.flush()

        for p in tqdm.tqdm(instances, desc="CPLEX", unit="inst"):
            if p.name in open_set:
                continue
            runtime = run_cplex(p)
            writer.writerow([p.name, f"{runtime:.2f}"])
            csv_f.flush()

    solved = sum(1 for _ in CSV_PATH.open()) - 1          # minus header
    print(f"✓ {solved} rows written to {CSV_PATH}")
    print(f"✓ logs saved in {LOG_DIR}")