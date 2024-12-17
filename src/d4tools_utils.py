from pathlib import Path
import subprocess
from typing import List
import sys
import traceback


def run_command(command: List[str]) -> subprocess.CompletedProcess:
    print(f"Executing command: {' '.join(command)}")
    try:
        results = subprocess.run(
            command,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except Exception as err:
        print("An error occured")
        print(err)
        traceback.print_exc()
        sys.exit(1)

    return results


def calculate_coverage(
    d4tools_command: str, d4_file: Path, bed_path: Path, out_path: Path
) -> List[str]:
    command = [
        d4tools_command,
        "stat",
        "--stat",
        "mean",
        "--region",
        str(bed_path),
        str(d4_file),
    ]
    results = run_command(command)

    print(f"Writing results to {out_path}")
    with out_path.open("w") as out_fh:
        print(results.stdout, file=out_fh)

    return results.stdout.splitlines()
    # return [row for row in results.stdout.split("\t") if row != ""]


def calculate_perc_at_thres(
    d4tools_command: str, d4_file: Path, bed_path: Path, thresholds: List[int], out_path
) -> List[str]:
    thresholds_str = ",".join([str(t) for t in thresholds])
    command = [
        d4tools_command,
        "stat",
        "--stat",
        f"perc_cov={thresholds_str}",
        "--region",
        str(bed_path),
        str(d4_file),
    ]
    results = run_command(command)

    print(f"Writing results to {out_path}")
    with out_path.open("w") as out_fh:
        print(results.stdout, file=out_fh)

    return results.stdout.splitlines()
    # return [row for row in results.stdout.split("\t") if row != ""]
