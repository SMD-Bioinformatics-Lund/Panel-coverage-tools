from pathlib import Path
import subprocess
from typing import Dict, List
import sys
import traceback


class Coverage:
    def __init__(self, chr: str, start: int, end: int, cov: float, perc_at_thres: Dict[int, float]):
        self.chr = chr
        self.start = start
        self.end = end
        self.cov = cov
        self.perc_at_thres = perc_at_thres

    def get_loc(self) -> str:
        return f"{self.chr}_{self.start}_{self.end}"


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


def get_loc_to_cov_dict(coverage_results: List[str]) -> Dict[str, float]:
    cov_dict = {}
    for cov_line in coverage_results:
        chr, start, end, cov = cov_line.split()
        loc = f"{chr}_{start}_{end}"
        cov_dict[loc] = float(cov)
    return cov_dict


def get_loc_to_perc_at_thres_dict(
    perc_at_thres_results: List[str], thresholds: List[int]
) -> Dict[str, Dict[int, float]]:
    cov_dict = {}
    for cov_line in perc_at_thres_results:
        fields = cov_line.split("\t")
        chr = fields[0]
        start = fields[1]
        end = fields[2]

        cov_at_thres_dict = {}
        for i, thres in enumerate(thresholds):
            cov_at_thres_dict[thres] = fields[3 + i]

        loc = f"{chr}_{start}_{end}"
        cov_dict[loc] = cov_at_thres_dict
    return cov_dict


def get_complete_coverage_dict(
    coverage_results: List[str], perc_at_thres_results: List[str], thresholds: List[int]
) -> Dict[str, Coverage]:
    loc_to_cov_dict = get_loc_to_cov_dict(coverage_results)
    perc_at_thres_dict = get_loc_to_perc_at_thres_dict(perc_at_thres_results, thresholds)

    cov_entries = {}

    for key in loc_to_cov_dict:
        cov = loc_to_cov_dict[key]
        perc_at_thres = perc_at_thres_dict[key]
        chr, start, end = key.split("_")
        cov_entry = Coverage(chr, int(start), int(end), cov, perc_at_thres)
        cov_entries[key] = cov_entry

    return cov_entries


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
