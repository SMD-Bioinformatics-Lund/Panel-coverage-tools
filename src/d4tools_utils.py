from pathlib import Path
import subprocess
from typing import Dict, List
import sys
import traceback


class Coverage:
    def __init__(
        self,
        chr: str,
        start: int,
        end: int,
        cov: float,
        perc_at_thres: Dict[int, float],
        location_str: str,
    ):
        self.chr = chr
        self.start = start
        self.end = end
        self.cov = cov
        self.perc_at_thres = perc_at_thres
        self.location_str = location_str

    def get_loc(self) -> str:
        return f"{self.chr}_{self.start}_{self.end}"

    def __str__(self) -> str:
        return f"{self.chr}_{self.start}_{self.end} {self.cov} {self.perc_at_thres}"


def collect_d4_coverages(
    all_exon_bed_rows: List[str],
    output_dir: Path,
    label: str,
    d4tools_command: str,
    d4_file: Path,
    thresholds: List[int],
) -> Dict[str, Coverage]:
    # Unique rows - only calculate coverage once
    regions_bed = output_dir / f"{label}.bed"
    cov_out = output_dir / f"{label}_coverage.tsv"
    exons_thres_out = output_dir / f"{label}_cov_at_thres.tsv"

    with regions_bed.open("w") as out_fh:
        for bed_row in set(all_exon_bed_rows):
            print(bed_row, file=out_fh)
    cov_results = calculate_coverage(d4tools_command, d4_file, regions_bed, cov_out)

    cov_at_thres_results = calculate_perc_at_thres(
        d4tools_command, d4_file, regions_bed, thresholds, exons_thres_out
    )

    cov_results_parsed = [row.replace("_", "-") for row in cov_results]
    cov_at_thres_parsed = [row.replace("_", "-") for row in cov_at_thres_results]

    cov_dict = get_complete_coverage_dict(cov_results_parsed, cov_at_thres_parsed, thresholds)

    return cov_dict


def run_command(command: List[str]) -> List[str]:
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

    return results.stdout.splitlines()


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
    result_lines = run_command(command)

    return result_lines
    # return [row for row in results.stdout.split("\t") if row != ""]


def get_loc_to_cov_dict(coverage_results: List[str]) -> Dict[str, float]:
    cov_dict: Dict[str, float] = {}
    for cov_line in coverage_results:
        chr, start, end, cov = cov_line.split()
        loc = f"{chr}_{start}_{end}"
        cov_dict[loc] = float(cov)
    return cov_dict


def get_loc_to_perc_at_thres_dict(
    perc_at_thres_results: List[str], thresholds: List[int]
) -> Dict[str, Dict[int, float]]:
    cov_dict: Dict[str, Dict[int, float]] = {}
    for cov_line in perc_at_thres_results:
        fields = cov_line.split("\t")
        chr = fields[0]
        start = fields[1]
        end = fields[2]

        cov_at_thres_dict: Dict[int, float] = {}
        for i, thres in enumerate(thresholds):
            cov_at_thres_dict[thres] = float(fields[3 + i])

        loc = f"{chr}_{start}_{end}"
        cov_dict[loc] = cov_at_thres_dict
    return cov_dict


def get_complete_coverage_dict(
    coverage_results: List[str], perc_at_thres_results: List[str], thresholds: List[int]
) -> Dict[str, Coverage]:
    loc_to_cov_dict = get_loc_to_cov_dict(coverage_results)
    perc_at_thres_dict = get_loc_to_perc_at_thres_dict(perc_at_thres_results, thresholds)

    cov_entries: Dict[str, Coverage] = {}

    for key in loc_to_cov_dict:
        cov = loc_to_cov_dict[key]
        perc_at_thres = perc_at_thres_dict[key]
        chr, start, end = key.split("_")
        location = f"{chr}:{start}-{end}"
        cov_entry = Coverage(chr, int(start), int(end), cov, perc_at_thres, location)
        cov_entries[key] = cov_entry

    return cov_entries


def calculate_perc_at_thres(
    d4tools_command: str, d4_file: Path, bed_path: Path, thresholds: List[int], out_path: Path
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
    result_lines = run_command(command)

    return result_lines
