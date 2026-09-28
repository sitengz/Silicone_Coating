#!/usr/bin/env python3
"""Block-average the apparent mechanical surface stress of a free coating film."""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path

ATM_ANGSTROM_TO_MN_PER_M = 0.0101325
REQUIRED = ("time_fs", "temp_K", "pxx_atm", "pyy_atm", "pzz_atm",
            "lx_A", "ly_A", "lz_A")


def load_info(path: Path) -> dict:
    info = json.loads(path.read_text(encoding="utf-8"))
    if info.get("geometry") != "film":
        raise ValueError("surface analysis requires a film .info file")
    if not info.get("case_name") or not isinstance(info["case_name"], str):
        raise ValueError(".info file has no valid case_name")
    if not info.get("simulation_template", {}).get("surface_measurement", {}).get("enabled"):
        raise ValueError(".info file does not describe a generated surface-measurement job")
    return info


def read_stress(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    columns: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                candidate = line[1:].strip().split()
                if all(name in candidate for name in REQUIRED):
                    columns = candidate
                continue
            if columns is None:
                raise ValueError(f"{path}:{line_number}: missing named column header")
            values = line.split()
            if len(values) != len(columns):
                raise ValueError(f"{path}:{line_number}: expected {len(columns)} columns")
            try:
                row = dict(zip(columns, map(float, values)))
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-numeric value") from error
            if any(not math.isfinite(value) for value in row.values()):
                raise ValueError(f"{path}:{line_number}: non-finite value")
            if min(row["lx_A"], row["ly_A"], row["lz_A"]) <= 0:
                raise ValueError(f"{path}:{line_number}: non-positive box length")
            if rows and row["time_fs"] <= rows[-1]["time_fs"]:
                raise ValueError(f"{path}:{line_number}: time must strictly increase")
            rows.append(row)
    if len(rows) < 2:
        raise ValueError(f"{path}: at least two stress records are needed")
    return rows


def gamma(row: dict[str, float]) -> float:
    lateral = (row["pxx_atm"] + row["pyy_atm"]) / 2.0
    return (row["lz_A"] / 2.0) * (row["pzz_atm"] - lateral) * ATM_ANGSTROM_TO_MN_PER_M


def make_blocks(rows: list[dict[str, float]], block_ns: float) -> list[dict]:
    width_fs = block_ns * 1.0e6
    origin = rows[0]["time_fs"]
    groups: dict[int, list[dict[str, float]]] = {}
    for row in rows:
        index = int((row["time_fs"] - origin) // width_fs)
        groups.setdefault(index, []).append(row)
    result = []
    for index, group in sorted(groups.items()):
        result.append({
            "start_ns": (origin + index * width_fs) / 1.0e6,
            "end_ns": (origin + (index + 1) * width_fs) / 1.0e6,
            "samples": len(group),
            "mean_mN_per_m": statistics.fmean(gamma(row) for row in group),
        })
    return result


def analyze(info_path: Path, stress_path: Path, output_dir: Path, block_ns: float) -> dict:
    if not math.isfinite(block_ns) or block_ns <= 0:
        raise ValueError("--block-ns must be positive and finite")
    info = load_info(info_path)
    rows = read_stress(stress_path)
    case = info["case_name"]
    values = [gamma(row) for row in rows]
    blocks = make_blocks(rows, block_ns)
    largest_block = max(block["samples"] for block in blocks)
    full_blocks = [block for block in blocks if block["samples"] >= 0.9 * largest_block]
    block_means = [block["mean_mN_per_m"] for block in full_blocks]
    block_se = (statistics.stdev(block_means) / math.sqrt(len(block_means))
                if len(block_means) > 1 else None)
    output_dir.mkdir(parents=True, exist_ok=True)
    series_path = output_dir / f"surface_series.{case}.dat"
    blocks_path = output_dir / f"surface_blocks.{case}.dat"
    summary_path = output_dir / f"surface_summary.{case}.json"
    with series_path.open("w", encoding="utf-8") as handle:
        handle.write("# time_ns temp_K pxx_atm pyy_atm pzz_atm lx_A ly_A lz_A gamma_apparent_mN_per_m\n")
        for row, value in zip(rows, values):
            handle.write(f"{row['time_fs']/1e6:.9g} {row['temp_K']:.9g} "
                         f"{row['pxx_atm']:.9g} {row['pyy_atm']:.9g} {row['pzz_atm']:.9g} "
                         f"{row['lx_A']:.9g} {row['ly_A']:.9g} {row['lz_A']:.9g} {value:.9g}\n")
    with blocks_path.open("w", encoding="utf-8") as handle:
        handle.write("# start_ns end_ns samples mean_gamma_apparent_mN_per_m\n")
        for block in blocks:
            handle.write(f"{block['start_ns']:.9g} {block['end_ns']:.9g} "
                         f"{block['samples']} {block['mean_mN_per_m']:.9g}\n")
    summary = {
        "case_name": case,
        "observable": "apparent mechanical surface stress of two free surfaces",
        "formula": "Lz*(Pzz-(Pxx+Pyy)/2)/2",
        "units": "mN/m",
        "source_info": str(info_path),
        "source_stress": str(stress_path),
        "samples": len(rows),
        "time_ns": [rows[0]["time_fs"] / 1e6, rows[-1]["time_fs"] / 1e6],
        "block_ns": block_ns,
        "blocks": len(blocks),
        "blocks_used_for_standard_error": len(full_blocks),
        "mean_sample_mN_per_m": statistics.fmean(values),
        "mean_full_block_mN_per_m": statistics.fmean(block_means) if block_means else None,
        "full_block_standard_error_mN_per_m": block_se,
        "warning": "Not surface free energy. Check drift, slab integrity, two free surfaces, "
                   "and matched bulk residual stress before interpretation.",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("info", type=Path, help="Matching film .info file")
    parser.add_argument("--file", type=Path, help="Stress file; default: film production file beside .info")
    parser.add_argument("--phase", choices=("prod", "equil"), default="prod")
    parser.add_argument("--block-ns", type=float, default=5.0)
    parser.add_argument("--output-dir", type=Path, help="Default: a case-named folder in the current directory")
    args = parser.parse_args()
    info = load_info(args.info)
    suffix = ".surface_eq.dat" if args.phase == "equil" else ".surface.dat"
    path = args.file or args.info.parent / f"stress.{info['case_name']}{suffix}"
    target = args.output_dir or Path.cwd() / info["case_name"]
    summary = analyze(args.info, path, target, args.block_ns)
    print(f"case: {summary['case_name']}")
    print(f"samples: {summary['samples']}; blocks: {summary['blocks']}")
    print(f"apparent surface stress: {summary['mean_sample_mN_per_m']:.6g} mN/m")
    print(f"results: {target}")
    print(summary["warning"])


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Error: {error}") from None
