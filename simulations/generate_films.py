#!/usr/bin/env python3
"""Generate V22/V35 films only after matched bulk NPT data files exist."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
Z_BOUNDS = re.compile(
    r"^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)\s+"
    r"([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)\s+zlo\s+zhi\s*$"
)


def config_output(path: Path) -> Path:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        key, separator, value = line.partition("=")
        if separator and key.strip() == "output":
            output = Path(value.strip().strip("\"'"))
            if output.name.startswith("data."):
                return output
    raise ValueError(f"{path} needs an output = .../data.<case> entry")


def bulk_thickness(path: Path) -> float:
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            match = Z_BOUNDS.match(raw)
            if match:
                thickness = float(match.group(2)) - float(match.group(1))
                if thickness <= 0:
                    raise ValueError(f"{path} has invalid z bounds")
                return thickness
    raise ValueError(f"{path} has no zlo zhi box bounds")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formulation", choices=("V22", "V35"), default="V22")
    parser.add_argument("--generator", type=Path, required=True,
                        help="Path to the freshly compiled matching generator")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    generator = args.generator.resolve()
    help_text = subprocess.check_output([str(generator), "--help"], text=True)
    if f"{args.formulation} four-component" not in help_text:
        raise ValueError(f"{generator} is not a {args.formulation} generator")
    series_root = ROOT / "simulations"
    if args.formulation == "V35":
        series_root /= "V35"
    configs = sorted(series_root.glob("*_*wt/model.conf"))
    if len(configs) != 14:
        raise ValueError(f"Expected 14 configs, found {len(configs)}")
    missing = []
    for config in configs:
        output = config_output(config)
        if output.is_absolute():
            raise ValueError(f"Expected repository-relative output in {config}")
        case = output.name[5:]
        bulk = ROOT / output.parent / case / f"data.{case}.npt_eq"
        if not bulk.is_file():
            missing.append(str(bulk.relative_to(ROOT)))
            continue
        thickness = bulk_thickness(bulk)
        film_output = output.parent / f"data.{case}_film"
        command = [str(generator), "--config", str(config.relative_to(ROOT)),
                   "--thickness", f"{thickness:.8f}",
                   "--output", str(film_output)]
        print(f"{config.parent.name}: bulk Lz={thickness:.8f} A -> {film_output}", flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)
    if missing:
        print(f"Skipped {len(missing)} cases without completed bulk NPT data:")
        for path in missing:
            print(f"  {path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"Error: {error}") from None
