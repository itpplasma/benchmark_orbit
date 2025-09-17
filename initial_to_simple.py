#!/usr/bin/env python3
"""Build start.dat from VMEC trajectory data."""
from __future__ import annotations

from pathlib import Path


def find_source_files(root: Path) -> list[Path]:
    pattern = "trajectory_data_vmec_"
    sources = sorted((root / "run_firm3d" / "trapped").glob(f"{pattern}*"))
    if not sources:
        raise FileNotFoundError(
            "no trajectory files found matching run_firm3d/trapped/" f"{pattern}*"
        )
    return sources


def parse_numeric(text: str, *, path: Path, line_no: int) -> float:
    token = text.replace("D", "e").replace("d", "e")
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(
            f"failed to parse numeric value {text!r} in {path}:{line_no}"
        ) from exc


def format_fortran(value: float) -> str:
    return f"{value:.18e}".replace("e", "d")


def first_coordinate(source: Path) -> tuple[str, str, str]:
    with source.open(encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"expected at least 4 columns in {source}:{line_no}, got {len(fields)}"
                )

            s_raw, theta_raw, phi_raw = fields[1:4]
            s_val = parse_numeric(s_raw, path=source, line_no=line_no)
            theta_val = parse_numeric(theta_raw, path=source, line_no=line_no)
            phi_val = parse_numeric(phi_raw, path=source, line_no=line_no)

            return (
                format_fortran(s_val),
                format_fortran(theta_val),
                format_fortran(phi_val),
            )
    raise ValueError(f"no data rows found in {source}")


def build_start_dat(root: Path) -> None:
    sources = find_source_files(root)

    rows: list[str] = []
    for source in sources:
        s_val, theta_val, phi_val = first_coordinate(source)
        rows.append(f"{s_val} {theta_val} {phi_val} 1.0d0 0.0d0")

    if not rows:
        raise ValueError("no data rows found in trajectory sources")

    start_path = root / "start.dat"
    start_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve().parent
    build_start_dat(root)


if __name__ == "__main__":
    main()
