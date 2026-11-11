#!/usr/bin/env python3
"""POLYD_Ecac storage-energy density calculator and visualizer."""

from __future__ import annotations

import argparse
import csv
import json
import math
import queue
import re
import sys
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Sequence, Tuple

EPS0 = 8.854e-12
VERSION = "1.1.0"
TOOL_ID = "polyd_ecac"
TOOL_NAME = "POLYD_Ecac"
CREDIT = "power by 陈东洋"
CITATIONS = (
    ("10.1002/adma.202402239", "https://www.x-mol.com/paperRedirect/1771775839390064640"),
    ("10.1016/j.cclet.2025.111436", "https://www.x-mol.com/paperRedirect/1932955793206231040"),
)
MIN_POINTS = 5
COLUMNS = [
    "file_tail_folder", "txt_filename", "source_relative_path",
    "charge energy density", "discharge energy density",
    "charge-discharge efficiency", "Emax", "Qmax", "er1", "er2",
    "loss energy", "error",
]
NUMERIC_COLUMNS = COLUMNS[3:-1]
SKIP_DIR_NAMES = {"output", "logs", "rubbish", "__pycache__"}
RADIANT_LISTING_HEADER = (
    "this is a list of files that are associated in any way with the dataset "
    "and/or the tasks in the dataset."
)
MV_M_UC_CM2_TO_J_CM3 = 0.01
COMPOSITION_COLORS = {
    "PEI_N_1%": "#0072B2",
    "PEI_N_5%": "#E69F00",
    "PEI_N_10%": "#009E73",
    "PEI": "#CC79A7",
    "PEIjian": "#6C6C6C",
}
SAMPLE_MARKERS = ("o", "s", "^", "D", "v", "P", "X", "h", "*")


@dataclass
class RunConfig:
    input_folder: Path
    output_folder: Path
    metadata_type: str = "polyk"


@dataclass
class CurveData:
    x: List[float]
    y: List[float]
    malformed_numeric_lines: List[int]


def is_digit_filename(path: Path) -> bool:
    return path.stem.isdigit()


def is_radiant_file_listing(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            prefix = handle.read(256)
    except OSError:
        return False
    encoding = "utf-16" if prefix.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8-sig"
    header = " ".join(prefix.decode(encoding, errors="replace").split()).casefold()
    return header.startswith(RADIANT_LISTING_HEADER)


def discover_txt_files(root: Path, mode: Optional[str] = None) -> Tuple[List[Path], List[Tuple[Path, str]]]:
    files: List[Path] = []
    skipped: List[Tuple[Path, str]] = []
    for path in root.rglob("*.txt"):
        relative_parts = {part.lower() for part in path.relative_to(root).parts[:-1]}
        if relative_parts.intersection(SKIP_DIR_NAMES):
            skipped.append((path, "excluded_output_or_log_folder"))
        elif not path.is_file():
            skipped.append((path, "not_a_regular_file"))
        elif path.stat().st_size == 0:
            skipped.append((path, "empty_file"))
        elif mode == "radiant" and is_radiant_file_listing(path):
            skipped.append((path, "radiant_dataset_file_listing"))
        else:
            files.append(path)
    files.sort(key=lambda p: (str(p.parent).lower(), 0 if is_digit_filename(p) else 1,
                               int(p.stem) if is_digit_filename(p) else 2**63 - 1, p.name.lower()))
    return files, skipped


def iter_txt_files(root: Path) -> List[Path]:
    """Compatibility helper returning only non-empty candidate curves."""
    return discover_txt_files(root)[0]


def parse_line_to_xy(line: str, mode: str) -> Optional[Tuple[float, float]]:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None
    parts = stripped.split()
    candidates = [(0, 1)] if mode == "radiant" or len(parts) < 3 else [(1, 2), (0, 1)]
    for ix, iy in candidates:
        if len(parts) <= max(ix, iy):
            continue
        try:
            x, y = float(parts[ix]), float(parts[iy])
        except ValueError:
            continue
        if math.isfinite(x) and math.isfinite(y):
            return x, y
    return None


def looks_like_numeric_data(line: str) -> bool:
    return bool(re.match(r"^\s*[+-]?(?:\d|\.\d)", line))


def read_xy_flexible(path: Path, mode: str) -> CurveData:
    for encoding in ("utf-16", "utf-8-sig", "utf-8", "mbcs"):
        try:
            text = path.read_text(encoding=encoding, errors="strict")
        except (OSError, UnicodeError):
            continue
        xs: List[float] = []
        ys: List[float] = []
        malformed: List[int] = []
        for line_number, line in enumerate(text.splitlines(), start=1):
            xy = parse_line_to_xy(line, mode)
            if xy is None:
                if xs and looks_like_numeric_data(line):
                    malformed.append(line_number)
                continue
            xs.append(xy[0])
            ys.append(xy[1])
        if len(xs) >= 2:
            return CurveData(xs, ys, malformed)
    raise ValueError("cannot parse two finite X/Y numeric columns")


def trapz(y_values: Sequence[float], x_values: Sequence[float]) -> float:
    if len(y_values) < 2 or len(y_values) != len(x_values):
        return math.nan
    return sum(((y_values[i] + y_values[i + 1]) * 0.5) * (x_values[i + 1] - x_values[i])
               for i in range(len(y_values) - 1))


def compute_eps_r(energy_j_cm3: float, emax_mv_m: float) -> float:
    if not math.isfinite(energy_j_cm3) or not math.isfinite(emax_mv_m) or emax_mv_m == 0:
        return math.nan
    return (2.0 * energy_j_cm3 * 1_000_000.0) / (EPS0 * (emax_mv_m * 1_000_000.0) ** 2)


def source_folder(path: Path, input_root: Optional[Path]) -> str:
    if input_root:
        try:
            return path.parent.relative_to(input_root).as_posix()
        except ValueError:
            pass
    return path.parent.as_posix()


def source_path(path: Path, input_root: Optional[Path]) -> str:
    if input_root:
        try:
            return path.relative_to(input_root).as_posix()
        except ValueError:
            pass
    return path.as_posix()


def tail_folder(path: Path) -> str:
    parts = path.parent.parts
    return "/".join(parts[-3:] if len(parts) >= 3 else parts)


def fmt(value: float) -> str:
    if value is None or not math.isfinite(value):
        return "NaN"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def empty_row(path: Path, error: str, input_root: Optional[Path] = None) -> dict:
    row = {column: "NaN" for column in COLUMNS}
    row.update({"file_tail_folder": source_folder(path, input_root), "txt_filename": path.name,
                "source_relative_path": source_path(path, input_root), "error": error})
    return row


def interpolate_zero(x: Sequence[float], y: Sequence[float], start: int) -> Optional[Tuple[int, List[float], List[float]]]:
    for index in range(start, len(x) - 1):
        x1, x2 = x[index], x[index + 1]
        if x1 == 0:
            return index, list(x[:index + 1]), list(y[:index + 1])
        if x1 * x2 <= 0:
            fraction = x1 / (x1 - x2) if x1 != x2 else 0.0
            y_zero = y[index] + fraction * (y[index + 1] - y[index])
            return index + 1, list(x[:index + 1]) + [0.0], list(y[:index + 1]) + [y_zero]
    return None


def process_curve(path: Path, data: CurveData, mode: str, input_root: Optional[Path] = None) -> dict:
    if len(data.x) < MIN_POINTS:
        return empty_row(path, f"too_few_points: need at least {MIN_POINTS}", input_root)
    if data.malformed_numeric_lines:
        return empty_row(path, "malformed_numeric_rows: " + ",".join(map(str, data.malformed_numeric_lines[:8])), input_root)
    if mode == "polyk":
        x = [(value - data.x[0]) / 1_000_000.0 for value in data.x]
        y = [(value - data.y[0]) * 100.0 for value in data.y]
    else:
        # Radiant columns are kV/cm and uC/cm^2.
        x = [(value - data.x[0]) * 100000.0 / 1_000_000.0 for value in data.x]
        y = [value - data.y[0] for value in data.y]
    positive_peak_index = max(range(len(x)), key=lambda index: x[index])
    negative_peak_index = min(range(len(x)), key=lambda index: x[index])
    # Preserve the established positive-field branch when it exists. A fully
    # negative acquisition is supported by reflecting both axes.
    if x[positive_peak_index] > 0:
        peak_index, polarity = positive_peak_index, 1.0
    else:
        peak_index, polarity = negative_peak_index, -1.0
    oriented_x = [value * polarity for value in x]
    oriented_y = [value * polarity for value in y]
    if peak_index < 1 or oriented_x[peak_index] <= 0:
        return empty_row(path, "segment_failed: no nonzero field peak", input_root)
    zero_segment = interpolate_zero(oriented_x, oriented_y, peak_index)
    if zero_segment is None:
        return empty_row(path, "incomplete_discharge: field does not return to zero", input_root)
    zero_index, segment_x, segment_y = zero_segment
    if zero_index - peak_index < 1:
        return empty_row(path, "segment_failed: too few discharge points", input_root)
    charge_x, charge_y = oriented_x[:peak_index + 1], oriented_y[:peak_index + 1]
    discharge_x, discharge_y = segment_x[peak_index:], segment_y[peak_index:]
    # Integrate E dP. Both modes now express E in MV/m and P in uC/cm^2.
    charge = abs(trapz(charge_x, charge_y)) * MV_M_UC_CM2_TO_J_CM3
    discharge = abs(trapz(discharge_x, discharge_y)) * MV_M_UC_CM2_TO_J_CM3
    emax, q_at_emax = oriented_x[peak_index], oriented_y[peak_index]
    efficiency = discharge / charge if charge > 0 else math.nan
    loss = charge - discharge
    if not all(math.isfinite(value) for value in (charge, discharge, emax, q_at_emax, efficiency, loss)):
        return empty_row(path, "invalid_numeric_result", input_root)
    if discharge > charge * 1.005:
        return empty_row(path, "invalid_energy_balance: discharge exceeds charge", input_root)
    return {
        "file_tail_folder": source_folder(path, input_root), "txt_filename": path.name,
        "source_relative_path": source_path(path, input_root),
        "charge energy density": fmt(charge), "discharge energy density": fmt(discharge),
        "charge-discharge efficiency": fmt(efficiency), "Emax": fmt(emax), "Qmax": fmt(q_at_emax),
        "er1": fmt(compute_eps_r(charge, emax)), "er2": fmt(compute_eps_r(discharge, emax)),
        "loss energy": fmt(loss), "error": "",
    }


def process_polyk(path: Path, x_raw: Sequence[float], y_raw: Sequence[float]) -> dict:
    return process_curve(path, CurveData(list(x_raw), list(y_raw), []), "polyk")


def process_radiant(path: Path, x_raw: Sequence[float], y_raw: Sequence[float]) -> dict:
    return process_curve(path, CurveData(list(x_raw), list(y_raw), []), "radiant")


def process_one(path: Path, mode: str, input_root: Optional[Path] = None) -> dict:
    try:
        return process_curve(path, read_xy_flexible(path, mode), mode, input_root)
    except Exception as exc:
        return empty_row(path, f"read_failed: {exc}", input_root)


def write_csv(rows: Iterable[dict], output_csv: Path) -> None:
    with output_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in COLUMNS})


def allocate_output_folder(requested: Path) -> Path:
    requested.mkdir(parents=True, exist_ok=True)
    if not any(requested.iterdir()):
        return requested
    stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    candidate, suffix = requested / stamp, 1
    while candidate.exists():
        candidate, suffix = requested / f"{stamp}_{suffix}", suffix + 1
    candidate.mkdir()
    return candidate


def write_skipped_csv(entries: Sequence[Tuple[Path, str]], root: Path, output_dir: Path) -> None:
    with (output_dir / "skipped_files.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=("source_relative_path", "reason"))
        writer.writeheader()
        for path, reason in entries:
            writer.writerow({"source_relative_path": source_path(path, root), "reason": reason})


def run_calculation(config: RunConfig, logger: Optional[Callable[[str], None]] = None) -> dict:
    if config.metadata_type not in {"polyk", "radiant"}:
        raise ValueError("metadataType must be polyk or radiant")
    if not config.input_folder.is_dir():
        raise FileNotFoundError(f"Input folder not found: {config.input_folder}")
    output_dir = allocate_output_folder(config.output_folder)
    files, skipped = discover_txt_files(config.input_folder, config.metadata_type)
    rows = []
    for index, path in enumerate(files, start=1):
        rows.append(process_one(path, config.metadata_type, config.input_folder))
        if logger and (index == 1 or index % 25 == 0 or index == len(files)):
            logger(f"Processed {index}/{len(files)} files")
    output_csv = output_dir / "segment_results.csv"
    write_csv(rows, output_csv)
    write_skipped_csv(skipped, config.input_folder, output_dir)
    error_count = sum(1 for row in rows if row.get("error"))
    valid_count = len(rows) - error_count
    result = {
        "ok": bool(valid_count), "status": "ok" if valid_count else "no_valid_curves",
        "tool": TOOL_ID, "version": VERSION, "metadata_type": config.metadata_type,
        "input_folder": str(config.input_folder), "requested_output_folder": str(config.output_folder),
        "output_folder": str(output_dir), "output_csv": str(output_csv), "file_count": len(rows),
        "valid_count": valid_count, "error_count": error_count, "skipped_count": len(skipped),
    }
    if logger:
        logger(f"Saved: {output_csv}")
        logger(f"Valid={valid_count}, errors={error_count}, skipped={len(skipped)}")
    return result


def load_input_json(path: Path, output_arg: Path) -> RunConfig:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    input_value = data.get("inputFolder") or data.get("input_folder")
    if not input_value:
        raise ValueError("input JSON must include inputFolder")
    input_folder = Path(input_value)
    if not input_folder.is_absolute():
        input_folder = (path.parent / input_folder).resolve()
    return RunConfig(input_folder, output_arg.resolve(), data.get("metadataType") or data.get("metadata_type") or "polyk")


def write_json(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def command_run(args: argparse.Namespace) -> int:
    logs: List[str] = []
    try:
        config = load_input_json(Path(args.input).resolve(), Path(args.output))
        result = run_calculation(config, logger=logs.append)
        output_dir = Path(result["output_folder"])
        write_json(output_dir / "result.json", result)
        (output_dir / "run.log").write_text("\n".join(logs) + "\n", encoding="utf-8")
        print(json.dumps(result))
        return 0 if result["ok"] else 2
    except Exception as exc:
        print(json.dumps({"ok": False, "tool": TOOL_ID, "version": VERSION, "error": str(exc)}), file=sys.stderr)
        return 1


def command_figure(args: argparse.Namespace) -> int:
    try:
        input_csv = Path(args.input).resolve()
        with input_csv.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        requested = Path(args.output).resolve()
        output_dir = allocate_output_folder(requested)
        result = export_publication_figures(rows, output_dir)
        result.update({"ok": True, "tool": TOOL_ID, "version": VERSION, "input_csv": str(input_csv)})
        print(json.dumps(result))
        return 0
    except Exception as exc:
        print(json.dumps({"ok": False, "tool": TOOL_ID, "version": VERSION, "error": str(exc)}), file=sys.stderr)
        return 1


def command_health_check() -> int:
    frozen = getattr(sys, "frozen", False)
    try:
        import matplotlib
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        visualization = {"available": True, "matplotlib_version": matplotlib.__version__, "error": ""}
    except Exception as exc:
        visualization = {"available": False, "matplotlib_version": "", "error": str(exc)}
    print(json.dumps({"ok": True, "tool": TOOL_ID, "name": TOOL_NAME, "version": VERSION,
                      "python": sys.version.split()[0], "frozen": frozen,
                      "executable": sys.executable, "credit": CREDIT,
                      "script": str(Path(__file__).resolve()),
                      "visualization": visualization}))
    return 0


def show_help() -> None:
    print(f"""POLYD_Ecac storage-energy density tool ({VERSION})

Usage:
  POLYD_Ecac-{VERSION}-win-x64.exe
  POLYD_Ecac-{VERSION}-win-x64.exe --health-check
  POLYD_Ecac-{VERSION}-win-x64.exe --licenses
  POLYD_Ecac-{VERSION}-win-x64.exe run --input input.json --output output/run
  POLYD_Ecac-{VERSION}-win-x64.exe figure --input output/run/segment_results.csv --output output/publication

Input JSON paths are resolved relative to the JSON file.
Existing output folders receive a timestamped run subfolder to preserve prior results.
""")


def to_float(value: object) -> Optional[float]:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def legend_label(group: str, depth: int) -> str:
    parts = [part for part in group.replace("\\", "/").split("/") if part]
    return "/".join(parts[-max(1, min(6, depth)):]) if parts else group


def group_colors(groups: Sequence[str]) -> dict:
    return {group: composition_color(publication_group_metadata(group)[1]) for group in sorted(set(groups))}


PLOT_LABELS = {
    "Emax": "Maximum electric field (MV m⁻¹)",
    "Qmax": "Polarization at maximum field (µC cm⁻²)",
    "charge energy density": "Charged energy density (J cm⁻³)",
    "discharge energy density": "Discharged energy density (J cm⁻³)",
    "charge-discharge efficiency": "Charge–discharge efficiency (%)",
    "er1": "Relative permittivity (charge)",
    "er2": "Relative permittivity (discharge)",
    "loss energy": "Energy loss (J cm⁻³)",
}
DEFAULT_Y_COLUMNS = ("discharge energy density", "charge-discharge efficiency")
EFFICIENCY_COLUMN = "charge-discharge efficiency"
PLOT_UNITS = {
    "Emax": "MV m⁻¹",
    "Qmax": "µC cm⁻²",
    "charge energy density": "J cm⁻³",
    "discharge energy density": "J cm⁻³",
    "loss energy": "J cm⁻³",
    EFFICIENCY_COLUMN: "%",
}
GROUP_PALETTE = ("#197D90", "#A44F2A", "#6753A1", "#567D30", "#B2476A", "#3166AC", "#826129")


def composition_color(composition: str) -> str:
    from hashlib import blake2s

    if composition in COMPOSITION_COLORS:
        return COMPOSITION_COLORS[composition]
    digest = blake2s(composition.casefold().encode("utf-8"), digest_size=2).digest()
    return GROUP_PALETTE[int.from_bytes(digest, "big") % len(GROUP_PALETTE)]


def build_plot_series(rows: Sequence[dict], groups: Sequence[str], x_column: str, y_columns: Sequence[str], colors: dict, legend_depth: int) -> List[dict]:
    from collections import Counter

    series: List[dict] = []
    group_rows = {group: [] for group in groups}
    for row in rows:
        group = row.get("file_tail_folder")
        if group in group_rows and not row.get("error"):
            group_rows[group].append(row)

    markers = group_markers(colors)
    labels = ["Root folder" if group == "." else legend_label(group, legend_depth) for group in groups]
    label_counts = Counter(labels)
    for group, short_label in zip(groups, labels):
        label = group if label_counts[short_label] > 1 else short_label
        for y_index, y_column in enumerate(y_columns):
            points = [(x_value, y_value) for row in group_rows[group]
                      for x_value, y_value in [(to_float(row.get(x_column)), to_float(row.get(y_column)))]
                      if x_value is not None and y_value is not None]
            points.sort(key=lambda pair: pair[0])
            if points:
                series.append({"group": group, "label": label, "metric": y_column,
                               "axis": "left" if y_index == 0 else "right",
                               "color": colors[group], "marker": markers[group], "points": points})
    return series


def group_markers(colors: dict) -> dict:
    markers = {}
    counts = {}
    for group, color in colors.items():
        count = counts.get(color, 0)
        markers[group] = SAMPLE_MARKERS[count % len(SAMPLE_MARKERS)]
        counts[color] = count + 1
    return markers


def make_figure(series: Sequence[dict], x_column: str, y_columns: Sequence[str], layout: str = "split", *, screen: bool = False):
    from textwrap import fill

    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, PercentFormatter, ScalarFormatter

    present = {item["metric"] for item in series}
    metrics = [column for column in y_columns if column in present]
    if not metrics:
        raise ValueError("No numeric data for the selected axes and groups.")
    overlay = layout == "dual" and len(metrics) == 2
    n_panels = 1 if overlay else len(metrics)
    height = max(4.72, 1.65 * n_panels) if not screen else max(4.5, 1.35 * n_panels)
    figure = Figure(figsize=(7.6 if screen else 15 / 2.54, height),
                    dpi=100, layout="constrained", facecolor="white")
    axes = list(figure.subplots(n_panels, 1, sharex=True, squeeze=False)[:, 0])
    metric_axes = axes + [axes[0].twinx()] if overlay else axes
    axes_by_metric = dict(zip(metrics, metric_axes))
    label_size = 10 if screen else 8
    tick_size = 9 if screen else 7
    for item in series:
        axis = axes_by_metric[item["metric"]]
        x_values, y_values = zip(*item["points"])
        dash = "--" if overlay and axis is metric_axes[-1] else "-"
        axis.plot(x_values, y_values, color=item["color"], linewidth=1.25, alpha=0.88,
                  linestyle=dash if len(x_values) > 1 else "None", marker=item["marker"],
                  markersize=4, markeredgewidth=0.4, markevery=max(1, len(x_values) // 12))
    for column, axis in axes_by_metric.items():
        label = PLOT_LABELS.get(column, column)
        unit = PLOT_UNITS.get(column)
        if screen and not overlay:
            name = label.removesuffix(f" ({unit})") if unit else label
            axis.set_title(name, loc="left", fontsize=label_size, color="#283446", pad=4)
            axis.set_ylabel(unit or "", fontsize=label_size, color="#283446")
        else:
            if unit:
                label = f"{label.removesuffix(f' ({unit})')}\n({unit})"
            axis.set_ylabel(label, fontsize=label_size, color="#283446")
        axis.tick_params(axis="both", labelsize=tick_size, width=0.75, length=3.5, colors="#475569")
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        axis.yaxis.set_major_formatter(PercentFormatter(1.0) if column == EFFICIENCY_COLUMN else ScalarFormatter(useOffset=False))
        for spine in axis.spines.values():
            spine.set_linewidth(0.75)
            spine.set_color("#798595")
        axis.spines["top"].set_visible(False)
        if overlay and axis is metric_axes[-1]:
            axis.spines["left"].set_visible(False)
            axis.spines["bottom"].set_visible(False)
        else:
            axis.spines["right"].set_visible(False)
            if screen:
                axis.grid(axis="y", color="#e9eef5", linewidth=0.65)
                axis.set_axisbelow(True)
    axes[-1].set_xlabel(PLOT_LABELS.get(x_column, x_column), fontsize=label_size, color="#283446")
    axes[-1].xaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=3))
    axes[-1].xaxis.set_major_formatter(PercentFormatter(1.0) if x_column == EFFICIENCY_COLUMN else ScalarFormatter(useOffset=False))
    if not overlay:
        figure.align_ylabels(axes)

    group_items = {item["group"]: item for item in series}
    if not screen or len(group_items) <= 8:
        handles = [
            Line2D([0], [0], color=item["color"], marker=item["marker"], linewidth=1.2,
                   markersize=4, label=fill(item["label"], 36))
            for item in group_items.values()
        ]
        if overlay:
            handles.extend([
                Line2D([0], [0], color="#475569", label="Left axis", linewidth=1.1),
                Line2D([0], [0], color="#475569", label="Right axis", linewidth=1.1, linestyle="--"),
            ])
        figure.legend(handles=handles, loc="outside upper center", ncols=min(len(handles), 2),
                      fontsize=tick_size, frameon=False, handlelength=1.8, columnspacing=1.3)
    return figure


def export_figures(target: Path, series: Sequence[dict], x_column: str, y_columns: Sequence[str], layout: str = "split") -> List[Path]:
    from matplotlib import rc_context

    if not series:
        raise RuntimeError("No plottable data for the selected axes and groups.")
    target.parent.mkdir(parents=True, exist_ok=True)
    groups = list(dict.fromkeys(item["group"] for item in series))
    pages = [groups[index:index + 6] for index in range(0, len(groups), 6)]
    output_paths: List[Path] = []
    with rc_context(publication_rcparams()):
        for index, page_groups in enumerate(pages, start=1):
            page_series = [item for item in series if item["group"] in page_groups]
            page_target = target if len(pages) == 1 else target.with_name(f"{target.stem}_{index:02d}{target.suffix}")
            figure = make_figure(page_series, x_column, y_columns, layout)
            figure.savefig(page_target, dpi=300, facecolor="white")
            output_paths.append(page_target)
    return output_paths


def publication_group_metadata(group: str) -> Tuple[str, str, str]:
    """Extract temperature, composition, and sample ID from a full CSV group path."""
    parts = [part for part in group.replace("\\", "/").split("/") if part]
    temperature = "unknown"
    # Temperatures are directories before the composition/sample leaf pair.
    for part in parts[:-2]:
        if re.fullmatch(r"\d{2,3}", part):
            temperature = part
            break
    composition = parts[-2] if len(parts) >= 2 else group
    sample = parts[-1] if parts else "unknown"
    return temperature, composition, sample


def publication_marker(sample: str) -> str:
    try:
        return SAMPLE_MARKERS[(int(sample) - 1) % len(SAMPLE_MARKERS)]
    except ValueError:
        return SAMPLE_MARKERS[sum(ord(character) for character in sample) % len(SAMPLE_MARKERS)]


def publication_sort_key(value: str) -> Tuple[int, str]:
    try:
        return int(value), value
    except ValueError:
        return 10**9, value


def publication_limits(values: Sequence[float], include_zero: bool = False) -> Tuple[float, float]:
    low, high = min(values), max(values)
    if include_zero:
        low = 0.0
    span = max(high - low, abs(high) * 0.04, 0.002)
    return low - (0.0 if include_zero else span * 0.08), high + span * 0.08


def publication_rcparams() -> dict:
    return {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }


def publication_rows(rows: Sequence[dict]) -> List[dict]:
    prepared: List[dict] = []
    for row in rows:
        if row.get("error"):
            continue
        field = to_float(row.get("Emax"))
        discharge = to_float(row.get("discharge energy density"))
        efficiency = to_float(row.get("charge-discharge efficiency"))
        group = row.get("file_tail_folder", "")
        if None in (field, discharge, efficiency) or not group:
            continue
        temperature, composition, sample = publication_group_metadata(group)
        prepared.append({
            **row,
            "temperature_C": temperature,
            "composition": composition,
            "sample_id": sample,
            "sample_marker": publication_marker(sample),
            "Emax_value": field,
            "discharge_value": discharge,
            "efficiency_percent": efficiency * 100.0,
        })
    return prepared


def write_publication_source_data(rows: Sequence[dict], path: Path) -> None:
    fields = list(COLUMNS) + ["temperature_C", "composition", "sample_id", "sample_marker", "efficiency_percent"]
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            source = {key: row.get(key, "") for key in fields}
            source["efficiency_percent"] = f"{row['efficiency_percent']:.6f}".rstrip("0").rstrip(".")
            writer.writerow(source)


def export_publication_figures(rows: Sequence[dict], output_dir: Path) -> dict:
    """Export raw-curve, temperature-aware energy/efficiency comparison figures."""
    from matplotlib import rc_context

    with rc_context(publication_rcparams()):
        return render_publication_figures(rows, output_dir)


def render_publication_figures(rows: Sequence[dict], output_dir: Path) -> dict:
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    prepared = publication_rows(rows)
    if not prepared:
        raise RuntimeError("No valid Emax, discharged-energy-density, and efficiency rows are available")
    temperatures = sorted({row["temperature_C"] for row in prepared}, key=publication_sort_key)
    if len(temperatures) != 2 or set(temperatures) != {"150", "200"}:
        raise RuntimeError("Publication export currently requires both 150 and 200 temperature groups")

    output_dir.mkdir(parents=True, exist_ok=True)
    x_limits = publication_limits([row["Emax_value"] for row in prepared], include_zero=True)
    energy_limits = publication_limits([row["discharge_value"] for row in prepared], include_zero=True)
    efficiency_limits = publication_limits([row["efficiency_percent"] for row in prepared])
    composition_order = {name: index for index, name in enumerate(COMPOSITION_COLORS)}
    compositions = sorted({row["composition"] for row in prepared},
                          key=lambda value: (composition_order.get(value, len(composition_order)), value))
    legend_columns = min(5, len(compositions))
    legend_rows = (len(compositions) + legend_columns - 1) // legend_columns
    figure = Figure(figsize=(183 / 25.4, 150 / 25.4), dpi=100, facecolor="white")
    figure.subplots_adjust(left=0.11, right=0.98, bottom=0.105, top=max(0.68, 0.90 - 0.035 * (legend_rows - 1)),
                           hspace=0.34, wspace=0.36)
    axes = figure.subplots(2, 2, sharex=True, sharey="col")
    metrics = (("discharge_value", PLOT_LABELS[DEFAULT_Y_COLUMNS[0]], energy_limits),
               ("efficiency_percent", PLOT_LABELS[DEFAULT_Y_COLUMNS[1]], efficiency_limits))
    labels = ("a", "b", "c", "d")
    for row_index, temperature in enumerate(temperatures):
        temperature_rows = [row for row in prepared if row["temperature_C"] == temperature]
        groups = sorted({row["file_tail_folder"] for row in temperature_rows})
        for column_index, (metric, y_label, y_limits) in enumerate(metrics):
            axis = axes[row_index][column_index]
            for group in groups:
                group_rows = sorted((row for row in temperature_rows if row["file_tail_folder"] == group),
                                    key=lambda row: row["Emax_value"])
                first = group_rows[0]
                axis.plot([row["Emax_value"] for row in group_rows], [row[metric] for row in group_rows],
                           color=composition_color(first["composition"]), linewidth=1.05,
                           alpha=0.84, marker=first["sample_marker"], markersize=3.1, markeredgewidth=0.35,
                          markevery=max(1, len(group_rows) // 6), zorder=2)
            axis.set_xlim(x_limits)
            axis.set_ylim(y_limits)
            axis.set_ylabel(y_label, fontsize=7)
            if row_index == len(temperatures) - 1:
                axis.set_xlabel(PLOT_LABELS["Emax"], fontsize=7)
            axis.tick_params(axis="both", labelsize=6.5, width=0.8, length=3)
            axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
            axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.text(-0.13, 1.03, labels[row_index * 2 + column_index], transform=axis.transAxes,
                      fontsize=8, fontweight="bold", va="bottom")
            if column_index == 0:
                axis.set_title(f"{temperature} °C · {len(groups)} curves",
                               fontsize=7, loc="left", pad=8)
    handles = [Line2D([0], [0], color=composition_color(composition), linewidth=1.5,
                      label=composition) for composition in compositions]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.55, 0.985), ncols=legend_columns,
                  fontsize=6.5, handlelength=2.1, columnspacing=1.3)
    base = output_dir / "temperature_energy_efficiency"
    outputs = []
    for suffix, kwargs in ((".svg", {}), (".pdf", {}), (".png", {"dpi": 600})):
        path = base.with_suffix(suffix)
        figure.savefig(path, facecolor="white", **kwargs)
        outputs.append(path)
    source_data = output_dir / "temperature_energy_efficiency_source_data.csv"
    write_publication_source_data(prepared, source_data)
    manifest = {
        "figure": "temperature_energy_efficiency",
        "version": VERSION,
        "archetype": "quantitative grid",
        "temperatures_C": temperatures,
        "raw_curve_count": len({row["file_tail_folder"] for row in prepared}),
        "point_count": len(prepared),
        "panels": {"a": "150 C discharged energy density", "b": "150 C efficiency",
                   "c": "200 C discharged energy density", "d": "200 C efficiency"},
        "statistics": "Raw curves only; no averaging, error bars, interpolation, or extrapolation.",
        "outputs": [str(path) for path in outputs],
        "source_data": str(source_data),
    }
    write_json(output_dir / "temperature_energy_efficiency_manifest.json", manifest)
    return {"output_dir": str(output_dir), "outputs": [str(path) for path in outputs],
            "source_data": str(source_data), "manifest": str(output_dir / "temperature_energy_efficiency_manifest.json"),
            "raw_curve_count": manifest["raw_curve_count"], "point_count": manifest["point_count"]}


def show_citations(window, tk, ttk) -> None:
    import webbrowser

    dialog = tk.Toplevel(window, name="citations")
    dialog.title("Cite POLYD_Ecac")
    dialog.transient(window)
    dialog.resizable(False, False)
    dialog.configure(background="#ffffff")
    content = ttk.Frame(dialog, style="Surface.TFrame", padding=(22, 18))
    content.pack(fill="both", expand=True)
    ttk.Label(content, text="Cite POLYD_Ecac", style="Section.TLabel").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 12))
    feedback = tk.StringVar(master=dialog)

    def copy_doi(doi: str) -> None:
        dialog.clipboard_clear()
        dialog.clipboard_append(doi)
        feedback.set(f"Copied {doi}")

    for row, (doi, url) in enumerate(CITATIONS, start=1):
        doi_field = ttk.Entry(content, width=34)
        doi_field.insert(0, doi)
        doi_field.configure(state="readonly")
        doi_field.grid(row=row, column=0, sticky="ew", pady=4)
        ttk.Button(content, text="Copy DOI", command=lambda value=doi: copy_doi(value)).grid(row=row, column=1, padx=(8, 0))
        ttk.Button(content, text="Open link", command=lambda link=url: webbrowser.open_new_tab(link)).grid(row=row, column=2, padx=(8, 0))
    ttk.Label(content, textvariable=feedback, style="Muted.TLabel").grid(row=3, column=0, columnspan=2, sticky="w", pady=(12, 0))
    close = ttk.Button(content, text="Close", command=dialog.destroy)
    close.grid(row=3, column=2, sticky="e", pady=(12, 0))
    dialog.bind("<Escape>", lambda _event: dialog.destroy())
    dialog.update_idletasks()
    left = window.winfo_rootx() + max(0, (window.winfo_width() - dialog.winfo_reqwidth()) // 2)
    top = window.winfo_rooty() + max(0, (window.winfo_height() - dialog.winfo_reqheight()) // 3)
    dialog.geometry(f"+{left}+{top}")
    if not dialog.winfo_viewable():
        dialog.wait_visibility()
    dialog.grab_set()
    close.focus_set()
    dialog.wait_window()


def configure_workspace_style(window, ttk) -> None:
    window.configure(background="#eef2f7")
    window.option_add("*Font", "{Segoe UI} 10")
    style = ttk.Style(window)
    style.theme_use("clam")
    style.configure(".", font=("Segoe UI", 10), foreground="#223148")
    style.configure("TFrame", background="#eef2f7")
    style.configure("Surface.TFrame", background="#ffffff")
    style.configure("TLabel", background="#ffffff", foreground="#26364c")
    style.configure("Title.TLabel", background="#ffffff", foreground="#152944", font=("Segoe UI", 20, "bold"))
    style.configure("Section.TLabel", background="#ffffff", foreground="#152944", font=("Segoe UI", 12, "bold"))
    style.configure("Muted.TLabel", background="#ffffff", foreground="#556981")
    style.configure("TButton", padding=(12, 7), background="#f6f8fb", bordercolor="#c6d2e2",
                    lightcolor="#f6f8fb", darkcolor="#f6f8fb", focuscolor="#254f9b")
    style.map("TButton", background=[("disabled", "#edf1f6"), ("pressed", "#d6e3f5"), ("active", "#e6eef9")],
              foreground=[("disabled", "#78889d")])
    style.configure("Primary.TButton", background="#2158b4", foreground="#ffffff",
                    bordercolor="#2158b4", lightcolor="#2158b4", darkcolor="#2158b4")
    style.map("Primary.TButton", background=[("disabled", "#bac9df"), ("pressed", "#153f88"), ("active", "#194b9e")],
              foreground=[("disabled", "#344e75"), ("!disabled", "#ffffff")])
    style.configure("TMenubutton", padding=(10, 7), background="#f6f8fb", bordercolor="#c6d2e2")
    style.map("TMenubutton", background=[("active", "#e6eef9")])
    style.configure("TEntry", padding=(5, 5), fieldbackground="#ffffff", bordercolor="#c6d2e2")
    style.configure("TCombobox", padding=(5, 5), fieldbackground="#ffffff", bordercolor="#c6d2e2")
    style.map("TCombobox", fieldbackground=[("readonly", "#ffffff")], selectforeground=[("readonly", "#223148")],
              selectbackground=[("readonly", "#ffffff")])
    style.configure("Horizontal.TProgressbar", thickness=7, background="#2158b4", troughcolor="#e1e9f4", borderwidth=0)


def launch_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    class PlotNavigation(NavigationToolbar2Tk):
        toolitems = tuple(item for item in NavigationToolbar2Tk.toolitems
                          if item[0] in {"Home", "Back", "Forward", "Pan", "Zoom"})

    root_dir = Path(sys.executable if getattr(sys, "frozen", False) else __file__).resolve().parent
    window = tk.Tk()
    window.title(f"{TOOL_NAME} {VERSION}")
    window.geometry("1200x760")
    window.minsize(850, 600)
    configure_workspace_style(window, ttk)

    input_var = tk.StringVar()
    output_var = tk.StringVar(value=str(Path.home() / "POLYD_Ecac" / "output" / "gui_run"))
    mode_var = tk.StringVar(value="Polyk")
    x_var = tk.StringVar(value=PLOT_LABELS["Emax"])
    x_metrics = {PLOT_LABELS.get(column, column): column for column in NUMERIC_COLUMNS}
    layout_var = tk.StringVar(value="Split panels")
    legend_depth_var = tk.IntVar(value=3)
    selected_metrics = {column: tk.BooleanVar(value=column in DEFAULT_Y_COLUMNS) for column in NUMERIC_COLUMNS}
    status_var = tk.StringVar(value="Ready")
    summary_var = tk.StringVar()
    rows: List[dict] = []
    group_values: List[str] = []
    colors: dict = {}
    worker_queue: queue.Queue = queue.Queue()
    calculation_running = False
    pending_redraw = None
    canvas = None
    toolbar = None

    header = ttk.Frame(window, style="Surface.TFrame", padding=(22, 12))
    header.pack(fill="x")
    ttk.Label(header, text=TOOL_NAME, style="Title.TLabel").pack(side="left")
    ttk.Label(header, text=VERSION, style="Muted.TLabel").pack(side="left", padx=(9, 0), pady=(6, 0))
    attribution = ttk.Frame(header, style="Surface.TFrame")
    attribution.pack(side="right")
    ttk.Label(attribution, text=CREDIT, style="Muted.TLabel").pack(side="left", padx=(0, 10))
    ttk.Button(attribution, text="Cite", command=lambda: show_citations(window, tk, ttk)).pack(side="left")

    actions = ttk.Frame(window, style="Surface.TFrame", padding=(22, 0, 22, 12))
    actions.pack(fill="x")
    run_button = ttk.Button(actions, text="Run calculation", style="Primary.TButton", command=lambda: run_gui())
    run_button.pack(side="left")
    open_button = ttk.Button(actions, text="Open CSV", command=lambda: load_csv_dialog())
    open_button.pack(side="left", padx=8)
    export_button = ttk.Menubutton(actions, text="Export")
    export_button.pack(side="left")
    export_menu = tk.Menu(export_button, tearoff=False)
    export_menu.add_command(label="Current view (PNG, SVG or PDF)…", command=lambda: export_current())
    export_menu.add_command(label="Publication set (150 / 200 °C)…", command=lambda: export_publication_set())
    export_button.configure(menu=export_menu)
    export_button.state(["disabled"])

    footer = ttk.Frame(window, style="Surface.TFrame", padding=(20, 8))
    footer.pack(side="bottom", fill="x")
    ttk.Label(footer, textvariable=status_var, style="Muted.TLabel").pack(side="left", fill="x", expand=True)
    progress = ttk.Progressbar(footer, mode="indeterminate", length=90)
    progress.pack(side="right", padx=(10, 0))
    ttk.Label(footer, textvariable=summary_var, style="Muted.TLabel").pack(side="right", padx=(10, 0))

    body = ttk.Frame(window, padding=(14, 12))
    body.pack(fill="both", expand=True)
    body.columnconfigure(1, weight=1)
    body.rowconfigure(0, weight=1)
    sidebar = ttk.Frame(body, style="Surface.TFrame", padding=(13, 12))
    sidebar.grid(row=0, column=0, sticky="nsew", padx=(0, 12))
    sidebar.columnconfigure(0, weight=1)
    sidebar.rowconfigure(8, weight=1)
    chart = ttk.Frame(body, style="Surface.TFrame", padding=(12, 10))
    chart.grid(row=0, column=1, sticky="nsew")
    chart.columnconfigure(0, weight=1)
    chart.rowconfigure(2, weight=1)

    ttk.Label(sidebar, text="Data and groups", style="Section.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 8))

    def browse(variable: tk.StringVar, prompt: str) -> None:
        candidate = Path(variable.get().strip()) if variable.get().strip() else root_dir
        initial = candidate if candidate.is_dir() else candidate.parent
        selected = filedialog.askdirectory(parent=window, title=prompt, initialdir=str(initial))
        if selected:
            variable.set(selected)

    def add_path(row: int, label: str, variable: tk.StringVar) -> None:
        ttk.Label(sidebar, text=label).grid(row=row, column=0, sticky="w", pady=(4, 2))
        fields = ttk.Frame(sidebar, style="Surface.TFrame")
        fields.grid(row=row + 1, column=0, sticky="ew", pady=(0, 6))
        ttk.Entry(fields, textvariable=variable, width=19).pack(side="left", fill="x", expand=True)
        ttk.Button(fields, text="Browse", command=lambda: browse(variable, f"Select {label.lower()}")).pack(side="left", padx=(5, 0))

    add_path(1, "Input folder", input_var)
    add_path(3, "Output folder", output_var)
    ttk.Label(sidebar, text="Data format").grid(row=5, column=0, sticky="w", pady=(4, 2))
    mode = ttk.Combobox(sidebar, textvariable=mode_var, values=("Polyk", "Radiant"), state="readonly", width=16)
    mode.grid(row=6, column=0, sticky="ew", pady=(0, 11))

    groups_header = ttk.Frame(sidebar, style="Surface.TFrame")
    groups_header.grid(row=7, column=0, sticky="ew", pady=(0, 5))
    ttk.Label(groups_header, text="Sample groups", style="Section.TLabel").pack(side="left", padx=(0, 4))
    ttk.Button(groups_header, text="None", command=lambda: select_groups(False)).pack(side="right")
    ttk.Button(groups_header, text="All", command=lambda: select_groups(True)).pack(side="right", padx=(0, 4))

    group_panel = ttk.Frame(sidebar, style="Surface.TFrame")
    group_panel.grid(row=8, column=0, sticky="nsew")
    group_panel.columnconfigure(0, weight=1)
    group_panel.rowconfigure(0, weight=1)
    group_list = tk.Listbox(group_panel, name="groups", height=5, width=23, selectmode="extended", exportselection=False,
                           background="#f8fafd", foreground="#253b57", selectbackground="#dce8fb", selectforeground="#142d53",
                           highlightthickness=1, highlightbackground="#c6d2e2", highlightcolor="#2158b4",
                           borderwidth=0, relief="flat", font=("Segoe UI", 10))
    group_list.grid(row=0, column=0, sticky="nsew")
    scroll_y = ttk.Scrollbar(group_panel, orient="vertical", command=group_list.yview)
    scroll_x = ttk.Scrollbar(group_panel, orient="horizontal", command=group_list.xview)
    scroll_y.grid(row=0, column=1, sticky="ns")
    scroll_x.grid(row=1, column=0, sticky="ew")
    group_list.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)

    ttk.Label(chart, text="Storage performance", style="Section.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 8))
    options = ttk.Frame(chart, style="Surface.TFrame")
    options.grid(row=1, column=0, sticky="ew", pady=(0, 6))
    options.columnconfigure(1, weight=1)
    ttk.Label(options, text="X axis").grid(row=0, column=0, padx=(0, 5), sticky="w")
    x_combo = ttk.Combobox(options, name="x_axis", textvariable=x_var, values=tuple(x_metrics), state="readonly", width=15)
    x_combo.grid(row=0, column=1, sticky="ew", padx=(0, 8))
    layout_combo = ttk.Combobox(options, name="plot_layout", textvariable=layout_var,
                                values=("Split panels", "Dual Y axis"), state="readonly", width=11)
    layout_combo.grid(row=0, column=2, sticky="w", padx=(0, 8))
    metrics_button = ttk.Menubutton(options, text="Y metrics")
    metrics_button.grid(row=0, column=3, sticky="e")
    metrics_menu = tk.Menu(metrics_button, tearoff=False)
    for column in NUMERIC_COLUMNS:
        metrics_menu.add_checkbutton(label=PLOT_LABELS.get(column, column), variable=selected_metrics[column], command=lambda: request_redraw())
    depth_menu = tk.Menu(metrics_menu, tearoff=False)
    for depth in range(1, 7):
        depth_menu.add_radiobutton(label=str(depth), variable=legend_depth_var, value=depth, command=lambda: request_redraw())
    metrics_menu.add_separator()
    metrics_menu.add_cascade(label="Legend path depth", menu=depth_menu)
    metrics_button.configure(menu=metrics_menu)
    plot_frame = ttk.Frame(chart, name="plot", style="Surface.TFrame")
    plot_frame.grid(row=2, column=0, sticky="nsew")
    navigation_frame = ttk.Frame(chart, style="Surface.TFrame")
    navigation_frame.grid(row=3, column=0, sticky="ew")

    def plot_state() -> Tuple[List[dict], List[str], str, str]:
        groups = [group_values[index] for index in group_list.curselection()]
        y_columns = [column for column in NUMERIC_COLUMNS if selected_metrics[column].get()]
        x_column = x_metrics[x_var.get()]
        layout = "dual" if layout_var.get() == "Dual Y axis" else "split"
        series = build_plot_series(rows, groups, x_column, y_columns, colors, legend_depth_var.get())
        if layout == "dual" and len({item["metric"] for item in series}) != 2:
            layout = "split"
        return series, y_columns, x_column, layout

    def select_groups(select: bool) -> None:
        if select:
            group_list.selection_set(0, "end")
        else:
            group_list.selection_clear(0, "end")
        request_redraw()

    def request_redraw(*_args) -> None:
        nonlocal pending_redraw
        if pending_redraw is None:
            pending_redraw = window.after_idle(redraw)

    def redraw() -> None:
        nonlocal canvas, toolbar, pending_redraw
        pending_redraw = None
        if toolbar is not None:
            toolbar.destroy()
            toolbar = None
        if canvas is not None:
            canvas.get_tk_widget().destroy()
            canvas = None
        for child in plot_frame.winfo_children():
            child.destroy()
        series, y_columns, x_column, layout = plot_state()
        export_menu.entryconfigure(0, state="normal" if series else "disabled")
        if not series:
            message = "Run a calculation or open a result CSV" if not rows else "Select groups and Y metrics with valid data"
            ttk.Label(plot_frame, text=message, style="Muted.TLabel").place(relx=0.5, rely=0.45, anchor="center")
            return
        if layout_var.get() == "Dual Y axis" and layout != "dual":
            layout_var.set("Split panels")
        figure = make_figure(series, x_column, y_columns, layout, screen=True)
        canvas = FigureCanvasTkAgg(figure, master=plot_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = PlotNavigation(canvas, navigation_frame, pack_toolbar=False)
        toolbar.configure(background="#ffffff")
        for child in toolbar.winfo_children():
            if isinstance(child, (tk.Button, tk.Checkbutton, tk.Label, tk.Frame)):
                child.configure(background="#ffffff")
        toolbar.update()
        toolbar.pack(side="left", fill="x", expand=True)
        canvas.draw()

    def load_rows(csv_path: Path) -> None:
        nonlocal rows, colors, group_values
        with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = set(reader.fieldnames or [])
            if "file_tail_folder" not in fields or not fields.intersection(NUMERIC_COLUMNS):
                raise ValueError("Choose a POLYD_Ecac result CSV with groups and numeric columns.")
            rows = list(reader)
        group_values = sorted({row["file_tail_folder"] for row in rows if row.get("file_tail_folder")})
        colors = group_colors(group_values)
        markers = group_markers(colors)
        glyphs = dict(zip(SAMPLE_MARKERS, "●■▲◆▼✚✕⬡★"))
        group_list.delete(0, "end")
        for index, group in enumerate(group_values):
            name = "Root folder" if group == "." else group
            group_list.insert("end", f"{glyphs[markers[group]]}  {name}")
            group_list.itemconfigure(index, foreground=colors[group], selectforeground=colors[group])
        group_list.selection_set(0, "end")
        export_button.state(["!disabled"] if rows and not calculation_running else ["disabled"])
        status_var.set(f"Loaded {csv_path.name}")
        summary_var.set(f"{sum(not row.get('error') for row in rows)} valid  ·  {sum(bool(row.get('error')) for row in rows)} flagged")
        request_redraw()

    def load_csv_dialog() -> None:
        selected = filedialog.askopenfilename(parent=window, filetypes=(("Result CSV", "*.csv"),))
        if selected:
            try:
                load_rows(Path(selected))
            except Exception as exc:
                messagebox.showerror(TOOL_NAME, str(exc), parent=window)

    def run_gui() -> None:
        nonlocal calculation_running
        if calculation_running:
            return
        input_path, output_path = input_var.get().strip(), output_var.get().strip()
        if not input_path or not output_path:
            messagebox.showerror(TOOL_NAME, "Select an input folder and an output folder.", parent=window)
            return
        config = RunConfig(Path(input_path).resolve(), Path(output_path).resolve(), mode_var.get().lower())
        calculation_running = True
        for control in (run_button, open_button, export_button):
            control.state(["disabled"])
        summary_var.set("")
        status_var.set("Calculating…")
        progress.start(14)

        def worker() -> None:
            run_logs: List[str] = []

            def log(message: str) -> None:
                run_logs.append(message)
                worker_queue.put(("progress", message))

            try:
                result = run_calculation(config, log)
                result_folder = Path(result["output_folder"])
                write_json(result_folder / "result.json", result)
                (result_folder / "run.log").write_text("\n".join(run_logs) + "\n", encoding="utf-8")
                worker_queue.put(("result", result))
            except Exception as exc:
                worker_queue.put(("error", str(exc)))

        def poll_worker() -> None:
            nonlocal calculation_running
            try:
                kind, payload = worker_queue.get_nowait()
            except queue.Empty:
                window.after(100, poll_worker)
                return
            if kind == "progress":
                status_var.set(payload)
                window.after(1, poll_worker)
                return
            calculation_running = False
            progress.stop()
            for control in (run_button, open_button):
                control.state(["!disabled"])
            export_button.state(["!disabled"] if rows else ["disabled"])
            if kind == "error":
                status_var.set("Calculation failed")
                messagebox.showerror(TOOL_NAME, payload, parent=window)
                return
            try:
                load_rows(Path(payload["output_csv"]))
                status_var.set("Calculation complete" if payload["ok"] else "No valid curves")
                summary_var.set(f"{payload['valid_count']} valid  ·  {payload['skipped_count']} skipped  ·  {payload['error_count']} flagged")
            except Exception as exc:
                status_var.set("Result could not be loaded")
                messagebox.showerror(TOOL_NAME, str(exc), parent=window)

        threading.Thread(target=worker, daemon=True).start()
        window.after(100, poll_worker)

    def export_current() -> None:
        series, y_columns, x_column, layout = plot_state()
        if not series:
            return
        folder = Path(output_var.get()).expanduser()
        target = filedialog.asksaveasfilename(parent=window, initialdir=str(folder if folder.is_dir() else root_dir),
                                             initialfile="polyd_ecac_plot.png", defaultextension=".png",
                                             filetypes=(("PNG image", "*.png"), ("SVG vector", "*.svg"), ("PDF vector", "*.pdf")))
        if target:
            try:
                outputs = export_figures(Path(target), series, x_column, y_columns, layout)
                status_var.set(f"Saved {len(outputs)} figure{'s' if len(outputs) != 1 else ''} to {outputs[0].parent.name}")
            except Exception as exc:
                messagebox.showerror(TOOL_NAME, str(exc), parent=window)

    def export_publication_set() -> None:
        folder = Path(output_var.get()).expanduser()
        selected = filedialog.askdirectory(parent=window, title="Select publication output folder",
                                           initialdir=str(folder if folder.is_dir() else root_dir))
        if selected:
            try:
                result = export_publication_figures(rows, allocate_output_folder(Path(selected) / "publication_figure"))
                status_var.set(f"Saved publication set to {Path(result['output_dir']).name}")
            except Exception as exc:
                messagebox.showerror(TOOL_NAME, str(exc), parent=window)

    group_list.bind("<<ListboxSelect>>", request_redraw)
    x_combo.bind("<<ComboboxSelected>>", request_redraw)
    layout_combo.bind("<<ComboboxSelected>>", request_redraw)
    window.bind("<Control-o>", lambda _event: open_button.invoke())
    window.bind("<Control-r>", lambda _event: run_button.invoke())
    request_redraw()
    window.mainloop()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv: return launch_gui()
    if argv[0] in {"--help", "-h", "help"}: show_help(); return 0
    if argv[0] == "--health-check": return command_health_check()
    if argv[0] == "--licenses":
        # Preserve the notice file's UTF-8 encoding when redirecting stdout.
        sys.stdout.buffer.write(Path(__file__).with_name("THIRD_PARTY_NOTICES.txt").read_bytes())
        return 0
    parser = argparse.ArgumentParser(add_help=False); commands = parser.add_subparsers(dest="command")
    run_parser = commands.add_parser("run"); run_parser.add_argument("--input", required=True); run_parser.add_argument("--output", required=True)
    figure_parser = commands.add_parser("figure"); figure_parser.add_argument("--input", required=True); figure_parser.add_argument("--output", required=True)
    parsed = parser.parse_args(argv)
    if parsed.command == "run": return command_run(parsed)
    if parsed.command == "figure": return command_figure(parsed)
    show_help(); return 1


if __name__ == "__main__":
    raise SystemExit(main())
