"""Numerical reports and publication-oriented SVG figures for cpptraj outputs.

Replicas remain the sampling unit: frame statistics are descriptive and replica
summaries weight each available replica equally. No frame-based confidence
intervals or pooled autocorrelation estimates are calculated.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


from .figure_names import metric_figure_stem
from .figure_style import PALETTES
from .labels import analysis_label, readable, observable_label, series_title, unit_label, coordinate_label
from .distributions import shared_edges, histogram, distribution_rows
PALETTE = PALETTES["lagoon"]
REPLICA_COLORS = {}


def _replica_color(name):
    return REPLICA_COLORS.get(name, PALETTE[0])
TIDY_COLUMNS = ("replica", "analysis", "section", "series", "kind", "coordinate", "coordinate_y", "value", "unit", "coordinate_unit", "source")


class TableFormatError(ValueError):
    """Input is not an unambiguous, rectangular numeric table."""


@dataclass
class NumericTable:
    data: np.ndarray
    columns: list[str]
    comments: list[str]
    row_labels: list[str] | None = None
    column_units: dict[str, str] | None = None


@dataclass
class Series:
    replica: str
    analysis: str
    section: str
    name: str
    kind: str
    x: np.ndarray
    values: np.ndarray
    unit: str
    xunit: str
    source: str
    discrete: bool = False
    circular: bool = False
    shared_basis: bool = False

    @property
    def key(self) -> tuple[str, ...]:
        return self.analysis, self.section, self.name, self.kind, self.unit, self.xunit, "circular" if self.circular else "categorical" if self.discrete else "linear"


@dataclass
class Matrix:
    replica: str
    analysis: str
    section: str
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    unit: str
    source: str
    matrix_type: str = ""
    xunit: str = "index"
    yunit: str = "index"
    discrete: bool = False
    betadetail: bool = False

    @property
    def key(self) -> tuple[str, ...]:
        return self.analysis, self.section, self.unit, self.matrix_type, self.xunit, self.yunit


def parse_numeric_table(path: str | Path) -> NumericTable:
    """Read whitespace/comma separated numeric cpptraj ASCII; reject ragged data.

    Comment headers are used only when their width agrees with the data. Fortran
    D exponents and explicit NaN values are supported. Non-numeric non-comment
    lines are rejected rather than silently discarding a possibly relevant row.
    """
    rows: list[list[float]] = []
    comments: list[str] = []
    candidates: list[list[str]] = []
    width = None
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(("#", "@", ";")):
                comments.append(line)
                if line.startswith("#"):
                    candidates.append(re.split(r"[\s,]+", line.lstrip("#").strip()))
                continue
            tokens = re.split(r"[\s,]+", line.split("#", 1)[0].strip())
            try:
                values = [float(token.replace("D", "E").replace("d", "e")) for token in tokens]
            except ValueError as exc:
                raise TableFormatError(f"{Path(path).name}:{line_number}: nonnumeric table row") from exc
            if width is None:
                width = len(values)
            if len(values) != width:
                raise TableFormatError(f"{Path(path).name}:{line_number}: ragged table ({len(values)} versus {width} columns)")
            rows.append(values)
    if not rows:
        raise TableFormatError(f"{Path(path).name}: no numeric rows")
    columns = next((candidate for candidate in reversed(candidates) if len(candidate) == width), [])
    if not columns:
        columns = [f"column_{index + 1}" for index in range(width or 0)]
    # A duplicate dataset header must not silently merge two columns.
    counts: dict[str, int] = defaultdict(int)
    unique = []
    for name in columns:
        counts[name] += 1
        unique.append(name if counts[name] == 1 else f"{name}_{counts[name]}")
    return NumericTable(np.asarray(rows, dtype=float), unique, comments)


def _parse_contact_table(path: Path) -> NumericTable:
    """Keep atom/contact identities alongside numeric H-bond/contact columns."""
    comments, header_candidates, rows = [], [], []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                header_candidates.append(line.lstrip("#").split())
                continue
            tokens = line.split()
            if rows and len(tokens) != len(rows[0]):
                raise TableFormatError(f"contact table line {line_number}: inconsistent number of fields")
            rows.append(tokens)
    if not rows:
        # Zero contacts is a valid scientific result; no fake zero occupancy row.
        return NumericTable(np.empty((0, 0)), [], comments, [])
    header = next((candidate for candidate in reversed(header_candidates) if len(candidate) == len(rows[0])), None)
    if header is None:
        raise TableFormatError("contact table header does not match the data")
    numeric, identity = [], []
    for index, label in enumerate(header):
        if label.lower() in {"acceptor", "donorh", "donor", "atom1", "atom2", "contact", "contact1", "contact2"}:
            identity.append(index)
            continue
        try:
            [float(row[index]) for row in rows]
            numeric.append(index)
        except ValueError:
            identity.append(index)
    if not numeric:
        raise TableFormatError("contact table contains no numeric observables")
    columns = [header[index] for index in numeric]
    labels = [";".join(f"{header[index]}={row[index]}" for index in identity) or f"contact_row={row_index}" for row_index, row in enumerate(rows, 1)]
    unit_lookup = {"frames": "count", "nframes": "count", "frac": "fraction", "fraction": "fraction", "avgdist": "angstrom", "avgdist.": "angstrom", "distance": "angstrom", "avgang": "degree", "avgang.": "degree", "angle": "degree", "stdev": "angstrom"}
    units = {label: unit_lookup.get(label.lower(), "") for label in columns}
    return NumericTable(np.asarray([[float(row[index]) for index in numeric] for row in rows]), columns, comments, labels, units)


def _safe_name(value: str) -> str:
    readable = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.") or "data"
    return readable[:100] + "_" + hashlib.sha1(value.encode()).hexdigest()[:7]


def _number(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        return f"{value:.10g}" if math.isfinite(value) else "nan"
    return str(value)


def _write_table(path: Path, columns: Iterable[str], rows: Iterable[Iterable[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_number(value) for value in row])


def descriptive_statistics(values: np.ndarray, block_size: int = 50) -> dict[str, Any]:
    """Describe finite observations; block SD uses complete contiguous blocks.

    A block containing a nonfinite sample is omitted in full. Frames are never
    compacted across missing samples to manufacture contiguous blocks.
    """
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    block_size = max(1, int(block_size))
    blocks = []
    for offset in range(0, len(array) - block_size + 1, block_size):
        block = array[offset:offset + block_size]
        if np.all(np.isfinite(block)):
            blocks.append(float(np.mean(block)))
    result: dict[str, Any] = {"n_total": len(array), "n_finite": len(finite), "mean": None, "sd": None, "median": None, "q05": None, "q95": None, "block_size": block_size, "n_complete_blocks": len(blocks), "block_mean_sd": None}
    if finite.size:
        result.update(mean=float(np.mean(finite)), median=float(np.median(finite)), q05=float(np.quantile(finite, 0.05)), q95=float(np.quantile(finite, 0.95)))
    if finite.size > 1:
        result["sd"] = float(np.std(finite, ddof=1))
    if len(blocks) > 1:
        result["block_mean_sd"] = float(np.std(blocks, ddof=1))
    return result


def circular_statistics(values: np.ndarray) -> dict[str, Any]:
    """Circular descriptions in degrees; mean is undefined at zero resultant."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    result = {"n_finite": len(finite), "circular_mean_degree": None, "resultant_length": None, "circular_sd_degree": None}
    if not len(finite):
        return result
    vector = np.mean(np.exp(1j * np.deg2rad(finite)))
    length = min(1.0, float(abs(vector)))
    result["resultant_length"] = length
    if length > 1e-12:
        result["circular_mean_degree"] = float(np.rad2deg(np.angle(vector)))
        result["circular_sd_degree"] = float(np.rad2deg(np.sqrt(-2 * np.log(length))))
    return result


def _nastruct_series(path: Path, replica: str, artifact: dict, config: dict) -> list[Series]:
    """Parse explicit current CPPTRAJ NAStruct headers, preserving pair identity.

    Format checked against Amber-MD/cpptraj src/Action_NAstruct.cpp Print().
    Unrecognized/malformed tables are rejected, not interpreted by position.
    """
    header = None
    groups: dict[tuple[str, str], list[list[float]]] = defaultdict(list)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                candidate = line.lstrip("#").split()
                if candidate[:3] in (["Frame", "Base1", "Base2"], ["Frame", "BP1", "BP2"]):
                    header = candidate
                continue
            if header is None:
                raise TableFormatError("NAStruct requires a recognized #Frame Base1 Base2 or #Frame BP1 BP2 header")
            if header[1] == "BP1":
                match = re.fullmatch(r"(\d+)\s+(\d+\s*-\s*\d+)\s+(\d+\s*-\s*\d+)\s+(.+)", line)
                if not match:
                    raise TableFormatError(f"NAStruct line {line_number}: malformed base-pair identifiers")
                tokens = [match.group(1), re.sub(r"\s+", "", match.group(2)), re.sub(r"\s+", "", match.group(3)), *match.group(4).split()]
            else:
                tokens = line.split()
            if len(tokens) != len(header):
                raise TableFormatError(f"NAStruct line {line_number}: header/data width mismatch")
            values = [float(tokens[0]), *[float("nan" if value == "----" else value) for value in tokens[3:]]]
            if header[1] == "Base1" and not all(re.fullmatch(r"\d+", token) for token in tokens[1:3]):
                raise TableFormatError("NAStruct base identities must be positive residue numbers")
            groups[(tokens[1], tokens[2])].append(values)
    if not groups or header is None:
        raise TableFormatError("NAStruct has no parameter rows (base pairs may be absent)")
    linear = {"Shear", "Stretch", "Stagger", "Major", "Minor", "DX", "DY", "DZ", "Shift", "Slide", "Rise", "Zp", "X-disp", "Y-disp"}
    angles = {"Buckle", "Propeller", "Opening", "RX", "RY", "RZ", "Tilt", "Roll", "Twist", "Incl.", "Tip"}
    recognized = linear | angles | {"BP", "HB"}
    if any(name not in recognized for name in header[3:]):
        raise TableFormatError("NAStruct contains an unrecognized parameter name")
    series = []
    for pair, rows in groups.items():
        data = np.asarray(sorted(rows, key=lambda row: row[0]), dtype=float)
        numeric = NumericTable(data, ["Frame", *header[3:]], [])
        x, xunit, _ = _coordinate(numeric, artifact, config, "timeseries")
        pair_label = "/".join(pair) if header[1] == "BP1" else "-".join(pair)
        for i, name in enumerate(header[3:], 1):
            unit = "angstrom" if name in linear else "degree" if name in angles else "fraction" if name == "BP" else "count"
            series.append(Series(replica, str(artifact.get("analysis", "nastruct")) + "/" + path.name.split(".", 1)[0], str(artifact.get("section", "nucleic")), f"{name}[{pair_label}]", "timeseries", x, data[:, i], unit, xunit, str(path), False, name in angles))
    return series


def concatenate_series(replica_series: list[Series]) -> tuple[np.ndarray, np.ndarray, list[float]]:
    """Join series by sample spacing for display, returning replica boundaries.

    Offsets intentionally remove any shared original start time. The resulting
    coordinate describes an appended display, not a continuous physical run.
    """
    all_x, all_y, boundaries = [], [], []
    end = None
    previous_step = 1.0
    for series in replica_series:
        if not len(series.x):
            continue
        differences = np.diff(series.x)
        positive = differences[np.isfinite(differences) & (differences > 0)]
        step = float(np.median(positive)) if len(positive) else previous_step
        offset = 0.0 if end is None else end + step
        if end is not None:
            boundaries.append((end + offset) / 2)
        shifted = series.x - series.x[0] + offset
        all_x.append(shifted)
        all_y.append(series.values)
        end = float(shifted[-1])
        previous_step = step
    if not all_x:
        return np.array([]), np.array([]), []
    return np.concatenate(all_x), np.concatenate(all_y), boundaries


def _coordinate(table: NumericTable, artifact: dict, config: dict, kind: str) -> tuple[np.ndarray, str, int | None]:
    metadata = artifact.get("metadata") or {}
    index = metadata.get("x_column")
    first = table.columns[0].lower().lstrip("#")
    recognizable = first in {"frame", "time", "res", "residue", "resnum", "resid", "atom", "index", "x", "lag", "bin"}
    if index is None and recognizable:
        index = 0
    if index is None:
        # Single values per row have an unambiguous sample/index coordinate.
        if table.data.shape[1] != 1:
            raise TableFormatError("no recognizable coordinate column; set artifact metadata.x_column")
        x = np.arange(1, len(table.data) + 1, dtype=float)
    else:
        index = int(index)
        if index < 0 or index >= table.data.shape[1]:
            raise TableFormatError("metadata.x_column is out of range")
        x = table.data[:, index].copy()
    if not np.all(np.isfinite(x)) or np.any(np.diff(x) <= 0):
        raise TableFormatError("coordinate must contain finite, strictly increasing values")
    xunit = metadata.get("x_unit")
    if kind == "timeseries" and first != "time" and not metadata.get("x_is_time"):
        frames = config.get("frames", {})
        dt = frames.get("dt_ps")
        if dt is not None:
            dt = float(dt)
            if not np.isfinite(dt) or dt <= 0:
                raise TableFormatError("frames.dt_ps must be positive and finite")
            start = int(frames.get("start", 1))
            stride = int(frames.get("stride", 1))
            x = ((start - 1) + (x - 1) * stride) * dt / 1000.0
            xunit = "time (ns)"
        else:
            xunit = "analyzed frame"
    elif not xunit:
        xunit = "residue" if first in {"res", "residue", "resnum", "resid"} else ("time (unspecified units)" if first == "time" else "index")
    return x, str(xunit), index


def _matrix_from_table(replica: str, artifact: dict, table: NumericTable, source: str) -> Matrix:
    metadata = artifact.get("metadata") or {}
    matrix_format = metadata.get("matrix_format", "auto").lower()
    data = table.data
    numeric_header_axes = None
    try:
        if len(table.columns) > 1 and not table.columns[0].startswith("column_"):
            numeric_header_axes = np.asarray([float(label) for label in table.columns[1:]])
    except ValueError:
        pass
    if matrix_format == "auto":
        first_three = [column.lower() for column in table.columns[:3]]
        if numeric_header_axes is not None:
            matrix_format = "indexed"
        elif data.shape[1] == 3 and (first_three[:2] in (["x", "y"], ["row", "column"]) or len(np.unique(data[:, 0])) * len(np.unique(data[:, 1])) == len(data)):
            matrix_format = "xyz"
        elif data.shape[0] == data.shape[1]:
            matrix_format = "dense"
        elif data.shape[1] == data.shape[0] + 1 and table.columns[0].lower() in {"frame", "index", "res", "residue", "x"}:
            matrix_format = "indexed"
        else:
            raise TableFormatError("ambiguous matrix format; set metadata.matrix_format to dense, indexed or xyz")
    if matrix_format == "xyz":
        if data.shape[1] != 3:
            raise TableFormatError("xyz matrix requires exactly three columns")
        x, y = np.unique(data[:, 0]), np.unique(data[:, 1])
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise TableFormatError("matrix coordinates must be finite")
        if len(x) * len(y) != len(data) or len(np.unique(data[:, :2], axis=0)) != len(data):
            raise TableFormatError("xyz matrix is incomplete or contains duplicate coordinates")
        values = np.full((len(y), len(x)), np.nan)
        values[np.searchsorted(y, data[:, 1]), np.searchsorted(x, data[:, 0])] = data[:, 2]
    elif matrix_format in {"dense", "indexed"}:
        values = data[:, 1:] if matrix_format == "indexed" else data
        atom_ids = metadata.get("atom_ids")
        default_x = atom_ids if atom_ids is not None and len(atom_ids) == values.shape[1] else numeric_header_axes if matrix_format == "indexed" and numeric_header_axes is not None else np.arange(1, values.shape[1] + 1)
        default_y = atom_ids if atom_ids is not None and len(atom_ids) == values.shape[0] else data[:, 0] if matrix_format == "indexed" else np.arange(1, values.shape[0] + 1)
        x = np.asarray(metadata.get("x_labels", default_x), dtype=float)
        y = np.asarray(metadata.get("y_labels", default_y), dtype=float)
        if len(x) != values.shape[1] or len(y) != values.shape[0]:
            raise TableFormatError("matrix axis labels do not match its shape")
    else:
        raise TableFormatError(f"unsupported matrix_format: {matrix_format}")
    analysis = str(artifact.get("analysis", "matrix"))
    matrix_type = str(metadata.get("matrix_type", "correlation" if analysis in {"dccm", "correlation"} or artifact.get("unit") == "correlation" else "covariance" if "covariance" in analysis or "covar" in analysis else ""))
    atom_ids = metadata.get("atom_ids")
    if atom_ids is not None and len(atom_ids) == values.shape[0] == values.shape[1]:
        if "x_labels" not in metadata and np.array_equal(x, np.arange(1, len(atom_ids) + 1)):
            x = np.asarray(atom_ids, dtype=float)
        if "y_labels" not in metadata and np.array_equal(y, np.arange(1, len(atom_ids) + 1)):
            y = np.asarray(atom_ids, dtype=float)
    atom_axis = atom_ids is not None and np.array_equal(x, atom_ids) and np.array_equal(y, atom_ids)
    default_axis = "topology atom number" if atom_axis else metadata.get("axis", "index")
    return Matrix(replica, analysis, str(artifact.get("section", "system")), x, y, values, str(artifact.get("unit", "")), source, matrix_type, str(metadata.get("x_unit", default_axis)), str(metadata.get("y_unit", default_axis)), bool(metadata.get("discrete", False)))


def _plotting(dpi: int, options=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .figure_style import configure
    global PALETTE
    options = options or {}
    PALETTE = PALETTES[options.get("palette", "lagoon")]
    plt._mdwb_style = configure(plt, options, dpi)
    plt._mdwb_figure_titles = {}
    return plt


def _save_figure(plt, figure, path: Path, files: list[str], output_dir: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    style = getattr(plt, "_mdwb_style", {})
    size = style.get("font_size_pt", 11)
    old_width, old_height = figure.get_size_inches()
    width = style.get("width_mm", 180) / 25.4
    figure.set_size_inches(width, width * old_height / old_width)
    for ax in figure.axes:
        ax.title.set_fontsize(size + 1)
        ax.xaxis.label.set_fontsize(size)
        ax.yaxis.label.set_fontsize(size)
        ax.tick_params(labelsize=size - 1)
        legend = ax.get_legend()
        if legend:
            for text in legend.get_texts():
                text.set_fontsize(size - 1)
    (output_dir / "reports").mkdir(parents=True, exist_ok=True)
    (output_dir / "reports" / "figure_style.json").write_text(json.dumps(style, indent=2) + "\n", encoding="utf-8")
    titles = [figure._suptitle.get_text()] if figure._suptitle is not None else []
    titles.extend(ax.get_title() for ax in figure.axes if ax.get_title())
    display_title = " · ".join(dict.fromkeys(title.replace("\n", " · ") for title in titles))
    plt._mdwb_figure_titles[path.relative_to(output_dir).as_posix()] = display_title
    figure.savefig(path, bbox_inches="tight", metadata={"Creator": "MD-bud", "Date": None, "Title": display_title})
    plt.close(figure)
    files.append(path.relative_to(output_dir).as_posix())


def _ylabel(series: Series) -> str:
    name = observable_label(series)
    return f"{name} ({unit_label(series.unit)})" if series.unit else name


def _running_average(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    count = np.cumsum(finite)
    return np.divide(np.cumsum(np.where(finite, values, 0)), count, out=np.full(len(values), np.nan), where=count > 0)


def _plot_values(series: Series) -> np.ndarray:
    if not series.circular:
        return series.values
    values = (series.values + 180) % 360 - 180
    values = values.copy()
    values[np.flatnonzero(np.abs(np.diff(values)) > 180) + 1] = np.nan
    return values


def _running_series_average(series: Series) -> np.ndarray:
    if not series.circular:
        return _running_average(series.values)
    valid = np.isfinite(series.values)
    vector = np.cumsum(np.where(valid, np.exp(1j * np.deg2rad(series.values)), 0))
    values = np.rad2deg(np.angle(vector))
    values[np.abs(vector) < 1e-12] = np.nan
    values[np.flatnonzero(np.abs(np.diff(values)) > 180) + 1] = np.nan
    return values


def _series_figures(plt, series_groups: dict, config: dict, output_dir: Path, figure_files: list[str]):
    alpha = float(config.get("plots", {}).get("alpha", 0.8))
    alpha = min(1.0, max(0.05, alpha))
    for key, group in series_groups.items():
        first = group[0]
        if first.kind not in {"timeseries", "profile"} or first.discrete:
            continue
        stem = metric_figure_stem(first.analysis, first.section, first.name, *key)
        title = series_title(first)
        for series in group:
            fig, ax = plt.subplots(figsize=(3.5, 2.65), layout="constrained")
            ax.plot(series.x, _plot_values(series), color=_replica_color(series.replica), alpha=alpha)
            ax.set(xlabel=coordinate_label(series.xunit), ylabel=_ylabel(series), title=series_title(series, readable(series.replica)))
            ax.grid(axis="y", color="#dddddd", linewidth=0.5, alpha=0.6)
            _save_figure(plt, fig, output_dir / "figures" / _safe_name(series.replica) / f"{stem}__{_safe_name(series.replica)}__{first.kind}.svg", figure_files, output_dir)
        if len(group) > 1:
            fig, ax = plt.subplots(figsize=(7, 3.15), layout="constrained")
            for i, series in enumerate(group):
                ax.plot(series.x, _plot_values(series), label=readable(series.replica), color=_replica_color(series.replica), alpha=alpha, linestyle=("-", "--", ":", "-.")[(i // len(PALETTE)) % 4])
            ax.set(xlabel=coordinate_label(first.xunit), ylabel=_ylabel(first), title=series_title(first, "Replica comparison"))
            ax.legend(frameon=False, ncols=min(4, len(group)), loc="best")
            ax.grid(axis="y", color="#dddddd", linewidth=0.5, alpha=0.6)
            _save_figure(plt, fig, output_dir / "figures" / "overlays" / f"{stem}__replica-comparison.svg", figure_files, output_dir)
        if first.kind == "timeseries":
            # Every line is drawn separately: no visual connection across replicas.
            if len(group) > 1:
                joined_x, _, boundaries = concatenate_series(group)
                fig, ax = plt.subplots(figsize=(7, 3.15), layout="constrained")
                offset = 0
                for i, series in enumerate(group):
                    x = joined_x[offset:offset + len(series.x)]
                    ax.plot(x, _plot_values(series), color=_replica_color(series.replica), alpha=alpha, label=readable(series.replica))
                    offset += len(series.x)
                for boundary in boundaries:
                    ax.axvline(boundary, color="#555555", ls="--", lw=0.7, alpha=0.7)
                ax.set(xlabel=f"Appended {first.xunit} (not a continuous trajectory)", ylabel=_ylabel(first), title=series_title(first, "Replicas shown in sequence"))
                ax.legend(frameon=False, ncols=min(4, len(group)))
                _save_figure(plt, fig, output_dir / "figures" / "concatenated" / f"{stem}__replicas-in-sequence.svg", figure_files, output_dir)
            _distribution_figures(plt, group, output_dir, figure_files, stem, alpha)
            edges = shared_edges(group)
            fig, axes = plt.subplots(1, 2, figsize=(7, 3.5), layout="constrained")
            for i, series in enumerate(group):
                finite = series.values[np.isfinite(series.values)]
                if series.circular:
                    finite = (finite + 180) % 360 - 180
                if len(finite):
                    # Per-replica densities prevent longer runs dominating a pool.
                    axes[0].hist(finite, bins=edges, density=True, histtype="step", color=_replica_color(series.replica), alpha=alpha, label=readable(series.replica))
                axes[1].plot(series.x, _running_series_average(series), color=_replica_color(series.replica), alpha=alpha, label=readable(series.replica))
            axes[0].set(xlabel=_ylabel(first), ylabel=_density_label(first), title="Probability density")
            axes[1].set(xlabel=coordinate_label(first.xunit), ylabel=_ylabel(first), title="Running circular mean" if first.circular else "Running mean")
            axes[1].legend(frameon=False, ncols=1)
            fig.suptitle(title)
            _save_figure(plt, fig, output_dir / "figures" / "diagnostics" / f"{stem}__distribution-and-running-mean.svg", figure_files, output_dir)


def _density_label(series):
    return f"Probability density (1 / {unit_label(series.unit)})" if series.unit else "Probability density"


def _distribution_figures(plt, group, output_dir, files, stem, alpha):
    edges = shared_edges(group)
    available = [(series, histogram(series, edges)) for series in group]
    available = [(series, data) for series, data in available if data["n_finite"]]
    selections = [([item], readable(item[0].replica), _safe_name(item[0].replica)) for item in available]
    if len(available) > 1:
        selections.append((available, "Replica comparison", "overlays"))
    for selection, context, folder in selections:
        fig, ax = plt.subplots(figsize=(7, 3.8), layout="constrained")
        for index, (series, data) in enumerate(selection):
            color = _replica_color(series.replica)
            ax.stairs(data["density"], edges, color=color, linewidth=1.7,
                      linestyle=("-", "--", ":", "-.")[(index // len(PALETTE)) % 4],
                      alpha=alpha, label=f"{readable(series.replica)} (n = {data['n_finite']:,})")
            ax.stairs(data["density"], edges, color=color, fill=True, alpha=0.10)
        ax.set(xlabel=_ylabel(group[0]), ylabel=_density_label(group[0]),
               title=series_title(group[0], f"Probability distribution · {context}"), ylim=(0, None))
        ax.legend(ncols=min(3, len(selection)), loc="best")
        ax.grid(axis="y", color="#dddddd", linewidth=0.5, alpha=0.6)
        _save_figure(plt, fig, output_dir / "figures" / "distributions" / folder / f"{stem}__{folder}__probability-distribution.svg", files, output_dir)


def _pca_figures(plt, series_list: list[Series], output_dir: Path, files: list[str]):
    """Project onto the first two reported PCs; no variance/FEL inference."""
    by_section: dict[str, dict[str, list[Series]]] = defaultdict(lambda: defaultdict(list))
    for series in series_list:
        if series.analysis == "pca" and series.kind == "timeseries":
            by_section[series.section][series.replica].append(series)
    for section, by_replica in by_section.items():
        pairs = []
        for replica, group in by_replica.items():
            def mode_number(series):
                match = re.search(r"(\d+)\]?$", series.name)
                return (int(match.group(1)) if match else 999999, series.name)
            ordered = sorted(group, key=mode_number)
            if len(ordered) < 2 or not np.array_equal(ordered[0].x, ordered[1].x):
                continue
            a, b = ordered[:2]
            valid = np.isfinite(a.values) & np.isfinite(b.values)
            if not np.any(valid):
                continue
            pairs.append((a, b, valid))
            fig, ax = plt.subplots(figsize=(3.5, 3.05), layout="constrained")
            dots = ax.scatter(a.values[valid], b.values[valid], c=a.x[valid], cmap="viridis", s=5, alpha=0.65, linewidths=0, rasterized=True)
            fig.colorbar(dots, ax=ax, label=a.xunit)
            ax.set(xlabel=_ylabel(a), ylabel=_ylabel(b), title=f"Principal component analysis\n{readable(section)} · {readable(replica)}")
            _save_figure(plt, fig, output_dir / "figures" / _safe_name(replica) / f"{_safe_name(section)}_pca_scatter.svg", files, output_dir)
        if len(pairs) > 1 and all(a.shared_basis and b.shared_basis for a, b, _ in pairs):
            fig, ax = plt.subplots(figsize=(3.5, 3.05), layout="constrained")
            for i, (a, b, valid) in enumerate(pairs):
                ax.scatter(a.values[valid], b.values[valid], color=_replica_color(a.replica), s=5, alpha=0.45, linewidths=0, label=a.replica, rasterized=True)
            ax.set(xlabel=_ylabel(pairs[0][0]), ylabel=_ylabel(pairs[0][1]), title=f"Principal component analysis\n{readable(section)} · Shared pooled basis")
            ax.legend(frameon=False, markerscale=1.7)
            _save_figure(plt, fig, output_dir / "figures" / "overlays" / f"{_safe_name(section)}_pca_scatter.svg", files, output_dir)


def _cluster_figures(plt, series_list: list[Series], output_dir: Path, files: list[str]):
    for series in series_list:
        if series.analysis != "cluster" or not series.discrete:
            continue
        finite = series.values[np.isfinite(series.values)]
        if not len(finite):
            continue
        states, counts = np.unique(finite, return_counts=True)
        fig, axes = plt.subplots(1, 2, figsize=(7, 2.7), layout="constrained")
        axes[0].step(series.x, series.values, color=PALETTE[0], where="mid", linewidth=0.8)
        axes[0].set(xlabel=series.xunit, ylabel="cluster assignment", title="Cluster assignments")
        axes[1].bar(np.arange(len(states)), counts / len(finite), color=PALETTE[0], alpha=0.8, width=0.75)
        axes[1].set(xticks=np.arange(len(states)), xticklabels=[_number(state) for state in states], xlabel="cluster", ylabel="fraction of analyzed frames", title="Within-replica populations")
        fig.suptitle(series_title(series, readable(series.replica)))
        _save_figure(plt, fig, output_dir / "figures" / _safe_name(series.replica) / f"{_safe_name(series.section)}_clusters.svg", files, output_dir)


def _contact_figures(plt, tables: list, output_dir: Path, files: list[str], figure_analyses=None):
    for replica, analysis, section, unit, source, table in tables:
        if table.row_labels is None or not len(table.data):
            continue
        fraction_index = next((i for i, column in enumerate(table.columns) if column.lower() in {"frac", "fraction"}), None)
        if fraction_index is None:
            continue
        values = table.data[:, fraction_index]
        valid = np.flatnonzero(np.isfinite(values))
        selected = valid[np.argsort(values[valid])[-20:]]
        if not len(selected):
            continue
        fig, ax = plt.subplots(figsize=(7, max(2.8, len(selected) * 0.2 + 1.1)), layout="constrained")
        ax.barh(np.arange(len(selected)), values[selected], color=PALETTE[0], alpha=0.85, height=0.75)
        labels = [table.row_labels[i].replace(";", "  ") for i in selected]
        ax.set(yticks=np.arange(len(selected)), yticklabels=labels, xlabel="fraction of analyzed frames", title=f"{analysis_label(analysis)}\n{readable(section)} · {readable(replica)}\nMost frequent contacts (up to 20)")
        ax.set_xlim(0, max(1, float(np.max(values[selected])) * 1.02))
        ax.grid(axis="x", color="#dddddd", linewidth=0.5, alpha=0.6)
        _save_figure(plt, fig, output_dir / "figures" / _safe_name(replica) / f"{_safe_name(analysis + '__' + section)}_contact_fractions.svg", files, output_dir)
        if figure_analyses is not None:
            figure_analyses[files[-1]] = analysis


def _matrix_figure(plt, matrix: Matrix, path: Path, files: list[str], output_dir: Path, title: str, limits: tuple[float | None, float | None], spread: bool = False):
    from matplotlib.colors import BoundaryNorm
    from matplotlib.colors import ListedColormap
    fig, ax = plt.subplots(figsize=(3.5, 3.15), layout="constrained")
    if matrix.discrete:
        cmap = ListedColormap(("#E6E6E6", "#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#333333"))
        norm = BoundaryNorm(np.arange(-0.5, 8.5), cmap.N)
        image = ax.imshow(matrix.values, origin="lower", aspect="auto", interpolation="nearest", cmap=cmap, norm=norm, rasterized=True)
        bar = fig.colorbar(image, ax=ax, ticks=np.arange(8))
        # Modern CPPTRAJ default uses extended/bridge (betadetail restores old labels).
        beta_labels = ("Parallel β", "Antiparallel β") if matrix.betadetail else ("Extended β", "β bridge")
        bar.ax.set_yticklabels(("None", *beta_labels, "3₁₀ helix", "α helix", "π helix", "Turn", "Bend"))
        bar.set_label("DSSP assignment")
    else:
        divergent = not spread and matrix.matrix_type in {"correlation", "covariance"}
        image = ax.imshow(matrix.values, origin="lower", aspect="equal" if matrix.values.shape[0] == matrix.values.shape[1] else "auto", interpolation="nearest", cmap="RdBu_r" if divergent else "viridis", vmin=limits[0], vmax=limits[1], rasterized=True)
        fig.colorbar(image, ax=ax, label=matrix.unit or ("correlation" if matrix.matrix_type == "correlation" else "value"))
    # Cells represent selected atoms/residues, not physical distances between
    # numeric identifiers. Equal-width cells avoid stretching gaps in masks.
    for axis, coordinates in ((ax.xaxis, matrix.x), (ax.yaxis, matrix.y)):
        positions = np.unique(np.linspace(0, len(coordinates) - 1, min(6, len(coordinates))).round().astype(int))
        axis.set_ticks(positions)
        axis.set_ticklabels([_number(coordinates[index]) for index in positions])
    ax.set(xlabel=matrix.xunit, ylabel=matrix.yunit, title=title)
    _save_figure(plt, fig, path, files, output_dir)


def _matrix_reports(plt, groups: dict, output_dir: Path, files: list[str], warnings: list[str], plot: bool):
    for key, matrices in groups.items():
        first = matrices[0]
        stem = _safe_name("__".join(key))
        finite_arrays = [matrix.values[np.isfinite(matrix.values)] for matrix in matrices]
        finite_arrays = [values for values in finite_arrays if len(values)]
        if first.matrix_type == "correlation":
            limits = (-1.0, 1.0)
        elif first.matrix_type == "covariance" and finite_arrays:
            bound = max(float(np.max(np.abs(values))) for values in finite_arrays)
            limits = (-bound if bound else -1, bound if bound else 1)
        elif finite_arrays:
            limits = (min(float(np.min(values)) for values in finite_arrays), max(float(np.max(values)) for values in finite_arrays))
        else:
            limits = (None, None)
        for matrix in matrices:
            if plot:
                _matrix_figure(plt, matrix, output_dir / "figures" / _safe_name(matrix.replica) / f"{stem}_matrix.svg", files, output_dir, f"{analysis_label(matrix.analysis)}\n{readable(matrix.section)} · {readable(matrix.replica)}", limits)
        if len(matrices) < 2 or first.discrete:
            continue
        if not all(matrix.values.shape == first.values.shape and np.array_equal(matrix.x, first.x) and np.array_equal(matrix.y, first.y) for matrix in matrices):
            warnings.append(f"{first.analysis}/{first.section}: matrix mean/SD omitted because replica axes or shapes differ")
            continue
        stack = np.stack([matrix.values for matrix in matrices])
        # Require each replica at each cell; no varying or implicit weights.
        valid = np.all(np.isfinite(stack), axis=0)
        mean = np.full(first.values.shape, np.nan)
        sd = np.full(first.values.shape, np.nan)
        mean[valid] = np.mean(stack[:, valid], axis=0)
        sd[valid] = np.std(stack[:, valid], axis=0, ddof=1)
        rows = ((x, y, mean[j, i], sd[j, i], len(matrices) if valid[j, i] else 0) for j, y in enumerate(first.y) for i, x in enumerate(first.x))
        _write_table(output_dir / "reports" / "matrices" / f"{stem}_replicate_summary.dat", ("x", "y", "equal_replica_mean", "between_replica_sd", "n_replicas"), rows)
        for name, values in (("mean", mean), ("between_replica_sd", sd)):
            if plot:
                aggregate = Matrix("replica summary", first.analysis, first.section, first.x, first.y, values, first.unit, "", first.matrix_type, first.xunit, first.yunit)
                _matrix_figure(plt, aggregate, output_dir / "figures" / "matrices" / f"{stem}_{name}.svg", files, output_dir, f"{analysis_label(first.analysis)}\n{readable(first.section)}\n{readable(name)} · {len(matrices)} replicas", (0, None) if name != "mean" else limits, spread=name != "mean")


def _series_tables(groups: dict, output_dir: Path, block_size: int, warnings: list[str]):
    stats_rows, replica_rows, circular_rows, circular_replica_rows, categorical_rows = [], [], [], [], []
    statistic_names = list(descriptive_statistics(np.array([]), block_size))
    for key, group in groups.items():
        first = group[0]
        stem = _safe_name("__".join(key))
        means = []
        circular_means = []
        for series in group:
            if series.discrete:
                finite = series.values[np.isfinite(series.values)]
                states, counts = np.unique(finite, return_counts=True)
                categorical_rows.extend((series.replica, series.analysis, series.section, series.name, state, count, len(finite), count / len(finite)) for state, count in zip(states, counts))
                continue
            if series.circular:
                statistics = circular_statistics(series.values)
                circular_rows.append((series.replica, series.analysis, series.section, series.name, *statistics.values()))
                if statistics["circular_mean_degree"] is not None:
                    circular_means.append(statistics["circular_mean_degree"])
                continue
            if series.kind != "timeseries":
                continue
            statistics = descriptive_statistics(series.values, block_size)
            stats_rows.append((series.replica, series.analysis, series.section, series.name, series.unit, *[statistics[name] for name in statistic_names]))
            if statistics["mean"] is not None:
                means.append(statistics["mean"])
        if means:
            replica_rows.append((first.analysis, first.section, first.name, first.unit, len(means), np.mean(means), np.std(means, ddof=1) if len(means) > 1 else None))
        if circular_means:
            statistics = circular_statistics(np.asarray(circular_means))
            if len(circular_means) < 2:
                statistics["circular_sd_degree"] = None
            circular_replica_rows.append((first.analysis, first.section, first.name, *statistics.values()))
        # Align on the actual coordinate; absent samples remain blank.
        coordinates = np.unique(np.concatenate([series.x for series in group]))
        if len(coordinates) * len(group) <= 2_000_000:
            maps = [dict(zip(series.x, series.values)) for series in group]
            _write_table(output_dir / "reports" / "combined" / f"{stem}.dat", (first.xunit, *[series.replica for series in group]), ((x, *[mapping.get(x) for mapping in maps]) for x in coordinates))
        else:
            warnings.append(f"{first.analysis}/{first.section}/{first.name}: wide table exceeds 2 million cells; use all_replicas.dat")
        if first.kind == "profile" and len(group) > 1:
            if all(np.array_equal(series.x, first.x) for series in group):
                stack = np.stack([series.values for series in group])
                valid = np.all(np.isfinite(stack), axis=0)
                rows = ((x, np.mean(stack[:, i]) if valid[i] else None, np.std(stack[:, i], ddof=1) if valid[i] else None, len(group) if valid[i] else 0) for i, x in enumerate(first.x))
                _write_table(output_dir / "reports" / "combined" / f"{stem}_replicate_summary.dat", (first.xunit, "equal_replica_mean", "between_replica_sd", "n_replicas"), rows)
            else:
                warnings.append(f"{first.analysis}/{first.section}: profile mean/SD omitted because residue/index axes differ")
    _write_table(output_dir / "reports" / "descriptive_statistics.dat", ("replica", "analysis", "section", "series", "unit", *statistic_names), stats_rows)
    _write_table(output_dir / "reports" / "replicate_statistics.dat", ("analysis", "section", "series", "unit", "n_replicas", "equal_replica_mean", "between_replica_mean_sd"), replica_rows)
    _write_table(output_dir / "reports" / "circular_statistics.dat", ("replica", "analysis", "section", "series", "n_finite", "circular_mean_degree", "resultant_length", "circular_sd_degree"), circular_rows)
    _write_table(output_dir / "reports" / "circular_replicate_statistics.dat", ("analysis", "section", "series", "n_replicas", "equal_replica_circular_mean_degree", "between_replica_resultant_length", "between_replica_circular_sd_degree"), circular_replica_rows)
    _write_table(output_dir / "reports" / "categorical_fractions.dat", ("replica", "analysis", "section", "series", "state", "count", "n_finite", "fraction"), categorical_rows)
    _write_table(output_dir / "reports" / "probability_distributions.dat",
                 ("replica", "analysis", "section", "series", "unit", "bin_left", "bin_right", "count",
                  "n_finite", "probability", "probability_density", "geometry"), distribution_rows(groups))


def _scalar_correlations(series_list: list[Series], output_dir: Path, warnings: list[str]):
    by_replica: dict[str, list[Series]] = defaultdict(list)
    for series in series_list:
        if series.kind == "timeseries" and not series.discrete and not series.circular:
            by_replica[series.replica].append(series)
    rows = []
    for replica, group in by_replica.items():
        if len(group) > 50:
            warnings.append(f"{replica}: scalar metric correlations omitted for >50 series")
            continue
        for i, a in enumerate(group):
            for b in group[i + 1:]:
                if a.xunit != b.xunit or not np.array_equal(a.x, b.x):
                    continue
                valid = np.isfinite(a.values) & np.isfinite(b.values)
                if np.count_nonzero(valid) < 3 or np.std(a.values[valid]) == 0 or np.std(b.values[valid]) == 0:
                    continue
                r = float(np.corrcoef(a.values[valid], b.values[valid])[0, 1])
                rows.append((replica, f"{a.analysis}/{a.section}/{a.name}", f"{b.analysis}/{b.section}/{b.name}", np.count_nonzero(valid), r))
    _write_table(output_dir / "reports" / "within_replica_metric_correlations.dat", ("replica", "metric_a", "metric_b", "paired_frames", "descriptive_pearson_r"), rows)


def _diagnostic_figures(plt, output, files, figure_analyses):
    directory = output / "reports" / "diagnostics"
    groups = defaultdict(list)
    with (directory / "acf_index.dat").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[(row["analysis"], row["section"], row["series"], row["lag_unit"])].append(row)
    for key, group in groups.items():
        fig, ax = plt.subplots(figsize=(7.2, 4.5))
        for i, row in enumerate(group):
            data = np.loadtxt(output / row["path"], skiprows=1, ndmin=2)
            ax.plot(data[:, 1], data[:, 2], label=row["replica"], color=_replica_color(row["replica"]), alpha=0.8)
        ax.axhline(0, color="#666666", linewidth=0.7)
        ax.set(xlabel=f"Separation between observations / lag ({key[3]})", ylabel="Autocorrelation", title=f"{analysis_label(key[0])} · Autocorrelation\n{readable(key[1])} · {readable(key[2])}")
        ax.legend(fontsize=8)
        token = metric_figure_stem(*key)
        _save_figure(plt, fig, output / "figures" / "diagnostics" / f"{token}__autocorrelation.svg", files, output)
        figure_analyses[files[-1]] = key[0]
    groups.clear()
    with (directory / "replica_distribution_distance.dat").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[(row["analysis"], row["section"], row["series"])].append(row)
    for key, rows in groups.items():
        names = list(dict.fromkeys(name for row in rows for name in (row["replica_a"], row["replica_b"])))
        matrix = np.zeros((len(names), len(names)))
        lookup = {name: i for i, name in enumerate(names)}
        for row in rows:
            i, j = lookup[row["replica_a"]], lookup[row["replica_b"]]
            matrix[i, j] = matrix[j, i] = float(row["js_distance_base2"])
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        plotted = ax.imshow(matrix, vmin=0, vmax=1, cmap="cividis", interpolation="nearest")
        step = max(1, math.ceil(len(names) / 24))
        positions = list(range(0, len(names), step))
        ax.set_xticks(positions, [names[i] for i in positions], rotation=45, ha="right")
        ax.set_yticks(positions, [names[i] for i in positions])
        ax.set_title(f"{analysis_label(key[0])}\n{readable(key[1])} · {readable(key[2])}\nBetween-replica distribution differences")
        fig.colorbar(plotted, ax=ax, label="Jensen–Shannon distance (base 2)")
        token = metric_figure_stem(*key)
        _save_figure(plt, fig, output / "figures" / "diagnostics" / f"{token}__between-replica-distribution-distance.svg", files, output)
        figure_analyses[files[-1]] = key[0]


def report_results(config: dict, manifest: dict, output_dir: str | Path) -> dict:
    """Read successful artifact files, export joint tables, and render SVGs.

    Missing/unrecognized artifacts are recorded without discarding successfully
    parsed replicas. ``plots.enabled=false`` allows dependency-light numeric
    reporting. cpptraj text summaries remain linked in ``artifact_index.dat``.
    """
    output_dir = Path(output_dir).resolve()
    (output_dir / "reports").mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    series_list: list[Series] = []
    matrices: list[Matrix] = []
    generic_tables = []
    artifact_index = []
    seen = set()
    replica_entries = list(manifest.get("replicas", []))
    actual_replica_names = {str(entry.get("name", "unnamed")) for entry in replica_entries}
    # Projections/assignments have already been split into real replicas by the
    # runner. Only pooled matrices and tables belong to these labeled entries.
    pooled_entries = [{"name": f"POOLED:{entry.get('name', 'analysis')}", "status": entry.get("status", "unknown"), "artifacts": [artifact for artifact in entry.get("artifacts", []) if artifact.get("kind") != "timeseries"]} for entry in manifest.get("pooled", [])]
    for replica_entry in [*replica_entries, *pooled_entries]:
        replica = str(replica_entry.get("name", "unnamed"))
        status = str(replica_entry.get("status", "unknown"))
        if status.lower() in {"failed", "error", "skipped", "pending"}:
            warnings.append(f"{replica}: replica status is {status}; only available artifacts are considered")
        for artifact in replica_entry.get("artifacts", []):
            source = Path(str(artifact.get("path", "")))
            if not source.is_absolute():
                source = output_dir / source
            source = source.resolve()
            analysis = str(artifact.get("analysis", source.stem))
            section = str(artifact.get("section", "system"))
            kind = str(artifact.get("kind", "table"))
            unit = str(artifact.get("unit", ""))
            parsed_status = "parsed"
            try:
                if (replica, str(source)) in seen:
                    raise TableFormatError("duplicate artifact entry")
                seen.add((replica, str(source)))
                if not source.is_file():
                    raise TableFormatError("artifact file is missing")
                metadata = artifact.get("metadata") or {}
                if metadata.get("format") in {"eigenvectors", "cluster_summary"}:
                    # These are multi-block analysis reports, not rectangular
                    # scalar tables. Preserve and index intentionally.
                    artifact_index.append((replica, status, analysis, section, kind, unit, "retained_native_format", str(source)))
                    continue
                if metadata.get("format") in {"hbond_occupancy", "contact_pairs"}:
                    table = _parse_contact_table(source)
                    generic_tables.append((replica, analysis, section, unit, str(source), table))
                    artifact_index.append((replica, status, analysis, section, kind, unit, parsed_status, str(source)))
                    continue
                if metadata.get("format") == "nastruct":
                    series_list.extend(_nastruct_series(source, replica, artifact, config))
                    artifact_index.append((replica, status, analysis, section, kind, unit, parsed_status, str(source)))
                    continue
                table = parse_numeric_table(source)
                if metadata.get("format") == "stacking_geometry":
                    from .stacking import geometry_series
                    expected = replica_entry.get("frame_count")
                    if expected is not None and len(table.data) != expected:
                        raise TableFormatError("Stacking geometry does not cover the planned replica frames")
                    series_list.extend(geometry_series(replica, artifact, table, config, str(source)))
                    artifact_index.append((replica, status, analysis, section, kind, unit, parsed_status, str(source)))
                    continue
                if kind == "matrix":
                    matrices.append(_matrix_from_table(replica, artifact, table, str(source)))
                elif kind in {"timeseries", "profile"}:
                    x, xunit, coordinate_index = _coordinate(table, artifact, config, kind)
                    columns = [(i, name) for i, name in enumerate(table.columns) if i != coordinate_index]
                    if not columns:
                        raise TableFormatError("table has a coordinate column but no measured columns")
                    discrete = bool(metadata.get("discrete", False) or metadata.get("categorical", False))
                    is_dssp = discrete and (metadata.get("palette") == "dssp" or "dssp" in analysis.lower() or metadata.get("discrete", False))
                    if is_dssp:
                        data = table.data[:, [i for i, _ in columns]].T
                        finite = data[np.isfinite(data)]
                        if np.any(finite != np.floor(finite)) or np.any((finite < 0) | (finite > 7)):
                            raise TableFormatError("discrete DSSP values must be integer assignments 0–7")
                        labels = metadata.get("residue_numbers")
                        if labels is None:
                            extracted = [re.search(r"(?:\[|_)(\d+)\]?$", name) for _, name in columns]
                            labels = [int(match.group(1)) for match in extracted] if all(extracted) else list(range(1, len(columns) + 1))
                        matrices.append(Matrix(replica, analysis, section, x, np.asarray(labels, dtype=float), data, "", str(source), "dssp", xunit, "residue" if metadata.get("residue_numbers") or all(re.search(r"(?:\[|_)(\d+)\]?$", name) for _, name in columns) else "selected residue index", True, bool(metadata.get("betadetail", False))))
                    for i, name in columns:
                        series_list.append(Series(replica, analysis, section, name, kind, x, table.data[:, i], unit, xunit, str(source), discrete, bool(metadata.get("circular", False)), bool(metadata.get("pooled", False))))
                else:
                    generic_tables.append((replica, analysis, section, unit, str(source), table))
            except (OSError, ValueError, IndexError) as exc:
                parsed_status = "unparsed"
                warnings.append(f"{replica}/{analysis}/{section}: {exc}; original artifact retained")
            artifact_index.append((replica, status, analysis, section, kind, unit, parsed_status, str(source)))
    from .stacking import write_stacking_report
    stacking_summary, stacking_unions = write_stacking_report(config, series_list, manifest, output_dir)
    series_list.extend(stacking_unions)
    if stacking_summary.get("omissions"):
        warnings.append(f"Stacking: {stacking_summary['omissions']} incomplete pair/target/residue replica results; see reports/stacking/omissions.dat")
    series_groups: dict[tuple, list[Series]] = defaultdict(list)
    for series in series_list:
        series_groups[series.key].append(series)
    matrix_groups: dict[tuple, list[Matrix]] = defaultdict(list)
    for matrix in matrices:
        matrix_groups[matrix.key].append(matrix)
    max_matrix_cells = int(config.get("reports", {}).get("max_matrix_cells", 1_000_000))

    def tidy_rows():
        for series in series_list:
            for x, value in zip(series.x, series.values):
                yield (series.replica, series.analysis, series.section, series.name, series.kind, x, None, value, series.unit, series.xunit, series.source)
        for matrix in matrices:
            if matrix.discrete:
                continue  # Already represented by the residue series above.
            if matrix.values.size > max_matrix_cells:
                warnings.append(f"{matrix.replica}/{matrix.analysis}: {matrix.values.size} matrix cells exceed reports.max_matrix_cells; source remains in artifact_index.dat")
                continue
            for j, y in enumerate(matrix.y):
                for i, x in enumerate(matrix.x):
                    yield (matrix.replica, matrix.analysis, matrix.section, "matrix", "matrix", x, y, matrix.values[j, i], matrix.unit, f"{matrix.xunit}; {matrix.yunit}", matrix.source)
        for replica, analysis, section, unit, source, table in generic_tables:
            for row_index, row in enumerate(table.data, 1):
                for column, value in zip(table.columns, row):
                    identity = table.row_labels[row_index - 1] if table.row_labels is not None else ""
                    series_name = f"{identity}|{column}" if identity else column
                    column_unit = table.column_units.get(column, unit) if table.column_units is not None else unit
                    yield (replica, analysis, section, series_name, "table", row_index, None, value, column_unit, "table row (no inferred physical coordinate)", source)

    _write_table(output_dir / "reports" / "all_replicas.dat", TIDY_COLUMNS, tidy_rows())
    _write_table(output_dir / "reports" / "artifact_index.dat", ("replica", "replica_status", "analysis", "section", "kind", "unit", "parse_status", "source"), artifact_index)
    _series_tables(series_groups, output_dir, int(config.get("stats", {}).get("block_size", 50)), warnings)
    from .diagnostics import write_diagnostics
    diagnostic_summary = write_diagnostics(config, series_list, output_dir)
    from .convergence import write_convergence
    convergence_summary = write_convergence(config, series_list, output_dir)
    if config.get("stats", {}).get("metric_correlations", False):
        _scalar_correlations(series_list, output_dir, warnings)
    figure_files: list[str] = []
    figure_analyses = {}
    plot = bool(config.get("plots", {}).get("enabled", True))
    plots_status = "disabled" if not plot else "complete"
    plt = None
    if plot:
        try:
            plt = _plotting(int(config.get("plots", {}).get("dpi", 600)), config.get("plots", {}))
        except ImportError:
            plot = False
            plots_status = "unavailable"
            warnings.append("Matplotlib is unavailable: numerical reports were written; install matplotlib and rerun report to create SVGs")
    if plot:
        REPLICA_COLORS.clear()
        REPLICA_COLORS.update({str(rep.get("name", "unnamed")): PALETTE[i % len(PALETTE)] for i, rep in enumerate(replica_entries)})
        if diagnostic_summary.get("enabled"):
            try:
                _diagnostic_figures(plt, output_dir, figure_files, figure_analyses)
            except (OSError, ValueError, RuntimeError, OverflowError) as exc:
                plots_status = "partial"
                plt.close("all")
                warnings.append(f"Diagnostic SVG rendering failed: {exc}; diagnostic tables remain available")
        if convergence_summary.get("enabled"):
            try:
                from .convergence import plot_convergence
                plot_convergence(plt, output_dir, _save_figure, figure_files, figure_analyses)
            except (OSError, ValueError, RuntimeError, OverflowError) as exc:
                plots_status = "partial"
                plt.close("all")
                warnings.append(f"Convergence SVG rendering failed: {exc}; numerical tables remain available")
        for key, group in series_groups.items():
            before = len(figure_files)
            try:
                _series_figures(plt, {key: group}, config, output_dir, figure_files)
            except (ValueError, RuntimeError, OverflowError) as exc:
                plots_status = "partial"
                plt.close("all")
                warnings.append(f"{group[0].analysis}/{group[0].section}/{group[0].name}: SVG rendering failed: {exc}; numerical reports remain available")
            figure_analyses.update({name: group[0].analysis for name in figure_files[before:]})
        for plot_function in (_pca_figures, _cluster_figures):
            before = len(figure_files)
            try:
                plot_function(plt, series_list, output_dir, figure_files)
            except (ValueError, RuntimeError, OverflowError) as exc:
                plots_status = "partial"
                plt.close("all")
                warnings.append(f"{plot_function.__name__}: SVG rendering failed: {exc}; numerical reports remain available")
            figure_analyses.update({name: "pca" if plot_function is _pca_figures else "cluster" for name in figure_files[before:]})
        try:
            _contact_figures(plt, generic_tables, output_dir, figure_files, figure_analyses)
        except (ValueError, RuntimeError, OverflowError) as exc:
            plots_status = "partial"
            plt.close("all")
            warnings.append(f"Contact-fraction SVG rendering failed: {exc}; numerical reports remain available")
    # Always finish matrix statistics, even when a plotting backend fails.
    _matrix_reports(None, matrix_groups, output_dir, figure_files, warnings, False)
    if plot:
        for key, group in matrix_groups.items():
            before = len(figure_files)
            try:
                _matrix_reports(plt, {key: group}, output_dir, figure_files, [], True)
            except (ValueError, RuntimeError, OverflowError) as exc:
                plots_status = "partial"
                plt.close("all")
                warnings.append(f"{group[0].analysis}/{group[0].section}: matrix SVG rendering failed: {exc}; numerical reports remain available")
            figure_analyses.update({name: group[0].analysis for name in figure_files[before:]})
    notes = ["Trajectory-frame SD, quantiles, histograms and running averages are descriptive, not uncertainty estimates.", "Between-replica SD describes the spread of replica means; each available replica receives equal weight. One replica has no between-replica SD.", "Angles flagged circular use circular means and resultant lengths, with means undefined when the resultant is zero. They are excluded from arithmetic means and ordinary Pearson correlations.", "Block-mean SD uses complete contiguous blocks of the requested analyzed-frame size; it is not an effective sample-size estimate or a confidence interval. Choose blocks longer than the observable correlation time.", "Independent replicas are appended only for visualization, with boundary markers; the appended axis is not a continuous physical trajectory.", "Matrix and profile means require matching axes and shapes. Matrix cell means require finite values from every included replica.", "Matrix index compatibility does not independently establish identical atom identities or alignment; consistent selection, topology and reference setup are required.", "Replica groups with different lengths are overlaid on their own coordinates. Missing observations are not interpolated.", "Unknown table formats are retained as original artifacts and indexed. Numeric generic tables are exported by row without assigning physical meaning.", "POOLED entries label shared analysis outputs, not additional independent replicas. Pooled PCA covariance weights frames; longer replicas contribute more to the common basis. PCA scatterplots use reported projection coordinates without inferring explained variance or free energy.", "DSSP integer states are categories. Their fractions are reported; arithmetic means of state IDs are not calculated."]
    if config.get("frames", {}).get("dt_ps") is None:
        notes.append("Saved-frame spacing was not supplied; time-series axes use analyzed frame rather than inferred time.")
    data_replica_names = {series.replica for series in series_list} | {matrix.replica for matrix in matrices} | {item[0] for item in generic_tables}
    summary = {"status": "complete_with_warnings" if warnings else "complete", "plots_status": plots_status, "replicas_requested": len(replica_entries), "replicas_with_numeric_data": len(data_replica_names & actual_replica_names), "pooled_groups_requested": len(pooled_entries), "series_count": len(series_list), "matrix_count": len(matrices), "figures": figure_files, "combined_data": str(output_dir / "reports" / "all_replicas.dat"), "warnings": warnings, "notes": notes}
    summary["diagnostics"] = diagnostic_summary
    summary["convergence"] = convergence_summary
    summary["stacking"] = stacking_summary
    summary["figure_style"] = getattr(plt, "_mdwb_style", None) if plt is not None else None
    summary["figure_analyses"] = figure_analyses
    summary["figure_titles"] = getattr(plt, "_mdwb_figure_titles", {}) if plt is not None else {}
    summary["viewer"] = "index.html"
    (output_dir / "reports" / "report_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# Analysis reporting notes", "", *[f"- {note}" for note in notes], "", "## Warnings", "", *([f"- {warning}" for warning in warnings] or ["None."]), "", "## Main outputs", "", "- all_replicas.dat: one tidy TSV with replica labels and source paths.", "- descriptive_statistics.dat: per-replica descriptive statistics and contiguous block diagnostics.", "- replicate_statistics.dat: equal-weight replica mean and between-replica SD.", "- circular_statistics.dat and circular_replicate_statistics.dat: angular descriptions in degrees.", "- categorical_fractions.dat: state fractions within each replica.", "- combined/: aligned per-observable tables; blank values represent missing coordinates.", "- probability_distributions.dat: common bin edges, counts, probabilities and densities per replica; scalar time series only.", "- artifact_index.dat: every requested artifact, its source, and parse status.", "- ../figures/: editable-text SVGs; dense matrix cells and PCA points are embedded rasters at the configured DPI."]
    (output_dir / "reports" / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    from .dashboard import write_dashboard
    from .interactive_report import write_preview_data
    write_preview_data(config, series_list, output_dir)
    write_dashboard(config, manifest, summary, output_dir)
    return summary
