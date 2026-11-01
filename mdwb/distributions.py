"""Comparable empirical histograms from full finite trajectory observations."""
import numpy as np


def shared_edges(group, count=32):
    if group[0].circular:
        return np.linspace(-180.0, 180.0, 37)
    low, high = np.inf, -np.inf
    for series in group:
        finite = series.values[np.isfinite(series.values)]
        if finite.size:
            low, high = min(low, float(finite.min())), max(high, float(finite.max()))
    if not np.isfinite(low):
        return np.linspace(0.0, 1.0, count + 1)
    if low == high:
        margin = max(0.5, abs(low) * 0.01)
        return np.array([low - margin, high + margin])
    return np.linspace(low, high, count + 1)


def histogram(series, edges):
    finite = series.values[np.isfinite(series.values)]
    if series.circular:
        finite = (finite + 180.0) % 360.0 - 180.0
    counts = np.histogram(finite, bins=edges)[0]
    probability = counts / finite.size if finite.size else np.zeros(counts.size)
    return {"count": counts, "probability": probability,
            "density": probability / np.diff(edges), "n_finite": int(finite.size)}


def distribution_rows(groups):
    """Long-form rows preserve replica weights, physical units and exact bin counts."""
    for group in groups.values():
        if group[0].kind != "timeseries" or group[0].discrete:
            continue
        edges = shared_edges(group)
        for series in group:
            data = histogram(series, edges)
            for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
                yield (series.replica, series.analysis, series.section, series.name, series.unit,
                       left, right, data["count"][index], data["n_finite"],
                       data["probability"][index] if data["n_finite"] else None,
                       data["density"][index] if data["n_finite"] else None,
                       "circular" if series.circular else "linear")
