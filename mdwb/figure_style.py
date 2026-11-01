"""Publication figure presets and explicit local-font resolution."""
PALETTES = {
    "lagoon": ("#087E8B", "#C86B3C", "#5975A4", "#8D6298", "#658B55", "#B14C66", "#A38337", "#4B5662"),
    "mineral": ("#396B89", "#B76B56", "#6B8C83", "#8B739B", "#BE994E", "#647080", "#AF6487", "#74854C"),
    "colorblind": ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#8B7400", "#333333"),
}


def configure(plt, options, dpi):
    from matplotlib import font_manager
    requested = options.get("font", "Arial")
    available = {font.name for font in font_manager.fontManager.ttflist}
    chosen = next((name for name in (requested, "Liberation Sans", "DejaVu Sans") if name in available), "DejaVu Sans")
    size = options.get("font_size", 11)
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [chosen], "font.size": size,
        "axes.labelsize": size, "axes.titlesize": size + 1, "axes.titleweight": "medium",
        "axes.titlepad": 12, "axes.labelpad": 7, "legend.fontsize": size - 1,
        "xtick.labelsize": size - 1, "ytick.labelsize": size - 1,
        "axes.linewidth": 0.8, "lines.linewidth": 1.55, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3.5, "ytick.major.size": 3.5, "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "text.color": "#202A32", "axes.labelcolor": "#202A32", "axes.edgecolor": "#45515A",
        "xtick.color": "#45515A", "ytick.color": "#45515A",
        "svg.fonttype": "path" if options.get("svg_text") == "paths" else "none",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.hashsalt": "cpptraj-workbench",
        "savefig.dpi": dpi, "savefig.pad_inches": 0.08,
    })
    return {"requested_font": requested, "resolved_font": chosen, "font_fallback": chosen != requested,
            "font_size_pt": size, "width_mm": options.get("width_mm", 180),
            "palette": options.get("palette", "lagoon"), "raster_dpi": dpi,
            "svg_text": options.get("svg_text", "editable"),
            "note": "Text remains editable unless paths is selected; raster DPI applies only to embedded raster layers."}
