"""Evidence-linked, deterministic interpretation pages for offline MD reports."""
from __future__ import annotations

import csv
from .figure_names import LAG_EXPLANATION
import hashlib
import html
import math
from pathlib import Path
from urllib.parse import quote


def _read(output, name):
    path = output / "reports" / name
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def _escape(value):
    return html.escape(str(value), quote=True)


def _table(rows, columns, limit=200):
    if not rows:
        return '<p>No applicable numerical records are available for this section.</p>'
    parts = ['<div class="scroll"><table><thead><tr>' + ''.join('<th>' + _escape(label) + '</th>' for _, label in columns) + '</tr></thead><tbody>']
    for row in rows[:limit]:
        cells = []
        for key, _ in columns:
            value = row.get(key, "")
            numeric = _number(value)
            displayed = f"{numeric:.5g}" if numeric is not None else value or "Unavailable"
            cells.append('<td>' + _escape(displayed) + '</td>')
        parts.append('<tr>' + ''.join(cells) + '</tr>')
    parts.append('</tbody></table></div>')
    if len(rows) > limit:
        parts.append(f'<p>Showing {limit} of {len(rows)} records. Download the linked data for all records.</p>')
    return ''.join(parts)


def _link(path, label):
    return '<a href="../' + quote(path, safe="/") + '">' + _escape(label) + '</a>'


GUIDANCE = {
    "stacking": "Inspect confirmed pi-system atoms and the centroid-distance/plane-angle criteria. Occupancies are geometric contacts, not interaction energies or lifetimes. Replica SD describes the spread of occupancy percentages; retain zero-occupancy pairs and review missing geometry.",
    "stacking_any": "A target/residue is occupied when any confirmed partner meets the stacking criteria. Simultaneous contacts count once; incomplete partner data prevent a whole-frame union estimate.",
    "rdf": "Inspect pair selections, box-volume normalization and the radial cutoff before comparing peaks; a heterogeneous solute environment need not approach g(r)=1, and peaks alone do not establish coordination numbers.",
    "watershell": "Inner and outer counts are cumulative, so the outer count includes the inner population; these geometric counts are not residence times, and complete one-residue solvent molecules give consistent CPU/CUDA definitions.",
    "rmsd": "Inspect plateaus, transitions and replica distributions. A flat RMSD can coexist with an unvisited state. Compare global and local fits to distinguish domain motion from internal deformation.",
    "rmsf": "Compare residue-wise profiles and between-replica spread on matching residue axes. Alignment and different state populations can change the apparent flexibility.",
    "rg": "Check for sustained compaction/expansion and different replica populations. Similar radius of gyration does not imply identical conformations.",
    "distance": "Inspect state transitions and dwell periods; verify atom/group definitions and minimum-image choices before interpreting separation.",
    "angle": "Check the ordered masks and vertex. Different means may conceal multimodal angular populations.",
    "dihedral": "Use circular summaries and distributions. Values near −180 and +180 degrees are neighbors; an arithmetic average is misleading.",
    "torsions": "Inspect circular populations and transitions for each torsion. A long trajectory in one rotamer can yield small fluctuations without representative sampling.",
    "hbond": "Compare occupancies across replicas and inspect geometric definitions. Occupancy is not a bond lifetime, and missing hydrogen/bond information can alter detection.",
    "contacts": "Compare contact fractions and residue maps across replicas. Contact cutoffs and the reference define what is called native; geometric contacts do not establish chemistry-specific interactions.",
    "dssp": "Compare secondary-structure state fractions per residue. State identifiers are categories and must not be averaged as continuous measurements.",
    "sasa": "Inspect exposure changes with the corresponding structures. Surface-area differences depend on selection and the LCPO model; they are not binding free energies.",
    "nastruct": "Check detected base-pair identities before comparing geometry. Missing/unrecognized pairs are not zero deformation; inspect each reported pair and parameter.",
    "pca": "Use the common pooled basis for comparisons. Longer replicas influence the basis more. Overlapping projections do not establish convergence of unplotted coordinates.",
    "cluster": "Compare state populations and transitions in the shared clustering. Results depend on atom selection, cluster count and sieve; these populations are not kinetic rate estimates.",
    "dccm": "Inspect corresponding matrix cells and their between-replica spread. Atom order and alignment must match; correlation is not causation or a direct allosteric pathway.",
    "pairwise_rmsd": "Inspect recurrent structural blocks and transitions. An isolated block may represent a state change, but the matrix alone does not establish state populations or kinetics.",
}


def write_review_pages(config, manifest, summary, output):
    output = Path(output)
    destination = output / "report_pages"
    destination.mkdir(exist_ok=True)
    stats = _read(output, "descriptive_statistics.dat")
    replicas = _read(output, "replicate_statistics.dat")
    circular = _read(output, "circular_replicate_statistics.dat")
    states = _read(output, "categorical_fractions.dat")
    artifacts = _read(output, "artifact_index.dat")
    enabled = summary.get("diagnostics", {}).get("enabled", False)
    sampling = _read(output, "diagnostics/sampling.dat") if enabled else []
    distances = _read(output, "diagnostics/replica_distribution_distance.dat") if enabled else []
    requested = len(manifest.get("replicas", []))
    findings = []

    def finding(subject, evidence, advice):
        findings.append({"subject": subject, "evidence": evidence, "advice": advice})

    if requested < 2:
        finding("Replica coverage", f"{requested} replica requested.", "Between-replica variability cannot be estimated with one replica; frames are not independent simulation replicas.")
    elif requested == 2:
        finding("Replica coverage", "Two replicas requested.", "The between-replica SD is based on two means and is sensitive to either run. Inspect both trajectories; agreement alone does not establish representative sampling.")
    for entry in manifest.get("replicas", []) + manifest.get("pooled", []):
        if entry.get("status") not in {"complete", "synthetic"}:
            finding(entry["name"], "Calculation status: " + entry.get("status", "unknown"), "Inspect the manifest and native logs. Missing analyses are not zero measurements.")
    for row in replicas:
        count = _number(row.get("n_replicas"))
        if count is not None and count < requested:
            finding('/'.join(row.get(k, "") for k in ("analysis", "section", "series")), f"Summary contains {int(count)} of {requested} replicas.", "Check which replicas supplied finite values before comparing averages.")
    for row in artifacts:
        if row.get("parse_status") == "unparsed":
            finding('/'.join(row.get(k, "") for k in ("replica", "analysis", "section")), "An artifact could not be interpreted numerically.", "Inspect the original artifact and parser warning; this output is not included in numerical summaries.")
    if enabled:
        for row in _read(output, "diagnostics/distribution_omissions.dat"):
            finding('/'.join(row.get(k, "") for k in ("analysis", "section", "series")), "Distribution comparison omitted: " + row.get("reason", "unknown"), "Inspect sample availability and diagnostic resource limits. No distribution agreement assessment was made for this group.")
    for row in sampling:
        subject = '/'.join(row.get(k, "") for k in ("replica", "analysis", "section", "series"))
        status = row.get("status", "unknown")
        if status != "estimate_requires_stationarity":
            finding(subject, "Sampling estimate unavailable: " + status, "Inspect raw values and sampling spacing; a missing estimate is not evidence of good or poor convergence.")
            continue
        ess = _number(row.get("effective_samples_estimate"))
        if ess is not None and ess < 50:
            finding(subject, f"Estimated effective samples: {ess:.4g} (<50 screening threshold).", "Inspect the ACF, transitions and trajectory length. The threshold is a review heuristic, not a validated convergence cutoff.")
        if row.get("lag_limit_reached", "").lower() in {"true", "1"}:
            finding(subject, "Autocorrelation summation reached its lag limit.", "The estimate may miss slow correlation. Inspect the ACF and consider longer sampling or a larger diagnostic lag limit.")
        if row.get("configured_block_shorter_than_suggestion", "").lower() in {"true", "1"}:
            finding(subject, "Configured block size is below ceil(5 × estimated g).", "Inspect block-size sensitivity; this recommendation is heuristic and does not ensure independent blocks.")
    if not enabled:
        finding("Sampling diagnostics", "Diagnostics were disabled or are unavailable in this report.", "Regenerate with diagnostics enabled to inspect ACF and effective-sample estimates. No sampling assessment was made here.")

    identity = [("analysis", "Analysis"), ("section", "Section"), ("series", "Observable")]
    replica_columns = identity + [("unit", "Unit"), ("n_replicas", "Replicas"), ("equal_replica_mean", "Mean of replica means"), ("between_replica_mean_sd", "SD of replica means")]
    within_columns = [("replica", "Replica"), *identity, ("unit", "Unit"), ("n_finite", "Finite samples"), ("mean", "Mean"), ("sd", "Within-replica SD"), ("block_mean_sd", "Block-mean SD")]

    def page(filename, title, body):
        document = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>' + _escape(title) + '</title><style>'
        document += 'body{font:16px/1.6 system-ui,sans-serif;background:#f3f6fa;color:#172b40;margin:0}main{max-width:1200px;margin:auto;padding:2rem}a{color:#006c92}nav{margin:1rem 0}section{background:white;padding:1.2rem;border:1px solid #ccd7e3;border-radius:12px;margin:1rem 0}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:.6rem;text-align:left;border-bottom:1px solid #dae1e9;vertical-align:top}th{background:#e8eff5}img{width:100%;max-height:600px;object-fit:contain}figure{margin:1rem 0}h1,h2{line-height:1.2}@media print{nav{display:none}body{background:white}section{break-inside:avoid}}'
        document += '</style></head><body><main><nav><a href="../index.html">Overview & figures</a> · <a href="sampling.html">Sampling review</a> · <a href="replicas.html">Replica comparison</a></nav><h1>' + _escape(title) + '</h1><p>' + _escape(config.get("study", {}).get("title", "MD analysis")) + '</p>' + body + '</main></body></html>'
        (destination / filename).write_text(document, encoding="utf-8")

    finding_table = _table(findings, [("subject", "Observable / scope"), ("evidence", "Recorded evidence"), ("advice", "What to inspect next")])
    body = '<section><h2>Items for scientific review</h2><p>These deterministic screening notes do not assign a convergence score. Absence of flags does not establish adequate exploration. Equilibration, stationarity and observable relevance require scientific judgment.</p>' + finding_table + '<p>' + _link("reports/review_findings.dat", "Download all review items") + '</p></section>'
    body += '<section><h2>Autocorrelation and changes over time</h2><p>' + LAG_EXPLANATION + '</p><p>ESS and SEM assume stationary sampling. Compare the second-half minus first-half mean with the original trace and its within-replica variability; this difference is descriptive, not a significance test. ACF estimates omit circular and categorical signals. Inspect their state populations and transitions separately.</p>'
    body += _table(sampling, [("replica", "Replica"), *identity, ("unit", "Unit"), ("status", "Estimate status"), ("n", "Samples"), ("effective_samples_estimate", "Estimated ESS"), ("sem_under_stationarity", "Conditional SEM"), ("second_minus_first_half_mean", "Second − first half mean")])
    if enabled:
        body += '<p>' + _link("reports/diagnostics/sampling.dat", "Full sampling table") + ' · ' + _link("reports/diagnostics/README.md", "Estimator definitions and limitations") + '</p>'
    body += '</section>'
    if summary.get("convergence", {}).get("enabled"):
        body += '<section><h2>Convergence over time</h2><p>Cumulative and trailing-window comparisons share an elapsed-time horizon. Self-derived pooled/tail references are descriptive targets, not independent gold standards. SD describes fluctuations, not uncertainty in the mean. Classical split R-hat assumes comparable independent replicas and is not rank-normalized; no automatic convergence threshold is applied.</p>'
        body += _table(_read(output, "convergence/within_replica.dat"), [("metric", "Metric"), ("replica", "Replica"), ("elapsed", "Elapsed"), ("coordinate_unit", "Axis unit"), ("scope", "Scope"), ("mean", "Mean"), ("sd", "SD"), ("mean_error", "Reference mean difference"), ("reference_js_distance", "Reference JS distance")])
        body += _table(_read(output, "convergence/between_replicas.dat"), [("metric", "Metric"), ("replica_a", "Replica A"), ("replica_b", "Replica B"), ("elapsed", "Elapsed"), ("scope", "Scope"), ("js_distance", "JS distance"), ("all_replica_classical_split_rhat", "All-replica classical split R-hat")])
        for name in ("within_replica", "between_replicas", "state_populations", "references", "reference_histograms", "omissions"):
            body += '<p>' + _link(f"reports/convergence/{name}.dat", name.replace('_', ' ')) + '</p>'
        body += '</section>'
    page("sampling.html", "Sampling and convergence review", body)

    body = '<section><h2>How to interpret the different standard deviations</h2><p><strong>Within-replica SD</strong> describes fluctuations of an observable along one trajectory. <strong>SD of replica means</strong> describes the spread of the independent simulation averages, with each available replica weighted equally. <strong>Block-mean SD</strong> describes complete contiguous block averages within a replica. None of these is automatically a confidence interval.</p><p>With one replica the between-replica SD is unavailable, not zero. Similar means can conceal different distributions. A large spread needs inspection of state populations, setup differences and sampling; there is no universal acceptable SD across observables and units.</p></section>'
    body += '<section><h2>Replica means and their spread</h2>' + _table(replicas, replica_columns) + '<p>' + _link("reports/replicate_statistics.dat", "Full replica summary") + '</p></section>'
    body += '<section><h2>Individual replica statistics</h2>' + _table(stats, within_columns) + '<p>' + _link("reports/descriptive_statistics.dat", "Full individual statistics") + '</p></section>'
    body += '<section><h2>Observed distribution differences</h2><p>Jensen–Shannon distance ranges from 0 to 1 using base-2 logarithms. Compare the actual distributions and shared bins. Values depend on binning and sample size; no pass/fail threshold is applied. The table is ordered by distance to help inspection, not statistical significance.</p>'
    body += _table(sorted(distances, key=lambda r: _number(r.get("js_distance_base2")) or 0, reverse=True), identity + [("replica_a", "Replica A"), ("replica_b", "Replica B"), ("n_a", "Samples A"), ("n_b", "Samples B"), ("js_distance_base2", "JS distance")])
    if enabled:
        body += '<p>' + _link("reports/diagnostics/replica_distribution_distance.dat", "Full distribution comparisons") + '</p>'
    body += '</section><section><h2>Circular and categorical observables</h2><p>Torsion means and spreads use circular definitions; secondary-structure states use fractions. Linear SDs are not substituted.</p>'
    body += _table(circular, identity + [("n_replicas", "Replicas"), ("equal_replica_circular_mean_degree", "Circular mean (degrees)"), ("between_replica_circular_sd_degree", "Circular spread (degrees)")]) + '</section>'
    if summary.get("stacking", {}).get("enabled"):
        body += '<section><h2>Pi-stacking occupancies</h2><p>Each replica contributes one occupancy fraction. The mean weights replicas equally; SD is the sample SD of replica occupancies (percentage points), unavailable for one replica. Target and residue occupancies count frames with any partner, without summing simultaneous contacts. Missing or invalid geometry is not counted as an unstacked frame; incomplete replicas are excluded from aggregate statistics and identified explicitly.</p>'
        body += _table(_read(output, "stacking/replica_summary.dat"), [("scope", "Scope"), ("display_name", "Pi system / pair"), ("n_replicas_complete", "Complete replicas"), ("n_replicas_requested", "Requested replicas"), ("mean_percent", "Mean occupancy (%)"), ("sd_percentage_points", "Replica SD (percentage points)"), ("status", "Status")])
        for name in ("per_replica", "replica_summary", "pairs", "rings", "criteria", "omissions"):
            body += '<p>' + _link(f"reports/stacking/{name}.dat", name.replace('_',' ')) + '</p>'
        body += '</section>'
    page("replicas.html", "Replica agreement and variability", body)

    analysis_pages = []
    names = sorted({row.get("analysis", "") for row in artifacts + stats + replicas + circular + states} - {""})
    for analysis in names:
        filename = "analysis_" + hashlib.sha256(analysis.encode()).hexdigest()[:16] + ".html"
        select = lambda rows: [row for row in rows if row.get("analysis") == analysis]
        body = '<section><h2>What to look for</h2><p>' + _escape(GUIDANCE.get(analysis, "Inspect the reported units, selections and original CPPTRAJ output before interpreting these values.")) + '</p></section>'
        body += '<section><h2>Numerical summaries</h2>' + _table(select(replicas), replica_columns) + _table(select(stats), within_columns)
        if select(circular):
            body += _table(select(circular), identity + [("n_replicas", "Replicas"), ("equal_replica_circular_mean_degree", "Circular mean (degrees)"), ("between_replica_circular_sd_degree", "Circular spread (degrees)")])
        if select(states):
            body += _table(select(states), [("replica", "Replica"), *identity, ("state", "State"), ("fraction", "Fraction")])
        body += '<p>Profiles and matrices retain their coordinate-wise summaries in the data files and figures; no scalar mean across residues or matrix cells is inferred here.</p></section>'
        body += '<section><h2>Artifact coverage</h2>' + _table(select(artifacts), [("replica", "Replica"), ("section", "Section"), ("kind", "Output type"), ("replica_status", "Replica status"), ("parse_status", "Parser status")]) + '<p>' + _link("reports/artifact_index.dat", "Full artifact index") + '</p></section>'
        body += '<section><h2>Figures</h2>'
        matched = [name for name in summary.get("figures", []) if summary.get("figure_analyses", {}).get(name) == analysis]
        for name in matched:
            path = Path(name)
            if path.is_absolute() or ".." in path.parts:
                continue
            relative = path.as_posix()
            body += '<figure><a href="../' + quote(relative, safe="/") + '"><img loading="lazy" src="../' + quote(relative, safe="/") + '" alt="' + _escape(path.stem) + '"></a><figcaption>' + _escape(path.stem) + '</figcaption></figure>'
        if not matched:
            body += '<p>No figures were indexed for this analysis in this reporting pass.</p>'
        body += '</section>'
        page(filename, analysis.upper() + " results", body)
        analysis_pages.append((analysis, "report_pages/" + filename))

    with (output / "reports" / "review_findings.dat").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["subject", "evidence", "advice"], delimiter="\t")
        writer.writeheader()
        writer.writerows(findings)
    intro = '<section class="card"><h2>Scientific results summary</h2><p>' + str(len(findings)) + ' items identified for review. These are observable-specific screening notes, not a convergence verdict.</p><p><a href="report_pages/sampling.html">Sampling and convergence review</a> · <a href="report_pages/replicas.html">Replica agreement and standard deviations</a></p><h3>Analysis pages</h3><ul>'
    for name, path in analysis_pages:
        intro += '<li><a href="' + quote(path, safe="/") + '">' + _escape(name.upper()) + '</a></li>'
    if not analysis_pages:
        intro += '<li>No parsed analysis entries available.</li>'
    intro += '</ul><h3>Numerical overview</h3><p>First eight linear-observable summaries; the replica-comparison page contains the complete table. SD here means the spread of replica means.</p>' + _table(replicas[:8], replica_columns)
    intro += '<p>Report narrative is generated from recorded outputs with explicit rules; it does not infer unmeasured mechanisms or certify equilibration.</p></section>'
    return intro
