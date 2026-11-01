"""Generate the self-contained report and supporting scientific review records."""
from pathlib import Path


def write_dashboard(config, manifest, summary, output):
    from .report_pages import write_review_pages
    from .interactive_report import write_interactive_report
    output = Path(output)
    # Retain the detailed review records and supporting pages in result bundles.
    write_review_pages(config, manifest, summary, output)
    write_interactive_report(config, manifest, summary, output)
