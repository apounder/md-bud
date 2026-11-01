"""Readable terminal summaries; no ANSI escapes in redirected cluster logs."""
import os
import shutil
import sys
import textwrap
from datetime import datetime


def panel(title, rows):
    width = max(50, min(96, shutil.get_terminal_size((88, 24)).columns - 2))
    color = sys.stdout.isatty() and "NO_COLOR" not in os.environ and os.environ.get("TERM") != "dumb"
    label = f"\033[1;36m{title}\033[0m" if color else title
    print("\n" + "=" * width + "\n" + label + "\n" + "-" * width)
    for key, value in rows:
        prefix = f"{key:<19} "
        for line in textwrap.wrap(str(value), width=width - 20, break_long_words=True) or [""]:
            print(prefix + line)
            prefix = " " * 20
    print("=" * width)


def event(message):
    print(f"[{datetime.now():%H:%M:%S}] {message}", flush=True)


def result_summary(summary, output):
    panel("REPORT READY", [("Report status", summary.get("status", "unknown")),
          ("Replicas with data", summary.get("replicas_with_numeric_data", "unknown")),
          ("Numerical series", summary.get("series_count", "unknown")),
          ("Publication figures", f"{len(summary.get('figures', []))} SVGs; {summary.get('plots_status', 'unknown')}"),
          ("Browser previews", "Embedded in index.html; copy to your computer for offline viewing"),
          ("Results", output), ("Review notes", f"{len(summary.get('warnings', []))} reporting warnings; inspect sampling review in HTML")])
