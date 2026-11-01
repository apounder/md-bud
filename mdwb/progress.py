"""Dependency-free terminal activity display; percentages only come from CPPTRAJ."""
import re
import sys
import threading
import time
from pathlib import Path


def batch_bar(done, total, failed=0):
    width = 24
    fill = int(width * done / total) if total else 0
    return f"Analyses finished [{'=' * fill}{'.' * (width - fill)}] {done}/{total} | failed: {failed} (batch count, not elapsed-work percentage)"


class CalculationProgress:
    def __init__(self, log_path, enabled=True):
        self.path = Path(log_path)
        self.enabled = enabled and sys.stderr.isatty()
        self.stop = threading.Event()
        self.thread = None

    def __enter__(self):
        if self.enabled:
            self.thread = threading.Thread(target=self._render, daemon=True)
            self.thread.start()
        return self

    def _render(self):
        start, tick, last_width = time.monotonic(), 0, 0
        while not self.stop.is_set():
            percent = None
            try:
                with self.path.open("rb") as handle:
                    handle.seek(0, 2)
                    handle.seek(max(0, handle.tell() - 65536))
                    tail = handle.read().decode("utf-8", errors="replace")
                matches = re.findall(r"(?<![\d.])(\d{1,3}(?:\.\d+)?)\s*%", tail)
                if matches and 0 <= float(matches[-1]) <= 100:
                    percent = float(matches[-1])
            except OSError:
                pass
            elapsed = int(time.monotonic() - start)
            if percent is None:
                spinner = "|/-\\"[tick % 4]
                detail = f"{spinner} working; no percentage reported"
            else:
                fill = int(percent / 5)
                detail = f"[{'=' * fill}{'.' * (20 - fill)}] cpptraj last reported {percent:g}%"
            line = f"  {detail} | {elapsed // 60:02}:{elapsed % 60:02} elapsed"
            try:
                sys.stderr.write("\r" + line.ljust(last_width))
                sys.stderr.flush()
            except (OSError, ValueError):
                return
            last_width = len(line)
            tick += 1
            self.stop.wait(0.25)
        try:
            sys.stderr.write("\r" + " " * last_width + "\r")
            sys.stderr.flush()
        except (OSError, ValueError):
            pass

    def __exit__(self, *exc):
        self.stop.set()
        if self.thread:
            self.thread.join(timeout=1)
