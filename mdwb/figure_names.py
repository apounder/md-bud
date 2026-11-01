"""Readable figure identities with a short suffix to prevent sanitized-name collisions."""
import hashlib
import json
import re
from types import SimpleNamespace

from .labels import analysis_label, observable_label


def metric_figure_stem(analysis, section, series, *identity):
    observable = observable_label(SimpleNamespace(name=series, section=section))
    parts = [analysis_label(analysis), str(section), observable]
    clean = [re.sub(r'[^a-z0-9]+', '-', p.lower()).strip('-') or 'metric' for p in parts]
    label = '__'.join(dict.fromkeys(clean))[:100].rstrip('-_')
    raw = json.dumps([analysis, section, series, *identity], ensure_ascii=True)
    suffix = hashlib.sha256(raw.encode()).hexdigest()[:8]
    return f'{label}__{suffix}'


LAG_EXPLANATION = (
    'Lag is the separation between two observations of the same metric within one replica. '
    'Lag 0 compares each observation with itself; lag 1 compares neighboring analyzed frames. '
    'If saved frames are 10 ps apart and the analysis stride is 5, one analyzed-frame lag is '
    '50 ps (0.05 ns), and lag 20 is 1 ns. The axis uses physical time when the saved-frame '
    'interval is known, otherwise frame units. Positive autocorrelation means values separated '
    'by that lag tend to vary together; values near zero suggest little linear correlation at '
    'that separation. Slow decay indicates fewer effectively independent samples. Lag is not '
    'elapsed simulation time or an equilibration cutoff, and an ACF near zero does not prove '
    'convergence. The maximum lag limits how far apart observations are compared; it does not '
    'discard trajectory frames.'
)
