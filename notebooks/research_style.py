"""Shared editorial boundaries and NEXUS plotting defaults for owned notebooks."""
import hashlib

from nbformat.v4 import new_code_cell, new_markdown_cell

NOTICE = """<img src="https://raw.githubusercontent.com/senanurcetin/enveda-casmi26/main/site/assets/nexus-midnight.svg" alt="NEXUS official research identity" width="360">

**NEXUS / MOLECULAR IDENTIFICATION RESEARCH** · Senanur Cetin · Noncommercial portfolio

**Research review — 23 October 2026.** Competition experiments and leaderboard tables below are historical author records from the campaign closed on 8 October. A new notebook execution recomputes only the cells shown; it does not reproduce the complete ranker or resubmit the historic predictions.

Later repository audits corrected sampled-query C1 reference exclusion and D2 eligibility handling. Historical validation tables were not recomputed under those fixes. Standard or tautomer-canonical InChIKey grouping is **not a Bemis–Murcko scaffold split**. The simplified shortcut notebook uses standard keys and is not a tautomer-leakage-resistant benchmark. The pool notebook computes coverage/decoy diagnostics; its ranker-score tables are quoted historical logs.

Treat directional observations as hypotheses under library/instrument shift and repeated experiment selection. No final private score, current rank or statistically established +0.002 fusion gain is claimed. [Model card and evidence](https://github.com/senanurcetin/enveda-casmi26/blob/main/MODEL_CARD.md).
"""

STYLE = """# NEXUS chart defaults; no remote assets or paid services needed for computation.
import matplotlib.pyplot as plt
from cycler import cycler
plt.rcParams.update({
    'figure.facecolor': '#0B0F19', 'axes.facecolor': '#061524',
    'savefig.facecolor': '#0B0F19', 'text.color': '#EAF4FF',
    'axes.labelcolor': '#B9C5CE', 'axes.edgecolor': '#52667D',
    'xtick.color': '#B9C5CE', 'ytick.color': '#B9C5CE',
    'grid.color': '#29384D', 'axes.grid': True, 'grid.alpha': 0.3,
    'axes.prop_cycle': cycler(color=['#00E0FF', '#8F9BFF', '#4EA8FF', '#DDB26F']),
    'figure.dpi': 120, 'font.size': 10, 'axes.spines.top': False,
    'axes.spines.right': False,
})
print('NEXUS research edition. Printed historical tables are not new leaderboard results.')
"""


def finalize(nb):
    nb.cells.insert(0, new_markdown_cell(NOTICE))
    nb.cells.insert(2, new_code_cell(STYLE))
    for index, cell in enumerate(nb.cells):
        cell.id = hashlib.sha256(f"{index}:{cell.cell_type}:{cell.source}".encode()).hexdigest()[:12]
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    return nb
