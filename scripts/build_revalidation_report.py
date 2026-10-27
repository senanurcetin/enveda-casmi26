"""Render a reviewer-readable summary only from verified completed study evidence."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.check_revalidation_publication import validate_result


def build():
    result = json.loads((ROOT / "reports/revalidation_study.json").read_text(encoding="utf-8"))
    validate_result(result)
    quality = result["candidate_source_quality"]
    structural = quality["structure_quarantine"]
    count = result["cohort"]["n_molecules"]
    lines = ["# Frozen-model retrospective revalidation", "",
             "Executed on the designated V-in population with rebuilt full-source artifacts. "
             "These are retrospective diagnostics, not new public/private competition scores or an untouched generalization benchmark.", "",
             f"- Evaluation denominator: **{count:,} molecules**, retained across all models, novelty classes and controls.",
             f"- Mass-conflicting identities quarantined: **{quality['quarantined_mass_conflict_groups']:,}**, "
             f"covering {quality['quarantined_source_rows']:,} provider rows.",
             f"- Unusable selected candidate identities quarantined: **{structural['quarantined_candidates']:,}**.",
             f"- Retained fingerprint-verified candidate structures: **{structural['retained_candidates']:,}**.",
             f"- Fresh D2 reference spectra: **{result['analog_index_provenance']['n_references']:,}**; no index cache used.",
             f"- Complete own-library exclusions add **{result['knockout_audit']['additional_same_library_rows_excluded']:,} source rows** "
             f"beyond the selected inputs for **{result['knockout_audit']['queries_with_additional_same_library_exclusions']:,} queries**.",
             f"- Queries whose standard truth key is absent from the candidate pool: **{result['queries_with_truth_absent_from_candidate_pool']:,}**; none removed from scoring.", "",
             "## Unperturbed MRR@25", "",
             "| Model | C1, 95% CI | C2, 95% CI | C3, 95% CI |", "|---|---:|---:|---:|"]
    for model, block in result["models"].items():
        cells = [f"{values['mrr']:.4f} [{values['ci_low']:.4f}, {values['ci_high']:.4f}]"
                 for values in (block["conditions"]["none"][cls] for cls in ("C1", "C2", "C3"))]
        lines.append("| " + model + " | " + " | ".join(cells) + " |")
    lines += ["", "Absolute intervals use 2,000 molecule bootstrap draws. C1 excludes all equivalent same-library "
              "references; C2 excludes equivalent references across libraries; C3 additionally removes candidate identities before generation.",
              "", "## Paired spectral controls", "",
              "Negative values mean the perturbed condition scored lower. Intervals are descriptive molecule-paired bootstrap "
              "intervals (4,000 draws), not model-promotion tests. Blank spectra retain precursor and acquisition information; "
              "intensity shuffling retains fragment positions.", "",
              "| Model | Control | C1 delta, 95% CI | C2 delta, 95% CI | C3 delta, 95% CI |",
              "|---|---|---:|---:|---:|"]
    for model, block in result["models"].items():
        for mode, control in block["controls"].items():
            cells = [f"{values['control_minus_reference']:+.4f} [{values['ci_low']:+.4f}, {values['ci_high']:+.4f}]"
                     for values in (control["classes"][cls] for cls in ("C1", "C2", "C3"))]
            lines.append("| " + model + " | " + mode + " | " + " | ".join(cells) + " |")
    lines += ["", "## Paired unperturbed model contrasts", "",
              "| Comparison | C1 delta, 95% CI | C2 delta, 95% CI | C3 delta, 95% CI |",
              "|---|---:|---:|---:|"]
    for comparison, block in result["paired_comparisons"].items():
        cells = [f"{values['delta_mrr']:+.4f} [{values['ci_low']:+.4f}, {values['ci_high']:+.4f}]"
                 for values in (block["classes"][cls] for cls in ("C1", "C2", "C3"))]
        lines.append("| " + comparison + " | " + " | ".join(cells) + " |")
    lines += ["", "## Interpretation", "",
              "In this selected population the frozen learned ranker improves C2/C3 while scoring below D1+D2 in C1. "
              "Removing fragments substantially reduces its performance. This tradeoff is class-dependent and does not establish universal model superiority.", "",
              "The simple models improve in C2 with blank spectra. With zero similarity scores their ranker falls back to "
              "mass-error tie-breaking; unperturbed residual similarities can instead promote decoys after truth references are excluded. "
              "This is an adverse diagnostic result, not a recommendation to discard spectra or a universal control pass.", "",
              "The machine-readable report also retains the prespecified fixed-grid mixture bands. "
              "Neither pointwise intervals nor grid bands correct historical model selection or domain shift.",
              "", "## Provenance and limits", "",
              "[Executed aggregate and input/model/code hashes](revalidation_study.json) · "
              "[Authenticated Kaggle readback](revalidation_publication.json) · "
              "[Protocol](../docs/studies/frozen-model-revalidation.md)", ""]
    lines += ["- " + value for value in result["limitations"]]
    (ROOT / "reports/revalidation_review.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Rendered verified retrospective study review")


if __name__ == "__main__":
    build()
