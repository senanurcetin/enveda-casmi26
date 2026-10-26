"""Derive descriptive loss accounting from the authenticated frozen study."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts.check_advanced_publication import validate_result
from scripts.check_kaggle_publications import source_hash
from casmi.validate.coverage_budget import coverage_budget


def derive(result, publication, source_sha256):
    validate_result(result)
    if (publication["aggregate_artifact_sha256"] != source_sha256 or publication["status"] != "pass"
            or publication["latest_run_status"] != "COMPLETE" or not publication["source_matches"]
            or publication["public"] is not True
            or publication["source_sha256"] != source_hash(json.loads(
                (ROOT / "kaggle/audits/advanced-internal/notebook.ipynb").read_text(encoding="utf-8")))):
        raise ValueError("Coverage accounting requires authenticated completed study evidence")
    budgets = {}
    for name, conditions in result["models"].items():
        budgets[name] = {}
        for cls in ("C1", "C2", "C3"):
            field = "pool_truth_present" if name in ("mass_only", "d1", "d1_d2") else "with_generation_truth_present"
            budgets[name][cls] = coverage_budget(conditions["none"][cls], result["candidate_coverage"][cls][field])
    return {"schema_version": 1, "status": "POSTHOC_DESCRIPTIVE_AGGREGATE_DECOMPOSITION",
            "source_artifact_sha256": source_sha256, "publication_source_sha256": publication["source_sha256"],
            "scope": "Unperturbed frozen internal test; no fitting, new test, selection or causal intervention",
            "identity": "1 = observed_mrr + within_coverage_ranking_loss + coverage_loss",
            "models": budgets,
            "limitations": [
                "Coverage is measured on the selected candidate lists, not a universal chemistry-recall estimate.",
                "The oracle ceiling assumes perfect ordering of existing candidates; it is not achieved model performance.",
                "Conditional metrics describe only covered molecules and must not replace the full test denominator.",
                "Per-query joint coverage/rank observations are not exported; conditional bootstrap intervals cannot be recovered from these aggregates.",
                "Loss accounting is algebraic and post hoc; changing generation can alter decoys and ranking, so gains are not causally identified.",
                "All prespecified variants are retained; no test-selected replacement model or promotion follows."]}


def figures(report):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.hashsalt": "nexus-casmi-coverage-v1"})
    fig, ax = plt.subplots(figsize=(12, 6), facecolor="#0B0F19")
    ax.set_facecolor("#0B0F19")
    values = [report["models"]["full"][cls] for cls in ("C1", "C2", "C3")]
    left = [0.0] * 3
    for field, label, color in (("observed_mrr", "Observed MRR@25", "#00E0FF"),
                                 ("within_coverage_ranking_loss", "Loss within existing candidates", "#B5A4FF"),
                                 ("coverage_loss", "Truth absent from candidates", "#FFCE85")):
        widths = [v[field] for v in values]
        ax.barh([2, 1, 0], widths, left=left, color=color, label=label, height=.5)
        for y, x, width in zip([2, 1, 0], left, widths):
            if width >= .06:
                ax.text(x + width / 2, y, f"{width:.3f}", ha="center", va="center", color="#0B0F19", fontsize=11, weight="bold")
        left = [x + width for x, width in zip(left, widths)]
    ax.set_yticks([2, 1, 0], ["C1", "C2", "C3"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of ideal MRR = 1 · observed selected test population", color="#C1CBD6", labelpad=12)
    ax.tick_params(colors="#C1CBD6")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, 1.25), ncol=3, frameon=False, labelcolor="#C1CBD6", fontsize=9)
    fig.text(.08, .94, "NEXUS / CASMI26", color="#00E0FF", fontsize=11, weight="bold")
    fig.text(.08, .87, "Where the fresh full model loses MRR", color="#F1F6FC", fontsize=20, weight="bold")
    c3 = report["models"]["full"]["C3"]
    fig.text(.08, .09, f"C3: {c3['truth_candidates_absent']}/250 truths absent; {c3['top25_hits']}/{c3['truth_candidates_present']} covered truths reach top 25; {c3['top1_hits']} rank first.", color="#C1CBD6", fontsize=11)
    fig.text(.08, .045, "Post hoc algebraic accounting. Oracle ordering ceiling is hypothetical; no causal gain or new validation is measured.", color="#91A9C1", fontsize=9)
    fig.subplots_adjust(left=.10, right=.97, top=.65, bottom=.25)
    folder = ROOT / "docs/assets"
    fig.savefig(folder / "candidate-loss-budget.png", dpi=200, facecolor=fig.get_facecolor())
    fig.savefig(folder / "candidate-loss-budget.svg", metadata={"Date": None})
    fig.savefig(folder / "candidate-loss-budget.pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def build(*, tables_only=False):
    raw = (ROOT / "reports/advanced_study.json").read_text(encoding="utf-8").encode()
    publication = json.loads((ROOT / "reports/advanced_publication.json").read_text(encoding="utf-8"))
    report = derive(json.loads(raw), publication, hashlib.sha256(raw).hexdigest())
    text = ["# Candidate coverage and ranking loss", "",
            "This post hoc descriptive analysis reuses the authenticated completed internal test. It fits no weights, "
            "opens no new test set and selects no model. The original scores, protocol and source snapshots remain frozen.", "",
            "![NEXUS candidate and ranking loss accounting](../docs/assets/candidate-loss-budget.png)", "",
            "[Vector SVG](../docs/assets/candidate-loss-budget.svg) · [PDF](../docs/assets/candidate-loss-budget.pdf)", "",
            "## Full model: exact loss accounting", "",
            "`1 = observed MRR + loss within covered candidates + loss from absent truth candidates`.", "",
            "| Class | MRR | Candidate oracle ceiling | Absent-truth loss | Within-coverage ranking loss |",
            "|---|---:|---:|---:|---:|"]
    for cls, v in report["models"]["full"].items():
        text.append(f"| {cls} | {v['observed_mrr']:.4f} | {v['coverage_ceiling']:.4f} | {v['coverage_loss']:.4f} | {v['within_coverage_ranking_loss']:.4f} |")
    v = report["models"]["full"]["C3"]
    text += ["", f"C3 lacks the correct candidate for **{v['truth_candidates_absent']}/250 molecules**. "
             f"Among the {v['truth_candidates_present']} covered molecules, {v['top25_hits']} reach the top 25, "
             f"{v['top1_hits']} rank first and {v['covered_outside_top25']} miss the top 25. "
             f"Perfect ordering of the same candidate lists could reach at most **{v['coverage_ceiling']:.4f} MRR**, "
             f"a hypothetical gain of {v['within_coverage_ranking_loss']:.4f}. Improving ranking alone cannot recover absent candidates.", "",
             f"Conditional C3 MRR on the covered subset is {v['conditional_mrr_given_coverage']:.4f}, and conditional hit@25 is "
             f"{v['conditional_hit25_given_coverage']:.1%}. These use only 26 molecules and must not replace the full-population "
             "MRR of 0.0373. No conditional uncertainty interval is inferred from aggregate-only evidence.", "",
             "## All prespecified C3 variants", "",
             "| Variant | Truth candidates | Top-25 hits | Top-1 hits | Full-population MRR | MRR given coverage |",
             "|---|---:|---:|---:|---:|---:|"]
    for name, classes in report["models"].items():
        v = classes["C3"]
        conditional = "Undefined" if v["conditional_mrr_given_coverage"] is None else f"{v['conditional_mrr_given_coverage']:.4f}"
        text.append(f"| {name} | {v['truth_candidates_present']} | {v['top25_hits']} | {v['top1_hits']} | {v['observed_mrr']:.4f} | {conditional} |")
    text += ["", "The larger C3 point estimate after removing generator features is already part of the prespecified "
             "secondary study. Its joint interval spans zero, so these descriptive counts do not establish superiority or justify promotion.", "",
             "## Interpretation limits", ""]
    text += ["- " + value for value in report["limitations"]]
    text += ["", "[Source study and intervals](advanced_review.md) · [Authenticated publication](advanced_publication.json) · "
             "[Derived JSON with source hashes](candidate_loss_budget.json)", ""]
    (ROOT / "reports/candidate_loss_budget.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (ROOT / "reports/candidate_loss_review.md").write_text("\n".join(text), encoding="utf-8")
    if not tables_only:
        figures(report)
    print("Verified and rendered aggregate-only candidate/ranking loss accounting")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-only", action="store_true")
    build(tables_only=parser.parse_args().tables_only)
