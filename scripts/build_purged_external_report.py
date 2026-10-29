"""Render authenticated family-purged study outcomes, including all adverse results."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.check_purged_external_publication import validate_result


def figures(result):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
        "svg.hashsalt": "nexus-casmi-purged-external-v1", "axes.spines.top": False,
        "axes.spines.right": False, "savefig.facecolor": "#0B0F19"})
    figure, axes = plt.subplots(1, 2, figsize=(14, 8), facecolor="#0B0F19")
    for axis in axes:
        axis.set_facecolor("#0B0F19"); axis.tick_params(colors="#C1CBD6")
        axis.xaxis.label.set_color("#C1CBD6")
        for side in ("left", "bottom"): axis.spines[side].set_color("#344456")
        axis.grid(axis="x", color="#263447", linewidth=.6); axis.set_axisbelow(True)
        axis.set_xlim(-.025, 1.025)
    colors = {"legacy": "#B5A4FF", "expanded": "#00E0FF"}
    coverage = result["candidate_coverage"]
    cases = [("internal", "C3", "Internal C3"), ("external", "C2", "External C2"), ("external", "C3", "External C3")]
    for index, (cohort, cls, label) in enumerate(reversed(cases)):
        for offset, policy in ((-.17, "legacy"), (.17, "expanded")):
            cell = coverage[cohort][policy][cls]; value = cell["covered"] / cell["n_molecules"]
            axes[0].barh(index + offset, value, height=.28, color=colors[policy], label=policy if index == 0 else None)
            axes[0].text(min(value + .025, .92), index + offset, f"{cell['covered']}/{cell['n_molecules']}", color="#E5EDF6", fontsize=9, va="center")
    axes[0].set_yticks(range(3), [v[2] for v in reversed(cases)])
    axes[0].set_xlabel("Truth present in candidate list · all selected molecules")
    axes[0].set_title("A / Candidate coverage", loc="left", color="#F1F6FC", pad=18)
    axes[0].legend(frameon=False, labelcolor="#C1CBD6", loc="lower right")
    names = [("mass_only", "Mass only"), ("d1_d2", "D1 + D2"), ("legacy", "Legacy budget"), ("expanded", "Expanded budget")]
    for index, (name, _) in enumerate(reversed(names)):
        for offset, cls, color in ((-.15, "C2", "#00E0FF"), (.15, "C3", "#B5A4FF")):
            cell = result["models"]["external"][name]["none"][cls]
            axes[1].hlines(index + offset, cell["ci_low"], cell["ci_high"], color=color, linewidth=2)
            axes[1].scatter(cell["mrr"], index + offset, color=color, s=32, label=cls if index == 0 else None, zorder=3)
    axes[1].set_yticks(range(4), [v[1] for v in reversed(names)])
    axes[1].set_xlabel("External MRR@25 · molecule bootstrap 95% intervals")
    axes[1].set_title("B / External ranking", loc="left", color="#F1F6FC", pad=18)
    axes[1].legend(frameon=False, labelcolor="#C1CBD6", loc="lower right")
    primary = result["primary"]
    figure.text(.07, .95, "NEXUS / CASMI26", color="#00E0FF", fontsize=11, weight="bold")
    figure.text(.07, .90, "Fresh training · family-purged external evaluation", color="#F1F6FC", fontsize=20, weight="bold")
    figure.text(.07, .86, f"300 train / 100 development / 150 internal test · {result['external_population']['n_molecules']} external molecules · eight fresh learned artifacts", color="#B7C4D3", fontsize=10)
    figure.text(.07, .10, f"Primary internal C3 coverage change: {primary['delta_coverage']:+.4f} [{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}]", color="#E5EDF6", fontsize=11)
    figure.text(.07, .06, "Selected external families removed from the full rebuilt training/reference corpus. Original-snapshot novelty feasibility did not pass.", color="#9DAEC0", fontsize=9)
    figure.text(.07, .035, "Shared independent candidate knowledge remains. Source annotations are not independently authenticated laboratory identities.", color="#9DAEC0", fontsize=9)
    figure.subplots_adjust(left=.15, right=.97, top=.77, bottom=.23, wspace=.6)
    folder = ROOT / "docs/assets"
    figure.savefig(folder / "purged-external-study.png", dpi=200)
    figure.savefig(folder / "purged-external-study.svg", metadata={"Date": None, "Creator": "NEXUS CASMI research"})
    figure.savefig(folder / "purged-external-study.pdf", metadata={"CreationDate": None, "ModDate": None, "Creator": "NEXUS CASMI research"})
    plt.close(figure)


def build(*, tables_only=False):
    raw = (ROOT / "reports/purged_external_study.json").read_bytes()
    result = json.loads(raw); validate_result(result)
    publication = json.loads((ROOT / "reports/purged_external_publication.json").read_text())
    if publication["status"] != "pass" or publication["aggregate_artifact_sha256"] != hashlib.sha256(raw).hexdigest():
        raise ValueError("Authenticated aggregate publication required before reporting")
    purge, primary, choice = result["purged_training"], result["primary"], result["training_seal"]["development_choice"]
    lines = ["![NEXUS](../site/assets/nexus-midnight.svg)", "", "# Family-purged external evaluation and C3 budgets", "",
        "This separately registered experiment evaluates freshly fitted models after removing selected external chemical families and molecular identities from the entire training/reference corpus. Original-snapshot source novelty feasibility remains unsuccessful: MassBank leaves 0/491 novel molecules and expanded GNPS leaves only 7/936. These findings are preserved; this study estimates generalization relative to a new purged training snapshot.", "",
        "[Protocol](../protocols/purged-external-v1.json) · [Execution details](../protocols/purged-external-v1-execution.json) · [Aggregate](purged_external_study.json) · [Authenticated publication](purged_external_publication.json)", "",
        "![NEXUS coverage and external ranking](../docs/assets/purged-external-study.png)", "",
        "[Vector SVG](../docs/assets/purged-external-study.svg) · [Standalone PDF](../docs/assets/purged-external-study.pdf)", "",
        "## Training/reference purge", "", "| Quantity | Count |", "|---|---:|",
        f"| Original training spectra | {purge['original_training_rows']:,} |",
        f"| Removed spectra | {purge['removed_training_rows']:,} |",
        f"| Retained training spectra | {purge['retained_training_rows']:,} |",
        f"| Removed declared training identities | {purge['removed_training_identities']:,} |",
        f"| Undefined-family labels quarantined | {purge['family_grouping_quarantine']['unassignable_family_labels_removed']:,} |",
        f"| Declared identities blocked by undefined family | {purge['family_grouping_quarantine']['unassignable_family_training_identities_removed']:,} |",
        f"| Selected external molecules | {purge['selected_external_molecules']} |",
        f"| Protected external family annotations | {purge['selected_external_families']} |",
        f"| Retained selected external families | {purge['retained_external_family_overlap']} |", "",
        "All selected external standard identities, official metric identities and exact stored spectral content have zero overlap with the purged training snapshot. Every selected external family is excluded from its full reference corpus, not only from supervised query labels. Earlier internal test-family exclusions are reconstructed before filtering and preserved. Independent COCONUT/PubChemLite candidate knowledge remains available under the declared C2 semantics.", "",
        "## Fresh fitting and development-only recommendation", "",
        "Internal query partitions contain 300 training, 100 development and 150 test molecules with no selected query/metric/family overlap. Two fresh D3 networks supply opposite-fold training features; six fixed epochs and development-only ranker stopping precede the eight-file model seal. Membership priors and audit markers are excluded from fitted features. Fixed NP-score priors are not freshly trained.", "",
        f"Development C3 coverage: legacy {choice['legacy_covered']}/100; expanded {choice['expanded_covered']}/100; gains {choice['gained']}; losses {choice['lost']}. Prespecified recommended policy: **{choice['selected_policy']}**. Both fixed budgets are reported on test regardless of that recommendation.", "",
        "## Primary C3 candidate-coverage comparison", "",
        f"Legacy: **{primary['legacy_covered']}/150**; expanded: **{primary['expanded_covered']}/150**. Paired change **{primary['delta_coverage']:+.4f}**, 95% percentile interval **[{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}]**, from 4,000 molecule resamples. Gains: {primary['gained']}; losses: {primary['lost']}. This is coverage, not achieved ranking accuracy."]
    for cohort in ("internal", "external"):
        lines += ["", f"## {cohort.title()} test metrics", "", "Every selected molecule remains in each denominator, including empty candidate lists and uncovered truths.", "",
            "| Model | Condition | Class | n | MRR@25 [95% CI] | Top 1 | Top 10 | Hit 25 |", "|---|---|---|---:|---:|---:|---:|---:|"]
        for name, conditions in result["models"][cohort].items():
            for mode, classes in conditions.items():
                for cls, v in classes.items():
                    lines.append(f"| {name} | {mode} | {cls} | {v['n']} | {v['mrr']:.4f} [{v['ci_low']:.4f}, {v['ci_high']:.4f}] | {v['top1']:.4f} | {v['top10']:.4f} | {v['recall_at_25']:.4f} |")
    lines += ["", "## Candidate coverage ceilings", "", "| Cohort | Policy | Class | Truth present | Oracle MRR ceiling | Mean / maximum candidates |", "|---|---|---|---:|---:|---:|"]
    for cohort, policies in result["candidate_coverage"].items():
        for policy, classes in policies.items():
            for cls, v in classes.items():
                lines.append(f"| {cohort} | {policy} | {cls} | {v['covered']}/{v['n_molecules']} | {v['oracle_mrr_ceiling']:.4f} | {v['mean_candidates']:.1f} / {v['maximum_candidates']} |")
    lines += ["", "Oracle ceilings assume perfect ordering of the existing lists. They are hypothetical bounds, not achieved accuracy. The three diagnostic baseline rules use original-pool candidates and cannot recover a removed C3 truth.", "",
        "## Fixed-family joint secondary uncertainty", "",
        "All 11 registered contrasts share a joint max-deviation bootstrap radius. Molecules are paired within each cohort and independently resampled across cohorts; absolute model intervals above do not replace these paired comparisons. No adverse comparison is omitted.", "",
        "| Comparison | Cohort | Difference [joint 95% CI] |", "|---|---|---:|"]
    for name, v in result["secondary_joint_intervals"]["outcomes"].items():
        lines.append(f"| {name} | {v['cohort']} | {v['delta']:+.4f} [{v['ci_low']:+.4f}, {v['ci_high']:+.4f}] |")
    lines += ["", "## Interpretation and limits", ""] + ["- " + text for text in result["limitations"]]
    lines += ["", "The historical 600/150/250 benchmark used a different training/reference corpus and population. Direct MRR comparison does not isolate a causal improvement. External spectra come from public source annotations, not an independent laboratory replication. The family-balanced source sample is not population prevalence. Unknown numeric collision energy is encoded zero, with no NCE/percent-to-eV conversion.", "",
        "No commercial service, billing, paid infrastructure, deployment, raw competition-record redistribution or competition submission is part of this project.", ""]
    (ROOT / "reports/purged_external_review.md").write_text("\n".join(lines), encoding="utf-8")
    if not tables_only: figures(result)
    print("Rendered authenticated family-purged external study with every declared outcome")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--tables-only", action="store_true")
    build(tables_only=parser.parse_args().tables_only)
