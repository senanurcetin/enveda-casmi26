"""Render verified internal-study tables and standalone NEXUS research figures."""
from __future__ import annotations

import json
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.check_advanced_publication import validate_result


def figures(result):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.hashsalt": "nexus-casmi-internal-v1", "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.facecolor": "#0B0F19"})
    figure, axes = plt.subplots(1, 2, figsize=(14, 8), facecolor="#0B0F19", gridspec_kw={"width_ratios": [1, 1.2]})
    colors = {"C1": "#FFCE85", "C2": "#00E0FF", "C3": "#B5A4FF"}
    for axis in axes:
        axis.set_facecolor("#0B0F19")
        axis.tick_params(colors="#C1CBD6", labelsize=10)
        axis.xaxis.label.set_color("#C1CBD6")
        axis.title.set_color("#F1F6FC")
        axis.spines["left"].set_color("#344456")
        axis.spines["bottom"].set_color("#344456")
        axis.grid(axis="x", color="#263447", linewidth=.6)
        axis.set_axisbelow(True)
    ax = axes[0]
    names = [("mass_only", "Mass only"), ("d1", "D1"), ("d1_d2", "D1 + D2"), ("full", "Fresh full ranker")]
    for index, (name, label) in enumerate(names):
        for offset, cls in zip((.20, 0, -.20), ("C1", "C2", "C3")):
            values = result["models"][name]["none"][cls]
            y = 3 - index + offset
            ax.hlines(y, values["ci_low"], values["ci_high"], color=colors[cls], linewidth=2)
            ax.scatter(values["mrr"], y, color=colors[cls], s=32, label=cls if index == 0 else None, zorder=3)
    ax.set_yticks(range(4), [value[1] for value in reversed(names)])
    ax.set_xlim(-.025, 1.025)
    ax.set_xlabel("Test MRR@25 · absolute molecule bootstrap 95% intervals")
    ax.set_title("A / Sealed test comparisons", loc="left", pad=20, fontsize=12, color="#F1F6FC")
    legend = ax.legend(loc="upper right", frameon=False, labelcolor="#C1CBD6")
    for text in legend.get_texts():
        text.set_color("#C1CBD6")
    ax = axes[1]
    comparisons = [("no_d1_minus_full:C2", "Remove D1 features"),
                   ("no_d2_minus_full:C2", "Remove D2 features"),
                   ("no_d3_minus_full:C2", "Remove D3 features"),
                   ("no_d4_minus_full:C2", "Remove D4 features"),
                   ("no_gen_features_minus_full:C2", "Remove generator features"),
                   ("no_np_score_minus_full:C2", "Remove NP-score feature"),
                   ("membership_only_minus_full:C2", "Membership-only diagnostic"),
                   ("full_blank_minus_none:C2", "Full model: blank fragments"),
                   ("full_intensity_shuffle_minus_none:C2", "Full model: shuffle intensity")]
    for index, (key, label) in enumerate(comparisons):
        values = result["secondary_joint_intervals"]["outcomes"][key]
        y = len(comparisons) - index - 1
        color = "#00E0FF" if index < 6 else "#B5A4FF"
        ax.hlines(y, values["ci_low"], values["ci_high"], color=color, linewidth=2)
        ax.scatter(values["delta_mrr"], y, color=color, s=32, zorder=3)
    ax.axvline(0, color="#8496AA", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(comparisons)), [value[1] for value in reversed(comparisons)])
    ax.set_xlabel("C2 ΔMRR versus fresh full model · joint 95% intervals")
    ax.set_title("B / Fixed-family secondary diagnostics", loc="left", pad=20, fontsize=11, color="#F1F6FC")
    primary = result["primary"]
    figure.text(.07, .95, "NEXUS / CASMI26", color="#00E0FF", fontsize=11, weight="bold")
    figure.text(.07, .905, "Fresh training · internal chemical-family holdout", color="#F1F6FC", fontsize=20, weight="bold")
    figure.text(.07, .865, "600 train / 150 development / 250 test queries · all learned fingerprint/ranker weights fitted from scratch", color="#B7C4D3", fontsize=10)
    figure.text(.07, .075, f"Primary C2 full − D1+D2: {primary['delta_mrr']:+.4f} [{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}]", color="#E5EDF6", fontsize=11)
    figure.text(.07, .045, "Internal query-family separation; retrieval references remain task-defined. Fixed candidate universe for feature ablations.", color="#9DAEC0", fontsize=9)
    figure.text(.07, .025, "Not external laboratory validation or a competition score. Source: verified reports/advanced_study.json.", color="#9DAEC0", fontsize=9)
    figure.subplots_adjust(left=.14, right=.97, top=.78, bottom=.20, wspace=.85)
    folder = ROOT / "docs/assets"
    figure.savefig(folder / "advanced-internal-study.png", dpi=200)
    figure.savefig(folder / "advanced-internal-study.svg", metadata={"Date": None, "Creator": "NEXUS CASMI research"})
    figure.savefig(folder / "advanced-internal-study.pdf", metadata={"CreationDate": None, "ModDate": None, "Creator": "NEXUS CASMI research"})
    plt.close(figure)


def build(*, tables_only=False):
    result = json.loads((ROOT / "reports/advanced_study.json").read_text(encoding="utf-8"))
    validate_result(result)
    protocol = result["protocol"]
    primary = result["primary"]
    lines = ["# Fresh-weight internal chemical-family benchmark", "",
             "Every fingerprint network and LambdaRank model was fitted from a fresh seeded initialization. "
             "The final test was opened only after all rankers, feature lists and learned artifact hashes were sealed. "
             "This is an internal query-family benchmark, not newly acquired external/laboratory validation or a competition score.", "",
             "![NEXUS internal benchmark and paired diagnostics](../docs/assets/advanced-internal-study.png)", "",
             "[Vector SVG](../docs/assets/advanced-internal-study.svg) · [Standalone PDF](../docs/assets/advanced-internal-study.pdf)", "",
             "## Predeclared population and fitting lineage", "",
             "| Partition | Molecules / families |", "|---|---:|"]
    for name in ("train", "dev", "test"):
        lines.append(f"| {name} | {result['partitions'][name]['n_molecules']} |")
    lines += ["", "Cross-partition selected metric/query/family overlap: **0**. One representative per generic "
              "Murcko/acyclic connectivity family defines a family-balanced population, not prevalence-weighted accuracy. "
              "The shared candidate universe and task-defined references remain available; held-out query separation is not reference-free discovery.", "",
              "| Fresh D3 fold | Training structures | Spectra | First / last training BCE |", "|---|---:|---:|---:|"]
    for fold in result["training_seal"]["fingerprint_training"]:
        lines.append(f"| {fold['fold']} | {fold['n_structures']} | {fold['n_spectra']} | "
                     f"{fold['training_bce_by_epoch'][0]:.4f} / {fold['training_bce_by_epoch'][-1]:.4f} |")
    lines += ["", "D3 has six fixed epochs, with no development/test checkpoint selection. Ranker-training queries "
              "use the opposite D3 fold. Training BCE describes fitted training data, not held-out fingerprint accuracy. "
              "Development data alone select ranker stopping iterations. Freshness concerns fingerprint and ranking weights; "
              "the RDKit NP-score is a fixed external chemical prior whose source-training overlap is not independently attested.", "",
              "## Test MRR@25", "", "| Model | C1, 95% CI | C2, 95% CI | C3, 95% CI |", "|---|---:|---:|---:|"]
    for name, conditions in result["models"].items():
        cells = [f"{v['mrr']:.4f} [{v['ci_low']:.4f}, {v['ci_high']:.4f}]"
                 for v in (conditions["none"][cls] for cls in ("C1", "C2", "C3"))]
        lines.append("| " + name + " | " + " | ".join(cells) + " |")
    mass = result["models"]["mass_only"]["none"]["C2"]["mrr"]
    baseline = result["models"]["d1_d2"]["none"]["C2"]["mrr"]
    lines += ["", f"The mass-only C2 baseline scores {mass:.4f}, compared with {baseline:.4f} for D1+D2. "
              "The primary comparison was fixed as D1+D2 before execution; its large effect should be interpreted alongside "
              "this simple-baseline weakness. No test-selected baseline replaces the declared comparison.", "",
              "## Single predeclared primary contrast", "", protocol["primary"] + ".", "",
              f"**ΔMRR {primary['delta_mrr']:+.4f}, paired 95% interval [{primary['ci_low']:+.4f}, {primary['ci_high']:+.4f}].** "
              "This estimates the difference for the selected internal population; no external-generalization or automatic model-promotion claim follows.", "",
              "## Joint secondary uncertainty", "",
              "All 27 prespecified feature/control × novelty contrasts share one approximate max-deviation bootstrap "
              "family. Negative effects mean the listed treatment scored below fresh full. All feature ablations are refitted; "
              "they share generated candidates. `no_gen_features` removes features, not the candidate generator.", "",
              "| Secondary comparison | C1 ΔMRR, joint CI | C2 ΔMRR, joint CI | C3 ΔMRR, joint CI |", "|---|---:|---:|---:|"]
    outcomes = result["secondary_joint_intervals"]["outcomes"]
    for name in [v + "_minus_full" for v in protocol["ablations"] if v != "full"] + ["full_blank_minus_none", "full_intensity_shuffle_minus_none"]:
        cells = [f"{v['delta_mrr']:+.4f} [{v['ci_low']:+.4f}, {v['ci_high']:+.4f}]"
                 for v in (outcomes[f"{name}:{cls}"] for cls in ("C1", "C2", "C3"))]
        lines.append("| " + name + " | " + " | ".join(cells) + " |")
    lines += ["", "## Candidate coverage", "", "| Class | Original-pool truth present | Pool + generation truth present |", "|---|---:|---:|"]
    for cls, counts in result["candidate_coverage"].items():
        lines.append(f"| {cls} | {counts['pool_truth_present']} / 250 | {counts['with_generation_truth_present']} / 250 |")
    lines += ["", "Every selected test query remains in the denominator, including empty candidate lists and recall misses. "
              "C3 original-pool truth coverage must be zero. Coverage is an upper bound on attainable ranking recall, not achieved MRR.", "",
              "## Descriptive C2 error slices", "",
              "Pointwise subgroup intervals are descriptive; the joint secondary guarantee above does not apply to "
              "subgroup selection. Cells below the predeclared minimum of 20 are suppressed. The test population contains "
              "141/250 pluskal_ms2 queries: family balancing does not make sources equally represented. Small timsTOF "
              "or library cells cannot establish cross-instrument or cross-library transfer.", "",
              "| Dimension | Slice | n | D1+D2 MRR | Full MRR | ΔMRR, descriptive CI |", "|---|---|---:|---:|---:|---:|"]
    for dimension, blocks in result["error_slices"].items():
        for name, block in blocks.items():
            if block["status"] == "SUPPRESSED_SMALL_CELL":
                lines.append(f"| {dimension} | {name} | {block['n']} | — | — | Suppressed |")
            else:
                v = block["classes"]["C2"]
                lines.append(f"| {dimension} | {name} | {block['n']} | {v['baseline_mrr']:.4f} | {v['full_mrr']:.4f} | "
                             f"{v['delta_mrr']:+.4f} [{v['ci_low']:+.4f}, {v['ci_high']:+.4f}] |")
    lines += ["", "## Reproducibility and limits", "",
              "[Protocol](../protocols/advanced-internal-v1.json) · [Model/input/code seals and all aggregate conditions](advanced_study.json) · "
              "[Authenticated publication readback](advanced_publication.json)", "",
              "This population, architecture and training procedure differ from historical v16. Comparing their raw scores "
              "does not isolate a model improvement. Fixed reference knowledge and the NP-score prior limit claims of wholly unseen chemistry.", ""]
    lines += ["- " + limitation for limitation in result["limitations"]]
    (ROOT / "reports/advanced_review.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not tables_only:
        figures(result)
    print("Rendered verified fresh-study review" + ("" if tables_only else " and standalone PNG/SVG/PDF figures"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-only", action="store_true", help="Verify/rebuild Markdown without importing plotting dependencies")
    build(tables_only=parser.parse_args().tables_only)
