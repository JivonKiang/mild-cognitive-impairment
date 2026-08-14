from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "03_RESULTS"
FIGURES = ROOT / "04_FIGURES"
SOURCE = FIGURES / "source_data"
PROTOCOL = ROOT / "00_PROTOCOL"
AUDIT = ROOT / "01_DATA_AUDIT"
COHORT_REGISTRY = PROTOCOL / "candidate_compatibility_registry.csv"
DATA_LEDGER = AUDIT / "data_ledger.csv"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7.2,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "axes.labelsize": 7.2,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    }
)

COLORS = {
    "ink": "#1E293B",
    "muted": "#64748B",
    "blue": "#2563EB",
    "teal": "#0F766E",
    "orange": "#D97706",
    "red": "#B91C1C",
    "gray": "#CBD5E1",
    "light": "#F1F5F9",
}


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontweight="bold", fontsize=10, color=COLORS["ink"])


def save_figure(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURES / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def load_fig1_cohort_records() -> tuple[pd.DataFrame, pd.DataFrame]:
    registry = pd.read_csv(COHORT_REGISTRY)
    ledger = pd.read_csv(DATA_LEDGER)
    required_registry = {
        "GSE63060": "ADMITTED_AS_DEVELOPMENT",
        "GSE63061": "ADMITTED_AS_DEVELOPMENT",
        "GSE140829": "ADMITTED",
        "GSE249477": "ADMITTED_AS_STRESS_TEST",
        "GSE282742": "ADMITTED_AS_SEPARATE_BRANCH",
        "GSE208623": "EXCLUDED_FROM_PRIMARY",
    }
    for accession, status in required_registry.items():
        rows = registry.loc[registry["accession"] == accession]
        if len(rows) != 1 or rows.iloc[0]["admission_status"] != status:
            raise ValueError(f"Fig. 1 registry contract failed for {accession}")
    required_ledger = {"GSE63060", "GSE63061", "GSE140829", "GSE249477", "GSE282742"}
    if not required_ledger.issubset(set(ledger["dataset"])):
        raise ValueError("Fig. 1 data ledger is missing a required cohort")
    return registry.set_index("accession", drop=False), ledger.set_index("dataset", drop=False)


def fig1_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    detail: str,
    facecolor: str,
    edgecolor: str,
) -> None:
    ax.add_patch(
        plt.Rectangle((x, y), width, height, facecolor=facecolor, edgecolor=edgecolor, lw=0.9)
    )
    ax.text(x + 0.018, y + height - 0.055, title, va="top", fontsize=7.1, fontweight="bold", color=COLORS["ink"])
    ax.text(x + 0.018, y + 0.045, detail, va="bottom", fontsize=5.7, color=COLORS["muted"])


def fig1_design() -> None:
    registry, ledger = load_fig1_cohort_records()
    development = ledger.loc[["GSE63060", "GSE63061"]]
    development_n = int(development["hc_n"].sum() + development["mci_n"].sum())
    external = ledger.loc["GSE140829"]
    external_n = int(external["hc_n"] + external["mci_n"])
    stress = ledger.loc["GSE249477"]
    stress_n = int(stress["hc_n"] + stress["mci_n"])
    progression = ledger.loc["GSE282742"]

    fig = plt.figure(figsize=(7.2, 4.85), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.50, 1.36, 0.90], hspace=0.06)
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.text(0.01, 0.85, "PUBLIC DATA / COHORT ADMISSION", fontsize=6.7, color=COLORS["muted"], fontweight="bold")
    ax.text(0.01, 0.43, "Leakage-aware transportability audit", fontsize=13.2, fontweight="bold", color=COLORS["ink"])
    ax.text(0.01, 0.08, "Primary estimand: peripheral whole-blood bulk-mRNA MCI versus HC. Counts are source-locked sample rows unless stated otherwise.", fontsize=6.1, color=COLORS["muted"])

    hero = fig.add_subplot(gs[1, 0])
    hero.set_xlim(0, 1)
    hero.set_ylim(0, 1)
    hero.axis("off")
    panel_label(hero, "a")
    hero.text(0.00, 0.98, "Primary evidence lane", fontsize=7.2, fontweight="bold", color=COLORS["ink"], va="top")
    fig1_box(
        hero, 0.01, 0.25, 0.255, 0.55,
        "Development family\nGSE63060 + GSE63061",
        f"{development_n} endpoint rows\nfold-contained selection\n{registry.loc['GSE63060', 'proposed_role']}",
        "#E8F0FE", COLORS["blue"],
    )
    fig1_box(
        hero, 0.375, 0.25, 0.245, 0.55,
        "Frozen MCI-versus-HC\npipeline",
        "endpoint + feature panel + preprocessing\nlocked before target diagnosis assessment",
        COLORS["light"], COLORS["gray"],
    )
    fig1_box(
        hero, 0.735, 0.50, 0.255, 0.35,
        "GSE140829\nprimary external audit",
        f"{external_n} HC/MCI rows\nsample-row audit; subject map unresolved",
        "#FDECEC", COLORS["red"],
    )
    fig1_box(
        hero, 0.735, 0.10, 0.255, 0.35,
        "GSE249477\ncross-platform stress test",
        f"{stress_n} HC/MCI rows\nRNA-seq; secondary transport stress test",
        "#E7F6F3", COLORS["teal"],
    )
    hero.annotate("", xy=(0.365, 0.53), xytext=(0.27, 0.53), arrowprops=dict(arrowstyle="->", lw=1.0, color=COLORS["muted"]))
    hero.annotate("", xy=(0.725, 0.66), xytext=(0.63, 0.53), arrowprops=dict(arrowstyle="->", lw=1.0, color=COLORS["muted"]))
    hero.annotate("", xy=(0.725, 0.26), xytext=(0.63, 0.53), arrowprops=dict(arrowstyle="->", lw=1.0, color=COLORS["muted"]))
    hero.text(0.00, 0.05, "Transport assessment is not positive validation when domain displacement prevents generalization.", fontsize=5.6, color=COLORS["red"])

    branches = fig.add_subplot(gs[2, 0])
    branches.set_xlim(0, 1)
    branches.set_ylim(0, 1)
    branches.axis("off")
    panel_label(branches, "b")
    branches.text(0.00, 0.95, "Separate, non-pooled evidence branches", fontsize=7.2, fontweight="bold", color=COLORS["ink"], va="top")
    fig1_box(
        branches, 0.01, 0.18, 0.30, 0.57,
        "GSE282742\nprogression association",
        f"{int(progression['unique_subjects'])} unique subjects; {int(progression['baseline_unique_subject_n'])} baseline records\nseparate progression endpoint",
        COLORS["light"], COLORS["gray"],
    )
    fig1_box(
        branches, 0.35, 0.18, 0.30, 0.57,
        "GSE208623\ncross-omics context",
        "blood-leukocyte methylation\nseparate association audit; not mRNA validation",
        COLORS["light"], COLORS["gray"],
    )
    fig1_box(
        branches, 0.69, 0.18, 0.30, 0.57,
        "Other incompatible sources",
        "PBMC, intervention, cell-specific, brain, or non-coding RNA\nexcluded from primary model performance",
        COLORS["light"], COLORS["gray"],
    )
    branches.text(0.50, 0.02, "Do not pool incompatible estimands", ha="center", fontsize=6.6, color=COLORS["red"], fontweight="bold")
    fig.suptitle("Fig. 1 | Cohort admission and evidence boundaries", x=0.02, ha="left", fontsize=8, fontweight="bold", color=COLORS["ink"])
    save_figure(fig, "Fig1_design_hierarchy")


def fig2_attack_and_transport() -> None:
    attack = pd.read_csv(RESULTS / "leakage_attack_replicates.csv")
    study = pd.read_csv(RESULTS / "study_heldout_metrics.csv")
    main = pd.read_csv(RESULTS / "metrics.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), gridspec_kw={"width_ratios": [1.0, 1.4, 1.15]}, constrained_layout=True)
    ax = axes[0]
    for idx, model in enumerate(["nested_clean", "full_scope_leakage_attack"]):
        vals = attack.loc[attack["model"] == model, "auc"]
        ax.boxplot(vals, positions=[idx], widths=0.46, patch_artist=True, boxprops=dict(facecolor=COLORS["blue"] if idx == 0 else COLORS["orange"], alpha=0.65, color="none"), medianprops=dict(color=COLORS["ink"], lw=1.1), whiskerprops=dict(color=COLORS["muted"]), capprops=dict(color=COLORS["muted"]), flierprops=dict(marker="", alpha=0))
        ax.scatter(np.full(len(vals), idx), vals, s=7, color=COLORS["ink"], alpha=0.4, zorder=3)
    ax.set_xticks([0, 1], ["clean\nfold-only", "attack\nfull-scope"])
    ax.set_ylim(0.74, 0.88)
    ax.set_ylabel("Internal AUC")
    ax.set_title("Selection scope changes\ninternal optimism", fontsize=7.5, loc="left")
    panel_label(ax, "a")
    ax.text(0.02, 0.03, "mean gap = 0.020; 20 repeated 5-fold runs", transform=ax.transAxes, fontsize=5.8, color=COLORS["muted"])
    ax = axes[1]
    rows = []
    for _, r in study[study["model"] == "transcriptome"].iterrows():
        rows.append((f"{r.train_domain} ->\n{r.test_domain}", r.auc, r.auc_ci_low, r.auc_ci_high, COLORS["teal"]))
    r = main[(main.endpoint == "stage_mci_vs_hc") & (main.model == "transcriptome") & (main.status == "EXTERNAL_VALIDATION")].iloc[0]
    rows.append(("GSE63060+61 ->\nGSE140829", r.auc, r.auc_ci_low, r.auc_ci_high, COLORS["red"]))
    for idx, (label, auc, low, high, color) in enumerate(rows):
        ax.errorbar(auc, idx, xerr=[[auc - low], [high - auc]], fmt="o", color=color, ecolor=color, elinewidth=1.2, capsize=2.5, ms=5)
        n = int(study[study["model"] == "transcriptome"].iloc[idx]["n"]) if idx < 2 else int(r["n"])
        ax.text(auc + 0.025, idx + 0.15, f"{auc:.3f} [{low:.3f}, {high:.3f}]; n={n}", fontsize=5.0, color=color, va="center")
    ax.axvline(0.5, color=COLORS["gray"], ls="--", lw=0.8)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.set_xlim(0.25, 0.95)
    ax.set_xlabel("AUC (bootstrap 95% CI)")
    ax.set_title("Study-held-out versus\nindependent transport", fontsize=7.5, loc="left")
    panel_label(ax, "b")
    ax = axes[2]
    metric_rows = main[(main.endpoint == "stage_mci_vs_hc") & (main.model == "transcriptome") & (main.status.isin(["NESTED_INTERNAL", "EXTERNAL_VALIDATION"]))].copy()
    labels = ["internal\npooled", "GSE140829\nexternal"]
    vals = [metric_rows.loc[metric_rows.status == "NESTED_INTERNAL", "auc"].iloc[0], metric_rows.loc[metric_rows.status == "EXTERNAL_VALIDATION", "auc"].iloc[0]]
    lows = [metric_rows.loc[metric_rows.status == "NESTED_INTERNAL", "auc_ci_low"].iloc[0], metric_rows.loc[metric_rows.status == "EXTERNAL_VALIDATION", "auc_ci_low"].iloc[0]]
    highs = [metric_rows.loc[metric_rows.status == "NESTED_INTERNAL", "auc_ci_high"].iloc[0], metric_rows.loc[metric_rows.status == "EXTERNAL_VALIDATION", "auc_ci_high"].iloc[0]]
    for idx, (value, low, high, color) in enumerate(zip(vals, lows, highs, [COLORS["blue"], COLORS["red"]])):
        ax.errorbar(idx, value, yerr=[[value - low], [high - value]], fmt="o", color=color, ecolor=color, elinewidth=2.4, capsize=3, ms=5)
    metric_ns = [int(metric_rows.loc[metric_rows.status == "NESTED_INTERNAL", "n"].iloc[0]), int(metric_rows.loc[metric_rows.status == "EXTERNAL_VALIDATION", "n"].iloc[0])]
    for idx, (value, low, high, color) in enumerate(zip(vals, lows, highs, [COLORS["blue"], COLORS["red"]])):
        ax.annotate(
            f"{value:.3f} [{low:.3f}, {high:.3f}]\nn={metric_ns[idx]}",
            xy=(idx, value),
            xytext=(-16, 8) if idx == 0 else (16, 8),
            textcoords="offset points",
            ha="right" if idx == 0 else "left",
            va="bottom",
            fontsize=5.0,
            color=color,
        )
    ax.set_xticks([0, 1], labels)
    ax.set_ylim(0.35, 0.95)
    ax.set_ylabel("AUC")
    ax.axhline(0.5, color=COLORS["gray"], ls="--", lw=0.8)
    ax.set_title("Primary endpoint\nexternal collapse", fontsize=7.5, loc="left")
    panel_label(ax, "c")
    ax.text(0.02, 0.03, "external OOD rate = 85.9%", transform=ax.transAxes, fontsize=6.5, color=COLORS["red"])
    save_figure(fig, "Fig2_leakage_transportability")


def fig3_external_calibration() -> None:
    pred = pd.read_csv(RESULTS / "predictions_anonymized.csv")
    main = pd.read_csv(RESULTS / "metrics.csv")
    curve = pd.read_csv(RESULTS / "selective_risk_curve.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), constrained_layout=True)
    models = ["age_plus_sex", "transcriptome"]
    colors = [COLORS["gray"], COLORS["red"]]
    ax = axes[0]
    for model, color in zip(models, colors):
        if model == "frozen_12_gene_score":
            continue
        sub = pred[
            (pred.cohort == "GSE140829")
            & (pred.endpoint == "stage_mci_vs_hc")
            & (pred.model == model)
        ]
        sub = sub.copy()
        sub["bin"] = pd.cut(sub.probability, bins=np.linspace(0, 1, 6), include_lowest=True)
        grouped = sub.groupby("bin", observed=True).agg(pred=("probability", "mean"), obs=("y", "mean"), n=("y", "size"))
        ax.plot(grouped.pred, grouped.obs, marker="o", ms=4, lw=1.2, color=color, label=model.replace("_", " "))
        for _, point in grouped.iterrows():
            ax.annotate(f"n={int(point.n)}", (point.pred, point.obs), xytext=(2, 3), textcoords="offset points", fontsize=4.5, color=color)
    ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color=COLORS["gray"])
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Mean predicted probability", ylabel="Observed MCI fraction")
    ax.set_title("Descriptive calibration\n5 display bins; ECE uses 10 bins", fontsize=7.2, loc="left")
    ax.legend(fontsize=5.7, loc="upper left")
    panel_label(ax, "a")
    ax = axes[1]
    sub = pred[
        (pred.cohort == "GSE140829")
        & (pred.endpoint == "stage_mci_vs_hc")
        & (pred.model == "transcriptome")
    ].copy()
    for truth, color, label in [(0, COLORS["blue"], "HC"), (1, COLORS["red"], "MCI")]:
        vals = sub.loc[sub.y == truth, "ood_score"].dropna()
        ax.hist(vals, bins=18, alpha=0.52, density=True, color=color, label=label)
    ax.axvline(main[(main.endpoint == "stage_mci_vs_hc") & (main.model == "transcriptome") & (main.status == "EXTERNAL_VALIDATION")]["ood_threshold_train_q95"].iloc[0], color=COLORS["ink"], ls="--", lw=0.9, label="training 95th percentile")
    ax.set_xlabel("Training-standardized distance")
    ax.set_ylabel("Density")
    ax.set_title("External domain shift\n(sample-row audit; n=383)", fontsize=7.5, loc="left")
    ax.legend(fontsize=5.3, loc="upper right")
    panel_label(ax, "b")
    ax = axes[2]
    threshold = main[(main.endpoint == "stage_mci_vs_hc") & (main.model == "transcriptome") & (main.status == "EXTERNAL_VALIDATION")]["selective_threshold_train_q20"].iloc[0]
    sub["confidence"] = np.maximum(sub.probability, 1 - sub.probability)
    curve = curve[
        (curve.cohort == "GSE140829")
        & (curve.endpoint == "stage_mci_vs_hc")
        & (curve.model == "transcriptome")
    ]
    ax.plot(
        curve.external_coverage,
        curve.external_selective_error,
        marker="o",
        ms=3.5,
        lw=1.2,
        color=COLORS["red"],
    )
    ax.axvline(float((sub.confidence >= threshold).mean()), color=COLORS["gray"], ls="--", lw=0.8, label="development q20 threshold")
    q20 = curve.loc[np.isclose(curve["development_confidence_quantile"], 0.2)]
    if not q20.empty:
        ax.text(0.03, 0.04, f"q20 accepted n={int(q20.accepted_n.iloc[0])}/{int(q20.external_n.iloc[0])}", transform=ax.transAxes, fontsize=5.5, color=COLORS["muted"])
    ax.set(xlim=(0, 1.02), ylim=(0, 1.02), xlabel="Coverage", ylabel="Selective error")
    ax.set_title("Development-threshold\nselective assessment", fontsize=7.5, loc="left")
    ax.legend(fontsize=5.2, loc="upper right")
    panel_label(ax, "c")
    save_figure(fig, "Fig3_external_uncertainty")


def fig4_mapping_and_sensitivity() -> None:
    main = pd.read_csv(RESULTS / "metrics.csv")
    stress = json.loads((RESULTS / "gse249477_stress_test.json").read_text(encoding="utf-8"))
    gse270 = json.loads((RESULTS / "gse270454_boundary_audit.json").read_text(encoding="utf-8"))
    internal = main[(main.model == "frozen_12_gene_score") & (main.status == "FROZEN_EXACT_INTERNAL") & (main.endpoint == "stage_mci_vs_hc")].iloc[0]
    partial = main[(main.model == "frozen_12_gene_score") & (main.status == "PARTIAL_EXTERNAL_SENSITIVITY") & (main.endpoint == "stage_mci_vs_hc")].iloc[0]
    panel = pd.DataFrame(
        {
            "dataset": ["Development", "GSE140829", "GSE249477", "GSE270454"],
            "mapped": [12, 10, 10, gse270["frozen_panel_mapping"]["mapped"]],
            "total": [12, 12, 12, gse270["frozen_panel_mapping"]["panel_size"]],
            "auc": [internal["auc"], partial["auc"], stress["auc_mci_vs_hc"], gse270["mci_vs_hc_row_level_auc"]],
            "low": [internal["auc_ci_low"], partial["auc_ci_low"], stress["auc_ci_low"], gse270["mci_vs_hc_auc_ci_low"]],
            "high": [internal["auc_ci_high"], partial["auc_ci_high"], stress["auc_ci_high"], gse270["mci_vs_hc_auc_ci_high"]],
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 2.85), gridspec_kw={"width_ratios": [1.0, 1.25]}, constrained_layout=True)
    ax = axes[0]
    x = np.arange(len(panel))
    ax.bar(x, panel.mapped / panel.total, color=[COLORS["teal"], COLORS["orange"], COLORS["orange"]], width=0.58)
    for idx, row in panel.iterrows():
        ax.text(idx, row.mapped / row.total + 0.04, f"{int(row.mapped)}/{int(row.total)}", ha="center", fontsize=7, fontweight="bold")
    ax.set_ylim(0, 1.2)
    ax.set_xticks(x, ["Development", "GSE140829", "GSE249477", "GSE270454"])
    ax.set_ylabel("Exact locked symbols mapped (fraction)")
    ax.set_title("Mapping is an evidence gate", fontsize=7.5, loc="left")
    ax.text(0.02, 0.03, "External bars are partial-score audits; no aliases substituted", transform=ax.transAxes, fontsize=5.1, color=COLORS["muted"])
    panel_label(ax, "a")
    ax = axes[1]
    colors = [COLORS["blue"], COLORS["orange"], COLORS["red"], COLORS["teal"]]
    for idx, row in panel.iterrows():
        ax.errorbar(idx, row.auc, yerr=[[row.auc - row.low], [row.high - row.auc]], fmt="o", color=colors[idx], ecolor=colors[idx], elinewidth=2.4, capsize=3, ms=5)
        n = [208, 383, 41, 35][idx]
        ax.text(idx, row.auc + 0.085, f"{row.auc:.3f} [{row.low:.3f}, {row.high:.3f}]\nn={n}", ha="center", fontsize=5.0, color=colors[idx])
    ax.axhline(0.5, color=COLORS["gray"], ls="--", lw=0.8)
    ax.set_ylim(0.15, 0.9)
    ax.set_xlim(-0.25, len(panel) - 0.60)
    ax.set_xticks(x, ["Development\nexact", "GSE140829\npartial", "GSE249477\npartial", "GSE270454\npartial"])
    ax.set_ylabel("MCI versus HC AUC")
    ax.set_title("Frozen score: internal versus\nexternal sensitivity audits", fontsize=7.5, loc="left")
    ax.text(0.02, 0.03, "GSE140829 = sample-row audit; GSE249477 = RNA-seq stress test; GSE270454 = verified small-cohort boundary audit", transform=ax.transAxes, fontsize=4.7, color=COLORS["muted"])
    panel_label(ax, "b")
    save_figure(fig, "Fig4_mapping_boundary")


def fig6_domain_separability() -> None:
    summary = pd.read_csv(RESULTS / "domain_separability_summary.csv").iloc[0]
    predictions = pd.read_csv(RESULTS / "domain_separability_predictions.csv")
    replicates = pd.read_csv(RESULTS / "domain_separability_replicates.csv")
    y = predictions["domain_label"].to_numpy(dtype=int)
    score = predictions["domain_probability"].to_numpy(dtype=float)
    fpr, tpr, _ = roc_curve(y, score)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), gridspec_kw={"width_ratios": [1.0, 1.25, 0.95]}, constrained_layout=True)
    ax = axes[0]
    ax.step(fpr, tpr, where="post", color=COLORS["ink"], lw=1.5)
    ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color=COLORS["gray"])
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="False-positive rate", ylabel="True-positive rate")
    ax.set_title(f"Domain classifier\nAUC = {float(summary.auc):.3f}", fontsize=7.5, loc="left")
    panel_label(ax, "a")
    ax = axes[1]
    for value, color, label in [(0, COLORS["blue"], "development\n(n=297)"), (1, COLORS["red"], "GSE140829\n(n=587)")]:
        vals = predictions.loc[predictions.domain_label == value, "domain_probability"]
        ax.hist(vals, bins=np.linspace(0, 1, 18), alpha=0.55, color=color, label=label)
    ax.set(xlim=(0, 1), xlabel="Predicted external-domain probability", ylabel="Sample-row count")
    ax.set_title("Domain separation\n(disease labels excluded from fit)", fontsize=7.2, loc="left")
    ax.legend(fontsize=5.4, loc="upper center")
    panel_label(ax, "b")
    ax = axes[2]
    ax.scatter(replicates["repeat"], replicates["auc"], s=18, color=COLORS["teal"], zorder=3)
    ax.axhline(0.5, color=COLORS["gray"], ls="--", lw=0.8)
    ax.set(xlim=(0.5, len(replicates) + 0.5), ylim=(0.45, 1.03), xlabel="Repeated CV run", ylabel="Domain AUC")
    ax.set_title("Replicate stability\n10 x 5-fold CV", fontsize=7.5, loc="left")
    panel_label(ax, "c")
    ax.text(0.03, 0.04, "disease labels: row alignment / post hoc strata only", transform=ax.transAxes, fontsize=4.8, color=COLORS["muted"])
    ax.text(0.03, 0.10, f"permutation p = {float(summary.permutation_p):.5f}", transform=ax.transAxes, fontsize=5.8, color=COLORS["muted"])
    save_figure(fig, "Fig6_domain_separability")


def fig7_cross_cohort_effect_heterogeneity() -> None:
    """Visualize aggregate effect concordance without introducing pooled inference."""
    pairwise = pd.read_csv(RESULTS / "cross_cohort_effect_pairwise_summary.csv")
    heterogeneity = pd.read_csv(RESULTS / "cross_cohort_effect_heterogeneity.csv")
    audit = json.loads((RESULTS / "cross_cohort_effect_audit.json").read_text(encoding="utf-8"))
    cohorts = ["GSE63060", "GSE63061", "GSE140829", "GSE249477"]
    short = {"GSE63060": "63060", "GSE63061": "63061", "GSE140829": "140829", "GSE249477": "249477"}
    correlation = np.eye(len(cohorts), dtype=float)
    concordance = np.ones((len(cohorts), len(cohorts)), dtype=float)
    count = np.full((len(cohorts), len(cohorts)), 11575, dtype=int)
    for _, row in pairwise.query("comparison_type == 'cohort_pair_descriptive'").iterrows():
        left = cohorts.index(row.cohort_a)
        right = cohorts.index(row.cohort_b)
        correlation[left, right] = correlation[right, left] = float(row.spearman_effect_correlation)
        concordance[left, right] = concordance[right, left] = float(row.sign_concordance_fraction)
        count[left, right] = count[right, left] = int(row.estimable_gene_count)

    fig = plt.figure(figsize=(7.2, 3.0), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25])
    ax = fig.add_subplot(grid[0, 0])
    image = ax.imshow(correlation, vmin=-0.25, vmax=0.75, cmap="RdBu_r")
    for row in range(len(cohorts)):
        for column in range(len(cohorts)):
            label = "1.00" if row == column else f"{correlation[row, column]:.2f}"
            ax.text(column, row, label, ha="center", va="center", fontsize=6.0, color="white" if abs(correlation[row, column]) > 0.48 else COLORS["ink"])
    ax.set_xticks(range(len(cohorts)), [short[cohort] for cohort in cohorts], rotation=45, ha="right")
    ax.set_yticks(range(len(cohorts)), [short[cohort] for cohort in cohorts])
    ax.set_title("MCI-HC effect correlation\n(Spearman)", fontsize=7.5, loc="left")
    ax.text(0.02, -0.38, "All are cohort-wise effects; no pooled estimate.", transform=ax.transAxes, fontsize=4.8, color=COLORS["muted"])
    panel_label(ax, "a")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.74, pad=0.03)
    colorbar.ax.tick_params(labelsize=5.5)

    ax = fig.add_subplot(grid[0, 1])
    image = ax.imshow(concordance, vmin=0.40, vmax=0.75, cmap="YlGnBu")
    for row in range(len(cohorts)):
        for column in range(len(cohorts)):
            if row == column:
                label = "1.00"
            else:
                label = f"{concordance[row, column]:.2f}\n(n={count[row, column] // 1000}k)"
            ax.text(column, row, label, ha="center", va="center", fontsize=5.3, color="white" if concordance[row, column] > 0.63 else COLORS["ink"])
    ax.set_xticks(range(len(cohorts)), [short[cohort] for cohort in cohorts], rotation=45, ha="right")
    ax.set_yticks(range(len(cohorts)), [short[cohort] for cohort in cohorts])
    ax.set_title("Direction concordance\n(nonzero effects)", fontsize=7.5, loc="left")
    ax.text(0.02, -0.38, "GSE249477 pairs: 10,923 estimable genes.", transform=ax.transAxes, fontsize=4.8, color=COLORS["muted"])
    panel_label(ax, "b")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.74, pad=0.03)
    colorbar.ax.tick_params(labelsize=5.5)

    ax = fig.add_subplot(grid[0, 2])
    i_squared = heterogeneity.loc[heterogeneity["status"] == "estimable_two_study_descriptive", "i_squared_percent"].dropna()
    bins = np.arange(0, 105, 10)
    ax.hist(i_squared, bins=bins, color=COLORS["teal"], edgecolor="white", linewidth=0.5)
    ax.axvline(float(i_squared.median()), color=COLORS["ink"], lw=1.0, ls="--")
    ax.set(xlim=(0, 100), xlabel="Two-study I-squared (%)", ylabel="Gene count")
    ax.set_title("AddNeuroMed within-family\nheterogeneity only", fontsize=7.5, loc="left")
    ax.text(0.58, 0.91, f"median = {float(i_squared.median()):.1f}%\n75th percentile = {float(i_squared.quantile(0.75)):.1f}%", transform=ax.transAxes, fontsize=5.4, va="top", color=COLORS["ink"])
    ax.text(0.03, -0.27, "Descriptive GSE63060/61 Q/I-squared; no P value.", transform=ax.transAxes, fontsize=4.8, color=COLORS["muted"])
    panel_label(ax, "c")
    fig.text(0.01, 0.005, "GSE140829 is a sample-row assessment with subject independence unresolved; it is excluded from pooled inference.", fontsize=5.1, color=COLORS["red"])
    if audit["gse140829_excluded_from_pooled_inference"] is not True:
        raise ValueError("Fig7 requires the GSE140829 pooled-inference exclusion")
    save_figure(fig, "Fig7_cross_cohort_effect_heterogeneity")


def figs1_gse208623_cross_omics_concordance() -> None:
    """Display aggregate, non-predictive methylation/transcriptome concordance."""
    cpg = pd.read_csv(RESULTS / "gse208623_cpg_associations.csv")
    genes = pd.read_csv(RESULTS / "gse208623_gene_concordance.csv")
    audit = json.loads((RESULTS / "gse208623_cross_omics_source_audit.json").read_text(encoding="utf-8"))
    concordance = audit["concordance"]
    fdr_significant = int((cpg["p_value_bh"] < 0.05).sum())

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), gridspec_kw={"width_ratios": [0.92, 1.35, 1.05]}, constrained_layout=True)
    ax = axes[0]
    categories = ["all\nCpGs", "FDR < 0.05\nCpGs", "unambiguous\ngenes"]
    values = [len(cpg), fdr_significant, int(audit["unambiguous_gene_count"])]
    colors = [COLORS["gray"], COLORS["orange"], COLORS["teal"]]
    log_floor = 1_000
    ax.bar(
        range(3),
        [value - log_floor for value in values],
        bottom=log_floor,
        color=colors,
        width=0.62,
    )
    for index, value in enumerate(values):
        ax.text(index, value * 1.035, f"{value:,}", ha="center", va="bottom", fontsize=6.2, fontweight="bold", color=COLORS["ink"])
    ax.set_yscale("log")
    ax.set_ylim(log_floor, 1_500_000)
    ax.set_xticks(range(3), categories)
    ax.set_ylabel("Count (log scale)")
    ax.set_title("Sex-adjusted MCI-HC\nCpG association accounting", fontsize=7.5, loc="left")
    ax.text(0.02, 0.07, "40 rows (20 MCI / 20 HC); AD excluded", transform=ax.transAxes, fontsize=4.7, color=COLORS["muted"])
    panel_label(ax, "a")

    ax = axes[1]
    ax.scatter(
        genes["transcriptome_hedges_g_mci_minus_hc"],
        genes["methylation_z"],
        s=3.2,
        alpha=0.20,
        color=COLORS["teal"],
        linewidths=0,
        rasterized=True,
    )
    ax.axhline(0.0, color=COLORS["gray"], lw=0.7, zorder=0)
    ax.axvline(0.0, color=COLORS["gray"], lw=0.7, zorder=0)
    ax.set_xlabel("Locked transcriptome Hedges' g (MCI - HC)")
    ax.set_ylabel("Aggregate methylation signed Z")
    ax.set_title("Gene-level cross-omics\nassociation", fontsize=7.5, loc="left")
    ax.text(0.03, 0.96, f"n = {int(concordance['overlap_gene_count']):,}\nSpearman rho = {float(concordance['spearman_correlation']):.3f}", transform=ax.transAxes, fontsize=6.0, va="top", color=COLORS["ink"])
    panel_label(ax, "b")

    ax = axes[2]
    observed_same = float(concordance["same_sign_fraction"])
    observed_opposite = float(concordance["opposite_sign_fraction"])
    null_mean = float(concordance["gene_sign_permutation_null_mean"])
    ax.bar([0, 1], [observed_same, observed_opposite], color=[COLORS["teal"], COLORS["orange"]], width=0.58)
    ax.axhline(null_mean, color=COLORS["ink"], lw=1.0, ls="--")
    for index, value in enumerate([observed_same, observed_opposite]):
        ax.text(index, value + 0.018, f"{value:.3f}", ha="center", fontsize=6.3, fontweight="bold")
    ax.set_ylim(0.0, 0.72)
    ax.set_xticks([0, 1], ["same sign", "opposite sign"])
    ax.set_ylabel("Fraction of sign-eligible genes")
    ax.set_title("Directional summary\nwith fixed null", fontsize=7.5, loc="left")
    ax.text(
        0.03,
        0.95,
        f"fixed-null mean = {null_mean:.3f}\n10,000 permutations; P = {float(concordance['gene_sign_permutation_p_value']):.5f}",
        transform=ax.transAxes,
        fontsize=5.0,
        va="top",
        color=COLORS["muted"],
    )
    panel_label(ax, "c")
    fig.suptitle("Fig. S1 | Source-locked cross-omics association audit", x=0.02, ha="left", fontsize=8, fontweight="bold", color=COLORS["ink"])
    save_figure(fig, "FigS1_gse208623_cross_omics_concordance")


def fig5_domain_invariant_sensitivity() -> None:
    rank = pd.read_csv(RESULTS / "domain_invariant_sensitivity.csv")
    gse249 = pd.read_csv(RESULTS / "gse249477_rank_transport_metrics.csv")
    tsp = pd.read_csv(RESULTS / "tsp_transport_metrics.csv")
    stable = pd.read_csv(RESULTS / "cross_study_stable_metrics.csv")
    coral = pd.read_csv(RESULTS / "transductive_coral_metrics.csv")
    coral_audit = json.loads((RESULTS / "transductive_coral_audit.json").read_text(encoding="utf-8"))
    stable_null = pd.read_csv(RESULTS / "cross_study_stable_permutation_null.csv")
    stable_audit = json.loads((RESULTS / "cross_study_stable_audit.json").read_text(encoding="utf-8"))
    conformal = pd.read_csv(RESULTS / "split_conformal_summary.csv").iloc[0]
    main = pd.read_csv(RESULTS / "metrics.csv")
    raw = main[
        (main.endpoint == "stage_mci_vs_hc")
        & (main.model == "transcriptome")
        & (main.status == "EXTERNAL_VALIDATION")
    ].iloc[0]
    fig = plt.figure(figsize=(8.3, 4.35), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.94])
    top = gs[0, 0].subgridspec(1, 5, width_ratios=[2.35, 0.72, 0.72, 0.72, 1.25])
    axes = [fig.add_subplot(top[0, idx]) for idx in range(5)]
    ax = axes[0]
    gse249_row = gse249.query("status == 'SECONDARY_INDEPENDENT_TRANSPORT_AUDIT'").iloc[0]
    tsp_140 = tsp.query("status == 'PRIMARY_EXTERNAL_TRANSPORT_AUDIT'").iloc[0]
    tsp_249 = tsp.query("status == 'SECONDARY_INDEPENDENT_CROSS_PLATFORM_TRANSPORT_AUDIT'").iloc[0]
    stable_140 = stable.query("cohort == 'GSE140829' and model == 'cross_study_stable_ridge'").iloc[0]
    stable_249 = stable.query("cohort == 'GSE249477' and model == 'cross_study_stable_ridge'").iloc[0]
    stable_60 = stable.query("cohort == 'GSE63060' and model == 'cross_study_stable_ridge'").iloc[0]
    stable_61 = stable.query("cohort == 'GSE63061' and model == 'cross_study_stable_ridge'").iloc[0]
    coral_140 = coral.query("cohort == 'GSE140829' and model == 'transductive_coral_pca'").iloc[0]
    coral_249 = coral.query("cohort == 'GSE249477' and model == 'transductive_coral_pca'").iloc[0]
    aucs = [
        float(raw.auc),
        float(rank.query("status == 'EXTERNAL_SENSITIVITY'").auc.iloc[0]),
        float(tsp_140.auc),
        float(stable_140.auc),
        float(coral_140.auc),
        float(gse249_row.auc),
        float(tsp_249.auc),
        float(stable_249.auc),
        float(coral_249.auc),
    ]
    lows = [
        float(raw.auc_ci_low),
        float(rank.query("status == 'EXTERNAL_SENSITIVITY'").auc_ci_low.iloc[0]),
        float(tsp_140.auc_ci_low),
        float(stable_140.auc_ci_low),
        float(coral_140.auc_ci_low),
        float(gse249_row.auc_ci_low),
        float(tsp_249.auc_ci_low),
        float(stable_249.auc_ci_low),
        float(coral_249.auc_ci_low),
    ]
    highs = [
        float(raw.auc_ci_high),
        float(rank.query("status == 'EXTERNAL_SENSITIVITY'").auc_ci_high.iloc[0]),
        float(tsp_140.auc_ci_high),
        float(stable_140.auc_ci_high),
        float(coral_140.auc_ci_high),
        float(gse249_row.auc_ci_high),
        float(tsp_249.auc_ci_high),
        float(stable_249.auc_ci_high),
        float(coral_249.auc_ci_high),
    ]
    colors = [COLORS["red"], COLORS["teal"], COLORS["blue"], COLORS["orange"], COLORS["ink"], COLORS["teal"], COLORS["blue"], COLORS["orange"], COLORS["ink"]]
    forest_labels = ["GSE140829 raw", "GSE140829 rank/z", "GSE140829 TSP", "GSE140829 stable", "GSE140829 CORAL*", "GSE249477 rank/z", "GSE249477 TSP", "GSE249477 stable", "GSE249477 CORAL*"]
    for idx, (auc, low, high, color) in enumerate(zip(aucs, lows, highs, colors)):
        ax.errorbar(auc, idx, xerr=[[auc - low], [high - auc]], fmt="o", color=color, ecolor=color, elinewidth=1.5, capsize=2.5, ms=4.5)
    ax.axvline(0.5, color=COLORS["gray"], ls="--", lw=0.8)
    ax.set_xlim(0.25, 0.80)
    ax.set_yticks(range(len(forest_labels)), forest_labels)
    ax.invert_yaxis()
    ax.set_xlabel("AUC (bootstrap 95% CI)")
    ax.set_title("No target-label-free method\nrestores transport", fontsize=7.5, loc="left")
    panel_label(ax, "a")
    # Keep discrimination and calibration on independent scales. A common
    # y-axis would make the three metrics look clinically comparable.
    metrics = ["AUC", "ECE", "Brier"]
    raw_vals = [float(raw.auc), float(raw.ece), float(raw.brier)]
    rank_row = rank.query("status == 'EXTERNAL_SENSITIVITY'").iloc[0]
    rank_vals = [float(rank_row.auc), float(rank_row.ece), float(rank_row.brier)]
    coral_vals = [float(coral_140.auc), float(coral_140.ece), float(coral_140.brier)]
    width = 0.24
    metric_limits = [(0.35, 0.75), (0.0, 0.40), (0.0, 0.40)]
    metric_axes = axes[1:4]
    for idx, (metric, raw_value, rank_value, coral_value, limits, metric_ax) in enumerate(
        zip(metrics, raw_vals, rank_vals, coral_vals, metric_limits, metric_axes)
    ):
        metric_ax.bar([-width], [raw_value], width, color=COLORS["red"], label="raw")
        metric_ax.bar([0], [rank_value], width, color=COLORS["teal"], label="rank/z")
        metric_ax.bar([width], [coral_value], width, color=COLORS["ink"], label="CORAL*")
        metric_ax.set_ylim(*limits)
        metric_ax.set_xticks([0], [metric])
        metric_ax.set_ylabel(metric, labelpad=1)
        metric_ax.set_title(metric, fontsize=7.0, loc="left")
        metric_ax.tick_params(axis="x", pad=2)
        if idx == 0:
            panel_label(metric_ax, "b")
            metric_ax.legend(fontsize=4.6, loc="upper left", handlelength=0.9)
    ax = axes[4]
    coverage = float(conformal.external_coverage_mean)
    ambiguous = float(conformal.ambiguous_rate_mean)
    singleton_error = float(conformal.singleton_error_mean)
    values = [coverage, ambiguous, singleton_error]
    labels = ["coverage", "ambiguous sets", "singleton error"]
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=[COLORS["orange"], COLORS["gray"], COLORS["red"]], height=0.56)
    ax.errorbar([coverage], [0], xerr=[[coverage - float(conformal.external_coverage_q025)], [float(conformal.external_coverage_q975) - coverage]], fmt="none", ecolor=COLORS["ink"], elinewidth=1.0, capsize=2.5, zorder=4)
    ax.axvline(0.9, color=COLORS["gray"], ls="--", lw=0.8)
    ax.set_xlim(0, 1.02)
    ax.set_yticks(y, labels)
    ax.set_xlabel("External rate")
    ax.invert_yaxis()
    ax.set_title("Conformal sets retain\ndomain-shift failure", fontsize=7.5, loc="left")
    ax.text(0.03, 0.03, "n=383; 20 repeats; whisker = 2.5th-97.5th percentile", transform=ax.transAxes, fontsize=4.7, color=COLORS["muted"])
    for bar, value in zip(bars, values):
        ax.text(value + 0.025, bar.get_y() + bar.get_height() / 2, f"{value:.2f}", va="center", fontsize=6.5, fontweight="bold")
    panel_label(ax, "c")
    bottom = gs[1, 0].subgridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.30)
    note = fig.add_subplot(bottom[0, 0])
    note.axis("off")
    note.add_patch(plt.Rectangle((0.0, 0.02), 1.0, 0.94, fc=COLORS["light"], ec=COLORS["gray"], lw=0.7))
    note.text(
        0.035,
        0.78,
        "Locked transport stress tests",
        fontsize=7.3,
        fontweight="bold",
        color=COLORS["ink"],
    )
    note.text(
        0.035,
        0.59,
        f"TSP: 140829 {float(tsp_140.auc):.3f}; 249477 {float(tsp_249.auc):.3f}. "
        f"Stable ridge: 140829 {float(stable_140.auc):.3f}; 249477 {float(stable_249.auc):.3f}.",
        fontsize=5.9,
        color=COLORS["blue"],
    )
    note.text(
        0.035,
        0.43,
        "TSP, stable ridge and CORAL scored target rows before endpoint assessment.",
        fontsize=5.45,
        color=COLORS["muted"],
    )
    note.text(
        0.035,
        0.27,
        f"CORAL*: 140829 {float(coral_140.auc):.3f}; 249477 {float(coral_249.auc):.3f}; target labels used only after lock.",
        fontsize=5.45,
        color=COLORS["red"],
    )
    note.text(
        0.035,
        0.11,
        "*batch-level transductive stress test, not single-sample deployment or clinical validation.",
        fontsize=5.35,
        color=COLORS["red"],
    )
    null_ax = fig.add_subplot(bottom[0, 1])
    null_values = stable_null["development_oof_auc"].astype(float).to_numpy()
    observed_auc = float(stable_audit["development_oof_auc"])
    empirical_p = float(stable_audit["development_label_permutation_p"])
    null_ax.hist(null_values, bins=np.linspace(0.30, 0.85, 12), color=COLORS["gray"], edgecolor="white", linewidth=0.55)
    null_ax.axvline(observed_auc, color=COLORS["red"], lw=1.8)
    null_ax.set_xlim(0.30, 0.85)
    null_ax.set_xlabel("Synchronized development OOF AUC")
    null_ax.set_ylabel("Permutations")
    null_ax.set_title("Observed development signal\nexceeds label null", fontsize=7.5, loc="left")
    null_ax.text(0.04, 0.88, f"observed = {observed_auc:.3f}\n100 permutations\nempirical p = {empirical_p:.4f}", transform=null_ax.transAxes, fontsize=5.9, va="top", color=COLORS["ink"])
    null_ax.text(0.04, 0.04, "development-family only;\nnot external validation", transform=null_ax.transAxes, fontsize=4.9, color=COLORS["red"])
    panel_label(null_ax, "d")
    save_figure(fig, "Fig5_domain_invariant_sensitivity")


def copy_source_data() -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    for name in [
        "metrics.csv",
        "study_heldout_metrics.csv",
        "leakage_attack_replicates.csv",
        "leakage_attack_summary.csv",
        "gse249477_stress_test.json",
        "gse249477_rank_transport_metrics.csv",
        "gse249477_rank_transport_audit.json",
        "tsp_transport_metrics.csv",
        "tsp_transport_selected_pairs.csv",
        "tsp_transport_audit.json",
        "cross_study_stable_metrics.csv",
        "cross_study_stable_predictions_anonymized.csv",
        "cross_study_stable_selected_genes.csv",
        "cross_study_stable_selection_events.csv",
        "cross_study_stable_permutation_null.csv",
        "cross_study_stable_audit.json",
        "transductive_coral_metrics.csv",
        "transductive_coral_predictions_anonymized.csv",
        "transductive_coral_audit.json",
        "cross_cohort_effect_pairwise_summary.csv",
        "cross_cohort_effect_heterogeneity.csv",
        "cross_cohort_effect_audit.json",
        "data_ledger.csv",
        "domain_invariant_sensitivity.csv",
        "split_conformal_summary.csv",
        "ood_distance_sensitivity.csv",
        "ood_distance_audit.json",
        "domain_separability_summary.csv",
        "domain_separability_replicates.csv",
        "domain_separability_strata.csv",
        "domain_separability_predictions.csv",
        "domain_separability_selected_features.csv",
        "domain_separability_audit.json",
        "panel_provenance_ledger.csv",
        "panel_provenance_audit.json",
        "gse270454_boundary_audit.json",
        "gse270454_join_audit.json",
        "gse270454_boundary_summary.csv",
        "gse108215_paired_summary.csv",
        "gse108215_paired_audit.json",
        "gse153881_label_conflict_audit.json",
        "predictions_anonymized.csv",
        "selective_risk_curve.csv",
        "gse208623_cpg_associations.csv",
        "gse208623_gene_concordance.csv",
        "gse208623_cross_omics_source_audit.json",
    ]:
        source = RESULTS / name if name != "data_ledger.csv" else ROOT / "01_DATA_AUDIT" / name
        if source.exists():
            shutil.copy2(source, SOURCE / name)
    shutil.copy2(COHORT_REGISTRY, SOURCE / COHORT_REGISTRY.name)


def write_panel_traceability() -> None:
    _, ledger = load_fig1_cohort_records()
    development = ledger.loc[["GSE63060", "GSE63061"]]
    development_n = int(development["hc_n"].sum() + development["mci_n"].sum())
    external = ledger.loc["GSE140829"]
    stress = ledger.loc["GSE249477"]
    progression = ledger.loc["GSE282742"]
    traceability = pd.DataFrame(
        [
            {
                "figure": "Fig1",
                "panel": "whole figure",
                "claim": "cohort admission separates the primary MCI-HC transport audit from incompatible endpoint, specimen, and assay branches",
                "analysis_unit": "sample rows; progression also reports unique subjects",
                "n_or_scope": f"AddNeuroMed primary endpoint={development_n}; GSE140829 MCI-HC={int(external['hc_n'] + external['mci_n'])}; GSE249477 MCI-HC={int(stress['hc_n'] + stress['mci_n'])}; GSE282742={int(progression['unique_subjects'])} unique subjects",
                "statistic": "descriptive counts",
                "source_data": "00_PROTOCOL/candidate_compatibility_registry.csv; 01_DATA_AUDIT/data_ledger.csv",
                "script": "02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "design context; no pooling of incompatible estimands; no validation, mechanism, or causal claim",
            },
            {
                "figure": "Fig2",
                "panel": "a",
                "claim": "label-aware selection changes internal optimism",
                "analysis_unit": "development sample rows; 20 repeated 5-fold runs",
                "n_or_scope": "n=208 primary endpoint",
                "statistic": "AUC distributions; mean optimism",
                "source_data": "03_RESULTS/leakage_attack_replicates.csv; 03_RESULTS/leakage_attack_summary.csv",
                "script": "02_ANALYSIS/run_leakage_attack.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "one leakage route; not all bias",
            },
            {
                "figure": "Fig2",
                "panel": "b-c",
                "claim": "same-family transfer exceeds external transport",
                "analysis_unit": "sample rows; GSE140829 subject independence unresolved",
                "n_or_scope": "study-held-out n=102/106; external n=383",
                "statistic": "AUC with bootstrap 95% CI",
                "source_data": "03_RESULTS/study_heldout_metrics.csv; 03_RESULTS/metrics.csv",
                "script": "02_ANALYSIS/run_study_heldout.py; 02_ANALYSIS/run_pairwise_external.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "transport audit; not population-wide generalization",
            },
            {
                "figure": "Fig3",
                "panel": "a-c",
                "claim": "external calibration and OOD uncertainty remain poor",
                "analysis_unit": "GSE140829 sample rows",
                "n_or_scope": "n=383; HC=249; MCI=134",
                "statistic": "5-bin display calibration; ECE=10 fixed bins; q95 OOD; selective risk",
                "source_data": "03_RESULTS/predictions_anonymized.csv; 03_RESULTS/metrics.csv; 03_RESULTS/selective_risk_curve.csv",
                "script": "02_ANALYSIS/run_pairwise_external.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "descriptive and diagnostic; no deployment claim",
            },
            {
                "figure": "Fig4",
                "panel": "a-b",
                "claim": "exact-symbol mapping bounds the frozen comparator",
                "analysis_unit": "sample rows; partial external score only",
                "n_or_scope": "development=12/12; GSE140829=10/12; GSE249477=10/12; GSE270454=10/12; GSE270454 join=45/45 unique",
                "statistic": "mapping fraction; partial-score AUC with bootstrap 95% CI",
                "source_data": "03_RESULTS/metrics.csv; 03_RESULTS/gse249477_stress_test.json; 03_RESULTS/gse270454_boundary_audit.json; 03_RESULTS/gse270454_join_audit.json",
                "script": "02_ANALYSIS/run_pairwise_external.py; 02_ANALYSIS/run_gse249477_stress_test.py; 02_ANALYSIS/audit_gse270454_join.py; 02_ANALYSIS/run_gse270454_locked12_joined.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "partial sensitivity; not complete-panel replication",
            },
            {
                "figure": "Fig5",
                "panel": "a",
                "claim": "rank/z, development-only relative-expression TSP, cross-study stable signal and target-unlabeled CORAL audits do not restore external transport",
                "analysis_unit": "GSE140829 sample rows; GSE249477 sample rows",
                "n_or_scope": "GSE140829 n=383; GSE249477 n=41 (HC=21; MCI=20)",
                "statistic": "AUC with bootstrap 95% CI",
                "source_data": "03_RESULTS/metrics.csv; 03_RESULTS/domain_invariant_sensitivity.csv; 03_RESULTS/gse249477_rank_transport_metrics.csv; 03_RESULTS/tsp_transport_metrics.csv; 03_RESULTS/cross_study_stable_metrics.csv; 03_RESULTS/cross_study_stable_audit.json; 03_RESULTS/transductive_coral_metrics.csv; 03_RESULTS/transductive_coral_audit.json",
                "script": "02_ANALYSIS/run_pairwise_external.py; 02_ANALYSIS/run_domain_invariant_sensitivity.py; 02_ANALYSIS/run_gse249477_rank_transport.py; 02_ANALYSIS/run_tsp_relative_expression_transport.py; 02_ANALYSIS/run_cross_study_stable_signal.py; 02_ANALYSIS/run_transductive_coral_audit.py; 02_ANALYSIS/transductive_alignment_core.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "GSE249477 is a secondary independent cross-platform transport audit. Cross-study stable ridge is development-family consistency, not independent validation. No target cohort supplies diagnoses to select, fit, prediction-calibrate, threshold or adapt model parameters. CORAL is a batch-level transductive stress test using target expression without target diagnoses before predictions lock; it is not single-sample deployment, target calibration or evidence that adaptation restores clinical transport. Target labels support only descriptive post-lock performance and calibration assessment.",
            },
            {
                "figure": "Fig5",
                "panel": "b-d and audit note",
                "claim": "ranking/TSP sensitivities do not restore calibrated coverage; the stable-signal development OOF AUC exceeds a cohort-preserving label-permutation null but remains a transport stress test",
                "analysis_unit": "GSE140829 sample rows; 20 conformal repeats; GSE249477 sample rows",
                "n_or_scope": "GSE140829 n=383; GSE249477 n=41",
                "statistic": "metric-specific AUC/ECE/Brier; coverage quantiles; 100 cohort-preserving label permutations of synchronized development OOF AUC; GSE249477 assessment-only Brier/ECE",
                "source_data": "03_RESULTS/domain_invariant_sensitivity.csv; 03_RESULTS/split_conformal_summary.csv; 03_RESULTS/gse249477_rank_transport_metrics.csv; 03_RESULTS/gse249477_rank_transport_audit.json; 03_RESULTS/tsp_transport_metrics.csv; 03_RESULTS/tsp_transport_audit.json; 03_RESULTS/cross_study_stable_metrics.csv; 03_RESULTS/cross_study_stable_permutation_null.csv",
                "script": "02_ANALYSIS/run_domain_invariant_sensitivity.py; 02_ANALYSIS/run_split_conformal.py; 02_ANALYSIS/run_gse249477_rank_transport.py; 02_ANALYSIS/run_tsp_relative_expression_transport.py; 02_ANALYSIS/run_cross_study_stable_signal.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "stress tests; GSE249477 labels are assessment-only and no target prediction calibration is fit. After prediction lock, target labels support only descriptive performance and calibration assessment. The stable-signal permutation null quantifies development-family discrimination only; it does not establish external transport, clinical utility, or deployment validity.",
            },
            {
                "figure": "Fig6",
                "panel": "a-c",
                "claim": "development and external expression domains are strongly separable",
                "analysis_unit": "sample rows; subject independence unresolved externally",
                "n_or_scope": "development=297; external=587; 10 x 5-fold CV; 5,000 permutations",
                "statistic": "domain AUC; balanced accuracy; permutation P",
                "source_data": "03_RESULTS/domain_separability_summary.csv; 03_RESULTS/domain_separability_replicates.csv; 03_RESULTS/domain_separability_predictions.csv; 03_RESULTS/domain_separability_strata.csv; 03_RESULTS/domain_separability_audit.json",
                "script": "02_ANALYSIS/run_domain_separability_audit.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "domain displacement; not biological cause or disease prediction",
            },
            {
                "figure": "Fig7",
                "panel": "a-c",
                "claim": "within-family MCI-HC gene effects are more concordant than the cross-platform GSE249477 comparisons",
                "analysis_unit": "cohort-wise aggregate gene effects; no sample identifiers exported",
                "n_or_scope": "11,575 fixed shared symbols; GSE249477 comparisons use 10,923 estimable genes",
                "statistic": "Spearman effect correlation; nonzero sign-concordance; AddNeuroMed-only descriptive Cochran Q/I-squared",
                "source_data": "03_RESULTS/cross_cohort_effect_pairwise_summary.csv; 03_RESULTS/cross_cohort_effect_heterogeneity.csv; 03_RESULTS/cross_cohort_effect_audit.json",
                "script": "02_ANALYSIS/run_cross_cohort_effect_heterogeneity.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "descriptive cohort-level effect audit; GSE140829 subject independence unresolved and excluded from pooled inference; no pooled four-cohort effect, heterogeneity P value, biological mechanism, or clinical utility claim",
            },
            {
                "figure": "FigS1",
                "panel": "a-c",
                "claim": "the source-locked GSE208623 methylation association shows a descriptive gene-level relationship with the frozen AddNeuroMed transcriptome effect reference",
                "analysis_unit": "GSE208623 sample rows for CpG association; aggregate gene statistics for cross-omics summary",
                "n_or_scope": "MCI=20; HC=20; AD=20 excluded before modeling; 738,242 CpGs; 10,422 overlap genes",
                "statistic": "sex-adjusted M-value OLS with BH FDR; descriptive signed-Z versus Hedges' g Spearman rho; fixed gene-sign permutation sensitivity",
                "source_data": "03_RESULTS/gse208623_cpg_associations.csv; 03_RESULTS/gse208623_gene_concordance.csv; 03_RESULTS/gse208623_cross_omics_source_audit.json",
                "script": "02_ANALYSIS/run_gse208623_cross_omics_audit.py; 02_ANALYSIS/make_figure_pack.py",
                "evidence_boundary": "source-locked, cohort-level orthogonal cross-omics association audit; not prediction, model validation, replication, mediation, cell-composition proof, mechanism, MR, or causal inference. Age and blood-cell composition are not adjusted.",
            },
        ]
    )
    traceability.to_csv(FIGURES / "panel_traceability.csv", index=False)


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig1_design()
    fig2_attack_and_transport()
    fig3_external_calibration()
    fig4_mapping_and_sensitivity()
    fig5_domain_invariant_sensitivity()
    fig6_domain_separability()
    fig7_cross_cohort_effect_heterogeneity()
    figs1_gse208623_cross_omics_concordance()
    copy_source_data()
    write_panel_traceability()
    manifest = []
    for path in sorted(FIGURES.rglob("*")):
        if path.is_file() and path.name != "figure_manifest.json":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest.append({"path": str(path.relative_to(ROOT)), "sha256": digest, "bytes": path.stat().st_size})
    (FIGURES / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figures": 7, "supplementary_figures": 1, "manifest_entries": len(manifest)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
