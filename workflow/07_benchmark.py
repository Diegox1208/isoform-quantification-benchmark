#!/usr/bin/env python3
#### created by Diego G. Salas
#
# Stage 07 compares every method with the simulated truth, per isoform and per gene.
# Error metrics use features with at least --min-truth true fragments. Fragments a method gives to
# features with no true fragment are reported apart, as a share of the true total.
# Isoform errors are also split by the number of isoforms in the gene and by the share of
# sequence an isoform does not share with its siblings.

import argparse
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

##### Getopts #####
#-----------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--info", required=True, help="samples.tsv")
parser.add_argument("--ref", required=True, help="reference folder from stage 02")
parser.add_argument("--sim", required=True, help="simulation folder from stage 04")
parser.add_argument("--quant", required=True, help="quantification folder from stage 06")
parser.add_argument("--prefix", required=True, help="output prefix, e.g. results/benchmark")
parser.add_argument("--min-truth", type=int, default=20)
args = parser.parse_args()

TX_METHODS = ["Salmon", "kallisto", "RSEM"]
GENE_METHODS = TX_METHODS + ["featureCounts"]
ISOFORM_BINS = [(1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 10**6, "5 or more")]
UNIQUE_BINS = [(0.0, 0.1, "under 10%"), (0.1, 0.5, "10-50%"), (0.5, 1.01, "50% or more")]

##### Functions #####
#-----------------------------------------------------------
class Timer:
    def __init__(self, block):
        self.block = block

    def __enter__(self):
        self.t = time.time()

    def __exit__(self, *exc):
        print(f"{self.block}\tdone in {round(time.time() - self.t)}s", flush=True)


def read_estimate(folder, method):
    if method == "Salmon":
        t = pd.read_csv(folder / "salmon" / "quant.sf", sep="\t")
        return t.set_index("Name")["NumReads"]
    if method == "kallisto":
        t = pd.read_csv(folder / "kallisto" / "abundance.tsv", sep="\t")
        return t.set_index("target_id")["est_counts"]
    t = pd.read_csv(folder / "rsem.isoforms.results", sep="\t")
    return t.set_index("transcript_id")["expected_count"]


def read_featurecounts(path, names):
    t = pd.read_csv(path, sep="\t", comment="#")
    cols = {c: next(n for n in names if f"/{n}/" in c) for c in t.columns[6:]}
    return t.rename(columns=cols).set_index("Geneid")[list(cols.values())]


def metrics(truth, estimate, min_truth):
    kept = truth >= min_truth
    rel = (estimate[kept] - truth[kept]) / truth[kept]
    return dict(features=int(kept.sum()),
                spearman=spearmanr(truth, estimate).statistic,
                median_abs_rel_error=rel.abs().median(),
                total_ratio=estimate.sum() / truth.sum(),
                false_share=estimate[truth == 0].sum() / truth.sum())


def bin_of(value, bins):
    for low, high, label in bins:
        if low <= value <= high if isinstance(low, int) else low <= value < high:
            return label
    return None

##### Data files #####
#-----------------------------------------------------------
with Timer("Data files"):
    info = pd.read_csv(args.info, sep="\t")
    catalog = pd.read_csv(Path(args.ref) / "catalog.tsv", sep="\t").set_index("transcript_id")
    catalog["isoform_bin"] = catalog["isoforms"].map(lambda v: bin_of(v, ISOFORM_BINS))
    catalog["unique_bin"] = catalog["unique_kmer_share"].map(lambda v: bin_of(v, UNIQUE_BINS))
    names = list(info["name"])

    truth, estimates = {}, {}
    for name in names:
        f = Path(args.sim) / name / "truth.tsv"
        if not f.is_file():
            print(f"/!\\ Error : {f} does not exist, run stage 04 first", file=sys.stderr)
            sys.exit(1)
        t = pd.read_csv(f, sep="\t").set_index("transcript_id")["fragments_simulated"]
        truth[name] = t.reindex(catalog.index, fill_value=0).astype(float)
        for method in TX_METHODS:
            est = read_estimate(Path(args.quant) / name, method)
            estimates[(name, method)] = est.reindex(catalog.index, fill_value=0).astype(float)
    fc = read_featurecounts(Path(args.quant) / "featureCounts_genes.txt", names)

##### Analysis #####
#-----------------------------------------------------------
with Timer("Analysis"):
    rows, long = [], []
    gene_of = catalog["gene_id"]
    for name in names:
        t = truth[name]
        t_gene = t.groupby(gene_of).sum()
        for method in GENE_METHODS:
            if method in TX_METHODS:
                e = estimates[(name, method)]
                rows.append(dict(name=name, level="isoform", method=method, **metrics(t, e, args.min_truth)))
                long.append(pd.DataFrame(dict(transcript_id=t.index, truth=t.to_numpy(), estimate=e.to_numpy(),
                                              name=name, method=method)))
                e_gene = e.groupby(gene_of).sum()
            else:
                e_gene = fc[name].reindex(t_gene.index, fill_value=0).astype(float)
            rows.append(dict(name=name, level="gene", method=method, **metrics(t_gene, e_gene, args.min_truth)))

    by_sample = pd.DataFrame(rows).merge(info[["name", "readLength", "seed"]], on="name")
    numeric = ["spearman", "median_abs_rel_error", "total_ratio", "false_share"]
    summary = by_sample.groupby(["level", "readLength", "method"], sort=False)[numeric].agg(["mean", "std"])
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index()

    long = pd.concat(long, ignore_index=True).merge(info[["name", "readLength", "seed"]], on="name")
    long = long.merge(catalog[["gene_id", "isoforms", "unique_kmer_share", "isoform_bin", "unique_bin"]],
                      left_on="transcript_id", right_index=True)
    kept = long[long["truth"] >= args.min_truth].assign(abs_rel_error=lambda d: (d["estimate"] - d["truth"]).abs() / d["truth"])

    strata = []
    for column, bins in (("isoform_bin", ISOFORM_BINS), ("unique_bin", UNIQUE_BINS)):
        per_sample = kept.groupby(["readLength", "seed", "method", column])["abs_rel_error"].agg(["median", "size"]).reset_index()
        s = per_sample.groupby(["readLength", "method", column]).agg(
            transcripts=("size", "mean"), median_abs_rel_error_mean=("median", "mean"),
            median_abs_rel_error_std=("median", "std")).reset_index()
        s = s.rename(columns={column: "stratum"}).assign(split_by=column.replace("_bin", ""))
        s["order"] = s["stratum"].map({label: i for i, (_, _, label) in enumerate(bins)})
        strata.append(s)
    strata = pd.concat(strata, ignore_index=True).sort_values(["split_by", "readLength", "method", "order"]).drop(columns="order")

##### Table #####
#-----------------------------------------------------------
with pd.option_context("display.width", 200, "display.max_columns", 20, "display.float_format", "{:.3f}".format):
    print(summary[["level", "readLength", "method"] + [f"{c}_mean" for c in numeric]].to_string(index=False))

##### Figure #####
#-----------------------------------------------------------
colours = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
markers = ["o", "s", "^", "D"]
ink, muted, grid, surface = "#1f1f1e", "#6b6b68", "#e6e6e3", "#fcfcfb"
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5), facecolor=surface, layout="constrained",
                             gridspec_kw=dict(width_ratios=[1, 1.2]))
lengths = sorted(summary["readLength"].unique())
mid = lengths[len(lengths) // 2]
iso = summary[summary["level"] == "isoform"]
ends = []
for method, c, mk in zip(TX_METHODS, colours, markers):
    s = iso[iso["method"] == method].set_index("readLength").loc[lengths]
    y, e = 100 * s["median_abs_rel_error_mean"], 100 * s["median_abs_rel_error_std"].fillna(0)
    a1.errorbar(lengths, y, yerr=e, color=c, marker=mk, markersize=8, linewidth=2, capsize=3,
                markeredgecolor=surface, markeredgewidth=1.5)
    ends.append([y.iloc[-1], y.iloc[-1], method])
ends.sort()
gap = 0.065 * a1.get_ylim()[1]
for i in range(1, len(ends)):
    ends[i][0] = max(ends[i][0], ends[i - 1][0] + gap)
step = lengths[-1] - lengths[-2] if len(lengths) > 1 else 50
for v, y_end, method in ends:
    leader = dict(arrowstyle="-", color=muted, linewidth=0.6) if abs(v - y_end) > gap / 3 else None
    a1.annotate(method, (lengths[-1], y_end), xytext=(lengths[-1] + 0.25 * step, v), textcoords="data",
                va="center", fontsize=9, color=ink, arrowprops=leader)
a1.set_xticks(lengths, [f"2 × {l} bp" for l in lengths])
a1.set_xlim(lengths[0] - 10, lengths[-1] + 1.6 * step)
a1.set_ylim(0, max(a1.get_ylim()[1], ends[-1][0] + gap))
a1.set_ylabel("median absolute error per isoform, %", color=ink)
a1.set_title("Isoform error by read length, mean of seeds", color=ink, loc="left")

st = strata[(strata["split_by"] == "isoform") & (strata["readLength"] == mid)]
labels = [label for _, _, label in ISOFORM_BINS]
x = np.arange(len(labels))
width = 0.26
values = []
for i, (method, c) in enumerate(zip(TX_METHODS, colours)):
    s = st[st["method"] == method].set_index("stratum").reindex(labels)
    vals = 100 * s["median_abs_rel_error_mean"].to_numpy(dtype=float)
    bars = a2.bar(x + (i - 1) * width, np.nan_to_num(vals), width - 0.03, color=c, label=method)
    for b, v in zip(bars, vals):
        if np.isfinite(v):
            values.append(v)
            a2.annotate(f"{v:.1f}" if v < 10 else f"{v:.0f}", (b.get_x() + b.get_width() / 2, v),
                        xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8, color=ink)
counts = st[st["method"] == TX_METHODS[0]].set_index("stratum").reindex(labels)["transcripts"]
a2.set_xticks(x, [f"{l} isoform{'' if l == '1' else 's'}\n{0 if pd.isna(n) else round(n)} transcripts"
                  for l, n in zip(labels, counts)])
a2.set_ylim(0, 1.15 * max([1] + values))
a2.set_ylabel("median absolute error per isoform, %", color=ink)
a2.set_title(f"Isoform error by isoforms per gene at 2 × {mid} bp", color=ink, loc="left")
a2.legend(frameon=False, fontsize=9, loc="upper left")
for a in (a1, a2):
    a.set_facecolor(surface)
    a.grid(axis="y", color=grid, linewidth=0.8)
    a.set_axisbelow(True)
    for side, sp in a.spines.items():
        sp.set_visible(side in ("left", "bottom"))
        sp.set_color("#d0d0cc")
    a.tick_params(colors=ink)

##### Saving data #####
#-----------------------------------------------------------
Path(args.prefix).parent.mkdir(parents=True, exist_ok=True)
long.to_csv(f"{args.prefix}_isoform_estimates.tsv.gz", sep="\t", index=False, float_format="%.2f")
by_sample.to_csv(f"{args.prefix}_by_sample.tsv", sep="\t", index=False, float_format="%.4f")
summary.to_csv(f"{args.prefix}_summary.tsv", sep="\t", index=False, float_format="%.4f")
strata.to_csv(f"{args.prefix}_strata.tsv", sep="\t", index=False, float_format="%.4f")
fig.savefig(f"{args.prefix}_error.png", dpi=120, facecolor=surface)
print(f"Saved {args.prefix}_isoform_estimates.tsv.gz, _by_sample.tsv, _summary.tsv, _strata.tsv and _error.png")
