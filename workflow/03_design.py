#!/usr/bin/env python3
#### created by Diego G. Salas
#
# Stage 03 fixes the ground truth, how many fragments each transcript should receive.
# Every sample shares this design, and the sample seeds only change where fragments fall and their errors.
# A random share of genes is expressed with log-normal weights, and each gene splits its fragments
# between its isoforms in Dirichlet proportions, so some isoforms dominate and some get almost nothing.

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

##### Getopts #####
#-----------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--catalog", required=True, help="catalog.tsv from stage 02")
parser.add_argument("--out", required=True, help="design table, one row per transcript")
parser.add_argument("--fragments", type=int, default=1_000_000, help="fragments per sample")
parser.add_argument("--expressed", type=float, default=0.80, help="share of genes that carry fragments")
parser.add_argument("--alpha", type=float, default=0.5, help="Dirichlet concentration of the isoform shares")
parser.add_argument("--seed", type=int, default=1)
args = parser.parse_args()

if not Path(args.catalog).is_file():
    print(f"/!\\ Error : {args.catalog} does not exist, run stage 02 first", file=sys.stderr)
    sys.exit(1)

##### Functions #####
#-----------------------------------------------------------
def allocate(weights, total):
    raw = weights / weights.sum() * total
    counts = np.floor(raw).astype(int)
    rest = total - counts.sum()
    counts[np.argsort(raw - counts)[::-1][:rest]] += 1
    return counts

##### Data files #####
#-----------------------------------------------------------
catalog = pd.read_csv(args.catalog, sep="\t")
rng = np.random.default_rng(args.seed)

##### Analysis #####
#-----------------------------------------------------------
genes = catalog["gene_id"].drop_duplicates().to_numpy()
expressed = rng.random(len(genes)) < args.expressed
gene_weight = dict(zip(genes, np.where(expressed, rng.lognormal(0, 1.5, len(genes)), 0.0)))

weights = np.zeros(len(catalog))
for gene_id, idx in catalog.groupby("gene_id", sort=False).indices.items():
    shares = rng.dirichlet(np.full(len(idx), args.alpha)) if len(idx) > 1 else np.ones(1)
    weights[idx] = gene_weight[gene_id] * shares * catalog["length"].to_numpy()[idx]

design = catalog[["transcript_id", "gene_id", "isoforms", "length", "unique_kmer_share"]].copy()
design["fragments"] = allocate(weights, args.fragments)

##### Table #####
#-----------------------------------------------------------
tx = design[design["fragments"] > 0]
print(f"{expressed.sum():,} of {len(genes):,} genes expressed, "
      f"{len(tx):,} of {len(design):,} transcripts carry {tx['fragments'].sum():,} fragments")

##### Saving data #####
#-----------------------------------------------------------
design.to_csv(args.out, sep="\t", index=False)
print(f"Saved {args.out}")
