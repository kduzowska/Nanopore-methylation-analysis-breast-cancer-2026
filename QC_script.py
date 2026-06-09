#!/usr/bin/env python3
"""
Script 0: BedMethyl QC
Author: KDuzowska
Date: 08/05/2026
Evaluates coverage and modification quality across all samples
within target gene regions.

Outputs:
    - QC_coverage_boxplot.png
    - QC_coverage_heatmap.png
    - QC_coverage_frac_heatmap.png
    - QC_summary.csv
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pyranges as pr

BEDMETHYL_DIR = "path/to/bedmethyl_files/"
REGIONS_FILE  = "path/to/genes_locations.txt"
OUTPUT_DIR    = "path/to/output/folder/bedmethyl_QC/"
MIN_COVERAGE  = 10
CHUNK_SIZE    = 100_000  # number of rows to read at a time — controls memory usage

os.makedirs(OUTPUT_DIR, exist_ok=True)

def sample_files(directory):
    """
    Find all bedmethyl.gz files in directory.
    Returns dict: sample_name -> filepath
    """
    pattern = os.path.join(directory, "*.bedmethyl.gz")
    files   = glob.glob(pattern)

    if not files:
        raise FileNotFoundError(f"No .bedmethyl.gz files found in {directory}")

    samples = {}
    for filepath in files:
        filename      = Path(filepath).name
        parts         = filename.split(".")
        name          = parts[1]
        samples[name] = filepath

    print(f"Found {len(samples)} samples")
    return samples


def load_gene_regions(regions_file):
    genes = pd.read_csv(regions_file, sep="\t", header=0)
    genes = genes.rename(columns={
        "chr":          "chrom",
        "start (-5kb)": "start",
        "stop (+5kb)":  "end",
        "Gene":         "gene"
    })
    genes["chrom"] = "chr" + genes["chrom"].astype(str)

    print(f"Loaded {len(genes)} gene regions across "
          f"{genes['chrom'].nunique()} chromosomes")
    return genes


def load_bedmethyl(filepath, sample_name, relevant_chroms):
    """
    Load a bedmethyl.gz file in chunks, filtering to relevant
    chromosomes on the fly to minimize memory usage.
    Chunked reading keeps only 100k rows in RAM at a time, discarding
    irrelevant chromosomes immediately before accumulating results.
    """
    cols_load = [
        "chrom", "chromStart", "name",
        "valid_coverage", "percent_modified",
        "count_modified", "count_fail"
    ]

    # pd.read_csv with chunksize returns an iterator
    chunk_iter = pd.read_csv(
        filepath,
        sep       = "\t",
        comment   = "#",
        header    = None,
        usecols   = [0, 1, 3, 9, 10, 11, 15],
        names     = cols_load,
        chunksize = CHUNK_SIZE
    )

    filtered = []
    for chunk in chunk_iter:
        chunk = chunk[chunk["chrom"].isin(relevant_chroms)]

        if len(chunk) > 0:
            filtered.append(chunk)

    if not filtered:
        print(f"  WARNING: no data retained for {sample_name} "
              f"— check chromosome name format")
        return pd.DataFrame(columns=cols_load + ["sample"])

    df            = pd.concat(filtered, ignore_index=True)
    df["sample"]  = sample_name
    return df

def filter_regions(df, genes):
    # PyRanges requires specific column names
    df_pr = pr.PyRanges(
        df.rename(columns={
            "chrom":      "Chromosome",
            "chromStart": "Start"
        }).assign(End=df["chromStart"] + 1)
    )

    genes_pr = pr.PyRanges(
        genes.rename(columns={
            "chrom": "Chromosome",
            "start": "Start",
            "end":   "End"
        })
    )

    result = df_pr.join(genes_pr)

    if result.df.empty:
        print("  WARNING: no overlaps found — check coordinate formats")

    return result.df

def load_all_samples(samples, genes):
    # Extract the set of chromosomes 
    relevant_chroms = set(genes["chrom"].unique())
    print(f"Filtering to chromosomes: {sorted(relevant_chroms)}")

    all_dfs = []
    for name, filepath in samples.items():
        print(f"Loading {name}...")
        df = load_bedmethyl(filepath, name, relevant_chroms)
        df = filter_regions(df, genes)
        all_dfs.append(df)
        print(f"  {name}: {len(df)} positions retained")

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nDone. Total rows across all samples: {len(combined)}")
    return combined

def calculate_coverage_summary(df, min_coverage):
    """
    Uses only 'name == m' rows to avoid double-counting —
    each CpG position appears twice in bedMethyl (h and m rows)
    but valid_coverage is identical for both.
    """
    df_m    = df[df["name"] == "m"]
    grouped = df_m.groupby(["sample", "gene"])

    summary = grouped["valid_coverage"].agg(
        mean_coverage   = "mean",
        median_coverage = "median",
        frac_above_min  = lambda x: (x >= min_coverage).mean()
    ).reset_index()

    return summary

def plot_coverage_boxplot(summary, output_dir):
    fig, ax = plt.subplots(figsize=(14, 5))

    sns.boxplot(
        data  = summary,
        x     = "sample",
        y     = "mean_coverage",
        ax    = ax,
        color = "steelblue"
    )

    ax.axhline(
        y         = MIN_COVERAGE,
        color     = "red",
        linestyle = "--",
        label     = f"Min coverage threshold ({MIN_COVERAGE}x)"
    )

    ax.set_title("Coverage distribution per sample across target genes",
                 fontsize=13)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Mean coverage (valid reads)")
    plt.xticks(rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()

    outpath = os.path.join(output_dir, "QC_coverage_boxplot.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved: {outpath}")

def plot_coverage_heatmap(summary, output_dir):
    pivot = summary.pivot(
        index   = "gene",
        columns = "sample",
        values  = "mean_coverage"
    )

    fig, ax = plt.subplots(figsize=(20, 80))

    sns.heatmap(
        pivot,
        cmap       = "YlOrRd",
        ax         = ax,
        linewidths = 0.3,
        cbar_kws   = {"label": "Mean coverage (valid reads)"}
    )

    ax.set_title("Per-gene mean coverage across all samples", fontsize=13)
    plt.tight_layout()

    outpath = os.path.join(output_dir, "QC_coverage_heatmap.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved: {outpath}")

def plot_fraction_heatmap(summary, output_dir):
    # Heatmap of fraction of CpGs above MIN_COVERAGE — genes × samples.
    pivot = summary.pivot(
        index   = "gene",
        columns = "sample",
        values  = "frac_above_min"
    )

    fig, ax = plt.subplots(figsize=(20, 80))

    sns.heatmap(
        pivot,
        cmap       = "RdYlGn",
        ax         = ax,
        linewidths = 0.3,
        vmin       = 0,   
        vmax       = 1,    
        cbar_kws   = {"label": f"Fraction of CpGs above {MIN_COVERAGE}x coverage"}
    )

    ax.set_title("Per-gene fraction of CpGs above coverage threshold",
                 fontsize=13)
    plt.tight_layout()

    outpath = os.path.join(output_dir, "QC_coverage_frac_heatmap.png")
    plt.savefig(outpath, dpi=300)
    plt.close()
    print(f"Saved: {outpath}")

def save_summary_table(summary, output_dir):
    outpath = os.path.join(output_dir, "QC_summary.csv")
    summary.to_csv(outpath, index=False)
    print(f"Saved: {outpath}")

# MAIN
samples    = sample_files(BEDMETHYL_DIR)
genes      = load_gene_regions(REGIONS_FILE)
quality_df = load_all_samples(samples, genes)
summary    = calculate_coverage_summary(quality_df, MIN_COVERAGE)

plot_coverage_boxplot(summary, OUTPUT_DIR)
plot_coverage_heatmap(summary, OUTPUT_DIR)
plot_fraction_heatmap(summary, OUTPUT_DIR)
save_summary_table(summary, OUTPUT_DIR)
