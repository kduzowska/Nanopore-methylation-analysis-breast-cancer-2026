#!/usr/bin/env python3
"""
Script 0: DMR with automatic segmentation QC
Author: KDuzowska
Date: 08/05/2026
Evaluates and filters the quality metrics in the DMR file

Outputs:
    - QC plots of cohen's h, score and num_sites
    - QC plot of size (bp) distribution
    - Excel file with statistics for  cohen's h, score and num_sites
    - a csv file cleaned, filtered (name = different, num_sites >= 5) and filtered with my regions txt file
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pyranges as pr

# Load the .txt DMR file
filepath = sys.argv[1]
REGIONS_FILE = "path/to/regions.txt"

dmr_file = pd.read_csv(filepath, sep = "\t")

# Distributions of Cohen's h, number of CpG sites per DMR and the scores (more positive = bigger difference in methylation between sample b (tested) and a (reference))
def basic_distrib(file, filepath):

    basename = Path(filepath).stem

    fig, axes = plt.subplots(1,3, figsize=(15,5))

    axes[0].hist(file["cohen_h"])
    axes[0].set_title("Distribution of cohen_h")

    axes[1].hist(file["num_sites"])
    axes[1].set_title("Distribution of CpG sites for each DMR")

    axes[2].hist(file["score"])
    axes[2].set_title("Distribution of scores")

    plt.tight_layout()

    plt.savefig(f"QC_distributions_{basename}.png", dpi = 300)

    
basic_distrib(dmr_file, filepath)

# Distribution of size (bp)
def size_distrib(file, filepath):

    basename = Path(filepath).stem

    size = file[["chrom_start", "chrom_end"]]

    diff = size["chrom_end"] - size["chrom_start"]
    plt.hist(diff)
    plt.savefig(f"size_distributions_{basename}.png", dpi = 300)
    
    return diff

size_distrib(dmr_file, filepath)

# Create an excel file with QC calculations
def excel_calc(file, filepath):
    basename = Path(filepath).stem

    with pd.ExcelWriter(f"QC_stats_{basename}.xlsx") as writer:
        file["cohen_h"].describe().to_excel(writer,sheet_name="cohen_h")
        file["num_sites"].describe().to_excel(writer,sheet_name="num_sites")
        file["score"].describe().to_excel(writer,sheet_name="score")

excel_calc(dmr_file, filepath)

# Separate columns into h and m columns for counts and percentages
def separate_hm(file):

    file = file.rename(columns={"#chrom": "chrom"})

    def parse_counts(x):
        h, m = np.nan, np.nan
        try:
            parts = str(x).split(",")
            for p in parts:
                if "h:" in p:
                    h = float(p.split(":")[1])
                if "m:" in p:
                    m = float(p.split(":")[1])
        except:
            pass
        return pd.Series([h, m])

    file[["a_counts_h", "a_counts_m"]] = file["a_counts"].apply(parse_counts)
    file[["b_counts_h", "b_counts_m"]] = file["b_counts"].apply(parse_counts)

    def parse_perc(x):
        h, m = np.nan, np.nan
        try:
            parts = str(x).split(",")
            for p in parts:
                if "h:" in p:
                    h = float(p.split(":")[1])
                if "m:" in p:
                    m = float(p.split(":")[1])
        except:
            pass
        return pd.Series([h, m])

    file[["a_percentages_h", "a_percentages_m"]] = file["a_percentages"].apply(parse_perc)
    file[["b_percentages_h", "b_percentages_m"]] = file["b_percentages"].apply(parse_perc)

    # drop original columns
    file = file.drop(columns=["a_counts", "b_counts", "a_percentages", "b_percentages"])

    return file

clean_dmr_file = separate_hm(dmr_file)

# Apply filters (cohen's h and filter out to my specified regions) and create a csv file for further analyses

def filtered_csv(file, filepath, regions_file):

    basename = Path(filepath).stem

    # Apply filters
    file = file[file["name"] == "different"]
    file = file[file["num_sites"] >= 5]

    # Compute deltas on fraction scale
    file["delta_5mC"]  = (file["b_percentages_m"] - file["a_percentages_m"]) / 100
    file["delta_5hmC"] = (file["b_percentages_h"] - file["a_percentages_h"]) / 100

    # Load regions
    regions = pd.read_csv(regions_file, sep="\t")
    regions = regions.rename(columns={
        "chr":          "chrom",
        "start (-5kb)": "reg_start",
        "stop (+5kb)":  "reg_end",
        "Gene":         "gene"
    })
    # Add chr prefix to match DMR file format
    regions["chrom"] = "chr" + regions["chrom"].astype(str)

    # Manual interval overlap — no PyRanges needed
    # For each DMR, check if it overlaps any gene region
    # Overlap condition: DMR start < region end AND DMR end > region start
    results = []
    for _, dmr in file.iterrows():
        matches = regions[
            (regions["chrom"]     == dmr["chrom"]) &
            (regions["reg_start"] <= dmr["chrom_end"]) &
            (regions["reg_end"]   >= dmr["chrom_start"])
        ]
        for _, region in matches.iterrows():
            row = dmr.copy()
            row["gene"] = region["gene"]
            results.append(row)

    if not results:
        print("WARNING: no overlaps found — check chromosome formats")
        return pd.DataFrame()

    result_file = pd.DataFrame(results)

    output_path = f"filtered_dmr_{basename}.csv"
    result_file.to_csv(output_path, index=False)
    print(f"Saved: {output_path} with {len(result_file)} rows")

    return result_file

filtered_csv(clean_dmr_file, filepath, REGIONS_FILE)



