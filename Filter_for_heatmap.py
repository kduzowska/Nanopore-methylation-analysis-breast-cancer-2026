"""
Script: Heatmap data preparation
Author: KDuzowska
Date: 12.05.2026
The goal of this script is to take .txt files, calculate delta 5mC and 5hmC and create an output csv file.

Next script (Create_heatmap.py): Will take all the .csv files and merge them together into one file and create a heatmap based on a list of chosen genes
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

filepath = sys.argv[1]
REGIONS_FILE = "path/to/target_regions.txt"

dmr_file = pd.read_csv(filepath, sep = "\t")

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

# Apply filters (filter out to my specified regions) and create a csv file for further analyses

def filtered_csv(file, filepath, regions_file):

    basename = Path(filepath).stem

    # Compute deltas
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

    output_path = f"Patient1_N_all_dmr_{basename}.csv"
    result_file.to_csv(output_path, index=False)
    print(f"Saved: {output_path} with {len(result_file)} rows")

    return result_file

filtered_csv(clean_dmr_file, filepath, REGIONS_FILE)
