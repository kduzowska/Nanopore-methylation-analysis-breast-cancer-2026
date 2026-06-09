"""
Script: Heatmap data preparation
Author: KDuzowska
Date: 12.05.2026
The goal of this script is to take .csv files after "Filter_for_heatmap.py" script and merge them together into one file and create a heatmap based on a list of chosen genes
"""

# libraries
import os
import glob
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path

# Load the .csv files
CSV_DIR = "path/to/csv/files/Data_for_heatmaps/"
OUTPUT_FOLDER = "path/to/output/folder/"
MARK = "delta_5mC"
GENES = ["CPT-1A", "SREBF1", "ELOVL5", "ELOVL1","PCCA", "CD36", "MGLL-MAGL", "PPARG", "ACSS3", "FASN", "LPIN1", "PBX1", "TBX3"]
#GENES = ["AKT1", "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CBFB", "CDH1",
#    "CHEK2", "ERBB3", "ESR1", "MAP3K1", "MED12", "NCOR1", "NF1", "PIK3CA",
#    "PTEN", "RAD50", "RAD51C", "RAD51D", "RB1", "STK11", "TBX3", "TP53",
#    "TSHR", "ANKRD30A", "ANKRD16", "ANKRD26", "POTEE", "POTEF", "POTEI",
#    "POTEJ", "POTE-8", "POTEM", "POTEG (POTE-14)", "POTEB (POTE-15)",
#    "POTEB2", "POTEB3", "POTEC (POTE-18)", "POTED (POTE-21)", "POTEH",
#    "PLIN1", "PLIN2", "PLIN4", "PLIN5", "LIPE-HSL", "LPL", "PNPLA2",
#    "PNPLA7", "MGLL-MAGL", "AQP7", "RBP4", "CD36", "FABP4", "FABP5",
#    "APOB", "AKR1C1", "ACSM1", "GPD1", "SDS", "PCCA", "AKR1C3",
#    "ELOVL7", "ELOVL1", "ELOVL2", "DHRS12", "ELOVL5", "ELOVL6",
#    "ACSS3", "ACSL1", "ACSL3", "ACACB", "PPARG", "FASN", "SCD",
#    "ACACA (ACC1)", "ADIPOQ", "LEP", "ENSA", "FLAD1", "PBX1",
#    "SREBF1", "SREBF2", "MYC", "CEBPD", "CPT-1A", "ACOX-1", "ACOX2",
#    "ETFA", "ETFB", "ETFDH", "ACADVL", "ACADM", "ACADS", "ACADL",
#    "ACADSB", "LPIN1", "SPHK1", "SPHK2", "CHKA", "CHKB", "ASAH1",
#    "FDFT1", "LEPROTL1", "BNIP3L", "IDI1", "INSIG1", "ACLY", "HMGCR"]

def load_csv_files(csv_dir, output_folder):
    pattern = os.path.join(csv_dir, "*.csv")
    files = glob.glob(pattern)

    all_dfs = []
    for filepath in files:
        filename=Path(filepath).stem
        parts = filename.split("_")

        patient = parts[0]
        category = parts[1]
        label = "_".join([parts[0], parts[1]])

        df = pd.read_csv(filepath)

        df["patient"] = patient
        df["category"] = category
        df["sample_label"] = label

        all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(f"{output_folder}combined_test.csv", index=False)
    return combined

data = load_csv_files(CSV_DIR, OUTPUT_FOLDER)

# Filter the dataframe so that all is left is: the genes of interest, delta_5mC and genes are rows and samples are columns

# Create a heatmap - samples grouped by their classification (SK, C1, C2, N, PT)
def build_matrix(combined, mark, genes):
    CATEGORY_ORDER = ["N", "C2", "C1", "PT"]
    # Filter to genes of interest only
    df = combined[combined["gene"].isin(genes)].copy()

    # For each gene, if there are multiple DMRs, keep the largest absolute delta
    # Create an absolute value column
    idx     = df.groupby(["sample_label", "gene"])[mark].apply(
    lambda x: x.abs().idxmax()
    ).reset_index(drop=True)  # ← flattens the MultiIndex to simple integers

    df_best = df.loc[idx.values]

    # Pivot to wide format
    matrix = df_best.pivot(index="gene",
                               columns="sample_label",
                               values= mark)
    
    matrix = matrix.dropna(how="all")
    
    def sort_key(col_name):
        category = col_name.split("_")[1]
        if category in CATEGORY_ORDER:
            return CATEGORY_ORDER.index(category)
        return 99  # unknown categories go to end

    sorted_cols = sorted(matrix.columns, key=sort_key)
    matrix      = matrix[sorted_cols]

    
    return matrix


matrix = build_matrix(data, MARK, GENES)
matrix = matrix.rename(index={"MGLL-MAGL": "MAGL"})

def build_clustermap(matrix, mark, output_dir):
    
    CATEGORY_COLORS = {
        "N":  "#2166ac",
        "C2": "#92c5de",
        "C1": "#f4a582",
        "PT": "#d6604d"
    }

    # Fill NaN with 0 for clustering — clustering needs finite values
    matrix_filled = matrix.fillna(0)
    
    # Create mask where original matrix was NaN
    # These cells will be shown as grey in the plot
    mask = matrix.isna()

    # Colors
    cmap = matplotlib.colormaps["vlag"].copy()
    cmap.set_bad("lightgrey")

    # Column annotation colors
    col_colors = matrix.columns.map(
        lambda col: CATEGORY_COLORS.get(col.split("_")[1], "grey")
    )

    # clustermap — pass filled matrix for clustering, mask for display
    g = sns.clustermap(
        matrix_filled,       # finite values — clustering works
        mask        = mask,  # hides cells that were originally NaN
        cmap        = cmap,
        center      = 0,
        vmin        = -1,
        vmax        = 1,
        col_colors  = col_colors,
        row_cluster = False,
        col_cluster = False,
        linewidths  = 0.5,
        linecolor   = "lightgrey",
        figsize     = (25, 30),
        cbar_pos = (0.02, 0.3, 0.05, 0.5),
        cbar_kws    = None
    )

    g.ax_heatmap.set_xlabel("")
    #g.ax_heatmap.set_ylabel("Gene",   fontsize=40)
    g.ax_heatmap.tick_params(axis="y", labelsize=70, rotation=0)
    plt.setp(
    g.ax_heatmap.get_xticklabels(),
    rotation=45,
    fontsize=50,
    ha="right",
    rotation_mode="anchor"
    )
    g.cax.tick_params(labelsize=40)
    #g.fig.suptitle(
    #    "Differential 5mC methylation at lipid metabolism genes\nin breast cancer field",
    #    fontsize = 16,
    #    y        = 1.02
    #)

    # Category legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color, label=cat)
        for cat, color in CATEGORY_COLORS.items()
    ]
    #g.ax_heatmap.legend(
    #    handles        = legend_elements,
    #    title          = "Sample category",
    #    bbox_to_anchor = (1.3, 1),
    #    loc            = "upper left",
    #    fontsize       = 10
    #)

    stem   = mark.replace("delta_", "")
    outpng = os.path.join(output_dir, f"clustermap_{stem}.png")
    outpdf = os.path.join(output_dir, f"clustermap_{stem}.pdf")
    plt.savefig(outpng, dpi=600, bbox_inches="tight")
    plt.savefig(outpdf,           bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpng}")

build_clustermap(matrix, MARK, OUTPUT_FOLDER)