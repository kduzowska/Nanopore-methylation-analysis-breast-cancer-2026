# Analysis of methylation differences in breast cancer field cancerization - POSTER ESHG 2026
## This repository details the analysis conducted to obtain results presented at the ESHG 2026 conference

The samples have been prepared and sequenced using Oxford Nanopore Technology (Minion XXX) with adaptive sampling and enrichment of targeted regions specified in regions.csv file attached.

The initial steps involve standard file processing via EPI2ME: basecalling (Dorado) and alignment to hg38.

1. Modkit 
Command:
$MODKIT" pileup "$bam" "$out" \
    --cpg \ <- report only counts from reference CpG dinucleotides!!
    --modified-bases 5mC 5hmC \
    --combine-strands \   <- base modification counts can be summed across strands!!
    --ref "$REF" \
    --bgzf \
    --header \
    --threads "$THREADS" \
    --bgzf-threads "$THREADS" \
    --sampling-threads "$THREADS" \
    --log-filepath "$log"

  echo "==> Indexing: $out"
  tabix -f -p bed "$out"

BAM files have been processed through modkit and created BEDMETHYL files with the following columns:
#chrom  chromStart      chromEnd        name    score   strand  thickStart      thickEnd        color   valid_coverage  percent_modified        count_modified  count_canonical count_other_mod count_delete    count_fail      count_diff      count_nocall

GLOBAL SAMPLE QC ANALYSIS - script: QC_script.py
Why/Results: A methods section and supplementary figure showing per-sample coverage distributions across target regions, global modification fractions, and QC metrics 

Input to the script: bedmethyl.gz files after modkit pileup

What the script does: 
- Summarizes all bedmethyl files in a given folder
- The user needs to adjust the initial filepaths and parameters: 
	- BEDMETHYL_DIR - where the bedmethyl.gz files are
	- REGIONS_FILE - where the txt file containing regions of interest (enriched regions) are
	- OUTPUT_DIR - where the output files should go
	- MIN_COVERAGE - standard: 10
	- CHUNK_SIZE - for memory usage optimization set a number of rows to be processed at a time (100 000 tested and recommended)
- The script filters the files to the regions of interest only - saves time and memory
- Calculates coverage summary, creates a heatmap (QC_coverage_heatmap.png) and a boxplot summary of coverage for the samples (QC_coverage_boxplot.png) and a heatmap that shows the fraction of CpGs that are above MIN_COVERAGE for each region and sample (QC_coverage_frac_heatmap.png).
- The QC_summary.csv file contains mean and median coverage as well as fraction above MIN_COVERAGE for each sample and region.

Output of the script: 
Coverage summaries:
    - QC_coverage_boxplot.png
    - QC_coverage_heatmap.png
    - QC_coverage_frac_heatmap.png
    - QC_summary.csv















