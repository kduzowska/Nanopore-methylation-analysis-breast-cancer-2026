# Analysis of methylation differences in breast cancer field cancerization - POSTER ESHG 2026
## This repository details the analysis conducted to obtain results presented at the ESHG 2026 conference

The samples have been prepared and sequenced using Oxford Nanopore Technology PromethION with adaptive sampling and enrichment of targeted regions specified in regions.csv file attached.

The initial steps involve standard file processing via EPI2ME: basecalling (Dorado) and alignment to hg38.

1. Modkit 
Command:
$MODKIT" pileup "$bam" "$out" \
    --cpg \ 
    --modified-bases 5mC 5hmC \
    --combine-strands \   
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

# GLOBAL SAMPLE QC ANALYSIS - script: QC_script.py
## What does it create: A methods section and supplementary figure showing per-sample coverage distributions across target regions, global modification fractions, and QC metrics 

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

#####################

After this step pairwise DMR fine-grained with automatic segmentation have been ran.
Example code and settings:
~/dist_modkit_v0.6.0_68b540b/modkit dmr pair \
  -b /path/to/dir/patientX_sampleA.bedmethyl.gz \
  -a /path/to/dir/patientX_sampleB.bedmethyl.gz \
  -o /path/to/dir/PatientX_sampleA_sampleB_dmr.bed \
  --ref /path/to/dir/hg38.analysisSet.fa \
  --base C \
  -t 14 \
  --io-threads 10 \
  -f \
  --header \
  --segment \
  /path/to/dir/PatientX_sampleA_sampleB_dmr.txt --fine-grained

The resulting .bed and .txt files have been saved to appropriate patient folders. 
For each .txt file a QC analysis and initial filtering has been conducted with the following script:

QC analysis and initial quality filtering DMR.txt files - script: QC_DMR_script.py OR Interactive_QC_DMR.ipynb

These two scripts are the same in terms of content but the Interactive_QC_DMR.ipynb script has been written in Jupyter Notebook in order to run each function step-by-step and inspect the file and the results, adjusting the parameters as preferred. This script has been used to adjust the filtering parameters and check if the script is correct. 
Next, the QC_DMR_script.py script - a .py version of the interactive .ipynb script, has been uploaded to a server hosting all the files of interest and ran there for each pariwise DMR analysis using a Shell script python_script_overlap.sh to iterate over all .txt files in the folder where the script is located in. 

Input to the script:
- DMR.txt files after running DMR with automatic segmentation
- Regions file with the regions of interest

What the script does:
- Calculates statistics and plots figures to visualize the distribution of cohen's h, score, num_sites and size (bp) distribution
- Filters the DMRs that were classified by the algorithm as "different" and have >= 5 num_sites

Output:
    - QC plots of cohen's h, score and num_sites
    - QC plot of size (bp) distribution
    - Excel file with statistics for  cohen's h, score and num_sites
    - a csv file cleaned, filtered (name = different, num_sites >= 5) and filtered with my regions txt file -> this file contains regions (genes) of interest that are differentially methylated between the compared samples











