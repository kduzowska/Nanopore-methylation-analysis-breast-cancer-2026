#!/bin/bash

# loop over all txt files in current directory
for file in *.txt
do
    echo "Processing $file"
    
    python3 QC_DMR_script.py "$file"
    
done