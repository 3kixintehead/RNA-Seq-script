#!/bin/bash

SECONDS=0

# change working directory
cd ~/RNASeq1

# List of FASTQ files for retrieval
fastq_files=("data/SRR8481467.fastq" "data/SRR8481466.fastq" "data/SRR8481464.fastq" "data/SRR8481465.fastq")

# Download genome indices
## wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# For loop to cycle fastqc, trimmomatic, hisat2 alignment and concatenation into counts file (.csv will need further edits for DESeq2 metadata)

for fastq_file in "${fastq_files[@]}"; do

    # Extract the base filename without path and extension
    base_name=$(basename "$fastq_file" .fastq)

    ### STEP 1: Run fastqc
    # run trimmomatic to trim reads with poor quality

    java -jar ~/RNASeq1/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 "$fastq_file" "data/${base_name}_trimmed.fastq" TRAILING:10 -phred33
    echo "Trimmomatic finished running."

    fastqc "data/${base_name}_trimmed.fastq" -o data/
    echo "fastqc finished running."

    ### STEP 2: Run HISAT2
    hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U "data/${base_name}_trimmed.fastq" -p 4 | samtools sort -o "HISAT2/${base_name}_trimmed.bam"
    echo "HISAT2 finished running."

    ### STEP 3: Run featureCounts - Quantifion
    featureCounts -S 2 -a HS.GRCh38.106/Homo_sapiens.GRCh38.106.gtf -o "quants/${base_name}_featurecounts.txt" HISAT2/"${base_name}_trimmed.bam"
    echo "featureCounts finished running."

    # Extract raw counts and append to a CSV file
    awk '{print $1, $7}' "quants/${base_name}_featurecounts.txt" | tail -n +3 > temp_${base_name}.txt
    temp_files+=("temp_${base_name}.txt")

done

# Use paste to merge all temp files column-wise into your final CSV
paste "${temp_files[@]}" | awk 'BEGIN{OFS=","}{$1=$1; print}' > GSE125554_zika_cts.csv

# Cleanup temp files
rm "${temp_files[@]}"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
