#!/bin/bash

mkdir -p '/home/dthuyy/Documents/Thesis/'
mkdir -p '/home/dthuyy/Documents/Thesis/align/'
mkdir -p '/home/dthuyy/Documents/Thesis/align/bam/'
mkdir -p '/home/dthuyy/Documents/Thesis/align/multiqc/'
mkdir -p '/home/dthuyy/Documents/Thesis/align/stats/'
mkdir -p '/home/dthuyy/Documents/Thesis/align/stats/mosdepth/'
mkdir -p '/home/dthuyy/Documents/Thesis/align/stats/qualimap/'
mkdir -p '/home/dthuyy/Documents/Thesis/raw'
mkdir -p '/home/dthuyy/Documents/Thesis/raw/qc_check/'
mkdir -p '/home/dthuyy/Documents/Thesis/ref'
mkdir -p '/home/dthuyy/Documents/Thesis/ref/genome/'
mkdir -p '/home/dthuyy/Documents/Thesis/ref/annotation/'
mkdir -p '/home/dthuyy/Documents/Thesis/sra'
mkdir -p '/home/dthuyy/Documents/Thesis/trim'
mkdir -p '/home/dthuyy/Documents/Thesis/trim/qc_check'
mkdir -p '/home/dthuyy/Documents/Thesis/variant_calling'


# SET UP PATH LINK 
export  p_align='/home/dthuyy/Documents/Thesis/align' \
        p_raw='/home/dthuyy/Documents/Thesis/raw' \
        p_ref='/home/dthuyy/Documents/Thesis/ref' \
        p_sra='/home/dthuyy/Documents/Thesis/sra' \
        p_tools='/home/dthuyy/Documents/Thesis/tools' \
        p_trim='/home/dthuyy/Documents/Thesis/trim' \
        p_var='/home/dthuyy/Documents/Thesis/variant_calling'

# RETRIEVE LINK FOR DOWNLOAD
links=(https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6410466/SRR6410466 \
    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6410467/SRR6410467 \
    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6410468/SRR6410468 \
    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6410469/SRR6410469 \
    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6410470/SRR6410470 
    )
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/ERR002/352/ERR2352264.sralite.1
# DOWNLOAD SRA
for link in ${links[@]}; do
    aria2c -x 16 -s 32 -j3 "$link" -d $p_sra
done

# ADD '.SRA' EXTENSION FOR DOWNLOADED FILE 
# Loop through all files in the directory
for file in $(ls $p_sra/SRR*); do
    # check if the file doesn't already have the .sra extension
    if [["$file" != *.sra]]; then
        mv "$file" "$file.sra"
    fi
# done

# CONVERT INTO ".fastq.gz" FILE
fastq-dump \
    --outdir $p_raw \
    --split-spot \
    --skip-technical \
    --split-files $p_sra/*.sra \
    && gzip -f $p_raw/*.fastq

# CHANGE FILE NAME 
mv $p_raw/SRR6410466_1.fastq.gz $p_raw/Samp_1_R1.fastq.gz
mv $p_raw/SRR6410466_2.fastq.gz $p_raw/Samp_1_R2.fastq.gz
mv $p_raw/SRR6410467_1.fastq.gz $p_raw/Samp_2_R1.fastq.gz
mv $p_raw/SRR6410467_2.fastq.gz $p_raw/Samp_2_R2.fastq.gz
mv $p_raw/SRR6410468_1.fastq.gz $p_raw/Samp_3_R1.fastq.gz
mv $p_raw/SRR6410468_2.fastq.gz $p_raw/Samp_3_R2.fastq.gz
mv $p_raw/SRR6410469_1.fastq.gz $p_raw/Samp_4_R1.fastq.gz
mv $p_raw/SRR6410469_2.fastq.gz $p_raw/Samp_4_R2.fastq.gz
mv $p_raw/SRR6410470_1.fastq.gz $p_raw/Samp_5_R1.fastq.gz
mv $p_raw/SRR6410470_2.fastq.gz $p_raw/Samp_5_R2.fastq.gz

# CHECK fastqc RAW DATA 
for file in $(ls $p_raw/*.fastq.gz); do
    fastqc $file -o $p_raw/qc_check
done 

# RAW DATA TRIMMING AND FILTERING 
## Loop over the files ending with ".gz"
for file in $(ls $p_raw/*_R1*.gz); do
    # Extract the base name without the extension
    base_name="${file%_*.*.*}"

    # Construct the input and output file names
    read_1=${base_name}_R1.fastq.gz
    read_2=${base_name}_R2.fastq.gz
    output_paired_1="${base_name}_R1_paired.fastq.gz"
    output_unpaired_1="${base_name}_R1_unpaired.fastq.gz"
    output_paired_2="${base_name}_R2_paired.fastq.gz"
    output_unpaired_2="${base_name}_R2_unpaired.fastq.gz"
    
    # Run Trimmomatic with the specified parameters
    java -jar $p_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads 4 \
        -trimlog "$p_trim/trim_test.log" \
        "$read_1" \
        "$read_2" \
        "$p_trim/$output_paired_1" \
        "$p_trim/$output_unpaired_1" \
        "$p_trim/$output_paired_2" \
        "$p_trim/$output_unpaired_2" \
        SLIDINGWINDOW:5:30 \
        MINLEN:30
    
    # Remove unpaired file
    rm -f "$p_trim/$output_unpaired_1" "$p_trim/$output_unpaired_2"

    # Check fastqc again
    fastqc $p_trim/$output_paired_1 -o $p_trim/qc_check
    fastqc $p_trim/$output_paired_2 -o $p_trim/qc_check
done


# MAPPING (~ALIGNMENT)
# Generate genome index (run once before alignment steps)
cd $p_trim

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir $p_ref/genome \
     --genomeFastaFiles $p_ref/genome/EV_sequence.fa \
     --sjdbGTFfile $p_ref/genome/edit_genome.gff3 \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbOverhang 100


for file in $(ls *_R1*.gz); do
    # Extract the base name without the extension
    base_name="${file%_*_*.*.*}"

    # Construct the input and output file names
    read_1=${base_name}_R1_paired.fastq.gz
    read_2=${base_name}_R2_paired.fastq.gz

    # Mapping with STAR
    STAR \
    --runThreadN 12 \
    --readFilesType Fastx \
    --readFilesCommand zcat \
    --genomeDir $p_ref/genome \
    --readFilesIn $read_1 $read_2 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmode Full \
    --outFileNamePrefix $p_align/${base_name}_ \
    > $p_align/${base_name}_STAR.log 2>&1
 
    # Check whether STAR run successful or not
    if [$? -ne 0]; then
        echo "STAR alignment failed for ${base_name}. Check $p_align/${base_name}_STAR.log for details."
    else
        echo "STAR alignment completed for ${base_name}."
    fi
done

# REMOVE DUPLICATE (Post-processing alignment data)
mv $p_align/*.bam $p_align/bam/
cd $p_align/bam/

# Define Picard jar path
picard_jar="/home/dthuyy/anaconda3/envs/python3/share/picard-3.2.0-0/picard.jar"

# Define read group information
RGID="id"
RGLB="library"
RGPL="platform"
RGPU="unit"
RGSM="sample"

# Add or replace read groups in all BAM files
for file in *.bam; do
    # Construct base name
    base_name="${file%Aligned*}"
    
    # Add or replace read groups
    java -jar $picard_jar AddOrReplaceReadGroups \
        I=${file} \
        O=${base_name}Aligned.sortedByCoord_RG.bam \
        RGID=${RGID} \
        RGLB=${RGLB} \
        RGPL=${RGPL} \
        RGPU=${RGPU} \
        RGSM=${RGSM}
    
    # Remove duplicates in the BAM file with read groups
    java -jar $picard_jar MarkDuplicates \
        I=${base_name}Aligned.sortedByCoord_RG.bam \
        O=${base_name}Aligned.sortedByCoord_RG_rmdup.bam \
        M=${base_name}Aligned.sortedByCoord_RG_rmdup.metrics2 \
        REMOVE_DUPLICATES=true
    
    # Optionally, you can remove the intermediate BAM files
    rm ${base_name}Aligned.sortedByCoord_RG.bam
done

# INDEX PREPROCESSED BAM 
samtools index $p_align/bam/*rmdup.bam 

# BAM QC
for bam_file in $p_align/bam/*_rmdup.bam; do
    # Extract the base name from the file
    base_name=$(basename "$bam_file" _rmdup.bam)

    # Run samtools stat
    samtools stat \
        --threads 1 \
        --reference "$p_ref/genome/EV_sequence.fa" \
        "$bam_file" \
    > "$p_align/stats/${base_name}.bam.stats"

    # Run mosdepth
    mosdepth \
        --threads 4 \
        --fasta "$p_ref/genome/EV_sequence.fa" \
        -n --fast-mode --by 500 \
        "$p_align/stats/mosdepth/${base_name}_final" \
        "$bam_file"
    
    # Run qualimap
    qualimap bamqc \
        -bam "$bam_file" \
        -p non-strand-specific \
        --collect-overlap-pairs \
        -outdir "$p_align/stats/qualimap/" \
        -nt 2
    
    # Run samtools coverage
    samtools coverage \
        "$bam_file" \
        > "$p_align/stats/${base_name}_cov.stats"

    # Run samtools flagstats
    samtools flagstats \
        "$bam_file" \
        > "$p_align/stats/${base_name}.flagstats"
done

# sample 2
STAR \
    --runThreadN 12 \
    --readFilesType Fastx \
    --readFilesCommand zcat \
    --genomeDir /home/dthuyy/Documents/Thesis/ref/genome \
    --readFilesIn /home/dthuyy/Documents/Thesis/trim/Samp_2_R1_paired.fastq.gz /home/dthuyy/Documents/Thesis/trim/Samp_2_R2_paired.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmode Full \
    --outFileNamePrefix /home/dthuyy/Documents/Thesis/align/Samp_2_ \
    --outSAMattrRGline ID:id LB:library PL:platform PU:unit SM:sample

java -jar $picard_jar MarkDuplicates \
    I=/home/dthuyy/Documents/Thesis/align/bam/Samp_2_Aligned.sortedByCoord.out.bam \
    O=/home/dthuyy/Documents/Thesis/align/bam/Samp_2_Aligned.sortedByCoord_rmdup.bam \
    M=${base_name}Aligned.sortedByCoord_rmdup.metrics2 \
    REMOVE_DUPLICATES=true
    
samtools index /home/dthuyy/Documents/Thesis/align/bam/Samp_2_Aligned.sortedByCoord_rmdup.bam


