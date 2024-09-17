#!/bin/bash

# VARIANT CALLING 
```Gatk```

# INDEX GENOME 
samtools faidx $p_ref/genome/EV_sequence.fa 
gatk CreateSequenceDictionary R=$p_ref/genome/EV_sequence.fa O=$p_ref/genome/EV_sequence.dict



gatk Mutect2 \
    -R $p_ref/genome/EV_sequence.fa \
    -I $p_align/bam/Samp_2_Aligned.sortedByCoord_rmdup.bam \
    -O $p_var/Samp_2_raw.vcf

gatk FilterMutectCalls \
    -V $p_var/Samp_2_raw.vcf \
    -O $p_var/Samp_2_filter.vcf \
    -R $p_ref/genome/EV_sequence.fa \
    --microbial-mode true \
    --stats $p_var/Samp_2_raw.vcf.stats


grep '^[^##]' /content/variant_calling_result/SRR28714678_filtered.vcf | head -n 20
bcftools view -f PASS -O v -o /content/variant_calling_result/SRR28714678_filtered_PASS.vcf /content/variant_calling_result/SRR28714678_filtered.vcf
grep '^[^##]' /content/variant_calling_result/SRR28714678_filtered_PASS.vcf | head -n 10    
echo "Before filter: $(grep '^[^##]' /content/variant_calling_result/SRR28714678_filtered.vcf | wc -l )"
echo "After filter: $(grep '^[^##]'  /content/variant_calling_result/SRR28714678_filtered_PASS.vcf | wc -l)"

# VARIANT ANNOTATION
```SnpEFF```   
