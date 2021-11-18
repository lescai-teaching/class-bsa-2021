###### allineamento
cd /home/lescai/DIDATTICA/simulations/chipseq/encode_data/EGR1/simulated_chr21/

bwa mem \
-t 8 \
-R "@RG\tID:chipsim\tSM:case\tPL:illumina\tLB:chipsim" \
/home/lescai/DIDATTICA/small_refs/chr21/reference/Homo_sapiens_assembly38_chr21.fasta \
chr21_case.fastq | samtools view -@ 8 -bhS -o chip_case.bam -
cd /home/lescai/DIDATTICA/simulations/chipseq/encode_data/EGR1/simulated_chr21/

bwa mem \
-t 8 \
-R "@RG\tID:chipsim\tSM:case2\tPL:illumina\tLB:chipsim" \
/home/lescai/DIDATTICA/small_refs/chr21/reference/Homo_sapiens_assembly38_chr21.fasta \
chr21_case2.fastq | samtools view -@ 8 -bhS -o chip_case2.bam -

bwa mem \
-t 8 \
-R "@RG\tID:chipsim\tSM:control\tPL:illumina\tLB:chipsim" \
/home/lescai/DIDATTICA/small_refs/chr21/reference/Homo_sapiens_assembly38_chr21.fasta \
chr21_control.fastq | samtools view -@ 8 -bhS -o chip_control.bam -

samtools sort -@ 8 chip_case.bam -o chip_case_sorted.bam
samtools sort -@ 8 chip_control.bam -o chip_control_sorted.bam
samtools sort -@ 8 chip_case2.bam -o chip_case2_sorted.bam

samtools index chip_case_sorted.bam
samtools index chip_control_sorted.bam
samtools index chip_case2_sorted.bam

###### mark duplicates

gatk MarkDuplicates \
-I chip_case_sorted.bam \
-M chip_case_metrics.txt \
-O chip_case_md.bam

gatk MarkDuplicates \
-I chip_case2_sorted.bam \
-M chip_case2_metrics.txt \
-O chip_case2_md.bam

gatk MarkDuplicates \
-I chip_control_sorted.bam \
-M chip_control_metrics.txt \
-O chip_control_md.bam


### peaks primo caso

macs2 callpeak -t chip_case_md.bam \
-c chip_control_md.bam \
-f BAM \
-n class_chip \
--nomodel --extsize 150 \
--outdir macs2 2> macs2/class_chip.log


### peaks secondo caso

macs2 callpeak -t chip_case2_md.bam \
-c chip_control_md.bam \
-f BAM \
-n class_chip2 \
--nomodel --extsize 150 \
--outdir macs2_case2 2> macs2_case2/class_chip2.log