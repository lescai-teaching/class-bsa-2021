
####### THIS PART IS EXECUTED ON CONTAINER VSCODE #######

### we create a folder for this exercise
cd /config/workspace 
mkdir -p chipseq_exercise
### and we move into it
cd chipseq_exercise

### then we create a folder for the reads
mkdir raw_data
### and we move into it
cd raw_data
### and we create symbolic links to the reads
ln -s /config/workspace/datasets_class/chipseq/reads/* .

### move back and create a folder for the alignments
cd /config/workspace/chipseq_exercise
mkdir -p alignments
### and we move into it
cd alignments


#################################################
## ALIGN READS TO GENOME ########################
#################################################

############## CASE 1 ##########################

bwa mem \
-t 2 \
-R "@RG\tID:chipsim\tSM:case\tPL:illumina\tLB:chipsim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/chipseq_exercise/raw_data/chr21_case.fastq.gz | samtools view -@ 2 -bhS -o chip_case.bam -

############## CASE 2 ##########################

bwa mem \
-t 2 \
-R "@RG\tID:chipsim\tSM:case2\tPL:illumina\tLB:chipsim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/chipseq_exercise/raw_data/chr21_case2.fastq.gz | samtools view -@ 2 -bhS -o chip_case2.bam -

############## CONTROL SAMPLE ###################

bwa mem \
-t 2 \
-R "@RG\tID:chipsim\tSM:control\tPL:illumina\tLB:chipsim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/chipseq_exercise/raw_data/chr21_control.fastq.gz | samtools view -@ 2 -bhS -o chip_control.bam -


####### As usual we sort the bam files

samtools sort -@ 2 chip_case.bam -o chip_case_sorted.bam
samtools sort -@ 2 chip_control.bam -o chip_control_sorted.bam
samtools sort -@ 2 chip_case2.bam -o chip_case2_sorted.bam

####### We index the bam files

samtools index chip_case_sorted.bam
samtools index chip_control_sorted.bam
samtools index chip_case2_sorted.bam


#################################################
## MARK DUPLICATES       ########################
#################################################

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


####### THIS PART IS EXECUTED ON CONTAINER RSTUDIO #######


#################################################
## CALL PEAKS CASE1 vs CONTROL ##################
#################################################

macs2 callpeak -t chip_case_md.bam \
-c chip_control_md.bam \
-f BAM \
-n class_chip \
--nomodel --extsize 150 \
--outdir macs2 2> macs2/class_chip.log


#################################################
## CALL PEAKS CASE2 vs CONTROL ##################
#################################################

macs2 callpeak -t chip_case2_md.bam \
-c chip_control_md.bam \
-f BAM \
-n class_chip2 \
--nomodel --extsize 150 \
--outdir macs2_case2 2> macs2_case2/class_chip2.log