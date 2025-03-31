#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=count_reads
#SBATCH --time=6:00:00
#SBATCH --mem=40G

# Check if featureCounts directory exists, otherwise create directory
if [ ! -d "featureCounts" ]; then
  mkdir featureCounts
  echo "### Created featureCounts directory ###"
fi

conda activate angsd

# Run featureCounts on the aligned data
# This is run in paired end mode and counting read pairs
# Summarization is performed at meta-feature (gene_id) level
featureCounts -a reference/gencode.vM36.basic.annotation.gtf \
-o featureCounts/read_counts_exon \
-p --countReadPairs \
-t exon \
-g gene_id \
-T 16 \
alignments/*.bam

exit
