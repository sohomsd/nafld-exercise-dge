#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=align_samples
#SBATCH --time=6:00:00
#SBATCH --mem=40G

cd ../

# Check if samples/trimmed directory exists, otherwise exit program
if [ ! -d "samples/trimmed" ]; then
  echo "samples/trimmed directory not found"
  exit 1
fi

# Check if genome directory exists, otherwise exit program
if [ ! -d "reference/GRCm39_STARindex" ]; then
  echo "reference/GRCm39_STARindex directory not found"
  exit 1
fi

# Check if alignments directory exists, otherwise create directory
if [ ! -d "alignments" ]; then
  mkdir alignments
  echo -e "### Created alignments directory ###\n"
fi

# Obtain the base names of the samples (without _1 and _2 from paired run identification)
sample_names=$(ls -l samples/trimmed/*.fq.gz | awk '{gsub("(samples/trimmed/|_[1-2]_val_[1-2].fq.gz)", "", $9); print $9}' | sort | uniq)

conda activate angsd

# Iterate over sample base names
for s in ${sample_names}; do

  # Check if alignment was already run for current sample, and continue to next sample if so
  if [ -f "alignments/${s}.Log.final.out" ]; then
    echo "### ${s} already aligned ###"
    continue 
  fi

  echo -e "### Running STAR on ${s} ###\n"

  # Run STAR for each paired end sample
  # Sort bam file by coordinate
  # Output all SAM attributes
  # Include unmapped reads within the BAM file
  STAR --runMode alignReads \
  --runThreadN 16 \
  --genomeDir reference/GRCm39_STARindex \
  --readFilesIn samples/trimmed/${s}_1_val_1.fq.gz samples/trimmed/${s}_2_val_2.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix alignments/${s}. \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within

  # Index sorted BAM file
  samtools index alignments/${s}.Aligned.sortedByCoord.out.bam

done


exit
