#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=trim_samples
#SBATCH --time=16:00:00
#SBATCH --mem=40G

# Check if samples directory exists, otherwise exit program
if [ ! -d "samples" ]; then
  echo "samples directory not found"
  exit 1
fi

# Check if qc/trimmed directory already exists, and create otherwise
if [ ! -d "qc/trimmed" ]; then
  mkdir -p "qc/trimmed"
  echo "Created qc/trimmed directory"
fi

# Obtain the base names of the samples (without _1 and _2 from paired run identification)
sample_names=$(ls -l samples/*.fastq.gz | awk '{gsub("_[1-2].fastq.gz", "", $9); print $9}' | sort | uniq)

# Create a temporary file to store MGI adapters for use in the fastqc command
echo -e "MGI_adapter1\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\nMGI_adapter2\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" > temp_adapters.txt

conda activate trim-galore

# Iterate through all sample base names
for s in ${sample_names}; do
  echo -e "Running TrimGalore on ${s}\n"

  # For a given sample, run trim galore in paired-end mode
  # Stringency of trimming is set to 3
  # MGI adapters are explicitly specified
  # Run fastqc on the trimmed samples, searching for MGI adapters in the temporary adapters file
  trim_galore -o samples/trimmed/ \
  --fastqc_args "--extract --delete -a temp_adapters.txt -o qc/trimmed/" \
  --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --stringency 3 \
  --paired \
  --cores 4 \
  ${s}_1.fastq.gz ${s}_2.fastq.gz

done

# Remove temporary files
rm temp_adapters.txt

exit
