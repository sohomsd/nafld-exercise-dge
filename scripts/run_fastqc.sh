#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --job-name=fastqc
#SBATCH --time=06:00:00
#SBATCH --mem=40G

cd ../

# Check if reference directory already exists
if [ ! -d "samples" ]; then
  echo "samples directory not found"
  exit 1
fi

# Check if qc directory already exists and create otherwise
if [ ! -d "qc" ]; then
  mkdir qc 
  echo -e "### Created 'qc' directory ###\n"
fi

# Create a temporary file to store MGI adapters for use in the fastqc command
echo -e "MGI_adapter1\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\nMGI_adapter2\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" > temp_adapters.txt

conda activate angsd

# Run fastqc on all of the samples, searching for adapters in the temporary adapters file
fastqc samples/*fastq.gz --extract --delete -t 6 -a temp_adapters.txt -o qc/

# Remove temporary files
rm temp_adapters.txt

exit
