#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=generate_reference
#SBATCH --time=06:00:00
#SBATCH --mem=40G

# Check if reference directory already exists and create otherwise
if [ ! -d "../reference" ]; then
  mkdir ../reference
  echo -e "### Created 'reference' directory ###\n"
fi

cd ../reference

# Check if genome fa file exists and download otherwise
if [ ! -f "GRCm39.primary_assembly.genome.fa" ]; then
  echo "### Downloading GENCODE genome assembly file ###"
  wget -nv 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz' \
  -O GRCm39.primary_assembly.genome.fa.gz

  gunzip GRCm39.primary_assembly.genome.fa.gz
fi

# Check if genome annotation gtf file exists and download otherwise
if [ ! -f "gencode.vM36.basic.annotation.gtf" ]; then
  echo "### Downloading GENCODE gtf file ###"
  wget -nv 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.basic.annotation.gtf.gz' \
  -O gencode.vM36.basic.annotation.gtf.gz

  gunzip gencode.vM36.basic.annotation.gtf.gz
fi

conda activate angsd

# Run STAR genome generate
# Max read length in the samples of interest is 150

STAR -- runMode genomeGenerate \
--runThreadN 16 \
--genomeDir GRCm39_STARindex \
--genomeFastaFiles GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile gencode.vM36.basic.annotation.gtf \
--sjdbOverhang 149

exit
