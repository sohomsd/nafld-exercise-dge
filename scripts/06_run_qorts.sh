#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=run_qorts
#SBATCH --time=1-00:00:00
#SBATCH --mem=40G

cd ../

# Check if qc directory already exists and create otherwise
if [ ! -d "qc/qorts" ]; then
  mkdir -p qc/qorts 
  echo -e "### Created 'qc/qorts' directory ###\n"
fi

# Obtain samples base names
sample_names=$(ls alignments/*.bam | awk -F '[/.]' '{print $2}')

conda activate qorts 

# Run qorts on all of the samples
for s in ${sample_names}; do
  java -Xmx4G -jar /athena/angsd/scratch/mef3005/share/envs/qorts/share/qorts-1.3.6-1/QoRTs.jar QC \
  --generatePdfReport \
  --maxReadLength 150 \
  --title ${s} \
  alignments/${s}*.bam \
  reference/gencode.vM36.basic.annotation.gtf \
  qc/qorts/${s}
done

# Run multiqc on qorts outputs
multiqc -o qc/ -n qorts_alignment_total_report -i "Alignment QORTS" qc/qorts/*/

exit
