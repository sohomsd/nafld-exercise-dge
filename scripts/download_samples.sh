#! /bin/bash

# Check if samples directory already exists and create otherwise
if [ ! -d "../samples" ]; then
  mkdir ../samples
  echo -e "### Created 'samples' directory ###\n"
fi

# Download temporary file containing ftp link and sample title
echo -e "### Download ENA Project Summary Table ###\n"
wget -nv 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1209934&result=read_run&fields=fastq_ftp,sample_title&format=tsv&download=true&limit=0' \
-O ftp_sample.tsv

# Format arguments to wget command so that each line has the sample_file_name followed by the ftp_url
wget_args=$(tail -n +2 ftp_sample.tsv | awk -F "\t|;" 'BEGIN {OFS = "\t"} {gsub("H-|-","", $2); print "../samples/" $2"_1.fastq.gz \"" $3 "\"\n" "../samples/" $2"_2.fastq.gz \"" $4 "\""}')
echo -e "\n### Processed Summary Table ###"
echo -e "${wget_args}\n"

# Download samples using xargs, reading from the formatted arguments
# Running wget in parallel
echo -e "### Downloading Samples ###\n"
echo ${wget_args} | xargs -P 0 -n2 wget -nv -O

# Remove temporary file containing ftp link and sample title
rm ftp_sample.tsv

echo -e "\n### Finished ###\n" 
exit
