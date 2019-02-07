#!/bin/bash

#For working with report after "Complete genome gene stats" usage
#Download gene features nucleotide fasta of complete genomes of interest
#NCBI homepage->Assembly database organism search->filter complete genomes
#->related search nucleotide refseq->select only chromosomes not plasmids->send
sequences=$1
gene_stats=$2

usage(){
  echo -e "\nUsage: $0 <sequences_file> <gene_stats_file>\n"
}

if [[ $# -ne 2 ]]
then
  usage
  exit
fi

cat $gene_stats | grep -v '#' | awk -F'\t' '{print $1}' > tmp
for locus in $(cat tmp)
do
  cat $sequences | sed -n '/.*'"$locus"'.*/,/>/p' | head -n -1 >> gene_db.fasta
done

cat gene_db.fasta | sed 's/>.*\[locus_tag=/>/g' | \
cut -d] -f1 > tmp && mv tmp gene_db.fasta
