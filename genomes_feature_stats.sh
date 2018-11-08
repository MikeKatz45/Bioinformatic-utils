#!/bin/bash

#Download all feature tables of interest from NCBI
#NCBI homepage->Assembly database organism search->filter complete genomes
#->filter refseq->download assemblies button->refseq feature table
genomes_folder=$1

usage(){
  echo -e "\nUsage: $0 <folder_with_NCBI_feature_tables>\n"
}

if [[ $# -ne 1 ]]
then
  usage
  exit
fi

for genome in $(ls "$genomes_folder")
do
  if [[ $genome != README.txt ]] && [[ $genome != md5checksums.txt ]]
  then
  zcat $genomes_folder/$genome | grep gene | \
  awk -F'\t' '{print $18}' >> header
  fi
done

genebasecount=0
for gene in $(cat header)
do
  genebasecount=$(expr "$genebasecount" + "$gene")
done

totalgenes=$(cat header | wc -l)
avggene=$(expr "$genebasecount" / "$totalgenes")


for genome in $(ls "$genomes_folder")
do
  if [[ $genome != README.txt ]] && [[ $genome != md5checksums.txt ]]
  then
  zcat $genomes_folder/$genome | grep protein_coding | awk -v avg="$avggene" \
  -F'\t' 'BEGIN{OFS="\t";} {if($18<avg)print $17,$15,$10,$18,$2;}' \
  >> gene_stats
  fi
done

cat gene_stats| awk -F'\t' '{print $1}' > header

for genome in $(ls "$genomes_folder")
do
  if [[ $genome != README.txt ]] && [[ $genome != md5checksums.txt ]]
  then
  zcat $genomes_folder/$genome | grep -f <( cat header) | \
  grep 'CDS.*with_protein' | awk -F'\t' '{print $14}' >> tmp
  fi
done

paste gene_stats tmp > header

for genome in $(ls "$genomes_folder")
do
  if [[ $genome != README.txt ]] && [[ $genome != md5checksums.txt ]]
  then
  zcat $genomes_folder/$genome | grep '^gene' | grep -v protein_coding | awk -v \
  avg="$avggene" -F'\t' 'BEGIN{OFS="\t";} {if($18<avg)print $17,$15,$10,$18,$2;}' \
  >> header
  fi
done

cat header | sort | uniq > tmp && mv tmp gene_stats

total_genomes=$(cat "$genomes_folder"/md5checksums.txt | wc -l)
tab="$(printf '\t')"

cat <<EOF > header
#Total number of analyzed genomes: $total_genomes
#Average gene lenght: $avggene
#Genes with shorter than average length listed below
#Locus_tag${tab}Gene_name${tab}Strand${tab}Gene_length${tab}Type${tab}Product
EOF

cat gene_stats >> header && mv header gene_stats
