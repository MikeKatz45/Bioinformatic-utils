#!/bin/bash

#Download all assembly reports of interest from NCBI
#NCBI homepage->Assembly database organism search->filter contig level
#->filter refseq->download assemblies button->refseq assembly reports
assemblies_folder=$1

usage(){
  echo -e "\nUsage: $0 <folder_with_NCBI_assembly_structure_reports>\n"
}

if [[ $# -ne 1 ]]
then
  usage
  exit
fi

total_assemblies=$(cat "$assemblies_folder"/md5checksums.txt | wc -l)
echo "#Total number of analyzed assemblies: $total_assemblies" > contigs_stats

echo -e "WGS_ID\tContigs\tShortest\tLongest\tCoverage\tDate\tPlatform\tAssembler\n" \
        >> contigs_stats

for assembly in $(ls "$assemblies_folder")
do

  if [[ $assembly != README.txt ]] && [[ $assembly != md5checksums.txt ]]
  then

    WGS_ID=$(cat "$assemblies_folder"/"$assembly" | grep WGS | cut -d' ' -f7)

    contigs=$(cat "$assemblies_folder"/"$assembly" | grep -v '#' | wc -l)

    short_contig=$(cat "$assemblies_folder"/"$assembly" | grep -v '#' | \
                  awk -F '\t' '{print $9}' | sort -n | head -1)

    long_contig=$(cat "$assemblies_folder"/"$assembly" | grep -v '#' | \
                  awk -F '\t' '{print $9}' | sort -n | tail -1)

    coverage=$(cat "$assemblies_folder"/"$assembly" | grep cov | cut -d' ' -f4)

    date=$(cat "$assemblies_folder"/"$assembly" | grep Dat | cut -d' ' -f13)

    platform=$(cat "$assemblies_folder"/"$assembly" | grep tec | cut -d' ' -f4-)

    assembler=$(cat "$assemblies_folder"/"$assembly" | grep met | cut -d' ' -f4-)
  fi

  echo -e "$WGS_ID\t$contigs\t$short_contig\t$long_contig\t$coverage\t$date\t$platform\t$assembler" \
          >> contigs_stats
done

cat contigs_stats | uniq > tmp && mv tmp contigs_stats
