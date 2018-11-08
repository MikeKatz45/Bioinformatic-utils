#!/bin/bash

file=$1
SS=$2
forward=$3
reverse=$4
kmer=$5
exp_cov=$6
cov_cutoff=$7

usage() {

  echo -e "\n\e[39mUsage: $(basename $0) [file_ids.txt] [velveth_Sequences_file] [forward.fastq.gz] [reverse.fastq.gz] [velvet_kmer] [velvet_coverage] [velvet_cutoff] \n"

}

if [  $# -le 6 ]

  then

    usage

    exit 1

fi

regex='^[0-9]+$'
if ! [[ $(basename $forward | cut -d{ -f2 | tr -d '.fastq.gz') =~ $regex ]]
  then
	mv $forward $(basename -s .fastq.gz $forward)
	forward=$(basename -s .fastq.gz $forward)
	mv $forward $(basename ${forward}_rm{0.fastq.gz)
	forward=$(basename ${forward}_rm{0.fastq.gz)
fi

if ! [[ $(basename $reverse | cut -d{ -f2 | tr -d '.fastq.gz') =~ $regex ]]
  then
        mv $reverse $(basename -s .fastq.gz $reverse)
        reverse=$(basename -s .fastq.gz $reverse)
        mv $reverse $(basename ${reverse}_rm{0.fastq.gz)
        reverse=$(basename ${reverse}_rm{0.fastq.gz)
fi


echo -e "\n\e[33mBegin: $(date) "

echo -e '\e[39m- Extracting read IDs from provided IDs file (1/4)...'
idstmp=`cat $file | grep \> | cut -d\> -f2 | sort -n | uniq`


echo '- Finding corresponding read IDs in Sequences file (2/4)...'
for i in ${idstmp}

do

  find . -maxdepth 1 -name "$SS" -type f -print0 | \

                xargs -0 -n1 -P$(nproc) grep -P "\t${i}\t" | \

                cut -f1 | cut -d\> -f2 | sort | uniq >> rm_ids.txt

done

echo '- IDs found and duplicates ignored (3/4)...'

past_frm=$(basename $forward | cut -d{ -f2 | tr -d '.fastq.gz')
past_rrm=$(basename $reverse | cut -d{ -f2 | tr -d '.fastq.gz')
end=$(cat rm_ids.txt | wc -l)

forw=$forward
reve=$reverse

forward=$(basename $forward | cut -d{ -f1)
reverse=$(basename $reverse | cut -d{ -f1)

echo '- Deleting reads and compressing edited files (4/4)...'
paste <(zcat $forw) <(zcat $reve) | paste - - - - | \

   awk -v FS="\t" -v OFS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | grep -v -f <( cat rm_ids.txt) | \

   tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes $(nproc) > ${forward}{$(expr $(echo $past_frm) + $(echo $end)).fastq.gz ) | \

   cut -f 5-8 | tr "\t" "\n" | pigz --best --processes $(nproc) > ${reverse}{$(expr $(echo $past_rrm) + $(echo $end)).fastq.gz

forward=${forward}{$(expr $(echo $past_frm) + $(echo $end)).fastq.gz
reverse=${reverse}{$(expr $(echo $past_rrm) + $(echo $end)).fastq.gz

echo -e "\e[33mDone: $(date) \n"
echo -e '\e[32mJob finished, generating report...'
count1=$(expr $(zcat $forw | wc -l) / 4)
count2=$(expr $(zcat $reve | wc -l) / 4)

echo -e "\e[39mReads present on input files: $count1 (Forward) $count2 (Reverse)"

start=$(cat "$file" | grep \> | wc -l)
echo  "Number of provided IDs: $start "

minus=$(expr $(echo $start) - $(echo $end))
echo "Number of ignored IDs due to duplications: $minus "
echo "Total number of deleted IDs: $end "

echo -e "\nRunning 10 velvet iterations...\n"

velveth fix_1 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_1 -exp_cov $6 -cov_cutoff $7

velveth fix_2 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_2 -exp_cov $6 -cov_cutoff $7

velveth fix_3 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_3 -exp_cov $6 -cov_cutoff $7

velveth fix_4 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_4 -exp_cov $6 -cov_cutoff $7

velveth fix_5 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_5 -exp_cov $6 -cov_cutoff $7

velveth fix_6 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_6 -exp_cov $6 -cov_cutoff $7

velveth fix_7 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_7 -exp_cov $6 -cov_cutoff $7

velveth fix_8 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_8 -exp_cov $6 -cov_cutoff $7

velveth fix_9 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_9 -exp_cov $6 -cov_cutoff $7

velveth fix_10 $5 -fastq.gz -shortPaired -separate $forward $reverse
velvetg fix_10 -exp_cov $6 -cov_cutoff $7

echo -e "\n\e[32m10 velvet iterations completed"
echo -e "\e[39mGenerating assemblies report..."

echo -e "Assembly#\tContigs#\tNucleotides#" > Assemblies.txt

fix1c=$(cat fix_1/contigs.fa | grep \> | wc -l)
fix1n=$(cat fix_1/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n1\t\t$fix1c\t\t$fix1n" >> Assemblies.txt

fix2c=$(cat fix_2/contigs.fa | grep \> | wc -l)
fix2n=$(cat fix_2/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n2\t\t$fix2c\t\t$fix2n" >> Assemblies.txt

fix3c=$(cat fix_3/contigs.fa | grep \> | wc -l)
fix3n=$(cat fix_3/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n3\t\t$fix3c\t\t$fix3n" >> Assemblies.txt

fix4c=$(cat fix_4/contigs.fa | grep \> | wc -l)
fix4n=$(cat fix_4/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n4\t\t$fix4c\t\t$fix4n" >> Assemblies.txt

fix5c=$(cat fix_5/contigs.fa | grep \> | wc -l)
fix5n=$(cat fix_5/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n5\t\t$fix5c\t\t$fix5n" >> Assemblies.txt

fix6c=$(cat fix_6/contigs.fa | grep \> | wc -l)
fix6n=$(cat fix_6/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n6\t\t$fix6c\t\t$fix6n" >> Assemblies.txt

fix7c=$(cat fix_7/contigs.fa | grep \> | wc -l)
fix7n=$(cat fix_7/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n7\t\t$fix7c\t\t$fix7n" >> Assemblies.txt

fix8c=$(cat fix_8/contigs.fa | grep \> | wc -l)
fix8n=$(cat fix_8/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n8\t\t$fix8c\t\t$fix8n" >> Assemblies.txt

fix9c=$(cat fix_9/contigs.fa | grep \> | wc -l)
fix9n=$(cat fix_9/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n9\t\t$fix9c\t\t$fix9n" >> Assemblies.txt

fix10c=$(cat fix_10/contigs.fa | grep \> | wc -l)
fix10n=$(cat fix_10/contigs.fa | grep -v \> | awk '{sum += length($1)}END{print sum}')
echo -e "\n10\t\t$fix10c\t\t$fix10n" >> Assemblies.txt

echo  "Creating amos file of best assembly..."
best=$(cat Assemblies.txt | grep -v Assem | awk -F '\t' '{print $3}' | grep '\S' | sort -n | head -1)
truebest=$(cat Assemblies.txt | grep "$best" | awk -F '\t' '{print $5}' | sort -n | tail -1)
amos=$(cat Assemblies.txt | grep "$truebest" | awk -F '\t' '{print $1}' | uniq)
velvetg fix_$amos -exp_cov $6 -cov_cutoff $7 -amos_file yes

echo -e "\e[32mDone\n"
