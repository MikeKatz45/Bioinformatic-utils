#!/bin/bash

#Dependencies: ripgrep (rg), awk, standard GNU commands (grep, sed, paste, tr, etc.), Velvet assembler suite, parallel gzip (pigz), gzip (zcat)

file=$1 #File with Tablet IDs from the reads to remove.
SS=$2 #Velvet sequences file created after using velveth
forward=$3 #Forward paired-end reads library
reverse=$4 #Reverse paired-end reads library
kmer=$5 #Kmer integer as estimated by Kmergenie or other kmer calculating tools (Used for velvet assembly)
exp_cov=$6 #Integer for this velvetg parameter
cov_cutoff=$7 #Integer for this velvetg parameter

usage() {

  echo -e "\n\e[39mUsage: $(basename $0) [file_ids.txt] [velveth_Sequences_file] [forward.fastq.gz] [reverse.fastq.gz] [velvet_kmer] [velvet_coverage] [velvet_cutoff] \n"

}

if [  $# -le 6 ]

  then

    usage

    exit 1

fi #Testing if script is started with enough arguments,

regex='^[0-9]+$'
if ! [[ $(basename $forward | cut -d{ -f2 | tr -d '.fastq.gz') =~ $regex ]]
  then
	mv $forward $(basename -s .fastq.gz $forward)
	forward=$(basename -s .fastq.gz $forward)
	mv $forward $(basename ${forward}_rm{0.fastq.gz)
	forward=$(basename ${forward}_rm{0.fastq.gz)
fi #Testing if paired-end libraries filenames are formatted (optional for more organized data handling when running multiple iterations of script)

if ! [[ $(basename $reverse | cut -d{ -f2 | tr -d '.fastq.gz') =~ $regex ]]
  then
        mv $reverse $(basename -s .fastq.gz $reverse)
        reverse=$(basename -s .fastq.gz $reverse)
        mv $reverse $(basename ${reverse}_rm{0.fastq.gz)
        reverse=$(basename ${reverse}_rm{0.fastq.gz)
fi #Same as above but for the reverse library.


echo -e "\n\e[33mBegin: $(date) "

echo -e '\e[39m- Extracting read IDs & parsing sequences file (1/4)...' 
ids=$(rg -j $(nproc) \> "$file" | sort -n | uniq | sed 's/>\(.*\)/\^\1;/g') Parsing tablet IDs into a variable
minus=$(echo "$ids" | tr ' ' '\n' | wc -l) #Counting number of unique tablet IDs proivided
rg -j $(nproc) \> $SS | cut -f1-2 | sed 's/>\(.*\)\t\(.*\)/\2 \1/g' > tmp #Parsing velvet sequences file for faster grep later


echo '- Finding corresponding read IDs in Sequences file (2/4)...'
for id in ${ids}
do
  rg -j $(nproc) $id tmp >> rm_ids
done #Finding read coordinates (Fastq header) with read ID (tablet) as query.

echo '- IDs found and duplicates ignored (3/4)...'

cut -d';' -f2 rm_ids | uniq > tmp && mv tmp rm_ids #Leaving only Fastq headers without read IDs (already ordered, hence no sort)

past_frm=$(basename $forward | cut -d{ -f2 | tr -d '.fastq.gz')
past_rrm=$(basename $reverse | cut -d{ -f2 | tr -d '.fastq.gz') #Filenames handling (optional)
end=$(cat rm_ids.txt | wc -l) #Counting total number of FastqHeaders obtained from initial read IDs

forw=$forward
reve=$reverse

forward=$(basename $forward | cut -d{ -f1) #Filename Handling (optional)
reverse=$(basename $reverse | cut -d{ -f1)

echo '- Deleting reads and compressing edited files (4/4)...'

paste <(zcat V03B_R1_trim_cor.fastq.gz) <(zcat V03B_R2_trim_cor.fastq.gz) | \ #Pairs both paired-end libraries (header\theader and so on) 
paste - - - - | rg -j $(nproc) -v -f rm_ids | \ #Puts the paired lines in a single line for each pair of reads (all foward and reverse 8 lines in a single line) & removes non-desired reads
awk -v FS="\t" -v OFS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | \ #Prints columns in order (all foward info first, then reverse)
tee >(cut -f1-4 | tr '\t' '\n' | pigz --best -p $(nproc) > ${forward}{$(expr $(echo $past_frm) + $(echo $end)).fastq.gz) \
| cut -f5-8 | tr '\t' '\n' | pigz --best -p $(nproc) > ${reverse}{$(expr $(echo $past_rrm) + $(echo $end)).fastq.gz #Substitutes tabs for new lines and reverts to original FASTQ format
#This block has had problems with pigz not running parallel and rg/grep being too slow so it might need tuning 

forward=${forward}{$(expr $(echo $past_frm) + $(echo $end)).fastq.gz
reverse=${reverse}{$(expr $(echo $past_rrm) + $(echo $end)).fastq.gz #Filename handling (optional)

echo -e "\e[33mDone: $(date) \n"
echo -e '\e[32mJob finished, generating report...'

count1=$(expr $(zcat $forw | wc -l) / 4)
count2=$(expr $(zcat $reve | wc -l) / 4)
echo -e "\e[39mReads present on input files: $count1 (Forward) $count2 (Reverse)"

start=$(cat "$file" | grep \> | wc -l)
echo  "Number of provided IDs: $start "
echo "Number of IDs left after duplicate removal: $minus "

minus=$(cat rm_ids | wc -l)
echo "Number of read headers that mapped to unique IDs (# of reads deleted from each file): $minus "

count3=$(($count1 - $minus))
count4=$(($count2 - $minus))
echo -e "\e[39mReads present on output files: $count3 (Forward) $count4 (Reverse)"

echo -e "\nRunning 10 velvet iterations...\n"

for i in $(seq 1 10)
do
  velveth fix_$i $5 -fastq.gz -shortPaired -separate $forward $reverse > /dev/null
  velvetg fix_$i -exp_cov $6 -cov_cutoff $7 > /dev/null
  echo "$i..."
done #To find best assembly out of 10 (they might vary up to 5 contigs)

echo -e "\n\e[32m10 velvet iterations completed"
echo -e "\e[39mGenerating assemblies report..."

for i in $(seq 1 10)
do
  echo $i >> Assemblies && cat fix_$i/contigs.fa | \
  tee >(rg -j $(nproc) \> | wc -l >> Assemblies) | \
  rg -j $(nproc) -v \> | awk '{sum += length($1)}END{print sum}' >> Assemblies
done #Creating a small assembly report to consult it after the script is done
cat Assemblies | paste - - - | \
sed '1s/^/Assembly#\t#ofContigs\t#ofNucleotides\n/' > tmp && mv tmp Assemblies

echo  "Creating amos file of best assembly..."
best=$(rg -v '#' Assemblies | awk -F'\t' '{print $2}' | sort -n | head -1)
truebest=$(rg "$best" Assemblies | awk -F'\t' '{print $3}' | sort -n | tail -1)
amos=$(rg "$best" Assemblies | rg "$truebest"| awk -F'\t' '{print $1}' | uniq)
velvetg fix_$amos -exp_cov $6 -cov_cutoff $7 -amos_file yes

echo -e "\e[32mDone\n"
