#!/bin/bash
#bump the assembler counter to the start of loop

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
rg -j $(nproc) \> $file | sort -n | uniq | sed 's/>\(.*\)/\t\1\t/g' > rm_ids
minus=$(cat rm_ids | wc -l)


echo '- Finding corresponding read IDs in Sequences file (2/4)...'
rg -j $(nproc) -f rm_ids $SS | cut -f1 | uniq | cut -d\> -f2 \
> tmp && mv tmp rm_ids

echo '- IDs found and duplicates ignored (3/4)...'

past_frm=$(basename $forward | cut -d{ -f2 | tr -d '.fastq.gz')
past_rrm=$(basename $reverse | cut -d{ -f2 | tr -d '.fastq.gz')
end=$(cat rm_ids.txt | wc -l)

forw=$forward
reve=$reverse

forward=$(basename $forward | cut -d{ -f1)
reverse=$(basename $reverse | cut -d{ -f1)

echo '- Deleting reads and compressing edited files (4/4)...'

paste <(zcat V03B_R1_trim_cor.fastq.gz) <(zcat V03B_R2_trim_cor.fastq.gz) | \
paste - - - - | rg -j $(nproc) -v -f rm_ids | \
awk -v FS="\t" -v OFS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | \
tee >(cut -f1-4 | tr '\t' '\n' | pigz --best --processes $(nproc) > ${forward}{$(expr $(echo $past_frm) + $(echo $end)).fastq.gz) \
| cut -f5-8 | tr '\t' '\n' | pigz --best --processes $(nproc) > ${reverse}{$(expr $(echo $past_rrm) + $(echo $end)).fastq.gz

forward=${forward}{$(expr $(echo $past_frm) + $(echo $end)).fastq.gz
reverse=${reverse}{$(expr $(echo $past_rrm) + $(echo $end)).fastq.gz

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
done

echo -e "\n\e[32m10 velvet iterations completed"
echo -e "\e[39mGenerating assemblies report..."

for i in $(seq 1 10)
do
  echo $i >> Assemblies && cat fix_$i/contigs.fa | \
  tee >(rg -j $(nproc) \> | wc -l >> Assemblies) | \
  rg -j $(nproc) -v \> | awk '{sum += length($1)}END{print sum}' >> Assemblies
done
cat Assemblies | paste - - - | \
sed '1s/^/Assembly#\t#ofContigs\t#ofNucleotides\n/' > tmp && mv tmp Assemblies

echo  "Creating amos file of best assembly..."
best=$(rg -v '#' Assemblies | awk -F'\t' '{print $2}' | sort -n | head -1)
truebest=$(rg "$best" Assemblies | awk -F'\t' '{print $3}' | sort -n | tail -1)
amos=$(rg "$best" Assemblies | rg "$truebest"| awk -F'\t' '{print $1}' | uniq)
velvetg fix_$amos -exp_cov $6 -cov_cutoff $7 -amos_file yes

echo -e "\e[32mDone\n"

