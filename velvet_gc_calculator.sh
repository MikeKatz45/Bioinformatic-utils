#!/bin/bash

#Thought for Velvet assembly output

totalcontigs=$(cat contigs.fa | grep \> | wc -l)
echo -e "#No. of contigs: $totalcontigs\n#GC counts in percentage:\n" \
        > gc_report

cat contigs.fa | grep \> > tmp
contigcounter=0

for header in $(cat tmp)
do
  contigcounter=$(expr "$contigcounter" + 1)

  formatheader=$(cat contigs.fa | grep "$header" | cut -d_ -f1-2 | cut -d\> -f2)

  totalnt=$(cat contigs.fa | sed -n '/'"$header"'/,/>/p' | grep -v \> | \
          awk '{sum += length($1)}'END'{print sum}')

  totalgc=$(cat contigs.fa | sed -n '/'"$header"'/,/>/p' | grep -v \> | \
          tr -d 'A' | tr -d 'T' | tr -d 'a' | tr -d 't' | \
          awk '{sum += length($1)}'END'{print sum}')

  gcpercent=$(bc -l <<< "scale=3; $totalgc*100/$totalnt")

  echo "$contigcounter- $formatheader: $gcpercent %" >> gc_report

done

rm tmp
