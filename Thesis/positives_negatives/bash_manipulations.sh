#!/bin/bash

#Example of manipulation for kidney renal cell carcinoma

#Download .tsv files from TCGA only with gene_symbol and mutation number
#Manually remove headers from files and add a line break at the end
cat *.tsv | tr '\t' ' ' > KIRC

#Find missing genes and manually add them with a value of zero for frequency
cat TCGA_geneset.txt | grep -v -f <(awk '{print $1}' KIRC) #prints missing
#For a lot of missing genes count the missing genes
#Make a file with line number equal to missing genes filled with zeros
for i in $(seq 1 number); do echo '0' >> tmp; done
#Paste missing genes with the zeros and append to cancer file
paste -d ' ' <(cat TCGA_geneset.txt | grep -v -f <(awk '{print $1}' KIRC)) <(cat tmp) >> KIRC

#Verify
diff <(cat TCGA_geneset.txt | sort) <(awk '{print $1}' KIRC | sort)

#Print frequencies, copy and paste them on corresponding cancer type
sort KIRC | awk '{print $2}' | tr '\n' ','
rm *.tsv && rm tmp
