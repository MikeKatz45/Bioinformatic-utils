#!/bin/bash

IDS=$1
SS=$2
forward=$3
reverse=$4

name_f=$(basename "$forward" .fastq.gz)
name_r=$(basename "$reverse" .fastq.gz)

usage() {
	echo -e "\nUso: $(basename $0) [file_ids.txt] [velveth_Sequences_file] [forwardreads] [reversereads]\n"
}

if [ $# -le 3 ]
	then
		usage
			exit 1
fi
echo -e "\n\e[33mInicio: $(date)"

echo -e '\e[39m- Generando archivos temporales (1/6)...'
	zcat $forward > tmp1
	zcat $reverse > tmp2
		cat $IDS | grep \> | tr -d '\>' > idstmp

echo '- Extrayendo lecturas del archivo de secuencias (2/6)...'
	for i in $(cat idstmp)
	 do cat $SS | grep $i | head -1 | cut -d\> -f2 | cut -f1
	  done > idstruetmp
rm idstmp

echo '- Identificando secuencias a eliminar en las librerías de lecturas (3/6)...'
	for i in $(cat idstruetmp)
	 do cat tmp1 | grep $i | head -1 |
	  cut -d' ' -f1 && cat tmp2 | grep $i |
	   head -1 | cut -d' ' -f1
	    done | uniq -D | uniq > Elim_ids.txt
echo '- Secuencias a eliminar verificadas en ambas librerías (4/6)...'
rm idstruetmp

echo '- Eliminando lecturas (5/6)...'
	for i in $(cat Elim_ids.txt)
	 do sed -i '/'"$i"'/,+3d' tmp1
	  done

echo 'Primera libreria (Forward) editada exitosamente.'
echo '- Eliminando lecturas (5/6)...'
	for i in $(cat Elim_ids.txt)
	 do sed -i '/'"$i"'/,+3d' tmp2
	  done

echo 'Segunda librería (Reverse) editada exitosamente.'

echo '- Compresión de archivos (6/6)...'
	cat tmp1 | pigz --best --processes $(nproc) > $name_f.elim$(cat Elim_ids.txt | wc -l).fastq.gz
	cat tmp2 | pigz --best --processes $(nproc) > $name_r.elim$(cat Elim_ids.txt | wc -l).fastq.gz
rm tmp1 tmp2

echo -e "\e[33mFinalizado: $(date) \n"
echo -e '\e[32mTarea completada, generando reporte...'

count1=$(expr $(zcat $forward | wc -l) / 4)
count2=$(expr $(zcat $reverse | wc -l) / 4)
	start=$(cat $IDS | grep \> | wc -l)
	 end=$(cat Elim_ids.txt | wc -l)
	  minus=$(expr $(echo $start) - $(echo $end))
count3=$(expr $(zcat $name_f.elim$(cat Elim_ids.txt | wc -l).fastq.gz | wc -l) / 4)
count4=$(expr $(zcat $name_r.elim$(cat Elim_ids.txt | wc -l).fastq.gz | wc -l) / 4)

echo -e '\e[39m- Resumen del proceso: '
echo  "Lecturas presentes en las librerias originales: $count1 (Forward) $count2 (Reverse)"
echo  "Número de lecturas proveídas para eliminar: $start "
echo  "Lecturas removidas por repetidos y/o no presentes en ambas librerías de lecturas: $minus "
echo  "Total de lecturas eliminadas: $end "
echo -e "Lecturas presentes en las nuevas librerias: $count3 (Forward) $count4 (Reverse)\n"
