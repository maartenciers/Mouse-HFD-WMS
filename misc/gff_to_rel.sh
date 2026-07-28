#!/bin/bash

indir=$1

for gff in ${indir}/*.gff; do

	echo $gff
      
	out=$(basename ${gff} .gff).rel
	
#	grep -P "\tCDS\t" "$gff" \
#		| sed 's/\t.*CDS//g' \
#		| sed 's/ID=//g' \
#		| sed 's/;.*;product=/\t/g' \
#		| sed 's/;protein_id=.*//g' \
#		| awk -F '\t' '{print $1 "\t" $7 "\t" $8 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' \
#		> ${indir}/${out}

        grep -P "\tCDS\t" "$gff" \
                > ${indir}/${out}

done
