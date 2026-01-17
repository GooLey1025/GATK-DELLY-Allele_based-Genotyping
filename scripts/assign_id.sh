#!/bin/bash
TYPE=$1
awk -v TYPE=${TYPE} 'BEGIN{OFS="\t"}
/^#/ {print; next}
{ 
key =TYPE"-"$1"-"$2
count[key]++
$3=key"-"count[key]
print
}' $2 