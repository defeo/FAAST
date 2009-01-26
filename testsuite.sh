#!/bin/bash

for (( i=1 ; $i<=$(wc -l < $1) ; i++ )) ;
do
    ./$3bench.sh $(head -n$i $1 | tail -n1) $2 ;
done
