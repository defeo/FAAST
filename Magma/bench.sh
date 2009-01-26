#!/bin/bash

magma=/users/crypto/defeo/beaufort/magma-2.11-2/magma;
output=/users/crypto/defeo/workspace/Artin-Schreier/tests;

(echo $4 "p =" $1 "d =" $2 "k =" $3 '@'$(hostname -s) 
$magma -b p:=$1 d:=$2 k:=$3 build_$4.mgm 
echo build_$4 done
) > $output/$4_$1\_$2\_$3.$(hostname -s).log
