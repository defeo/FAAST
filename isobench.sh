#!/bin/bash

bin=/users/crypto/defeo/cryptodata/workspace/Artin-Schreier/bin/grappe/test$4;
output=/users/crypto/defeo/cryptodata/workspace/Artin-Schreier/tests;

(echo isoNTL "p =" $1 "d =" $2 "k =" $3 '@'$(hostname -s) 
echo $1 $2 $3 | $bin
echo isoNTL done
) > $output/isoNTL$4_$1\_$2\_$3.$(hostname -s).log
