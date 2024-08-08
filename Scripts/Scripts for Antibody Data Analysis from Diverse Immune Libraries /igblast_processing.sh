#!/bin/bash


awk -F "\t" '{print $49 "\t" $50}' $1 | sort -k 1,2 | uniq -c | awk '{if (length($2)!=0) print $1 "\t" $2 "\t" $3 "\t" length($2) "\t" length($3)}' |  grep -v "\\*" | grep -v cdr3 | sort -n -r > $1"_processed"