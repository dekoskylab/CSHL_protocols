#!/bin/bash

FILE=$1
OUTFILE=$2

grep -v sequence $FILE | awk '{if(length($51)>5) print ">"$1 "\n" $2 }' > $OUTFILE"_nt"
