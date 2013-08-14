#!/bin/bash
FILES=/Users/anna/Research/LotsaGenomes/Chopped_10kb/*.fna
FOLDER=/Users/anna/Research/LotsaGenomes/Chopped_10kb
for f in $FILES
  do
    /Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 2 -i $f -d "${f%%.*}.db"
done
