#!/bin/bash
for infile in *.fna; do
  for k in 0{5..8}; do
    name=${infile%.*}
    outfile="Analyses/${name}_${k}_out.txt"
    ./kmerCount $k < $infile > $outfile
    echo "in=$infile  k=$k  out=$outfile"
  done
done
for infile in *.fna; do
  for k in 09 10 11; do
    name=${infile%.*}
    outfile="Analyses/${name}_${k}_out.txt"
    ./kmerCount $k < $infile > $outfile
    echo "in=$infile  k=$k  out=$outfile"
  done
done