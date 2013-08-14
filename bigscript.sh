#!/bin/bash

cd /Users/anna/Research/RAIphyCommandLine/
time (
for x in 1 3 5; do
    for n in 5 7 10 15 20; do
        for c in {0..99}; do
            super="/Users/anna/Research/LongContigsRecruit_${x}/${n}sp_10kbp_${c}/"
            dbin="${super}DB_input/"
            dbout="${super}db"
            seqs="${super}testseqs/"
            res="${super}res"
            ./raiphy -e .fna -m 2 -I $dbin -d $dbout
            ./raiphy -e .fna -m 0 -I $seqs -d $dbout -o $res
            rm -r $dbin
            rm -r $seqs
        done
    done
done
)
date
cd /Users/anna/Research/