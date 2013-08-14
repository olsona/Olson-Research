#!/bin/sh

#  graphviz.sh
#  
#
#  Created by Anna Olson on 7/5/13.
#

for i in 25 50 100 200 500; do
    time (
    dot -Tps2 /Users/anna/Research/ClusteringTest/out$i.gv -o /Users/anna/Research/ClusteringTest/graph$i.ps
    echo "Done with $i"
    )
done