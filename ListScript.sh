#!/bin/bash

cd ${1}
for l in *.fna; do
  echo -e ${l%'.fna'}'\t'$1'/'$l
done
cd ..
