#!/bin/bash

cd ${1}
for l in *.fa; do
  echo -e ${l%'.fa'}'\t'$1'/'$l
done
cd ..
