#!/bin/csh -f

set file=$1

egrep -v "^ATOM[[:space:]]+[0-9]+[[:space:]]+.?M[CN]" $file > tmp.pdb
cat tmp.pdb | sed 's/OC1/O1 /' | sed 's/OC2/O2 /'
rm tmp.pdb
