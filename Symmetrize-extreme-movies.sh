#!/bin/bash
# Usage: symmetrize-extreme-movies.sh extreme.pdb

in=$1

csplit -n 4 --silent $in '/TITLE/' '{*}'
rm -f xx0000

N=$(ls xx????|wc -l)
echo created $N xx-files

out=${in%%.pdb}-sym.pdb
{
for i in `seq 1 $N`; do
    nr=$(printf "%04d\n" $i)
    cat xx$nr
done 
for i in `seq 1 $N | tac`; do
    nr=$(printf "%04d\n" $i)
    cat xx$nr
done 
}> $out

rm -f xx????
