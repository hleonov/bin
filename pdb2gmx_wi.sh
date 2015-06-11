#!/bin/bash
#AUTHOR: Jan Henning Peters
#EMAIL: jpeters@gwdg.de
#DATE: 2. February 2010
#ARCH: bash
#REQUIRES: bash, his_from_whatif.py (in path)
#STATUSCOMMENTS: should work (but not extensively tested)
#DESCRIPTION: pdb2gmx_wi.sh runs pdb2gmx but uses the histidine-protonation states from WHATIF
#USAGE: Takes three input parameters: input file, output file, force-field
#EXAMPLE: pdb2gmx_wi.sh protein.pdb gmx.gro oplsaa

echo "getmol
$1

1
addhyd
tot 0
makmol
$1
whatif.pdb
tot 0

fullst
y" > wi_script
cp /usr/local/whatif-beta/dbdata/TOPOLOGY.H .

if [ -e whatif.pdb ]; then
  rm whatif.pdb
fi

DO_WHATIF.COM script wi_script &> whatif.log

rm wi_script

his_from_whatif.py whatif.pdb | pdb2gmx -f $1 -his -o $2 -ff $3 -ignh -water spce


