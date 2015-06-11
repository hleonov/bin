#!/bin/csh -f

cd $1
rm *.po[0-9]*
rm *.pe[0-9]*
rm *.o[0-9]*
rm *.e[0-9]*
rm step*
rm core*
rm bench*
rm err*npme*
rm \#*
cd ..
