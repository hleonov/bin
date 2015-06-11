#!/bin/bash

export GMXLIB=/home/hleonov/Programs/gromacs407/share/gromacs/top_amber99sb_gmx40

echo "setup.sh log" > setup.log
echo "Call: " $0 $* >> setup.log

#cp ../../mdp/*.mdp .

tar -xzf /home/hleonov/Programs/gromacs407/share/gromacs/top/ffamber99sb_mut.tgz

#echo "Protonating structure." | tee -a setup.log
#pdb2gmx -f $1 -o gmx.pdb -ff oplsaa &>> setup.log

echo "Converting structure to amber format." | tee -a setup.log
# ~/grubiscripts/python/make_amber.py -f $1 -o amber.gro &>> setup.log
python ~/bin/make_amber.py -f $1 -o amber.gro &>> setup.log
mv amber.gro amber.pdb

echo "Mutating structure." | tee -a setup.log
#echo "$2
#$3
#n" | mutate.py -f amber.gro -o mutated.gro -ff amber99sb &>> setup.log

echo "$2
$3
n" | python ~/bin/mutate_beta.py -f amber.pdb -o mut_$2$3.pdb &>> setup.log
# -script mutations.txt


echo "Creating topology." | tee -a setup.log
$gmx407/pdb2gmx -f mut_$2$3.pdb -o mut_$2$3\_gmx.pdb -ff amber99sb -water spce &>> setup.log

echo "defining Box." | tee -a setup.log
$gmx407/editconf -f mut_$2$3\_gmx.pdb -o mut_$2$3\_gmx.pdb -d 1.2 -bt dodecahedron  &>> setup.log

echo "Solvating in SPCE." | tee -a setup.log
$gmx407/genbox -cp mut_$2$3\_gmx.pdb -cs spc216.gro -o mut_$2$3\_water.pdb -p topol.top &>> setup.log

echo "Adding ions." | tee -a setup.log
$gmx407/grompp -f ~/Data/em.mdp -c mut_$2$3\_water.pdb -o water.tpr &>> setup.log

echo 13 | $gmx407/genion -s water.tpr -o mut_$2$3\_ions.pdb -p topol.top -pname NaJ -nname ClJ -conc 0.15 -neutral &>> setup.log

echo "Modifying topology for B-state." | tee -a setup.log
#make_bstate.py -p topol.top -ff amber99sb &>> setup.log
python ./make_bstate.py -p topol.top &>> setup.log
mv newtop.top.ndx newtop.top
mkdir lambda0
mkdir lambda1

echo "Setting up lambda=0 energy minimization." | tee -a setup.log
$gmx407/grompp -f ./mdp/em_lambda0.mdp -c mut_$2$3\_ions.pdb -p newtop.top -o lambda0/em.tpr -maxwarn 1 &>> setup.log

echo "Running lambda=0 energy minimization." | tee -a setup.log
cd lambda0
mpirun -c 4 mdrun -v -s em.tpr &>> em.log
cd ..

echo "Setting up lambda=1 energy minimization." | tee -a setup.log
$gmx407/grompp -f ./mdp/em_lambda1.mdp -c mut_$2$3\_ions.pdb -p newtop.top -o lambda1/em.tpr -maxwarn 1 &>> setup.log

echo "Running lambda=1 energy minimization." | tee -a setup.log
cd lambda1
mpirun -c 4 mdrun -v -s em.tpr &>> em.log
cd ..

mkdir lambda0/equil
mkdir lambda1/equil

echo "Setting up equilibration." | tee -a setup.log
### PR 
$gmx407/grompp -f ./mdp/equil_lambda0.mdp -c lambda0/confout.part0001.gro -p newtop.top -o lambda0/equil/equil.tpr -maxwarn 2 &>> setup.log
cd lambda0/equil
gsub -f equil.tpr -N lmb0_eq -r pr
cd ../..

$gmx407/grompp -f ./mdp/equil_lambda1.mdp -c lambda1/confout.part0001.gro -p newtop.top -o lambda1/equil/equil.tpr -maxwarn 2 &>> setup.log
cd lambda1/equil
gsub -f equil.tpr -N lmb1_eq -r pr
cd ../..

### Equilibration for A and B states - md without PR, 10ns?
cd lambda0/equil
$gmx407/grompp -f ../../mdp/crooks_equilibration_stateA.mdp -c pr.part0001.gro -p ../../newtop.top -o equil2.tpr -maxwarn 2
gsub -f equil2.tpr -N lmb0_eq2 -r eq
cd ../../

cd lambda1/equil
$gmx407/grompp -f ../../mdp/crooks_equilibration_stateB.mdp -c md.part0001.gro -p ../../newtop.top -o equil2.tpr -maxwarn 2
gsub -f equil2.tpr -N lmb1_eq2 -r eq
cd ../../

### crooks runs
cp ./lambda0/equil/eq.part0001.gro ./lambda0
cp ./lambda1/equil/eq.part0001.gro ./lambda1
python ~/bin/prepare_crooks_runs.py -d lambda0 -top ./newtop.top -mdp ./mdp/crooks_non_eq_50ps_lambda0.mdp -sw_time 200 -skip 50
python ~/bin/prepare_crooks_runs.py -d lambda1 -top ./newtop.top -mdp ./mdp/crooks_non_eq_50ps_lambda1.mdp -sw_time 200 -skip 50

# submit all the stuff
foreach S( lambda0/morphes lambda1/morphes )
 cd $S
 foreach fr (frame1??)
   cd $fr
   gsub -d 1 -N $fr -f topol.tpr
   cd ../
 end
 cd ../../
end

# and finally analyze the dgdl.xvg files ( or dhdl.xvg of you used gmx 4.5 )

python ~/software/pymacs_beta/scripts/analyze_crooks.py -pa runA/morphes/frame*/dgdl.xvg -pb runB/morphes/frame*/dgdl.xvg

# writes some plots and results.dat with dG and error estimation

