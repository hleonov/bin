/***************************************************************************************
keepbyz.c - A Script to remove water molecules added by genbox to the lipids.
Version 2. [changed Chris Neale's version that appears online in the Gromacs website]

This version accepts a solvent file in PDB format (instead of gro), and the Z limits
in integer format. It will output all solvent molecules with the water oxygen below
the minimum Z or above the maximum Z, thus removing water molecules from  the lipid itself.
Like the first version, it needs the number of atoms per water molecule in our solvent. (SOL = 3)

Usage: keepbyz_pdb new_waters.pdb 3 38 83 > water_to_keep.pdb
****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#define LINE_SIZE 1000

void showUsage(const char *c){
    printf("Usage: %s <solvent.gro> <numAtomsInWaterMolecule> <minZ> <maxZ>\n",c);
	printf("                                (e.g. 3 for TIP3P, 4 for TIP4P)\n");
}

int main(int argn, char *args[]){
    FILE *mf;
	char *c;
	char *begin;
	char linein[LINE_SIZE];
	float mx,my,mz, b1, b2;
	int count,keep;
	int atomInWater, maxZ, minZ;
	
	if(argn!=5){
	   showUsage(args[0]);
	   exit(1);
	}
	if((sscanf(args[2],"%d",&atomInWater)!=1) || (sscanf(args[3], "%d", &minZ) != 1) ||(sscanf(args[4], "%d", &maxZ) != 1) ){
	   printf("error: unable to determine number of atoms in the water model, or Z limits\n");
       showUsage(args[0]);
       exit(1);
	}
	mf=fopen(args[1],"r");
	if(mf==NULL){
	   printf("error: unable to open the membrane file %s\n",args[1]);
	   exit(1);
	}
	
	count=0;
	keep=0;
	while(fgets(linein,LINE_SIZE,mf)!=NULL){
	   if(count==atomInWater){
	       count=0;
		   keep=0;
	   }
	   count++;
	   c=&linein[strlen(linein)-38];           	   // would use -24 for gro without velocities or -48 with velocities. -38 for PDB
	   if(sscanf(c,"%f %f %f %f %f",&mx,&my,&mz, &b1, &b2)!=5) continue;
	   if(count==1){									//only check oxygen is OK, then keep the whole molecule.
          //Add selection terms here to set keep=1 when you want to keep that particular solvent
          if((mz>maxZ) || (mz<minZ)){
             keep=1;
			 //printf("%f\n", mz);
          }
       }
	   if(keep) printf("%s",linein);
   } 
  fclose(mf);
}
