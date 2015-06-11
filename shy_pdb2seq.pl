#!/usr/bin/perl
%aa=(ALA,A,CYS,C,ASP,D,GLU,E,PHE,F,GLY,G,HIS,H,ILE,I,LYS,K,LEU,L,
        MET,M,ASN,N,PRO,P,GLN,Q,ARG,R,SER,S,THR,T,VAL,V,TRP,W,TYR,Y);
foreach $file (@ARGV){
  open(IN,"$file");
  print ">$file\n";
  while(<IN>){
    if(/^ATOM.+ CA /){
      @fields=split;
      $i++;
      print"$aa{$fields[3]}";
      if($i == 10){
        print " ";
      }
      elsif($i == 20){
        print " ";
      }
      elsif($i == 30){
        print "  ";
      }
      elsif($i == 40){
        print " ";
      }
      elsif($i == 50){
        print "\n";
        $i = 0;
      }
    }
  }
  print "\n";
  close(IN);
}
