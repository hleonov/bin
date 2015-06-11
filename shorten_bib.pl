#!/usr/bin/env perl
# this script will take a ..bbl file that is produced by BibTex and shorten it
# to conatin just the 1st author and no title. a stanrda file will look like:
$* = 1;     # enable multi-line munging


# \bibitem{ArbelyE2004JMolBiol}
# E~Arbely, Z~Khattari, G~Brotons, M~Akkawi, T~Salditt, and I~T Arkin.
# \newblock {A Highly Unusual Palindromic Transmembrane Helical Hairpin Formed by
#   SARS Coronavirus E Protein.}
# \newblock {\em J Mol Biol}, 341(3):769--79, 2004.



#print "\\section*\{References\}\n";

while(<>){
  unless(/\\begin\{thebibliography\}/ || /\\end\{thebibliography\}/){
    if(/^\n/ && $cite){ # a new bibitem but not the 1st one
      $cite .= $_;
      &analyze;
      $cite = "";
    }
    else{
      $cite .= $_;
    }
  }
  else{
    print;
  }
}

sub analyze{
  $counter++;
  #print "[$counter] ";
  $cite =~ s/\n/ /g; # remove the spaces
  $cite =~ /(\\bibitem\{.+\})( .+) \\newblock \{.+\} (\\newblock \{.+\}.+)/;
  $bibitem = $1;
  $authors = $2;
  $journal = $3;
  
  print "$bibitem\n"; 
  # now analyze the authors
  @authors = split(/,/,$authors);
  if($#authors > 0){
    $authors = $authors[0].+" et al.\\ ";
  }
  print "$authors\n";
  # removing the 
  
  print "$journal\n";


}




