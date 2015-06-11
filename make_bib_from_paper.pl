#!/usr/bin/perl
# this script takes a LaTeX file and scans all of its \cite commands.
# It then takes all of the citation form the major new.bib file
# and makes a small bib file containing just the citations used in the paper.
$latex_file = $ARGV[0];
$latex_file =~ /^(.+)\.tex/; # this is the assumed ending for a LaTeX file
$just_file = $1;
$output = $just_file.".bib";
if($ARGV[1]){
  $citation_command = $ARGV[1];
}
else{
  $citation_command = "cite";
}
open(IN,$latex_file);
open(OUT,">$output");
# becuase there may be multiline or mac saved formats...
while(<IN>){
  $text .= $_;
}
close(IN);
@citations_gross = split(/\\$citation_command\{/,$text);
for($i = 1; $i <= $#citations_gross; $i++){
  @trim_citations = split(/\}/,$citations_gross[$i]);
  $trim_citations = $trim_citations[0];
  # so $trim_citations[0] is the citation/s 
  # first clean it upclean it up
  $trim_citations =~ s/\s//g;
  # however if there is more than one then do the following
  if($trim_citations =~ /,.+/){
    @all_trim_citations = split(/,/,$trim_citations);
    for ($j = 0; $j <= $#all_trim_citations; $j++){
      # check if the citation apeared erlier
      # if it didn't then add it to a string
      unless($citations{$all_trim_citations[$j]}){
        push(@citations,$all_trim_citations[$j]);
	$citations{$all_trim_citations[$j]} = 1;
      }
    }
  }
  # only 1 citation
  else{
    unless($citations{$trim_citations}){
      push(@citations,$trim_citations);
      $citations{$trim_citations} = 1;
    }
  }  
}
# now we have an array with all of te citation
# sift thru new.bib and creat a new bib file
foreach $i (@citations){
  open(IN,"/Users/isaiahar/Documents/new.bib");
  $print = "false";
  while(<IN>){
    if(/\@.+\{$i,/){
      $print = "true";
    }
    if($print eq "true"){
      print OUT;
      if($print eq "true" && /\};/){
        $print = "false";
	print OUT "\n";
      }
    }
  }
  close(IN);
}
$number_citations = $#citations + 1;
print "Made file \"$output\" with $number_citations citations using \"$citation_command\" as the command\n";



