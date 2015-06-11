#!/usr/bin/perl
#$* = 1;     # enable multi-line munging

# this script reads pub med data in medline format saved as text 
# and then copied
# the citatin name are the 1st author without any spaces
# and all capital.
# the same is true for the journal

# there are two arguments in the command line
# the 1st is the file name
# the 2nd is the optional file name to append to 

$in = shift || die "Usage: perl $0 infile outfile\n";
$out = shift || die "Usage: perl $0 infile outfile\n";;

open(IN,$in) || die "Input file \"$in\" couldn't be opened\n";
open(OUT,">>$out") || die "Output file \"$out\" couldn't be opened\n";


while($line = <IN>){
  if($line =~ m/^\s*PMID\s*-.+/ && $first_record){# new record
    &analyze;
    $citation = "";
    undef(@authors);
    $author = 0;

  }
  else{
    $citation .= $line; # append the contents of the citation
	$first_record=1;
    if($line =~ m/^\s*AU\s*-\s*(.+)\n/){
      $author++;
      $authors[$author] = $1;
    }
  }
}
&analyze;
# now parse
sub analyze {
  $output = "";
  $i++;
  # remove the new line
  $citation =~ s/\n//g;
  # title
  $citation =~ m/\s*TI\s*-\s*(.+)PG\s*-\s*/;
  $title = $1;
  # remove the extra tabbing spaces
  $title =~ s/\t/ /g;
  # now the source
  $citation =~ m/\s*SO\s*-\s*(.+)/;
  $source = $1;
  # now munge the source
  $source =~ m/^(.+)\s*(\d{4}).*;(.+):(\d+-\d+)/;
  $journal = $1;
  $year = $2;
  $volume = $3;
  $pages = $4;
  # increase dash on pages
  $pages =~ s/-/--/;
  # $output .= the name of the citation
  # 1st remove spaces form 1s auhotrs name
  $first_author = $authors[1];
  $first_author =~ s/\s//g;
  # 2nd remove spaces form journal
  $short_journal_name = $journal;
  $short_journal_name =~ s/\s//g;
  $name = $first_author.$year.$short_journal_name;
  $output .= "\@article\{$name,\n";
  # authors
  $output .= "author = \"";
  for($j = 1; $j <= $#authors; $j++){
    # munge the authors name
    # grab the last element of the authors name 
    # hopefully it is the initials
    @author_parts  = split(/\s+/,$authors[$j]);
    @initials = split(//,$author_parts[$#author_parts]);
    # now $output .= everything except the initials
    for($k = 0; $k < $#author_parts; $k++){
      if($k > 0){
        $output .= " $author_parts[$k]";
      }
      else{
        $output .= "$author_parts[$k]";
      }
    }
    $output .= ", ";
    # now $output .= the initials
    for($l = 0; $l <= $#initials; $l++){
      if($l > 0){
        $output .= " $initials[$l]";
      }
      else{
        $output .= "$initials[$l]";
      }
    }
    unless($j == $#authors){
      $output .= " and ";
    }
  }
  $output .= "\",\n";
  # done with the authors
  $output .= "journal = \"$journal\",\n";
  $output .= "pages = \"$pages\",\n";
  $output .= "title = \"\{$title\}\",\n";
  $output .= "volume = \"$volume\",\n";
  $output .= "year = $year\};\n\n";
  if($out){
    print OUT $output;
  }
  else{
    print $output;
  }
}









