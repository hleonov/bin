#!/usr/bin/perl
# requires you to save the files as medline
# we need to fix the files since <cr> is put instead of \n 

@files = split(/\s+/,`ls /Users/isaiahar/Desktop/*.fcgi`);
foreach $i(@files){
  &FixLineEndings($i);
  $input = $i;
  $output = $i."output";
  system("/Users/isaiahar/Bin/pub_med2bib.pl $input $output");
  system("cat /Users/isaiahar/Documents/new.bib $output > /Users/isaiahar/Documents/temporary.bib");
  system("mv /Users/isaiahar/Documents/temporary.bib /Users/isaiahar/Documents/new.bib");
  system("rm $output");
  system("rm $input");
}



sub FixLineEndings{ # this removes a carrige return and puts a new line instead
  my $out = $_[0]."_a_temporary_name";
  open(IN,$_[0]);
  open(OUT,">$out");
  while(<IN>){
    s/\r/\n/g;
    print OUT;
  }
  close(IN);
  close(OUT);
  system("mv $out $_[0]");
}




