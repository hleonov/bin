#!/usr/bin/perl -w
use strict;

my $head = undef;
my $body = "";

while (<>) {
  chomp;
  if (m/^>/) {
    something($head, $body) if ($body);
    $body = "";
    $head = $_;
    $head =~ s/\t/ /g;
  } else {
    $body .= $_;
  }
}
something($head, $body);

sub something {
  my $head = shift;
  my $body = shift;
  my $space = (24-length($head));
  print $head;
  for (my $i=1; $i<=$space; $i++) {
  	print " "; 
  }
  print "$body\n";
    
}
