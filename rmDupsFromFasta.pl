#!/usr/bin/perl

use strict;

use Bio::SeqIO;
use Digest::MD5 'md5_hex';

my %digests;
my $in = Bio::SeqIO->new(-fh => \*ARGV);
my $out = Bio::SeqIO->new;

while (my $seq = $in->next_seq) {
    next if $digests{md5_hex($seq->seq)}++;
    $out->write_seq($seq);
}



