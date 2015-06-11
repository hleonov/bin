#!/usr/bin/perl -w
use FileHandle;

$completed = 0;
$perc_file = "nothing";

$file = shift || die "Usage: perl $0 fasta-file id_perc [perc_file] [completed(1)]\n";
$id_perc = shift || die "Usage: perl $0 fasta-file id_perc [perc_file] [completed(1)]\n";
$perc_file = shift;
$completed = shift;

#$perc_file = "/data/secreted/secreted.nr.log";

$LEN_RANGE = 0.75;
readFasta($file);
$total = (keys %seq_hash) -1;
if (-e $perc_file) {
	print "$perc_file exists\n";
	readPercFile($perc_file);
	print "finished reading $perc_file\n";
} 

unless ($completed)
 {	
	open(LOG,">>$file.log") || die "Cannot open log file $file.log\n$!\n";
	print "opened $file.log - starting run\n";
	LOG->autoflush(1);
	for my $i (0 .. $total-1) {
		for my $j ($i+1 .. $total) {
			my $seq1 = $seq_hash{$i2seq{$i}};
			my $seq2 = $seq_hash{$i2seq{$j}};

			next if (defined $id_mat[$i][$j]);
				
			print LOG "$i $i2seq{$i}\t$j $i2seq{$j}\t";
			if (length($seq1) < $LEN_RANGE*length($seq2) ||
				 length($seq2) < $LEN_RANGE*length($seq1)) {
					 $id_mat[$i][$j] = 0;
					 print LOG "0\n";
			} else {
				$id_mat[$i][$j] = run_clustal($seq1,$seq2);		
			}
			
		}
	}
	close(LOG);
}
open(OUT,">$file.$id_perc.$LEN_RANGE.nr") || die "Cannot open outfile $file.$id_perc.nr\n$!\n";;
print OUT ">$i2seq{0}\n$seq_hash{$i2seq{0}}\n";
for my $j (1 .. $total) {
	my $flag=1;
	for $i (0 .. $j-1) {
		$flag=0 if ($id_mat[$i][$j] > $id_perc);
	}
	print OUT ">$i2seq{$j}\n$seq_hash{$i2seq{$j}}\n" if ($flag);
}

close(OUT);



sub run_clustal {
	my $fh = new FileHandle ">$file.msa" || print LOG "Error opening CL\n";  		
	print_seq($_[0],'seq1',$fh);
	print_seq($_[1],'seq2',$fh);
	$fh->close;
	my $id = 0;
	`rm $file.aln`;
	system("~/Programs/clustalw1.83/clustalw $file.msa");
#	if (system("~/Programs/clustalw1.83/clustalw $file.msa")) {
#		print LOG "Error\n";
	if (-e "$file.aln") { #alignment file was created
		$id = `grep '*' $file.aln | fold -1 | grep '*' | wc -l`;
		$id =~ s/\s//g;
		$id = 100*$id/(length $_[0]);
		print LOG "$id\n";
	} else {
		print LOG "Error\n";
		$id = "Error";
	}	
	return $id;
}
sub print_seq {	
	my $len = length($_[0]);
	my $fh = $_[2];
	print $fh ">$_[1]\n";
	for (my $i=0; $i<=$len; $i+=60) {
		if ($len-$i+1 > 0) {
			print $fh substr($_[0],$i,60)."\n";
		} else {
			print $fh substr($_[0],$i)."\n";
		}
	}
}
sub readFasta {
	open(IN,$_[0]) || die "Cannot open infile $_[0]\n$!\n";;
	my $seq = undef;
	my $name;
	my $i = 0;
	while (<IN>) {
		chomp;
		if (m/^>(.+)/) {
			if (defined $seq) {
				$seq_hash{$name} = $seq;
				$seq2i{$name} = $i;
				$i2seq{$i} = $name;								
				$i++;
			}
			$name = $1;
			$seq = "";						
		} else {
			$seq .= $_;
		}	
	}
	close(IN);
	
	$seq_hash{$name} = $seq;
	$seq2i{$name} = $i;
	$i2seq{$i} = $name;					
}
sub readPercFile {
	open(IN,$_[0]) || die "Cannot open infile $_[0]\n$!\n";
	while (<IN>) {
		chomp;   
		@parts = split(/\s+/,$_);
		#my $seq1 = $seq_hash{$i2seq{$parts[0]}};
		#my $seq2 = $seq_hash{$i2seq{$parts[2]}};
		#if (length($seq1) < $LEN_RANGE*length($seq2) ||
		#	 length($seq2) < $LEN_RANGE*length($seq1)) {
		#			 $id_mat[$parts[0]][$parts[2]] = 0;
		if ($parts[4] ne 'Error') {
			$id_mat[$parts[0]][$parts[2]] = $parts[4];
		} else {
			$id_mat[$parts[0]][$parts[2]] = undef;
		}
	#	print "$parts[0] $i2seq{$parts[0]} ".length($seq1)." $parts[2] $i2seq{$parts[2]} ".length($seq2)." $id_mat[$parts[0]][$parts[2]]\n";
	#	exit if ($parts[2] == 10);
	}
	close(IN);
}
