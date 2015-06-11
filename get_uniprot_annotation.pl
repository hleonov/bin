use HTTP::Request::Common qw(POST);
use URI::URL;
use LWP; 


open (IN,"data/cath_35_list_swiss_prot.txt");
open (OUT,">data/cath_35_list_swiss_prot_location.txt");
while(<IN>){
  chop;
  $entry = $_;
  ($pdb, $swissprot) = split;
  if($swissprot ne "000000") {
    $location = &location($swissprot);
    print OUT "$entry $location\n";
    print STDERR "$entry $location\n";
  }
  else {
    print OUT "$entry $location\n";
    print STDERR "$entry $location\n";
  }
}
close(IN);



sub location {
  $browser = LWP::UserAgent->new;
  $website = "http://www.uniprot.org/uniprot/".$_[0].".txt";
  $req = POST $website;
  $response = $browser->request($req);
  die "Error: " , $response->status_line, "\n" unless $response->is_success; 
  $string = $response->as_string;
  $string =~ /DR   GO; (.+)/;
  $location = $1;
  unless ($location) {
    $location = "000000"
  }
  return $location;
}


sub location_old {
  $browser = LWP::UserAgent->new;
  $website = "http://www.uniprot.org/uniprot/".$_[0];
  $req = POST $website;
  $response = $browser->request($req);
  die "Error: " , $response->status_line, "\n" unless $response->is_success; 
  $string = $response->as_string;
  $string =~ /Subcellular location\<\/acronym\>\<\/td\>\<td\>\<p\>\<a href=\"\/locations\/..-....\"\>(.{1,10})\<\/a\>/;
  $location = $1;
  unless ($location) {
    $location = "000000"
  }
  return $location;
}
