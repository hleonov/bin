#!/bin/csh -f

foreach dir (k*)
	cd $dir
	echo $dir >> ../stat.txt
	perl ~/bin/umb_stat.pl >> ../stat.txt
	cd ..
end
 
