#!/usr/bin/perl

# Open the xvg file for reading.
$file = @ARGV[0];
open(IN,"$file");


# Read the data line by line from the input file in a "while" loop,
# working upon each line as it is read.
while( $_ = <IN> ) {

     # pattern matching: search for lines that start with a digit,
     # and a pattern of 2 decimal numbers (with a possible dot in the middle).
     # Save the first decimal in $time, and the second decimal in $value.
     if (m/\s+(\d+\.?\d*)\s+(.\d+\.?\d*)/) { 
	$x_var = $1;
        $y_var = $2;
      	print "$x_var $y_var\n";
     }
	
     #exit if ($counter == 10); 
}

close(IN);
