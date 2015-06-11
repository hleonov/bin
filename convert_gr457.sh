cat $1 | sed 's/ -/-/g' |
awk '
  BEGIN{i=0}
  {if ($1~/IMPRO/) i=0;
   if (i) {
       print $1,9,$4,$3*4.1868/$2,$5 };
   if ($1~/DIHE/) i=1;
  }
' | sed 's/-/ /g' | awk '{printf "%-4s%-4s%-4s%-4s%3i%12.6f%12.6f  %i   ;add from JSL database\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > $name"bon.itp"


