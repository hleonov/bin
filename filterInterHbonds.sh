#rm [a-d].tmp
#cat output.tmp | tail +4 | awk '{if(1<=$10 && $10<=1397){print $0}}' | awk '{if(1>$12 || $12>1397){print $0}}' > a.tmp
#cat output.tmp | tail +4 | awk '{if(1398<=$10 && $10<=2865){print $0}}' | awk '{if(1398>$12 || $12>2865){print $0}}' > b.tmp
#cat output.tmp | tail +4 | awk '{if(2866<=$10 && $10<=4262){print $0}}' | awk '{if(2866>$12 || $12>4262){print $0}}' > c.tmp
#cat output.tmp | tail +4 | awk '{if(4263<=$10 && $10<=5730){print $0}}' | awk '{if(4263>$12 || $12>5730){print $0}}' > d.tmp

rm intra.tmp
cat output.out | tail +4 | awk '{if(1<=$10 && $10<=1397){print $0}}' | awk '{if(1>$12 || $12>1397){print $0}}' >> intra.tmp
cat output.out | tail +4 | awk '{if(1398<=$10 && $10<=2865){print $0}}' | awk '{if(1398>$12 || $12>2865){print $0}}' >> intra.tmp
cat output.out | tail +4 | awk '{if(2866<=$10 && $10<=4262){print $0}}' | awk '{if(2866>$12 || $12>4262){print $0}}' >> intra.tmp
cat output.out | tail +4 | awk '{if(4263<=$10 && $10<=5730){print $0}}' | awk '{if(4263>$12 || $12>5730){print $0}}' >> intra.tmp

