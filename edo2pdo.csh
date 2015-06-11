#!/bin/csh

foreach i (? 10 11)
if (-f $i/sam.edo ) then

rm -f mono?$i.pdo
touch monoa$i.pdo monob$i.pdo

cat > monoa$i.pdo << EOF
# UMBRELLA      3.0
# Component selection: 0 0 1
# nSkip 1
# Ref. Group 'Protein'
# Nr. of pull groups 1
EOF

echo $i | awk '{print "# Group 1 'POT'  Umb. Pos. "  (6.-$1*0.5) " Umb. Cons. 100" }' >> monoa$i.pdo

echo "#####" >> monoa${i}.pdo


cat > monob$i.pdo << EOF
# UMBRELLA      3.0
# Component selection: 0 0 1
# nSkip 1
# Ref. Group 'Protein'
# Nr. of pull groups 1
EOF

echo $i | awk '{print "# Group 1 'POT'  Umb. Pos. "  (6.-$1*0.5) " Umb. Cons. 100" }' >> monob$i.pdo

echo "#####" >> monob${i}.pdo


grep EDx $i/sam.edo | awk '((NR) %2) == 0 {print $0}' | awk -v x=$i '{print (0.2*t++) "  "  ($2-(6.-x*0.5))}' >> monoa${i}.pdo
grep EDx $i/sam.edo | awk '((NR+1) %2) == 0 {print $0}' | awk -v x=$i '{print (0.2*t++) "  "  ($2-(6.-x*0.5))}' >> monob${i}.pdo
endif
end
