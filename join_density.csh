foreach f (water_density*.xls)
   awk '{print $2}' $f > $f.col
end

paste *.col > density_results.txt
