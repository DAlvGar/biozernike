for g in */*hydro*dx; do 
	ligzernike moments -i $g -f ${g%.*}_fp.txt -o ${g%.*}_moments.txt -N 15 --skipflip --skipnormalize --maxCap 100 --multiplier 10; done
