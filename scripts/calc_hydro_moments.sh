for g in */*hydro*dx; do ligzernike moments -i $g -f ${g%.*}_fp.txt -o ${g%.*}_moments.txt -N 15 --maxCap 20; done
