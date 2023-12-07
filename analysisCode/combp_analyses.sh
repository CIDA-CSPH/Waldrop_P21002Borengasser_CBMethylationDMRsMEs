# COMBP ANALYSIS
## load the bedfiles, ewas results as saved in combp_setup.R

# generate bed, prep for combp

## sort bed files with bedtools
FILES="bedfiles/outcome/*.bed"
for file in $FILES ; do
output="${file}_sorted.bed"
sortBed -i $file > $output
done

## Add comb-p headers
FILES="bedfiles/outcome/*_sorted.bed"
for file in $FILES ; do
output="${file}_forCombp.bed"
echo -e "chrom\tstart\tend\trawp" | cat - $file > $output
done

# run combp
## Run through combp package folder -->
## Package can be found at: "local_path" -->
FILES="bedfiles/outcome/*_forCombp.bed"
for file in $FILES ; do
output="${file}.combp_results"
comb-p pipeline -c 4 --seed 1e-1 --dist 750 -p $output --anno hg19 $file
done