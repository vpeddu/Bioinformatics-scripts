# Moves and renames kallisto abundance files from their folders into one named abundance_tsv

for i in `find . -name *.tsv`; do  temp=`echo $i | cut -d / -f2`; newfilename="$temp"".tsv";echo $newfilename; echo $temp; cp $temp/abundance.tsv abundance_tsv/$newfilename; done
