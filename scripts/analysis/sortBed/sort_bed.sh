# Sort a file by first 2 columns
echo "Command: sort -k1,1 -k2,2n $1 > ${1%.*}.sorted.${1##*.}"
sort -k1,1 -k2,2n $1 > ${1%.*}.sorted.${1##*.}
