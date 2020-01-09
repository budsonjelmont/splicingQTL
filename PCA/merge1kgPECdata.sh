# reference: https://www.biostars.org/p/335605/

cd /sc/orga/projects/EPIASD/splicingQTL/PCA/1kg

# Get a list of all PLINK files
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;

sed -i 's/.bim//g' ForMerge.list ;

# Merge all projects into a single PLINK file

plink --merge-list ForMerge.list --out Merge ;

# Perform PCA

plink --bfile Merge --pca --pca-clusters


