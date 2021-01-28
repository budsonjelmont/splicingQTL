annotdir=/sc/arion/projects/EPIASD/splicingQTL/ext_data/other/gencode/hg19/v33lift37/
annotgtf=gencode.v33lift37.annotation.gtf

echo "python3 collapse_annotation.py $annotdir/$annotgtf $annotdir/${annotgtf%.*}.geneModel.gtf"
python3 collapse_annotation.py $annotdir/$annotgtf $annotdir/${annotgtf%.*}.geneModel.gtf
