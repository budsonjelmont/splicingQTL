script='make_GWAS_analysis_files-PGC2020_allSNPs.py'
cmd="ml bedtools; ml python/3.7.3; python3 /sc/arion/projects/EPIASD/splicingQTL/analysis/QTL_GWAS_coloc/scripts/${script}"

submitjob 2 -P acc_pintod02b -q express -c 1 -n 1 -m 25 -J GWAS_$script ${cmd}

