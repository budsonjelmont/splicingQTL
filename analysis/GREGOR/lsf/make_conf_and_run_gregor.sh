# Parameters
##################################

snptag=_clumped_LD0.6_1500kb

datadir=/sc/arion/projects/EPIASD/splicingQTL/analysis/annotationBEDs/Walker2018_data/sQTLs/
confdir=${datadir}GREGORconf/

snpfile=uniqueLeadSNPuniquephenos-chrpos${snptag}.txt

bedtag=RBP_CLIPDB

bedlistpath=/sc/arion/projects/EPIASD/splicingQTL/analysis/GREGOR/bedindexes/${bedtag}.bed.list.txt
refdir=/sc/hydra/projects/pintod02c/reference-databases/GREGOR_1kg_r2_0.7/
r2thresh=0.7
ldwindow=1000000

restag=${snptag}_r2${r2thresh}_ldwin${ldwindow}

outdir=${datadir}GREGORres${restag}/${bedtag}/
minneighbornum=500
bedfilesorted=false
pop=EUR
topnbed=2
jobnum=10
batchtype=local

##################################
# Declare parameter map
declare -A parammap
parammap=([INDEX_SNP_FILE]=${datadir}${snpfile} [BED_FILE_INDEX]=${bedlistpath} [REF_DIR]=${refdir} [R2THRESHOLD]=${r2thresh} [LDWINDOWSIZE]=${ldwindow} [OUT_DIR]=${outdir} [MIN_NEIGHBOR_NUM]=${minneighbornum} [POPULATION]=${pop} [TOPNBEDFILES]=${topnbed} [JOBNUMBER]=${jobnum} [BATCHTYPE]=${batchtype} [BEDFILE_IS_SORTED]=${bedfilesorted})

# Create conf file
mkdir -p ${confdir}

confpath=${confdir}GREGOR${restag}.conf
rm -f ${confpath}

for param in ${!parammap[@]}
do
  echo "${param} = ${parammap[$param]}" >> ${confpath}
done

# Run GREGOR on newly created .conf file
command="ml perl;ml gregor/1.4.0;perl /hpc/packages/minerva-centos7/gregor/1.4.0/GREGOR/script/GREGOR.pl --conf ${confpath}"
echo $command
bsub -J GREGOR${restag} -P acc_pintod02b -q premium -n 1 -W 12:00 -R rusage[mem=50000] -o /sc/arion/projects/EPIASD/splicingQTL/lsf_output/GREGORenrichment${restag}_output-%J.txt ${command}
