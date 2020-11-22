beddirbase=/sc/hydra/projects/pintod02c/datasets-external/Roadmap_Epigenomics/ChromHMM/
beddir=E081E082_FB_merge #{E073_DLPFC,E081_fetalBrain_M,E082_fetalBrain_F} 
listoutdir=/sc/arion/projects/EPIASD/splicingQTL/analysis/GREGOR/bedindexes/
listoutfile=${beddir}.bed.list.txt

ls ${beddirbase}${beddir}/* > ${listoutdir}${listoutfile}

