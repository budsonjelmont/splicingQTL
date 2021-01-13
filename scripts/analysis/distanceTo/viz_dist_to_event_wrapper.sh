# To plot distance from QTL SNP to intron using default axis limits:

Rscript viz_dist_to_event.R /sc/arion/projects/EPIASD/splicingQTL/output/sqtls/minCovars+seqPC9/4genoPCs/35HCPs/qtls/qtls.txt dist slope npval --logit --xlab "Distance to intron (kb)" --ylab "Effect size"  --colorlab "log10(pval)" --colorlo mediumpurple1 --colorhi purple4
