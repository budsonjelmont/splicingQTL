library(factoextra)

infile = '/sc/orga/work/belmoj01/sqtl/rename_workdir/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.gz.qqnorm.AllCombined'
#infile = '/sc/orga/work/belmoj01/sqtl/rename_workdir/out-extra3-100kb-covar_clusters_ilen100kb_reads50_ratio0.01_perind.counts.idsync.gz.phen_AllCombined'

df = read.table(infile)

pca.nsc = prcomp(df[,5:ncol(df)],center=FALSE,scale=FALSE)
pca.sc = prcomp(df[,5:ncol(df)],center=TRUE,scale=TRUE)

# Get eigenvalues
eigenvals.nsc = get_eigenvalue(pca.nsc)
eigenvals.sc = get_eigenvalue(pca.sc)

# Cumsum of variance explained 
vars.nsc = apply(pca.nsc$x, 2, var)  
props.nsc = vars.nsc / sum(vars.nsc)
varexp.nsc = cumsum(props.nsc)

vars.sc = apply(pca.sc$x, 2, var)
props.sc = vars.sc / sum(vars.sc)
varexp.sc = cumsum(props.sc)

# Make scree plot
scree.nsc = fviz_eig(pca.nsc,ncp=20, orientation='vertical', choice='eigenvalue')
scree.sc = fviz_eig(pca.sc,ncp=20, orientation='vertical', choice='eigenvalue')

ggsave('screeplot_unscaled.png',scree.nsc)
ggsave('screeplot_scaled.png',scree.sc)
