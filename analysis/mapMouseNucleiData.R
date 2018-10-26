# since we're using the mouse primary motor cortex nuclei data, then we need to manually assign clusters to these data
# my approach is just to use a simple correlation based classifier - for each nucleus, ask which of the N mouse clusters the nucleus is most correlated with
# then we assign the nucleus that cluster identity

library(nnet)
library(magrittr)

source('R/mapCellsToClusters.R')

# cluster average expression based on mouse AIBS single cell visp and alm data
tasic_expr_cluster= readRDS('data-raw/allenMouse/aibs_mouse_expr_cluster.rda') # this has been pre-computed and stored within the repo

# neuroexpresso markers is the list of markers used in the neuroexpresso paper
source('data-raw/markers/NeuroExpressoMarkers.R')

# define tasic clusters based on clusters observed in mouse primary visual cortex data, seems like excluding low quality clusters is good?
# this excludes the low quality and high intonic clusters (so we don't want to map to those, because they're probably ok for nuclei-based data?)
# but keeps the doublet and batch clusters, so if any cells map to those, we probably don't want them
tasic_2018_high_qual_clusters = tasic_expr_cluster$cluster[!grepl('Low Qual|High', tasic_expr_cluster$cluster)] %>% as.character()

neuroExpressoMarkers = unlist(neuroExpressoMarkers)

# let's get class and subclass info by getting this from the mouse metadata files
allMouseMeta = bind_rows(read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_samples-columns.csv'), 
                         read_csv('data-raw/allenALMMouse/mouse_ALM_2018-06-14_samples-columns.csv'))
mouse_cluster_class_info = allMouseMeta %>% distinct(class, subclass, cluster)

# map nuclei to clusters defined by mouse VISp and ALM data
mapped_cells = mapCellsToClusters(mouseNucleiBothCPMLog%>% t(), 
                   dissoc_cell_clusters_expr = tasic_expr_cluster, 
                   cell_type_markers = neuroExpressoMarkers, # we'll just use the markers defined in the neuroexpresso paper
                   use_cluster_names = tasic_2018_high_qual_clusters)

mapped_cells$cluster = factor(mapped_cells$cluster, levels = unique(mouse_cluster_class_info$cluster))

mapped_cells$cluster_correlation %>% hist()

# since some cells show overall low correlation to their mapped cluster, I think we should just remove them
mapped_cells[mapped_cells$cluster_correlation < .3, 'cluster'] = "Low Quality"

# add class and subclass info onto cluster labels
mapped_cells = left_join(mapped_cells, mouse_cluster_class_info)

# bind mouse nuclei metadata file with our cluster mappings
mouseNucleiMeta = cbind(mouseNucleiMeta, mapped_cells)

# write cluster info to a file
saveRDS(mouseNucleiMeta,file= 'data-raw/allenMouseMOpNuclei/mouseNucleiMeta.rds')


# the below code maps mouse nuclei to human nuclei, if that's wanted
# 
# humanNeuroExpressoMarkers = mouse2human(neuroExpressoMarkers)
# 
# colnames(human_cluster_expr) = plyr::mapvalues(colnames(human_cluster_expr), mouse_human_homologous_genes$humanGene, mouse_human_homologous_genes$mouseGene)
# mapped_cells_human_clusters = mapCellsToClusters(mouseNucleiBothCPMLog%>% t(), 
#                                   dissoc_cell_clusters_expr = human_cluster_expr, cell_type_markers = neuroExpressoMarkers, 
#                                   use_cluster_names = humanMeta$cluster %>% unique())



