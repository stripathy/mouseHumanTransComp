library(nnet)
library(magrittr)

# maps individual cells, in cells_to_map_df to transcriptomic clusters defined by dissoc_cell_clusters_expr
mapCellsToClusters = function(cells_to_map_df, # cells to be mapped
                              dissoc_cell_clusters_expr = tasic_expr_cluster,  # average expression of clusters to be mapped to
                              cell_type_markers, # which marker genes to use?
                              use_cluster_names = tasic_2018_inh_clusters # which clusters should be used?
){
  
  use_genes = Reduce(intersect, list(cell_type_markers, colnames(cells_to_map_df), colnames(dissoc_cell_clusters_expr)))
  
  expressed_genes = use_genes[apply(cells_to_map_df[, use_genes], 2, function(x) any(x>1))]
  
  use_genes = expressed_genes
  d = cells_to_map_df[, use_genes] %>% t()
  
  colnames(d) = rownames(cells_to_map_df)
  
  aibs_temp_df_t = dissoc_cell_clusters_expr[dissoc_cell_clusters_expr$cluster %in% use_cluster_names, use_genes] %>% t()
  colnames(aibs_temp_df_t) = dissoc_cell_clusters_expr$cluster[dissoc_cell_clusters_expr$cluster %in% use_cluster_names]
  
  m = cor(aibs_temp_df_t, d, method = 'spearman') 
  mappings_aibs = apply(m, 2, function(x) rownames(m)[which.is.max(x)])
  best_corrs = apply(m, 2, function(x) x[which.is.max(x)])
  
  mappings_aibs = data.frame(mappings_aibs, best_corrs)
  colnames(mappings_aibs) = c('cluster', 'cluster_correlation')
  
  return(mappings_aibs)
  
}