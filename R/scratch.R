getAvgMouseExpr = function(cell_type_name,
                           mouse_expression = mouseNucleiBothCPMLog, 
                           homologene_mapping = mouse_human_homologous_genes){
  
  mouse_cell_names = mouseMeta %>% filter(cluster == cell_type_name) %>% pull(sample_name)
  
  if (length(mouse_cell_names) < 2){
    return(NULL)
  }
  
  mouse_expr = mouse_expression[intersect(homologene_mapping$mouseGene, rownames(mouse_expression)),mouse_cell_names] %>% rowMeans()
  # human_expr = human_expression[intersect(homologene_mapping$humanGene, rownames(human_expression)),human_cell_names] %>% rowMeans()
  return(mouse_expr)
  
  
  mouse_expr = mouse_expr %>% as.data.frame() %>% tibble::rownames_to_column(var = 'mouseGene') 
  colnames(mouse_expr)[2] = cell_type_name
  
  human_expr = human_expr %>% as.data.frame() %>% tibble::rownames_to_column(var = 'humanGene') 
  colnames(human_expr)[2] = 'humanExpr'
  
  agg = merge(homologene_mapping, mouse_expr)
  agg = merge(agg, human_expr)
  
  return(agg)
}

getAvgHumanExpr = function(cell_type_name,
                           human_expression = humanBothCPMLog, 
                           homologene_mapping = mouse_human_homologous_genes){
  
  human_cell_names = humanMeta %>% filter(cluster == cell_type_name) %>% pull(sample_name)
  
  if (length(human_cell_names) < 2){
    return(NULL)
  }
  
  human_expr = human_expression[intersect(homologene_mapping$humanGene, rownames(human_expression)),human_cell_names] %>% rowMeans()
  # human_expr = human_expression[intersect(homologene_mapping$humanGene, rownames(human_expression)),human_cell_names] %>% rowMeans()
  return(human_expr)
}



valid_classes = c('GABAergic', 'Glutamatergic', 'Endothelial', 'Non-Neuronal', 'Non-neuronal')


valid_mouse_clusters = mouseNucleiMeta %>% group_by(cluster) %>% tally() %>% filter(n > 5) %>% pull(cluster)
mouse_clusters_avg_expr = lapply(valid_mouse_clusters, function(cluster_name) getAvgMouseExpr(cluster_name) ) %>% bind_cols() %>% as.data.frame()
colnames(mouse_clusters_avg_expr) = valid_mouse_clusters
rownames(mouse_clusters_avg_expr) = mouseCPMLog[intersect(mouse_human_homologous_genes$mouseGene, rownames(mouseCPMLog)),] %>% rownames()

valid_human_clusters = humanMeta %>% filter(class %in% valid_classes) %>% group_by(cluster) %>% tally() %>% filter(n > 5) %>% pull(cluster)
human_clusters_avg_expr = lapply(valid_human_clusters, function(cluster_name) getAvgHumanExpr(cluster_name) ) %>% bind_cols() %>% as.data.frame()
colnames(human_clusters_avg_expr) = valid_human_clusters
rownames(human_clusters_avg_expr) = humanBothCPMLog[intersect(mouse_human_homologous_genes$humanGene, rownames(humanBothCPMLog)),] %>% rownames()

human_cluster_expr = human_clusters_avg_expr %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'cluster')

mouse_human_avg_expr = merge(mouse_clusters_avg_expr %>% tibble::rownames_to_column(var = 'mouseGene') , mouse_human_homologous_genes)
mouse_human_avg_expr = merge(mouse_human_avg_expr, human_clusters_avg_expr %>% tibble::rownames_to_column(var = 'humanGene'), by = 'humanGene')

cor_mat = cor(mouse_human_avg_expr[, c(valid_mouse_clusters, valid_human_clusters)], method = 'pearson')    
