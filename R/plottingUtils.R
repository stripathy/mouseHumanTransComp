
# given a list of genes, plot mouse and human gene expression as boxplots
plotHumanMouseBoxplots = function(human_gene_names, mouse_samples = mouse_l23_samples, human_samples = human_l23_samples,
                                  mouse_expr = mouseCPMLog,
                                  human_expr = humanCPMLog, 
                                  mouseGenesDesc = mouseNucleiGenes, 
                                  homologous_genes = mouse_human_homologous_genes){
  
  use_genes = homologous_genes %>% filter(humanGene %in% human_gene_names) %>% pull(humanGene)
  
  use_genes_mouse_entrez = homologous_genes %>% filter(humanGene %in% human_gene_names) %>% pull(mouseID)
  use_genes_mouse = mouseGenesDesc %>% filter(entrez_id %in% use_genes_mouse_entrez) %>% pull(gene)
  
  mouse_subset_expr = mouse_expr[use_genes_mouse, mouse_samples] 
  rownames(mouse_subset_expr) = plyr::mapvalues(rownames(mouse_subset_expr), rownames(mouse_subset_expr), use_genes)
  mouse_subset_expr_gathered = mouse_subset_expr %>% t() %>% as.data.frame() %>% gather()
  mouse_subset_expr_gathered$species = 'mouse'
  
  use_genes = rownames(human_expr)[rownames(human_expr) %in% human_gene_names]
  
  m = human_expr[use_genes, human_samples] 
  m_gathered = m %>% t() %>% as.data.frame() %>% gather()
  m_gathered$species = 'human'
  
  mouse_human_subset_comb = bind_rows(m_gathered, mouse_subset_expr_gathered)
  
  boxplot = mouse_human_subset_comb %>% ggplot(aes(y = value, x = key, color = species)) + 
    geom_boxplot(position = position_dodge(0.9), outlier.size = .3) +
    stat_summary(fun.y=mean, geom="point", position = position_dodge(0.9), shape = 1, size = 3) + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab('') + ylab('Gene expression (log2 CPM+1)')
  
  return(boxplot)
  
}


# given sample names from mouse and human cells, average expression per group and then merge by gene homology
compareMouseHumanExpr = function(mouse_cell_names, human_cell_names, 
                                 mouse_expression = mouseNucleiBothCPMLog, human_expression = humanBothCPMLog, 
                                 homologene_mapping = mouse_human_homologous_genes){
  
  mouse_expr = mouse_expression[intersect(homologene_mapping$mouseGene, rownames(mouse_expression)),mouse_cell_names] %>% rowMeans()
  human_expr = human_expression[intersect(homologene_mapping$humanGene, rownames(human_expression)),human_cell_names] %>% rowMeans()
  
  mouse_expr = mouse_expr %>% as.data.frame() %>% tibble::rownames_to_column(var = 'mouseGene') 
  colnames(mouse_expr)[2] = 'mouseExpr'
  
  human_expr = human_expr %>% as.data.frame() %>% tibble::rownames_to_column(var = 'humanGene') 
  colnames(human_expr)[2] = 'humanExpr'
  
  agg = merge(homologene_mapping, mouse_expr)
  agg = merge(agg, human_expr)
  
  return(agg)
}
