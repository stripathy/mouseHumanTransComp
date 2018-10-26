library('biomaRt')
library(tidyr)


# calculate average expression values per cell type for set of conserved mouse and human genes


library(gplots)
heatmap.2(cor_mat[valid_mouse_clusters,valid_human_clusters], trace = 'none', margins = c(12, 12))

# figure out best correlated mouse clusters to large human l2/3 cluster
cor_mat[valid_mouse_clusters, 'Exc L2-3 LINC00507 FREM3'] %>% sort()




use_genes = mouse_human_homologous_genes %>% filter(humanGene %in% dopa_genes_human) %>% pull(humanGene)

use_genes_mouse = mouse_human_homologous_genes %>% filter(humanGene %in% dopa_genes_human) %>% pull(mouseGene)


mouse_subset_expr = mouseCPMLog[use_genes_mouse, mouse_l23_samples] 
rownames(mouse_subset_expr) = plyr::mapvalues(rownames(mouse_subset_expr), rownames(mouse_subset_expr), use_genes)
mouse_subset_expr_gathered = mouse_subset_expr %>% t() %>% as.data.frame() %>% gather()
mouse_subset_expr_gathered$species = 'mouse'

use_genes = rownames(humanCPMLog)[rownames(humanCPMLog) %in% htr_genes_human]

m = humanCPMLog[use_genes, human_l23_samples] 
m_gathered = m %>% t() %>% as.data.frame() %>% gather()
m_gathered$species = 'human'

mouse_human_subset_comb = bind_rows(m_gathered, mouse_subset_expr_gathered)

mouse_human_subset_comb %>% ggplot(aes(y = value, x = key, color = species)) + 
  geom_boxplot(position = position_dodge(0.9), outlier.size = .5) +
  stat_summary(fun.y=mean, geom="point", position = position_dodge(0.9), shape = 1, size = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

mouse_human_subset_comb %>% ggplot(aes(x = value, y = key, fill = species)) + geom_density_ridges(scale = 1) + facet_wrap(~species) # + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))

mouse_human_subset_comb %>% ggplot(aes(x = key, y = value, color = species)) + geom_dotpl



