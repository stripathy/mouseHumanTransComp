library(dplyr)
library(homologene)
library(ggplot2)
library(stringr)
library(devtools)
library(cowplot)

devtools::install_github("oganm/homologene")

library(homologene)

# reload saved intron-only expression data
humanCPMLog = readRDS(file = 'data-raw/allenHuman/humanIntron.rds')
mouseCPMLog = readRDS(file = 'data-raw/allenMouse/mouseIntron.rds')


# use ogan's homologene wrapper to find homologous genes btween mouse and human
mouse_human_homologous_genes = human2mouse(humanGenes$gene)

mouse_human_homologous_genes = mouse_human_homologous_genes %>% filter(mouseID %in% mouseGenes$gene_entrez_id)
intersect(mouseGenes$gene_entrez_id, mouse_human_homologous_genes$mouseID)

rownames(mouseCPMLog) = mouseGenes$gene_symbol
rownames(humanCPMLog) = humanGenes$gene


# given sample names from mouse and human cells, average expression per group and then merge by gene homology
compareMouseHumanExpr = function(mouse_cell_names, human_cell_names, 
                                 mouse_expression = mouseCPMLog, human_expression = humanCPMLog, 
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
# calculate average expression profiles for a small number of homologous cell types and plot 
# expression levels against eachother

human_micro_samples = humanMeta %>% filter(cluster == 'Micro L1-3 TYROBP') %>% pull(sample_name)
mouse_micro_samples = mouseMeta %>% filter(cluster == 'Microglia Siglech')  %>% pull(sample_name)

micro_agg = compareMouseHumanExpr(mouse_micro_samples, human_micro_samples)

# plot average expression profile of mouse genes vs human genes for microglia - this is Fig 8b from the Hodge preprint
ggplot(micro_agg, aes(x = mouseExpr, y = humanExpr)) + geom_point(alpha = .25) + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue') + 
  ggtitle('Microglia expression') + xlab('Mouse expr (log2 CPM+1)') + ylab('Human expr (log2 CPM+1)')

match(mouse_human_homologous_genes$humanGene, names(human_micro_expr))

## 
human_l3a_samples = humanMeta %>% filter(str_detect('Exc L3-4 RORB CARM1P1', cluster)) %>% pull(sample_name)
mouse_l3a_samples = mouseMeta %>% filter(str_detect('L5 IT VISp Hsd11b1 Endou', cluster))  %>% pull(sample_name)

l3a_agg_expr = compareMouseHumanExpr(mouse_l3a_samples, human_l3a_samples)

# - this is Fig 8a from the Hodge preprint
ggplot(l3a_agg, aes(x = mouseExpr, y = humanExpr)) + geom_point(alpha = .25) + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue')


human_l5_samples = humanMeta %>% filter(str_detect('Exc L4-5 FEZF2 SCN4B', cluster)) %>% pull(sample_name)
mouse_l5_samples = mouseMeta %>% filter(str_detect(cluster, 'L5 PT VISp'))  %>% pull(sample_name)

l5_agg_expr = compareMouseHumanExpr(mouse_l5_samples, human_l5_samples)

l5_agg_expr %>% ggplot(aes(x = mouseExpr, y = humanExpr)) + geom_point(alpha = .25) + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue')





human_sst_chodl = humanMeta %>% filter(str_detect('Inh L3-6 SST NPY', cluster)) %>% pull(sample_name)
mouse_sst_chodl = mouseMeta %>% filter(str_detect(cluster, 'Sst Chodl'))  %>% pull(sample_name)

sst_chodl_agg_expr = compareMouseHumanExpr(mouse_sst_chodl, human_sst_chodl)

l5_agg_expr %>% ggplot(aes(x = mouseExpr, y = humanExpr)) + geom_point(alpha = .25) + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue')



