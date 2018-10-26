library(dplyr)
library(homologene)
library(ggplot2)
library(stringr)
library(devtools)
library(cowplot)
library(homologene)
library('biomaRt')
library(tidyr)
library(readr)

source('R/plottingUtils.R') # has some basic utility functions

## let's focus the initial analysis on human MTG nuclei data, Both refers to both introns and exons
humanBothCPMLog = readRDS(file = 'data-raw/allenHuman/humanAll.rds')
humanGenes = read_csv('data-raw/allenHuman/human_MTG_2018-06-14_genes-rows.csv')
humanMeta = read_csv('data-raw/allenHuman/human_MTG_2018-06-14_samples-columns.csv')

# mouse MOp nuclei data 
mouseNucleiMeta = readRDS('data-raw/allenMouseMOpNuclei/mouseNucleiMeta.rds') # this is the metadata file plus my cluster annotations
mouseNucleiGenes = read_csv('data-raw/allenMouseMOpNuclei/mouse_MOp_nuclei_2018-10-04_genes-rows.csv')
mouseNucleiBothCPMLog = readRDS(file = 'data-raw/allenMouseMOpNuclei/mouseNucleiBothCPMLog.rds')

# mouseBothCPMLog = readRDS(file = 'data-raw/allenMouse/mouseAll.rds') # this is the file from mouse cell VISp data

# use ogan's homologene wrapper to find homologous genes btween mouse and human
mouse_human_homologous_genes = human2mouse(humanGenes$entrez_id)

mouse_human_homologous_genes = mouse_human_homologous_genes %>% filter(mouseID %in% mouseNucleiGenes$entrez_id)

# calculate average expression profiles for a small number of homologous cell types and plot 
# expression levels against eachother


### plot average expression profile of mouse genes vs human genes for microglia 
# this is Fig 8b from the Hodge preprint (except here we're plotting mouse nuclei and both exons and introns)
human_micro_samples = humanMeta %>% filter(cluster == 'Micro L1-3 TYROBP') %>% pull(sample_name)
mouse_micro_samples = mouseNucleiMeta %>% filter(cluster == 'Microglia Siglech')  %>% pull(sample_name)

micro_agg = compareMouseHumanExpr(mouse_micro_samples, human_micro_samples)
cor(micro_agg$mouseExpr, micro_agg$humanExpr, method = 'pearson')

ggplot(micro_agg, aes(x = mouseExpr, y = humanExpr)) + geom_point(alpha = .25) + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue') + 
  ggtitle('Microglia expression') + xlab('Mouse expr (log2 CPM+1)') + ylab('Human expr (log2 CPM+1)')

## plot average expression profile of mouse genes vs human genes for l3a/l5 cells
# - this is Fig 8a from the Hodge preprint
human_l3a_samples = humanMeta %>% filter(str_detect('Exc L3-4 RORB CARM1P1', cluster)) %>% pull(sample_name)
mouse_l3a_samples = mouseNucleiMeta %>% filter(str_detect('L5 IT VISp Hsd11b1 Endou', cluster))  %>% pull(sample_name)

l3a_agg_expr = compareMouseHumanExpr(mouse_l3a_samples, human_l3a_samples)
cor(l3a_agg_expr$mouseExpr, l3a_agg_expr$humanExpr, method = 'pearson')

ggplot(l3a_agg_expr, aes(x = mouseExpr, y = humanExpr)) + geom_point(alpha = .25) + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue')  + 
  ggtitle('Pyr L3c/L5a expression') + xlab('Mouse expr (log2 CPM+1)') + ylab('Human expr (log2 CPM+1)')

### 

## plot the scatter plot for layer 2/3 pyramidal cells
human_l23_samples = humanMeta %>% filter(str_detect('Exc L2-3 LINC00507 FREM3', cluster)) %>% pull(sample_name)
mouse_l23_samples = mouseNucleiMeta %>% filter(subclass == 'L2/3 IT')  %>% pull(sample_name)

l23_agg_expr = compareMouseHumanExpr(mouse_l23_samples, human_l23_samples, 
                                     mouse_expression = mouseNucleiBothCPMLog, human_expression = humanBothCPMLog)

library(ggrepel)

l23_agg_expr_filt = l23_agg_expr %>% filter(mouseExpr > .1 | humanExpr > .1)
# - this is Fig 8a from the Hodge preprint
ggplot(l23_agg_expr_filt, aes(x = mouseExpr, y = humanExpr)) + 
  geom_point(alpha = .1) + 
  geom_point(data = subset(l23_agg_expr_filt, str_detect(humanGene, 'HTR[0-9]')), color = 'red', alpha = 1) +
  geom_text_repel(data = subset(l23_agg_expr_filt, str_detect(humanGene, 'HTR[0-9]')), aes(label = humanGene), color = 'red', size = 5)+
  # geom_point(data = subset(l23_agg_expr_filt, str_detect(humanGene, 'CHRN|CHRM')), color = 'green', alpha = 1) +
  # geom_text_repel(data = subset(l23_agg_expr_filt, str_detect(humanGene, 'CHRN|CHRM')), aes(label = humanGene), color = 'green', size = 5)+
  # geom_text_repel(aes(label = ifelse(str_detect(humanGene, 'CHRM'), as.character(humanGene), '')), color = 'blue') + 
  geom_abline(slope = 1) + 
  geom_abline(slope = 1, intercept = log2(10), color = 'blue') + 
  geom_abline(slope = 1, intercept = -log2(10), color = 'blue') + 
  xlab('mouse expression (log2 CPM+1)') + 
  ylab('human expression (log2 CPM+1)') + ggtitle('layer 2/3 pyramidal cells')


# define gene ontology groups for various gene classes of interest
htr_go_groups = c('GO:0004993', 'GO:0022850')
achr_go_groups= c('GO:0015464', 'GO:0016907', 'GO:0022848')
adren_go_groups = c('GO:0004935')
dopa_go_groups = c('GO:0004952')


# we'll use biomart to get genes for specific Gene Ontology terms
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

gene.data <- getBM(attributes=c('hgnc_symbol', 'go_id'),
                   filters = 'go', values = htr_go_groups, mart = ensembl) %>% filter(go_id %in%  htr_go_groups)
htr_genes_human = gene.data$hgnc_symbol
htr_genes_human = htr_genes_human[str_detect(htr_genes_human, 'HTR')]
# plot boxplot for mouse and human serotonergic receptor genes
plotHumanMouseBoxplots(htr_genes_human, mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

gene.data <- getBM(attributes=c('hgnc_symbol', 'go_id'),
                   filters = 'go', values = achr_go_groups, mart = ensembl) %>% filter(go_id %in%  achr_go_groups)
achr_genes_human = gene.data$hgnc_symbol
achr_genes_human = achr_genes_human[str_detect(achr_genes_human, 'CHR')]
plotHumanMouseBoxplots(achr_genes_human, mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

gene.data <- getBM(attributes=c('hgnc_symbol', 'go_id'),
                   filters = 'go', values = adren_go_groups, mart = ensembl) %>% filter(go_id %in%  adren_go_groups)
adren_genes_human = gene.data$hgnc_symbol
adren_genes_human = adren_genes_human[str_detect(adren_genes_human, 'ADR')]
plotHumanMouseBoxplots(adren_genes_human, mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

gene.data <- getBM(attributes=c('hgnc_symbol', 'go_id'),
                   filters = 'go', values = dopa_go_groups, mart = ensembl) %>% filter(go_id %in%  dopa_go_groups)
dopa_genes_human = gene.data$hgnc_symbol
dopa_genes_human = dopa_genes_human[str_detect(dopa_genes_human, 'DRD')]
plotHumanMouseBoxplots(dopa_genes_human, mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()


# do same as above for VIP cells
human_vip = humanMeta %>% filter(str_detect(cluster, 'VIP')) %>% pull(sample_name)
mouse_vip = mouseNucleiMeta %>% filter(subclass == 'Vip')  %>% pull(sample_name)

plotHumanMouseBoxplots(htr_genes_human, mouse_samples = mouse_vip, human_samples = human_vip,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

plotHumanMouseBoxplots(achr_genes_human, mouse_samples = mouse_vip, human_samples = human_vip,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

plotHumanMouseBoxplots(adren_genes_human, mouse_samples = mouse_vip, human_samples = human_vip,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

plotHumanMouseBoxplots(dopa_genes_human, mouse_samples = mouse_vip, human_samples = human_vip,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog) %>% print()

# do same as above for layer 1 interneurons cells (LAMP5 or PAX6 in humans, Lamp5 in mouse)

human_lamp5 = humanMeta %>% filter(str_detect(cluster, 'LAMP5|PAX6')) %>% pull(sample_name)
mouse_lamp5 = mouseNucleiMeta %>% filter(subclass == 'Lamp5')  %>% pull(sample_name)

plotHumanMouseBoxplots(htr_genes_human, mouse_samples = mouse_lamp5, human_samples = human_lamp5,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog, 
                       mouseGenesDesc = mouseNucleiGenes) %>% print()

plotHumanMouseBoxplots(achr_genes_human, mouse_samples = mouse_lamp5, human_samples = human_lamp5,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog, 
                       mouseGenesDesc = mouseNucleiGenes) %>% print()

plotHumanMouseBoxplots(adren_genes_human, mouse_samples = mouse_lamp5, human_samples = human_lamp5,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog, 
                       mouseGenesDesc = mouseNucleiGenes) %>% print()

plotHumanMouseBoxplots(dopa_genes_human, mouse_samples = mouse_lamp5, human_samples = human_lamp5,
                       mouse_expr = mouseNucleiBothCPMLog, human_expr = humanBothCPMLog, 
                       mouseGenesDesc = mouseNucleiGenes) %>% print()


