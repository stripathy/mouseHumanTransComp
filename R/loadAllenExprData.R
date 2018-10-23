# load allen mouse and human transcriptomic data 
# copies a file written by oganm: https://github.com/oganm/metaNeighborTrial/blob/master/processAllenData.R

dir.create('data-raw/allenHuman',showWarnings = FALSE)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694416044',destfile = 'data-raw/allenHuman/humanExp.zip',method = 'wget')

dir.create('data-raw/allenMouse',showWarnings = FALSE)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985',destfile = 'data-raw/allenMouse/mouseExp.zip')


unzip('data-raw/allenHuman/humanExp.zip',exdir = 'data-raw/allenHuman')
unzip('data-raw/allenMouse/mouseExp.zip',exdir = 'data-raw/allenMouse')


library(readr)
library(edgeR)

humanIntrons = read_csv('data-raw/allenHuman/human_MTG_2018-06-14_intron-matrix.csv')
humanExons = read_csv('data-raw/allenHuman/human_MTG_2018-06-14_exon-matrix.csv')
humanGenes = read_csv('data-raw/allenHuman/human_MTG_2018-06-14_genes-rows.csv')
humanMeta = read_csv('data-raw/allenHuman/human_MTG_2018-06-14_samples-columns.csv')
assertthat::assert_that(all(humanExons$X1 == humanIntrons$X1) & all(humanExons$X1 == humanGenes$entrez_id))
assertthat::assert_that(all(colnames(humanExons) == colnames(humanIntrons)))
assertthat::assert_that(all(colnames(humanExons)[-1] ==humanMeta$sample_name))
# We used log2-transformed CPM of intronic plus exonic reads for both datasets.
humanSum = humanExons[,-1]+humanIntrons[,-1]
humanCPM = cpm(humanSum)
humanCPMLog = log2(humanCPM +1)
rownames(humanCPMLog) = humanGenes$gene

saveRDS(humanCPMLog,file= 'data-raw/allenHuman/humanAll.rds')

humanIntronCPM = cpm(humanIntrons)
humanIntronCPMLog = log2(humanCPM +1)
rownames(humanIntronCPMLog) = humanGenes$gene

saveRDS(humanIntronCPMLog,file= 'data-raw/allenHuman/humanIntron.rds')

mouseIntrons = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_intron-matrix.csv')
mouseExons = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_exon-matrix.csv')
mouseGenes = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_genes-rows.csv')
mouseMeta = read_csv('data-raw/allenMouse/mouse_VISp_2018-06-14_samples-columns.csv')
assertthat::assert_that(all(mouseExons$X1 == mouseIntrons$X1) & all(mouseExons$X1 == mouseGenes$entrez_id))
assertthat::assert_that(all(colnames(mouseExons) == colnames(mouseIntrons)))
assertthat::assert_that(all(colnames(mouseExons)[-1] == mouseMeta$sample_name))
# We used log2-transformed CPM of intronic plus exonic reads for both datasets.

mouseSum = mouseExons[,-1]+mouseIntrons[,-1]
mouseCPM = cpm(mouseSum)
mouseCPMLog = log2(mouseCPM +1)
rownames(mouseCPMLog) = mouseGenes$gene_symbol

saveRDS(mouseCPMLog,file= 'data-raw/allenMouse/mouseAll.rds')

mouseIntronCPM = cpm(mouseIntrons)
mouseIntronCPMLog = log2(mouseIntronCPM +1)
rownames(mouseIntronCPMLog) = mouseGenes$gene_symbol

saveRDS(mouseIntronCPMLog,file= 'data-raw/allenMouse/mouseIntron.rds')

rm(humanSum, humanCPM, humanCPMLog, humanIntron, humanIntronCPMLog, mouseSum, mouseCPM, mouseCPMLog, mouseIntronCPM, mouseIntronCPMLog)
