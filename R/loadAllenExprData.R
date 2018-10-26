# load allen mouse and human transcriptomic data 
# copies a file written by oganm: https://github.com/oganm/metaNeighborTrial/blob/master/processAllenData.R

dir.create('data-raw/allenHuman',showWarnings = FALSE)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694416044',destfile = 'data-raw/allenHuman/humanExp.zip',method = 'wget')

dir.create('data-raw/allenMouse',showWarnings = FALSE)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985',destfile = 'data-raw/allenMouse/mouseExp.zip')

dir.create('data-raw/allenALMMouse',showWarnings = FALSE)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/694413179',destfile = 'data-raw/allenALMMouse/mouseExp.zip')

# data from mouse primary motor cortex nuclei
# nuclei is good because it's a direct comparison to the human nuclei based data
dir.create('data-raw/allenMouseMOpNuclei', showWarnings = F)
download.file('http://celltypes.brain-map.org/api/v2/well_known_file_download/738608585',destfile = 'data-raw/allenMouseMOpNuclei/mouseExp.zip')

unzip('data-raw/allenHuman/humanExp.zip',exdir = 'data-raw/allenHuman')
unzip('data-raw/allenMouse/mouseExp.zip',exdir = 'data-raw/allenMouse')
unzip('data-raw/allenALMMouse/mouseExp.zip',exdir = 'data-raw/allenALMMouse')
unzip('data-raw/allenMouseMOpNuclei/mouseExp.zip',exdir = 'data-raw/allenMouseMOpNuclei')



library(readr)
library(edgeR)

### load human MTG nuclei data
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

### load mouse visual cortex dissociated cell data
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

rm(humanSum, humanCPM, humanCPMLog, humanIntron, humanIntronCPMLog, mouseSum, mouseCPM, mouseCPMLog, mouseIntronCPM, mouseIntronCPMLog, mouseIntrons)

### load mouse nuclei data
mouseNucleiIntrons = read_csv('data-raw/allenMouseMOpNuclei/mouse_MOp_nuclei_2018-10-04_intron-matrix.csv')
mouseNucleiExons = read_csv('data-raw/allenMouseMOpNuclei/mouse_MOp_nuclei_2018-10-04_exon-matrix.csv')
mouseNucleiGenes = read_csv('data-raw/allenMouseMOpNuclei/mouse_MOp_nuclei_2018-10-04_genes-rows.csv')
mouseNucleiMeta = read_csv('data-raw/allenMouseMOpNuclei/mouse_MOp_nuclei_2018-10-04_samples-columns.csv')
assertthat::assert_that(all(mouseNucleiExons$X1 == mouseNucleiIntrons$X1) & all(mouseNucleiExons$X1 == mouseNucleiGenes$gene)) 
# note that exon and intron gene names are gene symbols, not entrez ids

assertthat::assert_that(all(colnames(mouseNucleiExons) == colnames(mouseNucleiIntrons)))
assertthat::assert_that(all(colnames(mouseNucleiExons)[-1] == mouseNucleiMeta$seq_name)) # note that exon and intron 
# We used log2-transformed CPM of intronic plus exonic reads for both datasets.

mouseSum = mouseNucleiExons[,-1]+mouseNucleiIntrons[,-1]
mouseNucleiCPM = cpm(mouseSum)
mouseNucleiBothCPMLog = log2(mouseNucleiCPM +1)
rownames(mouseNucleiBothCPMLog) = mouseNucleiGenes$gene

colnames(mouseNucleiBothCPMLog) = mouseNucleiMeta$sample_name

saveRDS(mouseNucleiBothCPMLog,file= 'data-raw/allenMouseMOpNuclei/mouseNucleiBothCPMLog.rds')

rm(mouseNucleiExons, mouseNucleiIntrons, mouseNucleiIntrons, mouseSum, mouseNucleiCPM, mouseNucleiBothCPMLog)

