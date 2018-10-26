# look into inter-individual variability

library(limma)
library(Glimma)
library(edgeR)

samples = humanMeta[human_l23_samples, ] %>% 
  dplyr::select(donor, sex, sample_name) %>% 
  mutate(donor = as.factor(donor), sex = as.factor(sex), sample_name = as.factor(sample_name))

x = DGEList(humanSum[, human_l23_samples], 
            samples =  humanMeta[humanMeta$sample_name %in% human_l23_samples, ] %>% 
              mutate_if(., is.character, as.factor),
            , genes = humanGenes)

keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

cpm_l23 <- cpm(x)
lcpm_l23 = cpm(x, log=TRUE, prior.count=2)
keep.exprs <- rowSums(cpm_l23>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x$samples$sex

donor = x$samples$donor
sex = x$samples$sex
design <- model.matrix(~ sex)
colnames(design) <- gsub("sex", "", colnames(design))

v <- voom(x, design, plot=TRUE)

dupfit <- duplicateCorrelation(v, design, block=donor)

# colnames(design) <- gsub("donor", "", colnames(design))
fitRan <- lmFit(v,design,block=donor,correlation=dupfit$consensus)
fitRan <- eBayes(fitRan)
topTable(fitRan,coef=c(2,3)) ## for testing the effect of factor B, which has three levels.

contr.matrix <- makeContrasts(
  sex = M, 
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

m_v_f_l23 = topTreat(tfit, coef = 1, n = Inf)

# try to replicate some of the plots using group plots

library(gplots)
m_v_f_l23.topgenes <- m_v_f_l23$entrez_id[1:1000]
i <- which(v$genes$entrez_id %in% m_v_f_l23.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm_l23[i, ], scale="row",
          labRow=v$genes$gene[i], labCol=sex, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

rownames(humanMeta) = humanMeta$sample_name

use_genes = humanGenes %>% filter(entrez_id %in% mouse_human_homologous_genes$humanID) %>% pull(gene)

use_genes = c(use_genes, 'XIST')

use_genes = c('XIST', 'KDM5D', 'RPL26', 'BEX5',  'KCNQ1')

human_l23_samples = humanMeta %>% filter(str_detect(cluster, 'Exc')) %>% pull(sample_name)
human_l23_samples_split_1 = humanMeta[humanMeta$sample_name %in% human_l23_samples, ] %>% 
  group_by(donor) %>% sample_frac(size = .5) %>% pull(sample_name)

# human_l23_samples_split_1 = sample(human_l23_samples, length(human_l23_samples) * .5)
human_l23_samples_split_2 = setdiff(human_l23_samples, human_l23_samples_split_1)

# use_donors = top_gene_mat_1 %>% dplyr::select(donor) %>% group_by(donor) %>% tally() %>% filter(n > 1) %>% pull(donor)

top_gene_mat_1 = cbind(humanMeta[human_l23_samples_split_1, ] %>% mutate_if(., is.character, as.factor), 
                     humanBothCPMLog[use_genes, human_l23_samples_split_1] %>% t()) 

colnames(top_gene_mat_1) = make.names(colnames(top_gene_mat_1))

top_gene_mat_2 = cbind(humanMeta[human_l23_samples_split_2, ] %>% mutate_if(., is.character, as.factor), 
                       humanBothCPMLog[use_genes, human_l23_samples_split_2] %>% t()) 

colnames(top_gene_mat_2) = make.names(colnames(top_gene_mat_2))

use_genes = use_genes %>% make.names()
# 
# m1 = glmer(sex ~ (1|donor), data = top_gene_mat, family = 'binomial')
# m2 = glmer(sex ~ XIST + (1|donor), data = top_gene_mat, family = 'binomial')

library(parallel)
pvals_1 = mclapply(use_genes, function(gene){
  f1 = paste(gene, '~ (1|donor) + (1|cluster)')
  f2 = paste(gene, '~ sex + (1|donor) + (1|cluster)')
  m1 = lmer(formula = f1, data = top_gene_mat_1)
  m2 = lmer(formula = f2, data = top_gene_mat_1)
  return(tidy(anova(m1, m2))[2, 'p.value'] %>% as.numeric())
}, mc.cores = 40)

df_1 = cbind(gene = use_genes, pval = pvals_1 %>% unlist %>% as.numeric()) %>% as.data.frame
df_1$pval <- as.numeric(as.character(df_1$pval))


library(parallel)
pvals_2 = mclapply(use_genes, function(gene){
  f1 = paste(gene, '~ (1|donor) + (1|cluster)')
  f2 = paste(gene, '~ sex + (1|donor) + (1|cluster)')
  m1 = lmer(formula = f1, data = top_gene_mat_2)
  m2 = lmer(formula = f2, data = top_gene_mat_2)
  return(tidy(anova(m1, m2))[2, 'p.value'] %>% as.numeric())
}, mc.cores = 40)

df_2 = cbind(gene = use_genes, pval = pvals_2 %>% unlist %>% as.numeric()) %>% as.data.frame
df_2$pval <- as.numeric(as.character(df_2$pval))

m =merge(df_1 , df_2, by = 'gene') %>% mutate(fdr.x = p.adjust(pval.x, method = 'BH'), 
                                              fdr.y = p.adjust(pval.y, method = 'BH'))

m %>% arrange(-log10(pval.x) * log10(pval.y)) %>% head(50)

anova(m1, m2)

diff_e = top_gene_mat %>% group_by(donor, sex) %>% mutate(n = n()) %>% filter(n > 1) %>% summarize_at(vars(use_genes), funs(mean)) %>% 
  ungroup() %>% 
  gather(key = 'gene', value = 'expr', use_genes) %>%
  group_by(gene, sex) %>%
  summarise(expr = list(expr)) %>% 
  spread(sex, expr) %>% 
  group_by(gene) %>% 
  mutate(p_value = t.test(unlist(F), unlist(M))$p.value,
         t_value = t.test(unlist(F), unlist(M))$statistic, 
         fold_change = mean(unlist(F)) - mean(unlist(M)), 
         mean_M = mean(unlist(M)), mean_F = mean(unlist(F)))


use_genes = mouseGenes %>% filter(gene_entrez_id %in% mouse_human_homologous_genes$mouseID) %>% pull(gene_symbol)

use_genes = c(use_genes, 'XIST')

use_genes = c('XIST', 'KDM5D', 'RPL26', 'BEX5',  'KCNQ1')

mouse_cells_l23_samples = mouseMeta %>% filter(cluster == 'L2/3 IT VISp Agmat')  %>% pull(sample_name)

mouse_top_gene_mat = cbind(mouseMeta[mouseMeta$sample_name %in% mouse_cells_l23_samples, ], mouseBothCPMLog[use_genes, mouse_cells_l23_samples] %>% t())

diff_e_mouse = mouse_top_gene_mat %>% group_by(donor, sex) %>% mutate(n = n()) %>% filter(n > 5) %>% summarize_at(vars(use_genes), funs(mean)) %>% 
  ungroup() %>% 
  gather(key = 'gene', value = 'expr', use_genes) %>%
  group_by(gene, sex) %>%
  summarise(expr = list(expr)) %>% 
  spread(sex, expr) %>% 
  group_by(gene) %>% 
  mutate(p_value = t.test(unlist(F), unlist(M))$p.value,
         t_value = t.test(unlist(F), unlist(M))$statistic, 
         fold_change = mean(unlist(F)) - mean(unlist(M)), 
         mean_M = mean(unlist(M)), mean_F = mean(unlist(F)))

propnnz = function(vector){
  return(sum(vector > 0) / length(vector))
}

DotPlot <- function(
  object,
  genes.plot,
  cols.use = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE
) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(
    key = genes.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = col.max,
      min = col.min
    )) ->  data.to.plot
  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = genes.plot)
  )
  # data.to.plot$genes.plot <- factor(
  #   x = data.to.plot$genes.plot,
  #   levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  # )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  } else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}
