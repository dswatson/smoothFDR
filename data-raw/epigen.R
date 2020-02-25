# Set working directory
setwd('~/Documents/smoothFDR')

# Load libraries, register cores
library(data.table)
library(limma)
library(ChAMP)
library(dplyr)
library(devtools)
library(doMC)
registerDoMC(8)

# Import, normalize, transform data
anno <- readRDS('./other/450k_anno.rds')
dir <- system.file('extdata', package = 'ChAMPdata')
myLoad <- champ.load(dir, arraytype = '450K')
betas <- champ.norm(cores = 8)
mvals <- log2(betas / (1 - betas))

# Build model, compute stats, export
des <- model.matrix(~ Sample_Group, myLoad$pd)
fit <- eBayes(lmFit(mvals, des))
saveRDS(fit, './other/limma_fit.rds')

# Get top hits, annotate, arrange, filter, export
df <- topTable(fit, coef = 2, sort.by = 'none', number = Inf) %>%
  mutate(cpg = rownames(fit), z = qnorm(pt(t, df = fit$df.total))) %>%
  inner_join(anno, by = 'cpg') %>%
  arrange(start) %>% 
  filter(seqnames == 'chr13') %>%
  rename(p.value = P.Value, BH_q.value = adj.P.Val, 
         chr = seqnames, pos = start) %>%
  select(cpg, pos, chr, z, p.value, BH_q.value)
use_data(df, 'DNAm')
