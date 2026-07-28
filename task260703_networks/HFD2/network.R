
# Microbiome Analysis - Network ------------------------------------------------

## Setup -----------------------------------------------------------------------

# install.packages("psych", dependencies = TRUE)
# install.packages("igraph", dependencies = TRUE)
# install.packages("qgraph", dependencies = TRUE)
# install.packages("vegan", dependencies = TRUE)
# install.packages("MCL", dependencies = TRUE)
# install.packages("devtools", dependencies = TRUE); library("devtools")
# devtools::install_github("zdk123/SpiecEasi", dependencies = TRUE, ref = "v1.1.1")
#remotes::install_github("cran/NetSwan")

rm(list = ls())
graphics.off()
gc()

set.seed("123")

# Load libraries
library("ggrepel")
library("ggtext")
library("ggpubr")
library("ggcorrplot")
library("patchwork")
library("cowplot")
library("svglite")
library("openxlsx")
library("xlsx")
library("ape")
library("psych")
library("igraph")
library("qgraph")
library("vegan")
library("MCL")
library("SpiecEasi")
library("NetSwan")
library("phyloseq")
library("microbiome")
library("metagMisc")
library("tidyverse")

getwd()
setwd("~/YanomamiMouseHFD/")
getwd()

list.dirs()

path.wk.dir <- paste0(getwd(), "/task260703_networks/")
path.dt.dir <- paste0(getwd(), "/tables")

for (i in c(path.wk.dir)) {
  
  if (!dir.exists(i)) { dir.create(i) }
  
}; rm(i)

my.exp <- "HFD2"

path.fg.dir <- paste0(path.wk.dir, "/", my.exp ,"/figures")
path.tb.dir <- paste0(path.wk.dir, "/", my.exp,  "/tables")

for (i in c(path.fg.dir, path.tb.dir)) {
  
  if (!dir.exists(i)) { dir.create(i, recursive = T) }
  
}; rm(i)


## All -------------------------------------------------------------------------

# flattenCorrMatrix
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[row(cormat)[ut]],
    to = rownames(cormat)[col(cormat)[ut]],
    r  = cormat[ut],
    p = pmat[ut]
  )
}

## Helper functions
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}

### Data Wrangling -------------------------------------------------------------

# Metadata
tb.metad <- 
  read.delim(paste0(getwd(), "/misc/metadata_full.tsv")) %>% 
  filter(experiment == my.exp) %>% 
  mutate(timepoint = case_when(
    experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
    experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD",
  )) %>% 
  mutate(donor_ethnicity = ifelse(donor_ethnicity == "Western", "Belgian", donor_ethnicity)) %>% 
  mutate(mouse_n = ifelse(is_pooled == "N", str_glue("{experiment}_{donor_id}_{mouse_n}"), str_glue("{experiment}_{donor_id}"))) %>% 
  select(sample_id, assembly_id, donor_id:cage_id, mouse_n, timepoint) %>% 
  mutate(donor_ethnicity = factor(donor_ethnicity, c("Belgian", "Yanomami"))) %>% 
  mutate(timepoint = factor(timepoint, c("pre-HFD", "post-HFD"))) %>% 
  arrange(donor_ethnicity, timepoint, mouse_n) %>% 
  mutate(donor_ethnicity_timepoint = str_glue("{donor_ethnicity}: {timepoint}")) %>% 
  mutate(donor_ethnicity_timepoint = factor(donor_ethnicity_timepoint, unique(.$donor_ethnicity_timepoint))); head(tb.metad)


# Taxonomy + Counts

ls.counts <- list()

for (i in unique(tb.metad$assembly_id)) {
  
  ls.counts[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_counts.xlsx"), sheet = 1)
  
}; rm(i)

ls.taxa <- list()

for (i in unique(tb.metad$assembly_id)) {
  
  ls.taxa[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_taxonomy.xlsx"), sheet = 1)
  
}; rm(i)

ls.taxa.counts <- list()

for (i in unique(tb.metad$assembly_id)) {
  
  ls.taxa.counts[[i]] <- 
    merge(ls.taxa[[i]], ls.counts[[i]]) %>% 
    group_by(across(tax_id:Species)) %>%
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE))) %>% 
    data.frame(check.names = F) %>% 
    mutate("prev" = rowSums(.[tb.metad %>% filter(assembly_id == i) %>% pull(sample_id)] != 0)) %>% 
    filter(prev >= (tb.metad %>% filter(assembly_id == i) %>% pull(sample_id) %>% length())*.1) %>% 
    select(-prev)
  
}; rm(i)

tmp.taxa.counts <- 
  reduce(ls.taxa.counts, full_join) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)))

tmp.tb.metad <- tb.metad %>% rename(ID = sample_id)

tmp.tb.count <- tmp.taxa.counts %>% select(ID = tax_id, tmp.tb.metad$ID)

tmp.tb.taxon <- tmp.taxa.counts %>% select(ID = tax_id, Kingdom:Species)

row.names(tmp.tb.count) <- tmp.tb.count$ID
row.names(tmp.tb.taxon) <- tmp.tb.taxon$ID
row.names(tmp.tb.metad) <- tmp.tb.metad$ID

tmp.tb.count <- tmp.tb.count %>% select(-ID)
tmp.tb.taxon <- tmp.tb.taxon %>% select(-ID)

tmp.tb.count <- as.matrix(tmp.tb.count)
tmp.tb.taxon <- as.matrix(tmp.tb.taxon)

tmp.tb.count <- otu_table(tmp.tb.count, taxa_are_rows = TRUE)
tmp.tb.taxon <- tax_table(tmp.tb.taxon)
tmp.tb.metad <- sample_data(tmp.tb.metad)

ps.table <- phyloseq(tmp.tb.count, tmp.tb.taxon, tmp.tb.metad)
ps.table <- subset_taxa(ps.table, Kingdom != "d: Unassigned")
ps.table <- subset_taxa(ps.table, Phylum != "p: Chordata")
ps.table <- subset_taxa(ps.table, Phylum != "p: Streptophyta")

ps.table.tss <- phyloseq_standardize_otu_abundance(physeq = ps.table, method = 'total')

rm(list = ls(pattern = "^ls."))
rm(list = ls(pattern = "^tmp."))

# Glom tables

ls.glom.ps <- list(
  
  "Genus" = tax_glom(ps.table, taxrank="Genus")
  
)

ls.glom.tb <- list(
  
  "Genus" = merge(tax_table(ls.glom.ps[["Genus"]]), otu_table(ls.glom.ps[["Genus"]]), by = "row.names", all = T)[-1]
  
)


fun.get.filt.tax.glom <- function(my.glom.ps, min.prev = 0.1, min.max.abs = 0) {
  
  tmp.ps <- phyloseq_filter_prevalence(phyloseq_filter_sample_wise_abund_trim(my.glom.ps, minabund = min.max.abs, relabund = T), prev.trh = min.prev)
  
  return(merge(tax_table(tmp.ps), otu_table(tmp.ps), by = "row.names", all = T)[-1])
  
}


ls.glom.tb.filt <- list(

  "Genus" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001) %>% 
    select(1:6, 8:ncol(.))
  
)

fun.get.filt.tax.glom.ps <- function(my.glom.ps, min.prev = 0.1, min.max.abs = 0) {
  
  tmp.ps <- phyloseq_filter_prevalence(phyloseq_filter_sample_wise_abund_trim(my.glom.ps, minabund = min.max.abs, relabund = T), prev.trh = min.prev)
  
  return(tmp.ps)
  
}

ls.glom.ps.filt <- list(

  "Genus" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$"))
  
)


# Network -----------------------------------------------------------------


tb.glom.g <- ls.glom.tb.filt[["Genus"]]; tb.glom.g

top.phylum <- top.phylum <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota")

tb.glom.g.fix <- 
  tb.glom.g %>% 
  filter(!str_detect(Genus, " Un.\\)$")) %>% 
  mutate(Phylum2 = ifelse(Phylum %in% top.phylum, Phylum, "Others")) %>% 
  mutate(Phylum = case_when(
    Phylum2 == Phylum ~ Phylum,
    Phylum2 == "Others" & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
    Phylum2 == "Others" & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria"))

colors.tax <- c("#B06888FF", "#5080C8FF", "#E8B800FF", "#18C898FF", "#F88048FF", "#C8C8C8FF", "#303028FF")
names(colors.tax) <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota", "Other: Bacteria", "Other: Domains")

tb.net.g.abs <- tb.glom.g.fix 

my.group.1 <- (tb.metad$donor_ethnicity_timepoint %>% levels(.))[[1]]
my.group.2 <- (tb.metad$donor_ethnicity_timepoint %>% levels(.))[[2]]
my.group.3 <- (tb.metad$donor_ethnicity_timepoint %>% levels(.))[[3]]
my.group.4 <- (tb.metad$donor_ethnicity_timepoint %>% levels(.))[[4]]

my.samples.group.1 <- tb.metad %>% filter(donor_ethnicity_timepoint == my.group.1) %>% pull(sample_id)
my.samples.group.2 <- tb.metad %>% filter(donor_ethnicity_timepoint == my.group.2) %>% pull(sample_id)
my.samples.group.3 <- tb.metad %>% filter(donor_ethnicity_timepoint == my.group.3) %>% pull(sample_id)
my.samples.group.4 <- tb.metad %>% filter(donor_ethnicity_timepoint == my.group.4) %>% pull(sample_id)

tb.net.g.rel <- 
  tb.net.g.abs %>% 
  mutate(!!sym(my.group.1) := rowSums(select(., any_of(my.samples.group.1)))) %>% 
  mutate(!!sym(my.group.2) := rowSums(select(., any_of(my.samples.group.2)))) %>% 
  mutate(!!sym(my.group.3) := rowSums(select(., any_of(my.samples.group.3)))) %>% 
  mutate(!!sym(my.group.4) := rowSums(select(., any_of(my.samples.group.4)))) %>% 
  select(-any_of(tb.metad$sample_id)) %>% 
  mutate(across(where(is.numeric), ~./sum(.)))

names(tb.net.g.rel) <- paste0(names(tb.net.g.rel), ".rel")

tb.net.g <- 
  cbind(tb.net.g.abs, tb.net.g.rel) %>% 
  select(Phylum, Taxa = Genus, 
         any_of(my.samples.group.1), starts_with(my.group.1),
         any_of(my.samples.group.2), starts_with(my.group.2),
         any_of(my.samples.group.3), starts_with(my.group.3),
         any_of(my.samples.group.4), starts_with(my.group.4)) %>% 
  mutate(RowSum = (rowSums(select(., any_of(tb.metad$sample_id)))) / sum(rowSums(select(., any_of(tb.metad$sample_id)))),
         RowName = make.unique(Taxa)) %>% 
  arrange(-RowSum) %>% 
  column_to_rownames(., "RowName")

# Subset by treat
tb.count.1 <- tb.net.g %>% select(any_of(my.samples.group.1)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))
tb.count.2 <- tb.net.g %>% select(any_of(my.samples.group.2)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))
tb.count.3 <- tb.net.g %>% select(any_of(my.samples.group.3)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))
tb.count.4 <- tb.net.g %>% select(any_of(my.samples.group.4)) %>% t() %>% data.frame(check.names = F) %>% select(where(~ sum(.) != 0))

# Correlation tables
otu.corr.pearson.1 <- corr.test(tb.count.1, method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)
otu.corr.pearson.2 <- corr.test(tb.count.2, method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)
otu.corr.pearson.3 <- corr.test(tb.count.3, method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)
otu.corr.pearson.4 <- corr.test(tb.count.4, method = "pearson", use = "pairwise", adjust = "fdr", alpha = 0.05)

# SPARCC

get.sparcc.fun <- function(tb) {

  my.tb.sparcc <- tb

  sparcc.out <- sparcc(my.tb.sparcc)
  sparcc.out

  rownames(sparcc.out$Cor) <- colnames(my.tb.sparcc)
  colnames(sparcc.out$Cor) <- colnames(my.tb.sparcc)
  rownames(sparcc.out$Cov) <- colnames(my.tb.sparcc)
  colnames(sparcc.out$Cov) <- colnames(my.tb.sparcc)

  sparcc.out.portable <- sparcc.out$Cor %>% reorder_cormat %>% get_upper_tri() %>% reshape2::melt() %>% na.omit()

  # sparcc.out.plot.1 <- sparcc.out.portable %>% ggplot(aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # sparcc.out.plot.1

  sparcc.out.boot <- sparccboot(my.tb.sparcc, R = 10, ncpus = 10)
  sparcc.out.boot.p <- pval.sparccboot(sparcc.out.boot)

  sparcc.out.boot.cors <- sparcc.out.boot.p$cors
  sparcc.out.boot.pvals <- sparcc.out.boot.p$pvals

  sparcc.out.pcors <- diag(0.5, nrow = dim(sparcc.out$Cor)[1], ncol = dim(sparcc.out$Cor)[1])
  sparcc.out.pcors[upper.tri(sparcc.out.pcors, diag=FALSE)] <- sparcc.out.boot.cors
  sparcc.out.pcors <- sparcc.out.pcors + t(sparcc.out.pcors)

  sparcc.out.pval <- diag(0.5, nrow = dim(sparcc.out$Cor)[1], ncol = dim(sparcc.out$Cor)[1])
  sparcc.out.pval[upper.tri(sparcc.out.pval, diag=FALSE)] <- sparcc.out.boot.pvals
  sparcc.out.pval <- sparcc.out.pval + t(sparcc.out.pval)

  rownames(sparcc.out.pcors) <- colnames(my.tb.sparcc)
  colnames(sparcc.out.pcors) <- colnames(my.tb.sparcc)
  rownames(sparcc.out.pval) <- colnames(my.tb.sparcc)
  colnames(sparcc.out.pval) <- colnames(my.tb.sparcc)

  sparcc.reordered.all <- reorder_cor_and_p(sparcc.out.pcors, sparcc.out.pval)
  sparcc.reordered.cor <- sparcc.reordered.all$r
  sparcc.reordered.p <- sparcc.reordered.all$p

  sparcc.reordered.cor.processed <- sparcc.reordered.cor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(cor = value)
  sparcc.reordered.p.processed <- sparcc.reordered.p  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(p = value)

  sparcc.out.fdr <-
    left_join(sparcc.reordered.cor.processed, sparcc.reordered.p.processed, by = c("Var1", "Var2")) %>%
    filter(Var1 != Var2) %>%
    mutate(fdr = p.adjust(p, method = "fdr"))

  sparcc.out.fdr.sig <-
    sparcc.out.fdr %>%
    filter(p < 0.05 & abs(cor) > 0.5)

  sparcc.out.plot.2 <-
    sparcc.out.fdr %>%
    ggplot(aes(x = Var2, y = Var1, fill = cor)) +
    geom_tile() +
    geom_point(data = sparcc.out.fdr.sig, shape = 1) +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          aspect.ratio = 1)

  print(sparcc.out.plot.2)

  return(sparcc.out.fdr)

}

# otu.corr.sparcc.1 <- get.sparcc.fun(tb.count.1)
# otu.corr.sparcc.2 <- get.sparcc.fun(tb.count.2)
# otu.corr.sparcc.3 <- get.sparcc.fun(tb.count.3)
# otu.corr.sparcc.4 <- get.sparcc.fun(tb.count.4)

otu.corr.pearson.flatten.1 <- flattenCorrMatrix(otu.corr.pearson.1$r, otu.corr.pearson.1$p)
otu.corr.pearson.flatten.2 <- flattenCorrMatrix(otu.corr.pearson.2$r, otu.corr.pearson.2$p)
otu.corr.pearson.flatten.3 <- flattenCorrMatrix(otu.corr.pearson.3$r, otu.corr.pearson.3$p)
otu.corr.pearson.flatten.4 <- flattenCorrMatrix(otu.corr.pearson.4$r, otu.corr.pearson.4$p)

### Getting edges and nodes ----------------------------------------------------

# Net 1
otu.corr.edges.1 <- 
  otu.corr.pearson.flatten.1 %>% 
  filter(abs(r) >= 0.5 & p <= 0.05) %>%
  dplyr::rename(correlation = r, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.1 <- 
  merge(x = sort(unique(c(otu.corr.edges.1$from, otu.corr.edges.1$to))), y = tb.net.g, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa = x, any_of(my.samples.group.1), starts_with(my.group.1)) 

# Net 2
otu.corr.edges.2 <- 
  otu.corr.pearson.flatten.2 %>% 
  filter(abs(r) >= 0.5 & p <= 0.05) %>%
  dplyr::rename(correlation = r, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.2 <- 
  merge(x = sort(unique(c(otu.corr.edges.2$from, otu.corr.edges.2$to))), y = tb.net.g, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa = x, any_of(my.samples.group.2), starts_with(my.group.2)) 

# Net 3
otu.corr.edges.3 <- 
  otu.corr.pearson.flatten.3 %>% 
  filter(abs(r) >= 0.5 & p <= 0.05) %>%
  dplyr::rename(correlation = r, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.3 <- 
  merge(x = sort(unique(c(otu.corr.edges.3$from, otu.corr.edges.3$to))), y = tb.net.g, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa = x, any_of(my.samples.group.3), starts_with(my.group.3)) 

# Net 4
otu.corr.edges.4 <- 
  otu.corr.pearson.flatten.4 %>% 
  filter(abs(r) >= 0.5 & p <= 0.05) %>%
  dplyr::rename(correlation = r, p.value = p) %>% 
  mutate(type = if_else(correlation > 0, 2, 1))

otu.corr.nodes.4 <- 
  merge(x = sort(unique(c(otu.corr.edges.4$from, otu.corr.edges.4$to))), y = tb.net.g, by.x = 1, by.y = 0, all.x = T) %>% 
  select(Phylum, Taxa = x, any_of(my.samples.group.4), starts_with(my.group.4)) 

### Getting IDs ----------------------------------------------------------------

top.p <- names(colors.tax)

all.ids <- 
  reduce(list(otu.corr.nodes.1, otu.corr.nodes.2, otu.corr.nodes.3, otu.corr.nodes.4), full_join, by = c("Phylum", "Taxa")) %>% 
  mutate(Sum = rowSums(select(., matches(".rel")), na.rm = TRUE)) %>%
  arrange(desc(Sum)) %>% 
  mutate(IDT = seq_len(nrow(.))) %>%
  mutate(Phylum = factor(Phylum, levels = top.p)) %>%
  mutate(IDP = as.numeric(Phylum)) %>% 
  mutate(IDP2 = case_when(
    Phylum == "Other: Domains" ~ 2, 
    .default = 1)) %>% 
  mutate(ID = Taxa) %>% 
  select(ID, IDP, IDP2, Phylum, IDT, Taxa)

otu.corr.nodes.1 <- merge(all.ids, otu.corr.nodes.1, sort = F) %>% select(ID, IDP, IDP2, Phylum, IDT, Taxa, everything())
otu.corr.nodes.2 <- merge(all.ids, otu.corr.nodes.2, sort = F) %>% select(ID, IDP, IDP2, Phylum, IDT, Taxa, everything())
otu.corr.nodes.3 <- merge(all.ids, otu.corr.nodes.3, sort = F) %>% select(ID, IDP, IDP2, Phylum, IDT, Taxa, everything())
otu.corr.nodes.4 <- merge(all.ids, otu.corr.nodes.4, sort = F) %>% select(ID, IDP, IDP2, Phylum, IDT, Taxa, everything())

### Getting nets ---------------------------------------------------------------

net.1 <- graph_from_data_frame(d=otu.corr.edges.1, vertices=otu.corr.nodes.1, directed = F)
net.2 <- graph_from_data_frame(d=otu.corr.edges.2, vertices=otu.corr.nodes.2, directed = F)
net.3 <- graph_from_data_frame(d=otu.corr.edges.3, vertices=otu.corr.nodes.3, directed = F)
net.4 <- graph_from_data_frame(d=otu.corr.edges.4, vertices=otu.corr.nodes.4, directed = F)

### Getting net measures -------------------------------------------------------

# Adding degree, bet and clos

otu.corr.nodes.1 <- 
  merge(otu.corr.nodes.1, 
        reduce(list(data.frame(degree = degree(net.1, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.1, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.1, mode="all", normalized = T)) %>% rownames_to_column("Taxa"),
                    data.frame(hub_klein = hub.score(net.1)$vector) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>% 
  arrange(-hub_klein) %>% 
  mutate(tmp = 1:n()) %>% 
  mutate(sig = ifelse(tmp <= 5, "*", "")) %>% 
  mutate(sig2 = ifelse(tmp <= 5, tmp, "")) %>% 
  select(-tmp) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = paste0(my.group.1, ".rel"))

otu.corr.nodes.2 <- 
  merge(otu.corr.nodes.2, 
        reduce(list(data.frame(degree = degree(net.2, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.2, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.2, mode="all", normalized = T)) %>% rownames_to_column("Taxa"),
                    data.frame(hub_klein = hub.score(net.2)$vector) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>%   
  arrange(-hub_klein) %>% 
  mutate(tmp = 1:n()) %>% 
  mutate(sig = ifelse(tmp <= 5, "*", "")) %>% 
  mutate(sig2 = ifelse(tmp <= 5, tmp, "")) %>% 
  select(-tmp) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = paste0(my.group.2, ".rel"))

otu.corr.nodes.3 <- 
  merge(otu.corr.nodes.3, 
        reduce(list(data.frame(degree = degree(net.3, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.3, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.3, mode="all", normalized = T)) %>% rownames_to_column("Taxa"),
                    data.frame(hub_klein = hub.score(net.3)$vector) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>%   
  arrange(-hub_klein) %>% 
  mutate(tmp = 1:n()) %>% 
  mutate(sig = ifelse(tmp <= 5, "*", "")) %>% 
  mutate(sig2 = ifelse(tmp <= 5, tmp, "")) %>% 
  select(-tmp) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = paste0(my.group.3, ".rel"))

otu.corr.nodes.4 <- 
  merge(otu.corr.nodes.4, 
        reduce(list(data.frame(degree = degree(net.4, mode="all")) %>% rownames_to_column("Taxa"), 
                    data.frame(betweenness = betweenness(net.4, directed=T, weights=NA, normalized = T)) %>% rownames_to_column("Taxa"), 
                    data.frame(closeness = closeness(net.4, mode="all", normalized = T)) %>% rownames_to_column("Taxa"),
                    data.frame(hub_klein = hub.score(net.4)$vector) %>% rownames_to_column("Taxa")), 
               full_join, by = "Taxa")) %>%   
  arrange(-hub_klein) %>% 
  mutate(tmp = 1:n()) %>% 
  mutate(sig = ifelse(tmp <= 5, "*", "")) %>% 
  mutate(sig2 = ifelse(tmp <= 5, tmp, "")) %>% 
  select(-tmp) %>% 
  select(ID, everything()) %>% 
  dplyr::rename("rel" = paste0(my.group.4, ".rel"))


### Getting net layout ---------------------------------------------------------

# Net 1
net.1 <- graph_from_data_frame(d=otu.corr.edges.1, vertices=otu.corr.nodes.1, directed = F)

V(net.1)$name <- seq(1, length(V(net.1)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.1 <- get.edgelist(net.1) ; class(e.1) <- "numeric"
l.1 <- qgraph.layout.fruchtermanreingold(e.1, vcount=vcount(net.1), area=8*(vcount(net.1)^1.85), repulse.rad=(vcount(net.1)^2.6))
#marks.1 <- lapply(1:length(articulation.points(net.1)), function(x) which(V(net.1)$name == articulation.points(net.1)[[x]]))

# Net 2
net.2 <- graph_from_data_frame(d=otu.corr.edges.2, vertices=otu.corr.nodes.2, directed = F)

V(net.2)$name <- seq(1, length(V(net.2)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.2 <- get.edgelist(net.2) ; class(e.2) <- "numeric"
l.2 <- qgraph.layout.fruchtermanreingold(e.2, vcount=vcount(net.2), area=8*(vcount(net.2)^1.85), repulse.rad=(vcount(net.2)^2.6))
#marks.2 <- lapply(1:length(articulation.points(net.2)), function(x) which(V(net.2)$name == articulation.points(net.2)[[x]]))

# Net 2
net.3 <- graph_from_data_frame(d=otu.corr.edges.3, vertices=otu.corr.nodes.3, directed = F)

V(net.3)$name <- seq(1, length(V(net.3)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.3 <- get.edgelist(net.3) ; class(e.3) <- "numeric"
l.3 <- qgraph.layout.fruchtermanreingold(e.3, vcount=vcount(net.3), area=8*(vcount(net.3)^1.85), repulse.rad=(vcount(net.3)^2.6))
#marks.3 <- lapply(1:length(articulation.points(net.3)), function(x) which(V(net.3)$name == articulation.points(net.3)[[x]]))


# Net 2
net.4 <- graph_from_data_frame(d=otu.corr.edges.4, vertices=otu.corr.nodes.4, directed = F)

V(net.4)$name <- seq(1, length(V(net.4)$name)) # Customized layout to avoid nodes overlapping. Convert node label from names to numerical IDs
e.4 <- get.edgelist(net.4) ; class(e.4) <- "numeric"
l.4 <- qgraph.layout.fruchtermanreingold(e.4, vcount=vcount(net.4), area=8*(vcount(net.4)^1.85), repulse.rad=(vcount(net.4)^2.6))
#marks.4 <- lapply(1:length(articulation.points(net.4)), function(x) which(V(net.4)$name == articulation.points(net.4)[[x]]))



### Getting graphs -------------------------------------------------------------

# triangle vertex shape
mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 0.025
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- 2.5
  }
  
  symbols(
    x = coords[, 1], y = coords[, 2], bg = vertex.color,
    stars = cbind(vertex.size, vertex.size, vertex.size),
    add = TRUE, inches = FALSE
  )
}
# clips as a circle
add_shape("triangle",
          clip = shapes("circle")$clip,
          plot = mytriangle
)

# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 0.025
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- 2.5
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[, 1], coords[, 2], vertex.color, vertex.size, norays,
         FUN = function(x, y, bg, size, nor) {
           symbols(
             x = x, y = y, bg = bg,
             stars = matrix(c(size, size / 2), nrow = 1, ncol = nor * 2),
             add = TRUE, inches = FALSE
           )
         }
  )
}

# no clipping, edges will be below the vertices anyway
add_shape("star",
          clip = shape_noclip,
          plot = mystar, parameters = list(vertex.norays = 5)
)

# cols

cols.cor <- c("#D8263266", "#75B7D166")
names(cols.cor) <- c("Negative", "Positive")

shapes.p <- c("circle", "square")

# Main plot function

# SVG graphics device
# svg(paste0(path, "figures/Network/Network1.svg"),
#     width=25,
#     height=25,
#     pointsize=50)

#png(filename = paste0("task260224_taxProfile6/figures/network/network.png"), width = 5000 * 3, height = 5000 * .75, units = "px", res = 1200, type = "cairo")

par(mfrow = c(1, 4), mai = c(0, 0, 0.5, 0), bg = 'white')
# Net 1

plot.igraph(net.1, 
            vertex.size = 8,
            #vertex.size = 0,
            vertex.label = V(net.1)$sig2,
            vertex.label.cex = .8,
            vertex.label.color = "black",
            vertex.frame.color = "#00000000",
            vertex.color = colors.tax[V(net.1)$Phylum] %>% str_replace(., "FF$", "FF"),
            vertex.shape = shapes.p[V(net.1)$IDP2],
            edge.color = cols.cor[E(net.1)$type],
            edge.width = 1,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            #edge.curved=0.5,
            layout = l.1, 
            rescale = T,
            main = my.group.1, 
            margin = c(0, 0, 0, 0))

plot.igraph(net.2, 
            vertex.size = 8,
            #vertex.size = 0,
            vertex.label = V(net.2)$sig2,
            vertex.label.cex = .8,
            vertex.label.color = "black",
            vertex.frame.color = "#00000000",
            vertex.color = colors.tax[V(net.2)$Phylum] %>% str_replace(., "FF$", "FF"),
            vertex.shape = shapes.p[V(net.2)$IDP2],
            edge.color = cols.cor[E(net.2)$type],
            edge.width = 1,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            #edge.curved=0.5,
            layout = l.2, 
            rescale = T,
            main = my.group.2, 
            margin = c(0, 0, 0, 0))


plot.igraph(net.3, 
            vertex.size = 8,
            #vertex.size = 0,
            vertex.label = V(net.3)$sig2,
            vertex.label.cex = .8,
            vertex.label.color = "black",
            vertex.frame.color = "#00000000",
            vertex.color = colors.tax[V(net.3)$Phylum] %>% str_replace(., "FF$", "FF"),
            vertex.shape = shapes.p[V(net.3)$IDP2],
            edge.color = cols.cor[E(net.3)$type],
            edge.width = 1,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            #edge.curved=0.5,
            layout = l.3, 
            rescale = T,
            main = my.group.3, 
            margin = c(0, 0, 0, 0))

plot.igraph(net.4, 
            vertex.size = 8,
            #vertex.size = 0,
            vertex.label = V(net.4)$sig2,
            vertex.label.cex = .8,
            vertex.label.color = "black",
            vertex.frame.color = "#00000000",
            vertex.color = colors.tax[V(net.4)$Phylum] %>% str_replace(., "FF$", "FF"),
            vertex.shape = shapes.p[V(net.4)$IDP2],
            edge.color = cols.cor[E(net.4)$type],
            edge.width = 1,
            edge.arrow.mode = "-", 
            edge.arrow.size=0,
            #edge.curved=0.5,
            layout = l.4, 
            rescale = T,
            main = my.group.4, 
            margin = c(0, 0, 0, 0))


#dev.off() 

### Getting legends ------------------------------------------------------------

library(ggtext)
library(ggstar)

my.theme <-   
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_markdown(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_markdown(size = 10, colour = 'black', face = "bold"),
        strip.background = element_rect(fill = NA, color = NA), 
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        panel.border = element_rect(color = "black", linewidth = .25, fill = NA),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y.right = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black'),
        axis.text.y = element_markdown(size = 6, colour = 'black'),
        
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black'),
        plot.subtitle = element_markdown(size = 6, colour = 'black'),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

gg.leg.p.full <- 
  ggplot(top.p %>% data.frame(), aes(x = factor(., levels = unique(.)), y = 1:length(top.p))) +
  geom_point(aes(fill = factor(., levels = unique(.))), color = "black", shape = 22, size = 3) +
  scale_fill_manual(values = tax.colors2) +  
  labs(fill = "Phylum/Domain") + 
  theme_bw() +
  my.theme + 
  theme(legend.position = "right",
        axis.ticks.length.x = unit(0, "pt"),
        axis.ticks.length.y = unit(2.5, "pt"),
        axis.text.x = element_blank(),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_markdown(size = 7 + 2, face = "bold"),
        legend.text  = element_markdown(size = 7),
        legend.key.size = unit(0.65, "lines"))

gg.leg.p <- ggpubr::get_legend(gg.leg.p.full)

gg.leg.line.full <- 
  ggplot(data.frame(y = factor(rep(c("Positive", "Negative"), 2), levels = c("Positive", "Negative")), x = c(1, 2, 2, 1)), aes(x = x, y = y, group = y, color = y)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols.cor, name = "Correlation") +  
  theme_bw() +
  my.theme + 
  theme(legend.position = "right",
        axis.ticks.length.x = unit(0, "pt"),
        axis.ticks.length.y = unit(2.5, "pt"),
        axis.text.x = element_blank(),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_markdown(size = 7 + 2, face = "bold"),
        legend.text  = element_markdown(size = 7),
        legend.key.width = unit(1.25, "lines"),
        legend.key.size = unit(0.65, "lines"))

gg.leg.line <- ggpubr::get_legend(gg.leg.line.full)

gg.leg.degree.full <- 
  ggplot(data.frame(y = factor(c("Bacteria", "Other"),
                               levels = c("Bacteria", "Other"))), aes(x = 1, y = y, starshape = y)) +
  geom_star(fill = "black") + 
  theme_bw() +
  my.theme + 
  scale_starshape_manual(values = c(15, 13, 11, 1), name = "Domain") +
  scale_size_continuous(range = c(.5, 2), breaks=c(1, 2, 3), labels=c("", "", ""), name = "Degree") +
  theme(legend.position = "right",
        axis.ticks.length.x = unit(0, "pt"),
        axis.ticks.length.y = unit(2.5, "pt"),
        axis.text.x = element_blank(),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_markdown(size = 7 + 2, face = "bold"),
        legend.text  = element_markdown(size = 7),
        legend.key.size = unit(0.65, "lines"))

gg.leg.degree <- ggpubr::get_legend(gg.leg.degree.full)

ggpubr::as_ggplot(gg.leg.p) + ggpubr::as_ggplot(gg.leg.line) + ggpubr::as_ggplot(gg.leg.degree) + plot_layout(ncol = 1, heights = c(5, 2, 3))

png(filename = paste0("task260224_taxProfile6/figures/network/network_leg.png"), width = 5000 * .5, height = 5000 * 1, units = "px", res = 1200, type = "cairo")

ggpubr::as_ggplot(gg.leg.p) + ggpubr::as_ggplot(gg.leg.line) + ggpubr::as_ggplot(gg.leg.degree) + plot_layout(ncol = 1, heights = c(5, 2, 3))

dev.off() 

### Getting keystone -----------------------------------------------------------

tb.keystone <-
  bind_rows(
    otu.corr.nodes.1 %>% select(Phylum, Taxa, degree, betweenness, closeness, hub_klein, sig) %>% mutate(Group = levels(tb.metad$Group_tribe)[[1]]),
    otu.corr.nodes.2 %>% select(Phylum, Taxa, degree, betweenness, closeness, hub_klein, sig) %>% mutate(Group = levels(tb.metad$Group_tribe)[[2]]),
    otu.corr.nodes.3 %>% select(Phylum, Taxa, degree, betweenness, closeness, hub_klein, sig) %>% mutate(Group = levels(tb.metad$Group_tribe)[[3]]),
    otu.corr.nodes.4 %>% select(Phylum, Taxa, degree, betweenness, closeness, hub_klein, sig) %>% mutate(Group = levels(tb.metad$Group_tribe)[[4]])) %>% 
  mutate(Group = factor(Group, levels = c(levels(tb.metad$Group_tribe)[[1]], levels(tb.metad$Group_tribe)[[2]], levels(tb.metad$Group_tribe)[[3]], levels(tb.metad$Group_tribe)[[4]]))) %>% 
  mutate(my.label = ifelse(sig == "*", as.character(Taxa), NA_character_)) %>% 
  mutate(my.alpha = ifelse(sig == "*", 1, 0.35)) %>% 
  mutate(
    Domain = case_when(
      Phylum == "Other: Domains" ~ "Other",
      .default = "Bacteria") %>% 
      factor(., levels = c("Bacteria", "Other"))) 

ggplot(tb.keystone %>% mutate(Group = str_replace(Group, " \\(", "<br>(")), aes(y = hub_klein, x = Group)) +
  geom_text_repel(
    aes(point.size = 1.5, label = my.label), 
    size = 1.5, point.padding = 0, min.segment.length = 0.1, max.time = 1, color = "black", 
    max.iter = 1e5, position = position_jitter(width = .15,seed = 124), box.padding = 0.25, max.overlaps = Inf) +
  geom_point(aes(fill = Phylum, alpha = my.alpha, shape = Domain), size = 1.5, color = "black", stroke = 0.25,
             position = position_jitter(width = .15, seed = 124)) + 
  scale_y_continuous(limits = c(0, max(tb.keystone$hub_klein)*1.1), expand = expansion(mult = c(0.05, 0.1))) +
  scale_shape_manual(values = c(21, 22), name = "Domain") +
  ylab('Kleinberg’s hub score<br>') + xlab('') + 
  scale_fill_manual(name = "", values = tax.colors2) +
  scale_alpha_identity() +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 21, 
                                            #fill = cols.p, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0,
                                            alpha = 1),
                        ncol = 1)) +
  my.theme +
  theme(legend.position = "none", 
        legend.box.background = element_rect(fill = "white", colour = "black", linewidth = .75),
        legend.justification = c(0, .5), 
        legend.box.margin = margin(0, 0, 0, 0),
        legend.text  = element_markdown(size = 7),
        legend.key.size = unit(0.65, "lines"), 
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1), 
        axis.text.x = element_markdown(size = 7, colour = 'black', face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm')) -> a

a

ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_key.png"),
       plot = a,
       dpi = 1200,
       width = 5000 * 1.75,
       height = 5000 * .65,
       units = "px")

ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_key.svg"),
       plot = a,
       dpi = 1200,
       width = 5000 * 1.5,
       height = 5000 * .75,
       units = "px")

### CompM net measures ---------------------------------------------------------

fun.get.kruskal <- function(tb, gp_comp, gp_facet, value_var = "value") {
  
  tb %>%
    select(all_of(gp_comp), all_of(gp_facet), all_of(value_var)) %>% 
    group_by(pick(gp_facet)) %>%
    nest() %>%
    mutate(kruskal = map(.x = data, ~kruskal.test(get(value_var) ~ get(gp_comp), data = .x) %>% broom::tidy())) %>%
    unnest(kruskal) %>%
    ungroup() %>% 
    mutate(across(where(is.numeric), ~round(., 6))) %>% 
    mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>% 
    mutate(p.full = paste0("𝘱-adj: ", ifelse(adj.p.value < 0.001, "<0.001", sprintf("%1.3f", round(adj.p.value, 3))))) %>% 
    select(-data) %>% 
    data.frame(check.names = F) -> p
  
  return(p)
  
}

fun.get.post.hoc.lsd <- function(tb, tb.summary, tb.stat, gp_comp, gp_facet, value_var = "value", ls.invert = NULL) {
  
  tb.compm <- 
    data.frame(gp_comp = "", gp_facet = "", group = "")[-1,] %>% 
    rename(!!gp_comp := gp_comp, !!gp_facet:= gp_facet)
  
  for (i in unique(tb[[gp_facet]])) {
    
    tmp.tb <- tb %>% filter(!!as.name(gp_facet) == i)
    
    tmp.lsd <- with(tmp.tb, agricolae::kruskal(get(value_var), get(gp_comp), alpha = 0.05, p.adj = "fdr"))
    
    tmp.lsd <-
      cbind(gp_comp = row.names(tmp.lsd[["groups"]]), tmp.lsd[["groups"]]) %>%
      rename(!!gp_comp := gp_comp) %>% 
      remove_rownames() %>%
      mutate(!!gp_facet := i, group = str_trim(groups)) %>%
      select(all_of(gp_facet), all_of(gp_comp), group)
    
    if (tb.stat %>% filter(!!as.name(gp_facet) == i) %>% pull(p.value) > 0.05) {
      
      tmp.lsd <- mutate(tmp.lsd, group = NA_character_)
      
    }
    
    tb.compm <- bind_rows(tb.compm, tmp.lsd)
    
  }
  
  tb.compm <- merge(tb.summary, tb.compm, sort = F, all.x = T)
  
  return(tb.compm)
  
}

fun.upper.whisker <- function(x) {
  boxplot.stats(x)$stats[5] # The 5th element is the upper whisker
}

tb.keystone

tr.colors <- c("#FFB900FF", "#088158FF", "#5773CCFF",  "#8B4F82FF")
names(tr.colors) <- levels(tb.metad$Group_tribe)

my.theme <-     
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_markdown(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_markdown(size = 10, colour = 'black', face = "bold"),
        strip.background = element_rect(fill = NA, color = NA), 
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        panel.border = element_rect(color = "black", linewidth = .25, fill = NA),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y.right = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black'),
        axis.text.y = element_markdown(size = 6, colour = 'black'),
        
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black'),
        plot.subtitle = element_markdown(size = 6, colour = 'black'),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

net.measure.table <- merge(tb.keystone, tb.metad %>% select(Group = Group_tribe, Ethnicity) %>% distinct(), by = "Group", sort = F)

net.measure.table.longer <- 
  net.measure.table %>% 
  pivot_longer(cols = degree:closeness, names_to = "measure", values_to = "values") %>% 
  data.frame(check.names = F) %>% 
  mutate(measure = str_to_title(measure)) %>% 
  mutate(measure = factor(measure, levels = unique(measure)))

net.measure.table.summary <-
  net.measure.table.longer %>% 
  group_by(Group, measure) %>% 
  summarise(temp = list(c(psych::describe(values), quantile(values, na.rm = T), IQR = matrixStats::iqr(values, na.rm = T)))) %>% 
  unnest_wider(temp) %>% 
  select(-c(vars, trimmed, mad)) %>% 
  mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
  data.frame(check.names = F)

net.measure.kruskal <- 
  fun.get.kruskal(tb = net.measure.table.longer, gp_comp = "Group", gp_facet = "measure", value_var = "values") %>% 
  arrange(measure); head(net.measure.kruskal)

net.measure.posthoc <- 
  fun.get.post.hoc.lsd(tb = net.measure.table.longer, 
                       tb.summary = net.measure.table.summary, 
                       tb.stat = net.measure.kruskal, gp_comp = "Group", gp_facet = "measure", value_var = "values") %>% 
  arrange(measure, Group); head(net.measure.posthoc)

net.measure.gg <- 
  merge(net.measure.table.longer, net.measure.posthoc) %>% 
  merge(., net.measure.kruskal) %>% 
  mutate(p.full = paste0("KW, 𝘱-adj: ", ifelse(adj.p.value < 0.001, "<0.001", sprintf("%1.3f", round(adj.p.value, 3))))) %>% 
  arrange(measure, Group) %>% 
  
  ggplot(., aes(x = Group, y = values, fill = Group)) + 
  facet_wrap( ~ measure, scales = "free", nrow = 1, ncol = 5) +
  geom_jitter(shape = 16, alpha = 1, size = 1.5, width = .25, height = 0, color = "grey75") + 
  geom_boxplot(width = .45, alpha = .85, color = "black", outlier.shape = NA, linewidth = 0.25) +
  stat_summary(fun = "mean", color = "black", fill = "white", shape = 23, size = 0.35, stroke = 0.25) +
  stat_summary(fun = fun.upper.whisker,
               geom = "text",
               aes(label = group), # Label with rounded whisker value
               vjust = -.5, # Adjust vertical position to be above the whisker
               size = 2.5, 
               fontface = "bold",
               color = "black") +
  geom_text(mapping = aes(label = p.full, x = 0.5, y = Inf),
            vjust = 2, position = position_dodge(width = 0),
            size = 1.75, hjust = 0, color = "black", check_overlap = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.5)), breaks = scales::extended_breaks(n = 4)) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  scale_fill_manual(values = tr.colors, name = NULL) +
  ylab('') + xlab('') +
  guides(fill = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0, 
                                            alpha = 1),
                        ncol = 2)) +
  my.theme + 
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        legend.title = element_blank(),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black', lineheight = 1.25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

net.measure.gg

ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_measures.png"), plot = net.measure.gg, dpi = 1200, width = 5000 * 1.25, height = 5000 * 0.65, units = "px")
ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_measures.svg"), plot = net.measure.gg, dpi = 1200, width = 5000 * 1.25, height = 5000 * 0.65, units = "px")

### Parameters table -----------------------------------------------------------

nets <- c("net.1", "net.2", "net.3", "net.4")

net.name <- c(levels(tb.metad$Group_tribe)[[1]], levels(tb.metad$Group_tribe)[[2]], 
              levels(tb.metad$Group_tribe)[[3]], levels(tb.metad$Group_tribe)[[4]])

nodes <- c("otu.corr.nodes.1", "otu.corr.nodes.2", "otu.corr.nodes.3", "otu.corr.nodes.4")

net.attr.table <- data.frame("Attributes" = c("n. of nodes", 
                                              "n. of edges", 
                                              "positive edges",
                                              "negative edges",
                                              "n. of clusters",
                                              "size of clusters",
                                              "edge density",
                                              "avg. degree",
                                              "max. degree",
                                              "avg. closeness",
                                              "avg. betweenness",
                                              "max. betweenness",
                                              "keystone taxa",
                                              "modularity coef.",
                                              "articulation points",
                                              "clustering coefficient",
                                              "diameter"))


for (i in seq_along(nets)) {
  
  n <- net.name[i]
  
  net <- get(nets[i])
  
  node <- get(nodes[i])
  
  # Getting IDG back to name
  V(net)$name <- paste0(V(net)$IDT, ") ", V(net)$Taxa)
  
  # n. of vertices
  attr01 <- vcount(net)
  
  # n. of edges
  attr02 <- ecount(net)
  
  # n. of positive edges
  attr03 <- paste0(unname(table(E(net)$type)[2]), " (", round((unname(table(E(net)$type)[2])/ecount(net))*100, 2), "%)")
  
  # n. of negative edges
  attr04 <- paste0(unname(table(E(net)$type)[1]), " (", round((unname(table(E(net)$type)[1])/ecount(net))*100, 2), "%)")
  
  # n. of clusters
  # clusters(net)$membership                                  # Membership of vertices
  attr05 <- components(net)$no                                   # Number of clusters
  attr06 <- paste(components(net)$csize, collapse = ", ")        # Sizes of clusters
  
  # Hub detection
  attr07 <- round(edge_density(net, loops=F), 3) # The proportion of present edges from all possible edges in the network
  
  attr08 <- paste0(round(mean(degree(net, mode = "all")), 3),
                   " (±", round(sd(degree(net, mode = "all")), 3), ")") # Average degrees of nodes
  
  attr09 <- max(degree(net, mode="all")) # Node/s with the most degrees
  
  attr10 <- paste0(round(mean(closeness(net, mode = "all", weights = NA, normalized = T)), 3),
                   " (±", round(sd(closeness(net, mode = "all", weights = NA, normalized = T)), 3), ")") # Average closeness of nodes
  
  attr11 <- paste0(round(mean(betweenness(net, directed=T, weights=NA, normalized = T)), 3),
                   " (±", round(sd(betweenness(net, directed=T, weights=NA, normalized = T)), 3), ")") # Average betweenness of nodes
  
  attr12 <- round(max(betweenness(net, directed=T, weights=NA, normalized = T)), 3) # Node/s with the most degrees
  
  # Sort the species based on hubbiness score
  attr13 <- node %>% 
    mutate(hub_sig = ifelse(sig == "*", paste0(Taxa, "(", round(hub_klein, 2),")"), NA_character_)) %>% 
    filter(sig == "*") %>% pull(hub_sig) %>% paste(collapse = ", ")
  
  # Modularity
  attr14 <- round(modularity(net, membership(components(net))), 4)
  
  # Articulation.points
  attr15 <- length(articulation_points(net))
  
  # Clustering Coeficient
  attr16 <- round(transitivity(net, type="global"), 3)
  
  # Diameter
  attr17 <- diameter(net)
  
  attr.col <- t(rbind(lapply(ls(pattern = "^attr[0-9]"), get)))
  
  net.attr.table <- cbind(net.attr.table, attr.col)
  names(net.attr.table)[i+1] <- n 
  
}

net.attr.table

### Robustness -----------------------------------------------------------------

rob.net <- 
  
  rbind(
    
    data.frame(node_frac = 0, Betweenness = 0, Degree = 0, Cascading = 0, Random = 0, Group = levels(tb.metad$Group_tribe)[[1]]),
    swan_combinatory(net.1,10) %>% data.frame() %>% dplyr::rename(node_frac = 1, Betweenness = 2, Degree = 3, Cascading = 4, Random = 5) %>% 
      mutate(Group = levels(tb.metad$Group_tribe)[[1]]),
    
    data.frame(node_frac = 0, Betweenness = 0, Degree = 0, Cascading = 0, Random = 0, Group = levels(tb.metad$Group_tribe)[[2]]),
    swan_combinatory(net.2,10) %>% data.frame() %>% dplyr::rename(node_frac = 1, Betweenness = 2, Degree = 3, Cascading = 4, Random = 5) %>% 
      mutate(Group = levels(tb.metad$Group_tribe)[[2]]),
    
    data.frame(node_frac = 0, Betweenness = 0, Degree = 0, Cascading = 0, Random = 0, Group = levels(tb.metad$Group_tribe)[[3]]),
    swan_combinatory(net.3,10) %>% data.frame() %>% dplyr::rename(node_frac = 1, Betweenness = 2, Degree = 3, Cascading = 4, Random = 5) %>% 
      mutate(Group = levels(tb.metad$Group_tribe)[[3]]),
    
    data.frame(node_frac = 0, Betweenness = 0, Degree = 0, Cascading = 0, Random = 0, Group = levels(tb.metad$Group_tribe)[[4]]),
    swan_combinatory(net.4,10) %>% data.frame() %>% dplyr::rename(node_frac = 1, Betweenness = 2, Degree = 3, Cascading = 4, Random = 5) %>% 
      mutate(Group = levels(tb.metad$Group_tribe)[[4]])) %>% 
  
  pivot_longer(cols = -c(Group, node_frac), names_to = "attack", values_to = "loss") %>% 
  
  mutate(attack = factor(attack, levels = c("Degree", "Betweenness", "Cascading", "Random")))

ggplot(rob.net, aes(x = node_frac, y = loss)) + 
  facet_wrap(.~attack, scales = "fixed", nrow = 1) +
  geom_hline(yintercept = 0.8, size = 0.5, linetype = 2, color = "black") +
  #geom_point(fill = NA, aes(color = Group), size = 1, shape = 16) + 
  geom_line(aes(color = Group)) + 
  scale_color_manual(values = tr.colors) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.025), add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.025), add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = c(0, .2, .4, .6, .8, 1)) +
  #scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = c(0, .2, .4, .6, .8, 1)) +
  ylab('Connectivity loss (%)<br>') +
  xlab('<br>Fraction of nodes removal (%)') +
  guides(color = 
           guide_legend(order = 1, 
                        override.aes = list(shape = 22, 
                                            fill = tr.colors, 
                                            color = "black", 
                                            size = 3,
                                            linetype = 0, 
                                            alpha = 1),
                        ncol = 2)) +
  theme_bw() + 
  my.theme + 
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        panel.spacing.x = unit(5, "mm"),
        legend.title = element_blank(),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black', lineheight = 1.25),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm')) -> net.rob.gg

net.rob.gg

ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_rob.png"),
       plot = net.rob.gg,
       dpi = 1200,
       width = 5000 * 1.75,
       height = 5000 * 0.75,
       units = "px")

ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_rob.svg"),
       plot = net.rob.gg,
       dpi = 1200,
       width = 5000 * 1.75,
       height = 5000 * 0.75,
       units = "px")

### ZiPi plot ------------------------------------------------------------------

get.membership <- function(net) {
  
  new.net <- net
  
  cluster.net <- cluster_louvain(new.net)
  
  V(new.net)$community <- cluster.net$membership
  
  rain <- 
    colorRampPalette(
      paletteer::paletteer_d("awtools::bpalette")[c(1:14, 16)],
      alpha = T)(max(V(new.net)$community)) %>% str_replace_all(., "FF$", "80")
  
  
  V(new.net)$color <- rain[V(new.net)$community]
  
  E(new.net)$color <- apply(as.data.frame(get.edgelist(new.net)), 1, 
                            function(x) ifelse(V(new.net)$community[x[1]] == V(new.net)$community[x[2]], 
                                               rain[V(new.net)$community[x[1]]], '#00000010'))
  
  
  E(new.net)$lty <- apply(as.data.frame(get.edgelist(new.net)), 1, 
                          function(x) ifelse(V(new.net)$community[x[1]] == V(new.net)$community[x[2]], 
                                             1, 2))
  
  return(list(new.net, cluster.net, rain))
  
}

net.1.n <- get.membership(net.1)
net.2.n <- get.membership(net.2)
net.3.n <- get.membership(net.3)
net.4.n <- get.membership(net.4)

png(filename = paste0("task260224_taxProfile6/figures/network/network_membership.png"), width = 5000 * 3, height = 5000 * .75, units = "px", res = 1200, type = "cairo")

par(mfrow = c(1, 4), mai = c(0, 0, 0.5, 0), bg = 'white')

plot(net.1.n[[2]], net.1.n[[1]], 
     col = V(net.1.n[[1]])$color,
     vertex.size = 6,
     #vertex.size = 0,
     vertex.label = "",
     vertex.label.cex = 1 + log2(V(net.1.n[[1]])$degree)/100,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.shape = shapes.p[V(net.1.n[[1]])$IDP2],
     vertex.color = tax.colors2[V(net.1.n[[1]])$IDP],
     edge.color=E(net.1.n[[1]])$color, 
     edge.width = 1,
     edge.arrow.mode = "-", 
     edge.arrow.size=0,
     #edge.curved=0.5,
     layout = l.1, 
     rescale = T,
     main = levels(tb.metad$Group_tribe)[[1]], 
     mark.border=str_replace(net.1.n[[3]], "80$", "00"), 
     mark.col=str_replace(net.1.n[[3]], "80$", "20"),
     margin = c(0, 0, 0, 0))

plot(net.2.n[[2]], net.2.n[[1]], 
     col = V(net.2.n[[1]])$color,
     vertex.size = 6,
     #vertex.size = 0,
     vertex.label = "",
     vertex.label.cex = 1 + log2(V(net.2.n[[1]])$degree)/100,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.shape = shapes.p[V(net.2.n[[1]])$IDP2],
     vertex.color = tax.colors2[V(net.2.n[[1]])$IDP],
     edge.color=E(net.2.n[[1]])$color, 
     edge.width = 1,
     edge.arrow.mode = "-", 
     edge.arrow.size=0,
     #edge.curved=0.5,
     layout = l.2, 
     rescale = T,
     main = levels(tb.metad$Group_tribe)[[2]], 
     mark.border=str_replace(net.2.n[[3]], "80$", "00"), 
     mark.col=str_replace(net.2.n[[3]], "80$", "20"),
     margin = c(0, 0, 0, 0))

plot(net.3.n[[2]], net.3.n[[1]], 
     col = V(net.3.n[[1]])$color,
     vertex.size = 6,
     #vertex.size = 0,
     vertex.label = "",
     vertex.label.cex = 1 + log2(V(net.3.n[[1]])$degree)/100,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.shape = shapes.p[V(net.3.n[[1]])$IDP2],
     vertex.color = tax.colors2[V(net.3.n[[1]])$IDP],
     edge.color=E(net.3.n[[1]])$color, 
     edge.width = 1,
     edge.arrow.mode = "-", 
     edge.arrow.size=0,
     #edge.curved=0.5,
     layout = l.3, 
     rescale = T,
     main = levels(tb.metad$Group_tribe)[[3]], 
     mark.border=str_replace(net.3.n[[3]], "80$", "00"), 
     mark.col=str_replace(net.3.n[[3]], "80$", "20"),
     margin = c(0, 0, 0, 0))

plot(net.4.n[[2]], net.4.n[[1]], 
     col = V(net.4.n[[1]])$color,
     vertex.size = 6,
     #vertex.size = 0,
     vertex.label = "",
     vertex.label.cex = 1 + log2(V(net.4.n[[1]])$degree)/100,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     vertex.shape = shapes.p[V(net.4.n[[1]])$IDP2],
     vertex.color = tax.colors2[V(net.4.n[[1]])$IDP],
     edge.color=E(net.4.n[[1]])$color, 
     edge.width = 1,
     edge.arrow.mode = "-", 
     edge.arrow.size=0,
     #edge.curved=0.5,
     layout = l.4, 
     rescale = T,
     main = levels(tb.metad$Group_tribe)[[4]], 
     mark.border=str_replace(net.4.n[[3]], "80$", "00"), 
     mark.col=str_replace(net.4.n[[3]], "80$", "20"),
     margin = c(0, 0, 0, 0))

# Close the graphics device

dev.off() 

#Module hubs, these nodes have a high degree of internal connectivity within their communities, but few connections with other communities. They are the core nodes within the communities, maintaining the stability and consistency of the communities. These nodes are characterized by Zi > 2.5 and Pi < 0.62.
#Connectors, these nodes have a high degree of internal connectivity within their communities, as well as some external connections. They are important nodes within the communities, as well as bridge nodes between communities, facilitating communication and cooperation among communities. These nodes are characterized by Zi < 2.5 and Pi > 0.62.
#Network hubs, these nodes have a high degree of internal connectivity within their communities, as well as many external connections. They are key nodes in the network, connecting different communities, forming the backbone structure of the network, and having a significant impact on the network’s integrity and efficiency. These nodes are characterized by Zi > 2.5 and Pi > 0.62.
#Peripherals, these nodes have a low degree of internal connectivity within their communities, and almost no external connections. They are marginal nodes within the communities, contributing and influencing little to the communities.

#R1, Ultra-peripheral nodes, these nodes have a high degree of internal connectivity within their communities, but very few connections with other communities. They are isolated nodes within the communities, having almost no contact with other nodes. These nodes are characterized by Zi > 2.5 and Pi < 0.05.
#R2, Peripheral nodes, these nodes have a high degree of internal connectivity within their communities, but few connections with other communities. They are marginal nodes within the communities, contributing and influencing little to the communities. These nodes are characterized by Zi > 2.5 and 0.05 <= Pi < 0.62.
#R3, Non-hub connectors, these nodes have a high degree of internal connectivity within their communities, as well as some external connections. They are connecting nodes within the communities, but not bridge nodes between communities, having some role in communication and cooperation among communities. These nodes are characterized by Zi > 2.5 and Pi > 0.62.
#R4, Non-hub kinless nodes, these nodes have a low degree of internal connectivity within their communities, and almost no external connections. They are kinless nodes within the communities, having no contact with other nodes, nor a sense of belonging to the communities. These nodes are characterized by -2.5 <= Zi <= 2.5 and Pi < 0.05.
#R5, Provincial hubs, these nodes have a low degree of internal connectivity within their communities, but few connections with other communities. They are central nodes within the communities, having a significant impact on the stability and consistency of the communities. These nodes are characterized by -2.5 <= Zi <= 2.5 and 0.05 <= Pi < 0.62.
#R6, Connector hubs, these nodes have a low degree of internal connectivity within their communities, as well as many external connections. They are bridge nodes between communities, connecting different communities, having a significant role in communication and cooperation among communities. These nodes are characterized by -2.5 <= Zi <= 2.5 and Pi > 0.62.
#R7, Kinless hubs, these nodes are the key nodes in the network, connecting different communities and forming the backbone structure of the network, has an important impact on the integrity and efficiency of the network. These nodes are characterized by Zi < -2.5 and Pi > 0.62.

get.zipi <- function(net.w.memb, net.name) {
  
  my.net <- net.w.memb[[1]]
  memb.net <- net.w.memb[[2]]
  
  zp.seqdeg <- degree(my.net); zp.seqdeg
  zp.n.nodes <- length(zp.seqdeg); zp.n.nodes
  
  zp.Z <- zp.seqdeg
  zp.Z[] <- 0
  zp.P <- zp.Z
  
  zp.memb <- membership(memb.net)
  zp.seq <- seq(1:zp.n.nodes)
  
  for (i in 1:zp.n.nodes) {
    
    zp.L <- zp.memb==zp.memb[i]         
    zp.neigh <- neighbors(my.net, i)               
    zp.K.i <- sum(zp.L[zp.neigh])
    
    zp.sum <- 0
    zp.sum.sq <- 0	
    zp.sum.P <- 0
    zp.M.iv <- zp.seq[zp.L]
    
    for (j in 1:sum(zp.L)) {
      
      zp.neigh.j <- neighbors(my.net,zp.M.iv[j])
      zp.K.j <- sum(zp.L[zp.neigh.j])
      
      zp.sum <- zp.sum + zp.K.j
      zp.sum.sq <- zp.sum.sq + zp.K.j^2
      
    }
    
    zp.Z[i] <- (zp.K.i - zp.sum/sum(zp.L))/sqrt(zp.sum.sq/sum(zp.L)-(zp.sum/sum(zp.L))^2)
    
    if (zp.K.i - zp.sum/sum(zp.L) == 0) { zp.Z[i] <- 0 }
    
    for (k in 1:max(zp.memb)) {
      
      zp.L.p <- zp.memb == k
      zp.K.i.p <- sum(zp.L.p[zp.neigh])
      zp.sum.P <- zp.sum.P + (zp.K.i.p/zp.seqdeg[i])^2
      
    }
    
    zp.P[i] <- 1 - zp.sum.P
    
  }
  
  return(data.frame(name = V(net.w.memb[[1]])$name, 
                    Taxa = V(net.w.memb[[1]])$Taxa, 
                    IDP = V(net.w.memb[[1]])$IDP, 
                    IDP2 = V(net.w.memb[[1]])$IDP2, 
                    Phylum = V(net.w.memb[[1]])$Phylum, 
                    IDT = V(net.w.memb[[1]])$IDT, 
                    degree = unname(zp.seqdeg), 
                    module = unname(zp.memb), 
                    Pi = unname(zp.P), 
                    Zi = unname(zp.Z)) %>% 
           mutate(role = case_when(
             Zi > 2.5 & Pi < 0.62 ~ 'Module hubs',
             Zi < 2.5 & Pi > 0.62 ~ 'Connectors',
             Zi > 2.5 & Pi > 0.62 ~ 'Network hubs',
             .default = 'Peripherals')))
  
}

net.zp.1 <- get.zipi(net.1.n, levels(tb.metad$Group_tribe)[[1]]) %>% mutate(Group = levels(tb.metad$Group_tribe)[[1]]); table(net.zp.1$role)
net.zp.2 <- get.zipi(net.2.n, levels(tb.metad$Group_tribe)[[2]]) %>% mutate(Group = levels(tb.metad$Group_tribe)[[2]]); table(net.zp.2$role)
net.zp.3 <- get.zipi(net.3.n, levels(tb.metad$Group_tribe)[[3]]) %>% mutate(Group = levels(tb.metad$Group_tribe)[[3]]); table(net.zp.3$role)
net.zp.4 <- get.zipi(net.4.n, levels(tb.metad$Group_tribe)[[4]]) %>% mutate(Group = levels(tb.metad$Group_tribe)[[4]]); table(net.zp.4$role)

bind_rows(net.zp.1, net.zp.2, net.zp.3, net.zp.4) %>% 
  mutate(my.alpha = ifelse(role != "Peripherals", 1, 1)) %>% 
  
  ggplot(., aes(y = Zi, x = Pi)) +
  geom_star(aes(color = Phylum, starshape = as.factor(IDP2), alpha = my.alpha), size = 1.25) +
  facet_wrap( ~ Group, scales = "fixed", nrow = 2, ncol = 2) +
  scale_y_continuous(limits = c(-4, 4)) +
  scale_x_continuous(limits = c(0, 1)) +
  annotate("richtext", x = Inf, y = -Inf, size = 3, label = paste0("Connectors"), color = "gray60", fill = NA, label.color = NA, hjust = 1, vjust = 0) +
  annotate("richtext", x = Inf, y = Inf, size = 3, label = paste0("Network hubs"), color = "gray60", fill = NA, label.color = NA, hjust = 1, vjust = 1) +
  annotate("richtext", x = -Inf, y = Inf, size = 3, label = paste0("Module hubs"), color = "gray60", fill = NA, label.color = NA, hjust = 0, vjust = 1) +
  annotate("richtext", x = -Inf, y = -Inf, size = 3, label = paste0("Peripherals"), color = "gray60", fill = NA, label.color = NA, hjust = 0, vjust = 0) +
  geom_hline(yintercept = 2.5, linewidth = 0.25, linetype = "dotted", color = "gray40") +
  geom_vline(xintercept = 0.62, linewidth = 0.25, linetype = "dotted", color = "gray40") +
  scale_starshape_manual(values = c(15, 13, 11, 1), name = "Domain") +
  scale_alpha_identity() +
  ylab('Within-module connectivity (Z<sub>i</sub>)<br>') + xlab("<br>Among-module connectivity (P<sub>i</sub>)") + 
  scale_color_manual(name = "", values = tax.colors2) +
  my.theme +
  theme(legend.position = "none", 
        #panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1), 
        plot.margin = unit(c(1, 1, 1, 1), 'cm'))

bind_rows(net.zp.1, net.zp.2, net.zp.3, net.zp.4) %>% 
  mutate(my.alpha = ifelse(role != "Peripherals", 1, 1)) %>% 
  
  ggplot(., aes(y = Zi, x = Pi)) +
  geom_star(aes(fill = Phylum, starshape = as.factor(IDP2), alpha = my.alpha), color = "black", size = 1.5) +
  facet_wrap( ~ Group, scales = "fixed", nrow = 2, ncol = 2) +
  scale_y_continuous(limits = c(-4, 4)) +
  scale_x_continuous(limits = c(0, 1)) +
  annotate("richtext", x = Inf, y = -Inf, size = 3, label = paste0("Connectors"), color = "black", fill = NA, label.color = NA, hjust = 1, vjust = 0) +
  annotate("richtext", x = Inf, y = Inf, size = 3, label = paste0("Network hubs"), color = "black", fill = NA, label.color = NA, hjust = 1, vjust = 1) +
  annotate("richtext", x = -Inf, y = Inf, size = 3, label = paste0("Module hubs"), color = "black", fill = NA, label.color = NA, hjust = 0, vjust = 1) +
  annotate("richtext", x = -Inf, y = -Inf, size = 3, label = paste0("Peripherals"), color = "black", fill = NA, label.color = NA, hjust = 0, vjust = 0) +
  geom_hline(yintercept = 2.5, linewidth = 0.25, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0.62, linewidth = 0.25, linetype = "solid", color = "black") +
  scale_starshape_manual(values = c(15, 13, 11, 1), name = "Domain") +
  scale_alpha_identity() +
  ylab('Within-module connectivity (Z<sub>i</sub>)<br>') + xlab("<br>Among-module connectivity (P<sub>i</sub>)") + 
  scale_fill_manual(name = "", values = tax.colors2) +
  my.theme +
  theme(legend.position = "none", 
        #panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1), 
        plot.margin = unit(c(1, 1, 1, 1), 'cm'))


bind_rows(net.zp.1, net.zp.2, net.zp.3, net.zp.4) %>% 
  mutate(my.alpha = ifelse(role != "Peripherals", 0.6, 0.6)) %>% 
  
  ggplot(., aes(y = Zi, x = Pi)) +
  geom_star(aes(fill = Group, starshape = as.factor(IDP2), alpha = my.alpha), color = "black", size = 2) +
  facet_wrap( ~ Group, scales = "fixed", nrow = 2, ncol = 2) +
  scale_y_continuous(limits = c(-4, 4)) +
  scale_x_continuous(limits = c(0, 1)) +
  annotate("richtext", x = Inf, y = -Inf, size = 3, label = paste0("Connectors"), color = "black", fill = NA, label.color = NA, hjust = 1, vjust = 0) +
  annotate("richtext", x = Inf, y = Inf, size = 3, label = paste0("Network hubs"), color = "black", fill = NA, label.color = NA, hjust = 1, vjust = 1) +
  annotate("richtext", x = -Inf, y = Inf, size = 3, label = paste0("Module hubs"), color = "black", fill = NA, label.color = NA, hjust = 0, vjust = 1) +
  annotate("richtext", x = -Inf, y = -Inf, size = 3, label = paste0("Peripherals"), color = "black", fill = NA, label.color = NA, hjust = 0, vjust = 0) +
  geom_hline(yintercept = 2.5, linewidth = 0.25, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0.62, linewidth = 0.25, linetype = "solid", color = "black") +
  scale_starshape_manual(values = c(15, 13, 11, 1), name = "Domain") +
  scale_alpha_identity() +
  ylab('Within-module connectivity (Z<sub>i</sub>)<br>') + xlab("<br>Among-module connectivity (P<sub>i</sub>)") + 
  scale_fill_manual(name = "", values = tr.colors) +
  my.theme +
  theme(legend.position = "none", 
        #panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1), 
        plot.margin = unit(c(1, 1, 1, 1), 'cm'))


bind_rows(net.zp.1, net.zp.2, net.zp.3, net.zp.4) %>% 
  mutate(my.alpha = ifelse(role != "Peripherals", 1, 0.25)) %>% 
  
  ggplot(., aes(y = Zi, x = Pi)) +
  geom_star(aes(color = Group, starshape = as.factor(IDP2), alpha = my.alpha), size = 1.25) +
  #facet_wrap( ~ Group, scales = "fixed", nrow = 2, ncol = 2) +
  scale_y_continuous(limits = c(-4, 4)) +
  scale_x_continuous(limits = c(0, 1)) +
  annotate("richtext", x = Inf, y = -Inf, size = 3.5, label = paste0("Connectors"), color = "gray30", fill = NA, label.color = NA, hjust = 1, vjust = 0) +
  annotate("richtext", x = Inf, y = Inf, size = 3.5, label = paste0("Network hubs"), color = "gray30", fill = NA, label.color = NA, hjust = 1, vjust = 1) +
  annotate("richtext", x = -Inf, y = Inf, size = 3.5, label = paste0("Module hubs"), color = "gray30", fill = NA, label.color = NA, hjust = 0, vjust = 1) +
  annotate("richtext", x = -Inf, y = -Inf, size = 3.5, label = paste0("Peripherals"), color = "gray30", fill = NA, label.color = NA, hjust = 0, vjust = 0) +
  geom_hline(yintercept = 2.5, linewidth = 0.25, linetype = "solid", color = "gray40") +
  geom_vline(xintercept = 0.62, linewidth = 0.25, linetype = "solid", color = "gray40") +
  scale_starshape_manual(values = c(15, 13, 11, 1), name = "Domain") +
  scale_alpha_identity() +
  ylab('Within-module connectivity (Z<sub>i</sub>)<br>') + xlab("<br>Among-module connectivity (P<sub>i</sub>)") + 
  scale_color_manual(name = "", values = tr.colors) +
  my.theme +
  theme(legend.position = "none", 
        #panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1), 
        plot.margin = unit(c(1, 1, 1, 1), 'cm'))


bind_rows(net.zp.1, net.zp.2, net.zp.3, net.zp.4) %>% 
  mutate(my.alpha = ifelse(role != "Peripherals", .8, .2)) %>% 
  
  ggplot(., aes(y = Zi, x = Pi)) +
  geom_point(aes(fill = Group, shape = as.factor(IDP2), alpha = my.alpha), color = "black", size = 2, stroke = 0.25) +
  #facet_wrap( ~ Group, scales = "fixed", nrow = 2, ncol = 2) +
  scale_y_continuous(limits = c(-4, 4)) +
  scale_x_continuous(limits = c(0, 1)) +
  annotate("richtext", x = Inf, y = -Inf, size = 2, label = paste0("Connectors"), fontface = 'italic', color = "black", fill = NA, label.color = NA, hjust = 1, vjust = 0) +
  annotate("richtext", x = Inf, y = Inf, size = 2, label = paste0("Network hubs"), fontface = 'italic', color = "black", fill = NA, label.color = NA, hjust = 1, vjust = 1) +
  annotate("richtext", x = -Inf, y = Inf, size = 2, label = paste0("Module hubs"), fontface = 'italic', color = "black", fill = NA, label.color = NA, hjust = 0, vjust = 1) +
  annotate("richtext", x = -Inf, y = -Inf, size = 2, label = paste0("Peripherals"), fontface = 'italic', color = "black", fill = NA, label.color = NA, hjust = 0, vjust = 0) +
  geom_hline(yintercept = 2.5, linewidth = 0.25, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0.62, linewidth = 0.25, linetype = "dashed", color = "black") +
  scale_shape_manual(values = c(21, 22), name = "Domain") +
  scale_alpha_identity() +
  ylab('Within-module connectivity (Z<sub>i</sub>)<br>') + xlab("<br>Among-module connectivity (P<sub>i</sub>)") + 
  scale_fill_manual(name = "", values = tr.colors) +
  my.theme +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 1), 'cm')) -> gg.zipi; gg.zipi

ggsave(filename = paste0("task260224_taxProfile6/figures/network/network_zipi.png"),
       plot = gg.zipi,
       dpi = 1200,
       width = 5000 * .9,
       height = 5000 * .9,
       units = "px")

ls.tables <- 
  
  list(
    
    "Net (Summary)" = 
      net.attr.table %>% 
      mutate(Attributes = str_to_title(Attributes) %>% str_replace(., "Of", "of")) %>% 
      mutate(across(everything(), ~as.character(.) %>% str_replace("NA", NA_character_))) %>% 
      data.frame(check.names = F),
    
    "Net (Measures)" = 
      net.measure.table %>% 
      select(Net = Group, Phylum, Taxa, degree:hub_klein, keystone = sig) %>%
      mutate(across(where(is.factor), ~as.character(.))),
    
    "Net Measures (Kruskal)" = 
      net.measure.kruskal %>%
      select(-c(p.full)) %>% 
      dplyr::rename(Measure = measure) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      data.frame(check.names = F),
    
    "Net Measures (Summary)" =
      net.measure.posthoc %>% 
      select(Group, Measure = measure, everything()) %>% 
      mutate(across(where(is.factor), ~as.character(.))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "<.?sup>"))) %>% 
      mutate(across(where(is.character), ~str_remove_all(., "\\*"))),
    
    "Net (Robustness)" = 
      rob.net %>% 
      select(Group, everything()) %>% 
      pivot_wider(., names_from = "attack", values_from = "loss") %>% 
      arrange(Group, node_frac) %>% 
      data.frame(check.names = F),
    
    "Net (ZiPi)" = 
      bind_rows(net.zp.1, net.zp.2, net.zp.3, net.zp.4) %>% 
      select(-c(name, IDP, IDP2, IDT)) %>% 
      data.frame(check.names = F)
    
  )

openxlsx::write.xlsx(x = ls.tables, file = paste0("task260224_taxProfile6/tables/network_all.xlsx"), rowNames = F)
