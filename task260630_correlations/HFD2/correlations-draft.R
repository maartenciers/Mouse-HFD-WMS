

## Taxa gloms ------------------------------------------------------------------


ls.glom.ps <- list(
  
  "Kingdom" = tax_glom(ps.table, taxrank="Kingdom"),
  "Phylum" = tax_glom(ps.table, taxrank="Phylum"),
  "Class" = tax_glom(ps.table, taxrank="Class"),
  "Order" = tax_glom(ps.table, taxrank="Order"),
  "Family" = tax_glom(ps.table, taxrank="Family"),
  "Genus" = tax_glom(ps.table, taxrank="Genus")
  
)

ls.glom.tb <- list(
  
  "Kingdom" = merge(tax_table(ls.glom.ps[["Kingdom"]]), otu_table(ls.glom.ps[["Kingdom"]]), by = "row.names", all = T)[-1],
  "Phylum" = merge(tax_table(ls.glom.ps[["Phylum"]]), otu_table(ls.glom.ps[["Phylum"]]), by = "row.names", all = T)[-1],
  "Class" = merge(tax_table(ls.glom.ps[["Class"]]), otu_table(ls.glom.ps[["Class"]]), by = "row.names", all = T)[-1],
  "Order" = merge(tax_table(ls.glom.ps[["Order"]]), otu_table(ls.glom.ps[["Order"]]), by = "row.names", all = T)[-1],
  "Family" = merge(tax_table(ls.glom.ps[["Family"]]), otu_table(ls.glom.ps[["Family"]]), by = "row.names", all = T)[-1],
  "Genus" = merge(tax_table(ls.glom.ps[["Genus"]]), otu_table(ls.glom.ps[["Genus"]]), by = "row.names", all = T)[-1]
  
)


fun.get.filt.tax.glom <- function(my.glom.ps, min.prev = 0.1, min.max.abs = 0) {
  
  tmp.ps <- phyloseq_filter_prevalence(phyloseq_filter_sample_wise_abund_trim(my.glom.ps, minabund = min.max.abs, relabund = T), prev.trh = min.prev)
  
  return(merge(tax_table(tmp.ps), otu_table(tmp.ps), by = "row.names", all = T)[-1])
  
}


ls.glom.tb.filt <- list(
  
  "Phylum" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001) %>% 
    select(1:2, 8:ncol(.)) %>% group_by(across(c(1:2))) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
    data.frame(check.names = F),
  
  "Class" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001) %>% 
    select(1:3, 8:ncol(.)) %>% group_by(across(c(1:3))) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
    data.frame(check.names = F),
  
  "Order" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001) %>% 
    select(1:4, 8:ncol(.)) %>% group_by(across(c(1:4))) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
    data.frame(check.names = F),
  
  "Family" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001) %>% 
    select(1:5, 8:ncol(.)) %>% group_by(across(c(1:5))) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
    data.frame(check.names = F),
  
  "Genus" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001) %>% 
    select(1:6, 8:ncol(.))
  
)

fun.get.filt.tax.glom.ps <- function(my.glom.ps, min.prev = 0.1, min.max.abs = 0) {
  
  tmp.ps <- phyloseq_filter_prevalence(phyloseq_filter_sample_wise_abund_trim(my.glom.ps, minabund = min.max.abs, relabund = T), prev.trh = min.prev)
  
  return(tmp.ps)
  
}

ls.glom.ps.filt <- list(
  
  "Phylum" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$")) %>% tax_glom(., "Phylum"),
  "Class" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$")) %>% tax_glom(., "Class"),
  "Order" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$")) %>% tax_glom(., "Order"),
  "Family" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$")) %>% tax_glom(., "Family"),
  "Genus" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$"))
  
)

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

top.phylum <- c("Bacillota", "Bacteroidota", "Verrucomicrobiota", "Pseudomonadota", "Actinomycetota")


# MaAsLin2 - DA taxa -----------------------------------------------------------


# Comps:			       
# 1. T1: Bel vs. Yanom	 
# 2. T2: Bel vs. Bic		
# 3. Bel: T1 vs. T2		
# 4. Yanom: T1 vs. T2		

comps <- c("Belgian: pre-HFD|Yanomami: pre-HFD",
           "Belgian: post-HFD|Yanomami: post-HFD",
           "Belgian: pre-HFD|Belgian: post-HFD",
           "Yanomami: pre-HFD|Yanomami: post-HFD")

ls.msl <- list()

for (j in seq_along(comps)) {
  
  
  for (i in seq_along(ranks)) {
    
    tmp.comp.factor <- "donor_ethnicity_timepoint" 
    tmp.comp <- comps[[j]]
    tmp.meta <- tb.metad %>% filter(!!as.symbol(tmp.comp.factor) %in% unlist(str_split(tmp.comp, "\\|"))) %>% column_to_rownames("sample_id")
    tmp.tax.table <- ls.glom.tb.filt[[i]] %>% filter(!str_detect(!!as.symbol(ranks[[i]]), "Un.\\)$")) %>% mutate(tmpID = paste0("tmpID.", 1:n()))
    tmp.rel.tax.id <- tmp.tax.table %>% select(tmpID, any_of(ranks[[i]]))
    tmp.tax.table <- tmp.tax.table %>% column_to_rownames("tmpID") %>% select(any_of(rownames(tmp.meta)))
    
    tmp.msl.comp <- 
      Maaslin2(input_data = tmp.tax.table, 
               input_metadata = tmp.meta, 
               output = "/tmp/",
               normalization = "TSS",
               correction = "BH",
               fixed_effects = c(tmp.comp.factor),
               #random_effects = "cage_id",
               reference = paste0(tmp.comp.factor, ",", str_remove({tmp.comp}, '\\|.*')),
               min_abundance = 0,
               min_prevalence = 0,
               plot_heatmap = F, 
               plot_scatter = F); head(tmp.msl.comp)
    
    ls.msl[[ranks[[i]]]][[j]] <- 
      merge(tmp.rel.tax.id, tmp.msl.comp$results, by.x = "tmpID", by.y = "feature") %>% 
      mutate(enrich_group = ifelse(coef < 0, str_remove({tmp.comp}, '\\|.*'), str_remove({tmp.comp}, '.*\\|'))) %>% 
      mutate(comp = tmp.comp, 
             rank = ranks[[i]]) %>% 
      select(rank, feature = ranks[[i]], enrich_group, coef, stderr, pval, qval, N, N.not.zero, comp)
    
    
  }
  
}

rm(i, j)
rm(comps)
rm(list = ls(pattern = "[.]?tmp[.]?"))

# Filtering

ls.tax.ann <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]] 
  
  ls.tax.ann[[i]] <- 
    tax_table(ps.table) %>% 
    data.frame(check.names = F) %>% 
    remove_rownames() %>% 
    mutate(feature = .[[tmp.rank]]) %>% 
    select(Kingdom, Phylum, feature) %>% 
    distinct()
  
}; rm(i, tmp.rank)

tb.tax.ann <- 
  bind_rows(ls.tax.ann) %>% 
  mutate(Phylum2 = case_when(
    ! Phylum %in% names(colors.tax) & str_detect(Phylum, "Unassigned") ~ "Unassigned",
    ! Phylum %in% names(colors.tax) & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
    ! Phylum %in% names(colors.tax) & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
    Phylum %in% names(colors.tax) ~ Phylum)) %>%
  filter(Phylum2 != "Unassigned") %>% 
  mutate(Phylum2 = factor(Phylum2, names(colors.tax)))

df.msl.raw <- 
  bind_rows(ls.msl) %>% 
  mutate(comp = factor(comp, levels = unique(comp))) %>% 
  mutate(enrich_group = factor(enrich_group, levels = unique(c(levels(tb.metad$donor_ethnicity_timepoint))))) %>% 
  mutate(rank = factor(rank, levels = ranks)) %>% 
  merge(., tb.tax.ann, "feature", all.x = T) %>% 
  select(rank, Kingdom:Phylum2, feature, comp, enrich_group, coef:N.not.zero) %>% 
  arrange(rank, comp, enrich_group, -abs(coef))

df.msl.filt <-
  df.msl.raw %>% 
  filter(qval < 0.05 & abs(coef) > 1) %>% 
  filter(enrich_group != "Others") %>% 
  filter(!str_detect(feature, " Un.\\)$"))

df.msl <- list("Raw" = df.msl.raw, "Filtered" = df.msl.filt)

openxlsx::write.xlsx(x = df.msl, file = paste0(path.tb.dir, "/maaslin.xlsx"), rowNames = F)


# Correlations -----------------------------------------------------------------


## MaAsLin2 - Taxa-Meta --------------------------------------------------------

comps <- c("Belgian: post-HFD|Yanomami: post-HFD")

ls.msl.meta <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.meta.var <- names(tb.metad)[10:ncol(tb.metad)]
  
  for (j in seq_along(tmp.meta.var)) {
    
    tmp.var <- tmp.meta.var[[j]]
    
    tmp.comp.factor <- "donor_ethnicity_timepoint" 
    tmp.comp <- comps[[1]]
    
    tmp.meta <- 
      tb.metad %>% 
      filter(!!as.symbol(tmp.comp.factor) %in% unlist(str_split(tmp.comp, "\\|"))) %>% 
      column_to_rownames("sample_id") %>% 
      select(any_of(c(tmp.comp.factor, tmp.var)))
    
    tmp.tax.table <- 
      ls.glom.ps.filt[[tmp.rank]] %>% 
      subset_samples(., timepoint == "post-HFD") %>% 
      prune_taxa(taxa_sums(.) > 0, .) %>% 
      microbiome::transform(., "total")
    
    tmp.tax.table <- 
      merge(
        tax_table(tmp.tax.table), 
        otu_table(tmp.tax.table), 
        by = "row.names", all = T)[-1] %>% 
      select(Kingdom:!!as.symbol(tmp.rank), 8:ncol(.)) %>% 
      mutate(tmpID = paste0("tmpID.", 1:n()))
    
    tmp.prev.table <- 
      tmp.tax.table %>% 
      select(tmpID, any_of(rownames(tmp.meta))) %>% 
      pivot_longer(cols = -tmpID, names_to = "tmp", values_to = "prev") %>% 
      mutate(prev = ifelse(prev > 0, 1, 0)) %>% 
      merge(tmp.meta, ., by.x = 0, by.y = "tmp") %>% 
      group_by(tmpID, !!as.symbol(tmp.comp.factor)) %>% 
      summarize(prev = sum(prev, na.rm = T)) %>% 
      pivot_wider(names_from = tmp.comp.factor, values_from = "prev") %>% 
      data.frame(check.names = F)
    
    tmp.rel.tax.id <- tmp.tax.table %>% select(tmpID, any_of(ranks[[i]]))
    tmp.tax.table <- tmp.tax.table %>% column_to_rownames("tmpID") %>% select(any_of(rownames(tmp.meta)))
    
    tmp.msl.comp <- 
      Maaslin2(input_data = tmp.tax.table, 
               input_metadata = tmp.meta, 
               output = "/tmp/",
               normalization = "TSS",
               correction = "BH",
               fixed_effects = c(tmp.var),
               #random_effects = "cage_id",
               #reference = paste0(tmp.comp.factor, ",", str_remove({tmp.comp}, '\\|.*')),
               min_abundance = 0.01,
               min_prevalence = 0.5,
               plot_heatmap = F, 
               plot_scatter = F); head(tmp.msl.comp)
    
    ls.msl.meta[[ranks[[i]]]][[j]] <- 
      merge(tmp.rel.tax.id, tmp.msl.comp$results, by.x = "tmpID", by.y = "feature") %>% 
      merge(., tmp.prev.table) %>% 
      mutate(comp = tmp.comp, rank = ranks[[i]]) %>% 
      select(rank, feature = ranks[[i]], metadata, coef, stderr, pval, qval, 
             N, N.not.zero, any_of(unlist(str_split(tmp.comp, pattern = "\\|"))), comp) %>% 
      filter(!str_detect(feature, "Un.\\)$"))
    
  }
  
}

rm(i, j)
rm(comps)
rm(list = ls(pattern = "[.]?tmp[.]?"))

# Filtering

ls.tax.ann <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]] 
  
  ls.tax.ann[[i]] <- 
    tax_table(ps.table) %>% 
    data.frame(check.names = F) %>% 
    remove_rownames() %>% 
    mutate(feature = .[[tmp.rank]]) %>% 
    select(Kingdom, Phylum, feature) %>% 
    distinct()
  
}; rm(i, tmp.rank)

tb.tax.ann <- 
  bind_rows(ls.tax.ann) %>% 
  mutate(Phylum2 = case_when(
    ! Phylum %in% names(colors.tax) & str_detect(Phylum, "Unassigned") ~ "Unassigned",
    ! Phylum %in% names(colors.tax) & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
    ! Phylum %in% names(colors.tax) & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
    Phylum %in% names(colors.tax) ~ Phylum)) %>%
  filter(Phylum2 != "Unassigned") %>% 
  mutate(Phylum2 = factor(Phylum2, names(colors.tax)))

df.msl.meta.raw <- 
  bind_rows(ls.msl.meta) %>% 
  mutate(comp = factor(comp, levels = unique(comp))) %>% 
  mutate(rank = factor(rank, levels = ranks)) %>% 
  mutate(metadata = factor(metadata, levels = names(tb.metad)[10:ncol(tb.metad)])) %>% 
  merge(., tb.tax.ann, "feature", all.x = T) %>% 
  select(rank, Kingdom:Phylum2, feature, metadata, coef:comp) %>% 
  arrange(rank, metadata, -abs(coef), qval)

df.msl.meta.filt <- df.msl.meta.raw %>% filter(qval < 0.1)

df.msl.meta <- list("Raw" = df.msl.meta.raw, "Filtered" = df.msl.meta.filt)

openxlsx::write.xlsx(x = df.msl.meta, file = paste0(path.tb.dir, "/maaslin_meta.xlsx"), rowNames = F)

tmp.corr.flatted <-
  df.msl.meta.raw %>% 
  filter(rank == "Genus") %>% 
  mutate(sig.f = case_when(
    #abs(coef) >=  1 & pval <= 0.1 ~ "*", 
    qval <= 0.01 ~ "**", 
    qval <= 0.1 ~ "*", 
    .default = "")
  ) %>% 
  mutate(r.vjust = ifelse(qval <= 0.1, .75, .5)) %>% 
  mutate(r.tier = case_when(
    coef >=  1 ~ "[1, Inf)",
    coef >=  0.50 & coef <  1 ~ "[0.5, 1)",
    coef >=  0.00 & coef <  0.50 ~ "[0, 0.5)",
    coef >  -0.50 & coef <  0.00 ~ "(−0.5, 0)",
    coef >  -1 & coef <= -0.50 ~ "(−1, −0.5]",
    coef <= -1 ~ "(-Inf, −1]",
  )) %>% 
  mutate(r.color = ifelse(abs(coef) >= 1, "white", "black")) %>% 
  mutate(r.tier = factor(r.tier, c("[1, Inf)", "[0.5, 1)", "[0, 0.5)", "(−0.5, 0)", "(−1, −0.5]", "(-Inf, −1]")))

#tmp.taxa.has.sig <- tmp.corr.flatted %>% filter(sig.f != "") %>% pull(feature) %>% unique()

tmp.taxa.has.sig <-
  df.msl.filt %>%
  filter(rank == "Genus" & comp == "Belgian: post-HFD|Yanomami: post-HFD") %>%
  # mutate(sig = ifelse(qval < 0.05 & abs(coef) > 1, as.character(enrich_group), "Not Sig.")) %>%
  # arrange(sig, -abs(coef), -qval) %>%
  # group_by(sig) %>%
  # mutate(top = 1:n()) %>%
  # ungroup() %>%
  # filter(top <= 10) %>%
  pull(feature)

tmp.corr.flatted.filt <- 
  tmp.corr.flatted %>% 
  #filter(feature %in% tmp.taxa.has.sig) %>% 
  rename(metad = metadata) %>% 
  mutate(feature2 = str_glue("{feature} (B={`Belgian: post-HFD`};Y={`Yanomami: post-HFD`})"))

ggplot(tmp.corr.flatted.filt, aes(x = metad, y = feature2, fill = r.tier)) +
  geom_tile(colour = "white", size = .25) +
  geom_text(aes(label = sig.f, vjust = r.vjust, color = r.color), size = 3, hjust= 0.5, fontface='bold') +
  scale_color_identity() +
  scale_fill_manual(values = c("#67001FFF", "#B2182BFF", "#F4A582FF", "#D1E5F0FF", "#2166ACFF", "#053061FF")) + 
  
  #scale_fill_distiller(palette = "RdBu", na.value="#353535", breaks=c(-1, -0.5, 0, 0.5, 1), limits=c(-1,1)) +
  labs(fill = "Coef. de correlação (r)") + xlab("") + ylab("") +
  #guides(fill = guide_colourbar(barheight = 0.5, barwidth = 7.5, label.position = "bottom", title.position = "top", title.hjust = .5, draw.ulim = FALSE, draw.llim = FALSE)) + 
  my.theme.1 + 
  theme(legend.position = "right", 
        legend.key.height = unit(3, "mm"), 
        axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1))


## MaAsLin2 - Pathway-Meta -----------------------------------------------------


tb.work.fun.ko <- 
  openxlsx::read.xlsx(xlsxFile = str_glue("{getwd()}/task260619_functional-profiler/{my.exp}/tables/profiles/tb-work-ko.xlsx"), sheet = 1)

#KO
ls.info.ko <- list()

ls.info.ko[["def"]] <- read.delim("misc/ko-info/KOs_table_wide.tsv") %>% select(ko_id = L4, D1:D4) %>% mutate_if(is.character, str_trim)
ls.info.ko[["md"]] <- read.delim("misc/ko-info/modules_steps.tsv") %>% mutate_if(is.character, str_trim)

tb.aux.ann.ko <-
  ls.info.ko[["def"]] %>%
  select(ko_class = D1, ko_subclass = D2) %>%
  distinct() %>%
  filter(! ko_class %in% c("Not Included in Pathway or Brite",
                           "Human Diseases",
                           "Brite Hierarchies",
                           "Organismal Systems")) %>%
  mutate(ko_class = factor(ko_class, c("Metabolism",
                                       "Cellular Processes",
                                       "Genetic Information Processing",
                                       "Environmental Information Processing"))) %>%
  arrange(ko_class, ko_subclass)

tb.aux.ann.ko.def <-
  ls.info.ko[["def"]] %>%
  select(ko_id, ko_class = D1, ko_subclass = D2, ko_pathway = D3) %>%
  distinct() %>%
  filter(! ko_class %in% c("Not Included in Pathway or Brite",
                           "Human Diseases",
                           "Brite Hierarchies",
                           "Organismal Systems")) %>%
  mutate(ko_class = factor(ko_class, c("Metabolism",
                                       "Cellular Processes",
                                       "Genetic Information Processing",
                                       "Environmental Information Processing"))) %>%
  arrange(ko_id, ko_class, ko_subclass)

tb.aux.ann.ko.met <-
  ls.info.ko[["def"]] %>%
  select(ko_class = D1, ko_subclass = D2, ko_pathway = D3) %>%
  distinct() %>%
  filter(ko_class %in% c("Metabolism") & ! ko_subclass %in% c("Not included in regular maps")) %>%
  arrange(ko_class, ko_subclass, ko_pathway) %>%
  mutate(across(.cols = ko_class:ko_subclass, .fns = ~factor(., unique(.))))

# Prep

tb.ko.pw.counts <-
  tb.work.fun.ko %>%
  filter(!is.na(ko_pathway)) %>%
  select(ortholog_id, ko_pathway, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_pathway) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_pathway) %>%
  data.frame(check.names = F)

tmp.tb.metad <- tb.metad %>% rename(ID = sample_id)
tmp.tb.count <- tb.ko.pw.counts %>% select(ID = ko_pathway, tmp.tb.metad$ID)

tmp.tb.taxon <-
  tb.ko.pw.counts %>%
  mutate(Kingdom = ko_pathway,
         Phylum = ko_pathway,
         Class = ko_pathway,
         Order = ko_pathway,
         Family = ko_pathway,
         Genus = ko_pathway,
         Species = ko_pathway) %>%
  select(ID = ko_pathway, Kingdom:Species)

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

ps.table.ko.pw <- phyloseq(tmp.tb.count, tmp.tb.taxon, tmp.tb.metad)

rm(list = ls(pattern = "^tmp."))

# #
# 
# ls.ps <-
#   list(
#     "ko_pathway" = ps.table.ko.pw
#   )
# 
# ls.lf <- list()
# 
# ranks <- c("ko_pathway")
# comps <- c("Belgian: post-HFD|Yanomami: post-HFD")
# 
# for (j in seq_along(comps)) {
#   
#   for (i in seq_along(ranks)) {
#     
#     tmp.comp.factor <- "donor_ethnicity_timepoint"
#     
#     tmp.comp <- comps[[j]]
#     
#     tmp.ps.table <-
#       merge_phyloseq(subset_samples(ls.ps[[i]], donor_ethnicity_timepoint == unlist(str_split(tmp.comp, "\\|"))[[1]]),
#                      subset_samples(ls.ps[[i]], donor_ethnicity_timepoint == unlist(str_split(tmp.comp, "\\|"))[[2]])) %>%
#       prune_samples(sample_sums(.) > 0, .)
#     
#     tmp.lf.comp <-
#       tmp.ps.table %>%
#       run_lefse(taxa_rank = "Kingdom",
#                 wilcoxon_cutoff = Inf,
#                 group = tmp.comp.factor,
#                 kw_cutoff = Inf,
#                 multigrp_strat = TRUE,
#                 lda_cutoff = 0) %>%
#       marker_table() %>%
#       data.frame() %>%
#       remove_rownames() %>%
#       mutate(comp = tmp.comp, padj = p.adjust(pvalue, method = "fdr")) %>%
#       mutate(rank = ranks[[i]])
#     
#     ls.lf[[ranks[[i]]]][[j]] <- tmp.lf.comp
#     
#   }
#   
# }; rm(i, j)
# 
# df.lefse.raw <-
#   bind_rows(ls.lf) %>%
#   mutate(comp = factor(comp, levels = unique(comp)))
# 
# df.lefse.raw <-
#   merge(tb.aux.ann.ko.def %>% select(ko_class:ko_pathway) %>% distinct(), df.lefse.raw, by.x = "ko_pathway", by.y = "feature", all.y = T) %>%
#   select(ko_class, ko_subclass, feature = ko_pathway, everything()) %>%
#   arrange(comp, enrich_group, -ef_lda)
# 
# df.lefse.filt <- df.lefse.raw %>% filter(ef_lda > 2, padj < 0.05)
# 
# df.lefse <- list("Raw" = df.lefse.raw, "Filtered LDA>2" = df.lefse.filt)
# 
# openxlsx::write.xlsx(x = df.lefse, file = paste0(path.tb.dir, "/lefse-pw.xlsx"), rowNames = F)


comps <- c("Belgian: post-HFD|Yanomami: post-HFD")

ls.msl.meta <- list()

tmp.meta.var <- names(tb.metad)[10:ncol(tb.metad)]

for (j in seq_along(tmp.meta.var)) {
  
  tmp.var <- tmp.meta.var[[j]]
  
  tmp.comp.factor <- "donor_ethnicity_timepoint" 
  tmp.comp <- comps[[1]]
  
  tmp.meta <- 
    tb.metad %>% 
    filter(!!as.symbol(tmp.comp.factor) %in% unlist(str_split(tmp.comp, "\\|"))) %>% 
    column_to_rownames("sample_id") %>% 
    select(any_of(c(tmp.comp.factor, tmp.var)))
  
  tmp.tax.table <- 
    ps.table.ko.pw %>% 
    subset_samples(., timepoint == "post-HFD") %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    microbiome::transform(., "total")
  
  tmp.tax.table <- 
    merge(
      tax_table(tmp.tax.table), 
      otu_table(tmp.tax.table), 
      by = "row.names", all = T)[-1] %>% 
    select(Species, 8:ncol(.)) %>% 
    mutate(tmpID = paste0("tmpID.", 1:n()))
  
  tmp.prev.table <- 
    tmp.tax.table %>% 
    select(tmpID, any_of(rownames(tmp.meta))) %>% 
    pivot_longer(cols = -tmpID, names_to = "tmp", values_to = "prev") %>% 
    mutate(prev = ifelse(prev > 0, 1, 0)) %>% 
    merge(tmp.meta, ., by.x = 0, by.y = "tmp") %>% 
    group_by(tmpID, !!as.symbol(tmp.comp.factor)) %>% 
    summarize(prev = sum(prev, na.rm = T)) %>% 
    pivot_wider(names_from = tmp.comp.factor, values_from = "prev") %>% 
    data.frame(check.names = F)
  
  tmp.rel.tax.id <- tmp.tax.table %>% select(tmpID, Species)
  tmp.tax.table <- tmp.tax.table %>% column_to_rownames("tmpID") %>% select(any_of(rownames(tmp.meta)))
  
  tmp.msl.comp <- 
    Maaslin2(input_data = tmp.tax.table, 
             input_metadata = tmp.meta, 
             output = "/tmp/",
             normalization = "TSS",
             correction = "BH",
             fixed_effects = c(tmp.var),
             #random_effects = "cage_id",
             #reference = paste0(tmp.comp.factor, ",", str_remove({tmp.comp}, '\\|.*')),
             #min_abundance = 0.001,
             min_prevalence = 0.5,
             plot_heatmap = F, 
             plot_scatter = F); head(tmp.msl.comp)
  
  ls.msl.meta[[j]] <- 
    merge(tmp.rel.tax.id, tmp.msl.comp$results, by.x = "tmpID", by.y = "feature") %>% 
    merge(., tmp.prev.table) %>% 
    mutate(comp = tmp.comp, rank = "KO pw") %>% 
    select(rank, feature = Species, metadata, coef, stderr, pval, qval, 
           N, N.not.zero, any_of(unlist(str_split(tmp.comp, pattern = "\\|"))), comp)
  
}

rm(j)
rm(comps)
rm(list = ls(pattern = "[.]?tmp[.]?"))

# Filtering

df.msl.meta.raw <- 
  bind_rows(ls.msl.meta) %>% 
  mutate(comp = factor(comp, levels = unique(comp))) %>% 
  mutate(metadata = factor(metadata, levels = names(tb.metad)[10:ncol(tb.metad)])) %>% 
  merge(., tb.aux.ann.ko.def %>% select(ko_class:ko_pathway) %>% distinct(), by.x = "feature", by.y = "ko_pathway", all.x = T) %>% 
  select(ko_class, ko_subclass, feature, metadata, coef:comp) %>%
  arrange(metadata, -abs(coef), qval)

df.msl.meta.filt <- df.msl.meta.raw %>% filter(qval < 0.1)

df.msl.meta <- list("Raw" = df.msl.meta.raw, "Filtered" = df.msl.meta.filt)

openxlsx::write.xlsx(x = df.msl.meta, file = paste0(path.tb.dir, "/maaslin_pw-meta.xlsx"), rowNames = F)

tmp.corr.flatted <-
  df.msl.meta.raw %>% 
  mutate(sig.f = case_when(
    #abs(coef) >=  1 & pval <= 0.1 ~ "*", 
    abs(coef) >=  1 & pval <= 0.1 ~ "*", 
    .default = "")
  ) %>% 
  mutate(r.vjust = ifelse(qval <= 0.1, .75, .5)) %>% 
  mutate(r.tier = case_when(
    coef >=  1 ~ "[1, Inf)",
    coef >=  0.50 & coef <  1 ~ "[0.5, 1)",
    coef >=  0.00 & coef <  0.50 ~ "[0, 0.5)",
    coef >  -0.50 & coef <  0.00 ~ "(−0.5, 0)",
    coef >  -1 & coef <= -0.50 ~ "(−1, −0.5]",
    coef <= -1 ~ "(-Inf, −1]",
  )) %>% 
  mutate(r.color = ifelse(abs(coef) >= 1, "white", "black")) %>% 
  mutate(r.tier = factor(r.tier, c("[1, Inf)", "[0.5, 1)", "[0, 0.5)", "(−0.5, 0)", "(−1, −0.5]", "(-Inf, −1]")))

# tmp.taxa.has.sig <- tmp.corr.flatted %>% filter(sig.f != "") %>% pull(feature) %>% unique()

tmp.taxa.has.sig <-
  df.lefse.filt %>%
  filter(ko_subclass == "Xenobiotics biodegradation and metabolism" & comp == "Belgian: post-HFD|Yanomami: post-HFD") %>%
  # mutate(sig = ifelse(qval < 0.05 & abs(coef) > 1, as.character(enrich_group), "Not Sig.")) %>%
  # arrange(sig, -abs(coef), -qval) %>%
  # group_by(sig) %>%
  # mutate(top = 1:n()) %>%
  # ungroup() %>%
  # filter(top <= 10) %>%
  pull(feature)

tmp.corr.flatted.filt <- 
  tmp.corr.flatted %>% 
  #filter(feature %in% tmp.taxa.has.sig) %>% 
  rename(metad = metadata) %>% 
  mutate(feature2 = str_glue("{feature} (B={`Belgian: post-HFD`};Y={`Yanomami: post-HFD`})"))

ggplot(tmp.corr.flatted.filt %>% filter(ko_class == "Environmental Information Processing"), aes(x = metad, y = feature2, fill = r.tier)) +
  facet_grid(ko_subclass~., scales = "free", space = "free") + 
  geom_tile(colour = "white", size = .25) +
  geom_text(aes(label = sig.f, vjust = r.vjust, color = r.color), size = 3, hjust= 0.5, fontface='bold') +
  scale_color_identity() +
  scale_fill_manual(values = c("[1, Inf)" = "#67001FFF", 
                               "[0.5, 1)" = "#B2182BFF",
                               "[0, 0.5)" = "#F4A582FF", 
                               "(−0.5, 0)" = "#D1E5F0FF", 
                               "(−1, −0.5]" = "#2166ACFF", 
                               "(-Inf, −1]" = "#053061FF")) + 
  
  #scale_fill_distiller(palette = "RdBu", na.value="#353535", breaks=c(-1, -0.5, 0, 0.5, 1), limits=c(-1,1)) +
  labs(fill = "Coef. de correlação (r)") + xlab("") + ylab("") +
  #guides(fill = guide_colourbar(barheight = 0.5, barwidth = 7.5, label.position = "bottom", title.position = "top", title.hjust = .5, draw.ulim = FALSE, draw.llim = FALSE)) + 
  my.theme.1 + 
  theme(legend.position = "right", 
        legend.key.height = unit(3, "mm"), 
        axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1))




tmp.tax.table <- 
  ps.table.ko.pw %>% 
  subset_samples(., timepoint == "post-HFD") %>% 
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  microbiome::transform(., "total")

tmp.tax.meta.table <- 
  merge(
    tax_table(tmp.tax.table), 
    otu_table(tmp.tax.table), 
    by = "row.names", all = T)[-1] %>% 
  select(Species, 8:ncol(.)) %>% 
  filter(Species %in% c(tmp.corr.flatted.filt$feature %>% unique())) %>% 
  column_to_rownames("Species") %>% 
  t() %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  merge(., tb.metad)

tmp.tax.meta.table

ggscatter(tmp.tax.meta.table, x = "liver_pathology", y = "Pyruvate metabolism [PATH:ko00620]", add = "reg.line", color = "donor_ethnicity", merge = T) +
  stat_cor() +
  stat_regline_equation()

ggscatter(tmp.tax.meta.table,
          x = "colon_mucus", y = "Pyruvate metabolism [PATH:ko00620]",
          color = "donor_ethnicity", palette = "jco", # Color by groups "cyl"
          shape = "donor_ethnicity" # Change point shape by groups "cyl"
) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(label.x = min(tmp.tax.meta.table$`colon_mucus`, na.rm = T),
           label.y = min(tmp.tax.meta.table$`Pyruvate metabolism [PATH:ko00620]`, na.rm = T), vjust = -2)


# -------------------------------------------------------------------------


# Prep

tb.ko.pw.counts <-
  tb.work.fun.ko %>%
  select(ortholog_id, ko_class, ko_subclass, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_class, ko_subclass) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_class, ko_subclass) %>%
  data.frame(check.names = F)

tmp.tb.metad <- tb.metad %>% rename(ID = sample_id)
tmp.tb.count <- tb.ko.pw.counts %>% select(ID = ko_subclass, tmp.tb.metad$ID)

tmp.tb.taxon <-
  tb.ko.pw.counts %>%
  mutate(Kingdom = ko_subclass,
         Phylum = ko_subclass,
         Class = ko_subclass,
         Order = ko_subclass,
         Family = ko_subclass,
         Genus = ko_subclass,
         Species = ko_subclass) %>%
  select(ID = ko_subclass, Kingdom:Species)

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

ps.table.ko.pw <- phyloseq(tmp.tb.count, tmp.tb.taxon, tmp.tb.metad)

rm(list = ls(pattern = "^tmp."))

#

ls.ps <-
  list(
    "ko_pathway" = ps.table.ko.pw
  )


comps <- c("Belgian: post-HFD|Yanomami: post-HFD")

ls.msl.meta <- list()

tmp.meta.var <- names(tb.metad)[10:ncol(tb.metad)]

for (j in seq_along(tmp.meta.var)) {
  
  tmp.var <- tmp.meta.var[[j]]
  
  tmp.comp.factor <- "donor_ethnicity_timepoint" 
  tmp.comp <- comps[[1]]
  
  tmp.meta <- 
    tb.metad %>% 
    filter(!!as.symbol(tmp.comp.factor) %in% unlist(str_split(tmp.comp, "\\|"))) %>% 
    column_to_rownames("sample_id") %>% 
    select(any_of(c(tmp.comp.factor, tmp.var)))
  
  tmp.tax.table <- 
    ps.table.ko.pw %>% 
    subset_samples(., timepoint == "post-HFD") %>% 
    prune_taxa(taxa_sums(.) > 0, .)
  
  # %>% 
  #   microbiome::transform(., "total")
  
  tmp.tax.table <- 
    merge(
      tax_table(tmp.tax.table), 
      otu_table(tmp.tax.table), 
      by = "row.names", all = T)[-1] %>% 
    select(Species, 8:ncol(.)) %>% 
    mutate(tmpID = paste0("tmpID.", 1:n()))
  
  tmp.prev.table <- 
    tmp.tax.table %>% 
    select(tmpID, any_of(rownames(tmp.meta))) %>% 
    pivot_longer(cols = -tmpID, names_to = "tmp", values_to = "prev") %>% 
    mutate(prev = ifelse(prev > 0, 1, 0)) %>% 
    merge(tmp.meta, ., by.x = 0, by.y = "tmp") %>% 
    group_by(tmpID, !!as.symbol(tmp.comp.factor)) %>% 
    summarize(prev = sum(prev, na.rm = T)) %>% 
    pivot_wider(names_from = tmp.comp.factor, values_from = "prev") %>% 
    data.frame(check.names = F)
  
  tmp.rel.tax.id <- tmp.tax.table %>% select(tmpID, Species)
  tmp.tax.table <- tmp.tax.table %>% column_to_rownames("tmpID") %>% select(any_of(rownames(tmp.meta)))
  
  tmp.msl.comp <- 
    Maaslin2(input_data = tmp.tax.table, 
             input_metadata = tmp.meta, 
             output = "/tmp/",
             normalization = "TSS",
             correction = "BH",
             fixed_effects = c(tmp.var),
             #random_effects = "cage_id",
             #reference = paste0(tmp.comp.factor, ",", str_remove({tmp.comp}, '\\|.*')),
             #min_abundance = 0.005,
             min_abundance = 0,
             min_prevalence = 0.5,
             plot_heatmap = F, 
             plot_scatter = F); head(tmp.msl.comp)
  
  ls.msl.meta[[j]] <- 
    merge(tmp.rel.tax.id, tmp.msl.comp$results, by.x = "tmpID", by.y = "feature") %>% 
    merge(., tmp.prev.table) %>% 
    mutate(comp = tmp.comp, rank = "KO pw") %>% 
    select(rank, feature = Species, metadata, coef, stderr, pval, qval, 
           N, N.not.zero, any_of(unlist(str_split(tmp.comp, pattern = "\\|"))), comp)
  
}

rm(j)
rm(comps)
rm(list = ls(pattern = "[.]?tmp[.]?"))

# Filtering

df.msl.meta.raw <- 
  bind_rows(ls.msl.meta) %>% 
  mutate(comp = factor(comp, levels = unique(comp))) %>% 
  mutate(metadata = factor(metadata, levels = names(tb.metad)[10:ncol(tb.metad)])) %>% 
  merge(., tb.aux.ann.ko.def %>% select(ko_class:ko_subclass) %>% distinct(), by.x = "feature", by.y = "ko_subclass", all.x = T) %>% 
  select(ko_class, feature, metadata, coef:comp) %>%
  arrange(metadata, -abs(coef), qval)

#df.msl.meta.filt <- df.msl.meta.raw %>% filter(qval < 0.1)

#df.msl.meta <- list("Raw" = df.msl.meta.raw, "Filtered" = df.msl.meta.filt)

#openxlsx::write.xlsx(x = df.msl.meta, file = paste0(path.tb.dir, "/maaslin_pw-meta.xlsx"), rowNames = F)

tmp.corr.flatted <-
  df.msl.meta.raw %>% 
  mutate(sig.f = case_when(
    #abs(coef) >=  1 & pval <= 0.1 ~ "*", 
    abs(coef) >=  1 & pval <= 0.1 ~ "*", 
    .default = "")
  ) %>% 
  mutate(r.vjust = ifelse(qval <= 0.1, .75, .5)) %>% 
  mutate(r.tier = case_when(
    coef >=  1 ~ "[1, Inf)",
    coef >=  0.50 & coef <  1 ~ "[0.5, 1)",
    coef >=  0.00 & coef <  0.50 ~ "[0, 0.5)",
    coef >  -0.50 & coef <  0.00 ~ "(−0.5, 0)",
    coef >  -1 & coef <= -0.50 ~ "(−1, −0.5]",
    coef <= -1 ~ "(-Inf, −1]",
  )) %>% 
  mutate(r.color = ifelse(abs(coef) >= 1, "white", "black")) %>% 
  mutate(r.tier = factor(r.tier, c("[1, Inf)", "[0.5, 1)", "[0, 0.5)", "(−0.5, 0)", "(−1, −0.5]", "(-Inf, −1]")))

# tmp.taxa.has.sig <- tmp.corr.flatted %>% filter(sig.f != "") %>% pull(feature) %>% unique()

tmp.taxa.has.sig <-
  df.lefse.filt %>%
  filter(ko_subclass == "Xenobiotics biodegradation and metabolism" & comp == "Belgian: post-HFD|Yanomami: post-HFD") %>%
  # mutate(sig = ifelse(qval < 0.05 & abs(coef) > 1, as.character(enrich_group), "Not Sig.")) %>%
  # arrange(sig, -abs(coef), -qval) %>%
  # group_by(sig) %>%
  # mutate(top = 1:n()) %>%
  # ungroup() %>%
  # filter(top <= 10) %>%
  pull(feature)

tmp.corr.flatted.filt <- 
  tmp.corr.flatted %>% 
  #filter(feature %in% tmp.taxa.has.sig) %>% 
  rename(metad = metadata) %>% 
  mutate(feature2 = str_glue("{feature} (B={`Belgian: post-HFD`};Y={`Yanomami: post-HFD`})"))

ggplot(tmp.corr.flatted.filt, aes(x = metad, y = feature2, fill = r.tier)) +
  facet_grid(ko_class~., scales = "free", space = "free") + 
  geom_tile(colour = "white", size = .25) +
  geom_text(aes(label = sig.f, vjust = r.vjust, color = r.color), size = 3, hjust= 0.5, fontface='bold') +
  scale_color_identity() +
  scale_fill_manual(values = c("[1, Inf)" = "#67001FFF", 
                               "[0.5, 1)" = "#B2182BFF",
                               "[0, 0.5)" = "#F4A582FF", 
                               "(−0.5, 0)" = "#D1E5F0FF", 
                               "(−1, −0.5]" = "#2166ACFF", 
                               "(-Inf, −1]" = "#053061FF")) + 
  
  #scale_fill_distiller(palette = "RdBu", na.value="#353535", breaks=c(-1, -0.5, 0, 0.5, 1), limits=c(-1,1)) +
  labs(fill = "Coef. de correlação (r)") + xlab("") + ylab("") +
  #guides(fill = guide_colourbar(barheight = 0.5, barwidth = 7.5, label.position = "bottom", title.position = "top", title.hjust = .5, draw.ulim = FALSE, draw.llim = FALSE)) + 
  my.theme.1 + 
  theme(legend.position = "right", 
        legend.key.height = unit(3, "mm"), 
        axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1))




tmp.tax.table <- 
  ps.table.ko.pw %>% 
  subset_samples(., timepoint == "post-HFD") %>% 
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  microbiome::transform(., "total")

tmp.tax.meta.table <- 
  merge(
    tax_table(tmp.tax.table), 
    otu_table(tmp.tax.table), 
    by = "row.names", all = T)[-1] %>% 
  select(Species, 8:ncol(.)) %>% 
  filter(Species %in% c(tmp.corr.flatted.filt$feature %>% unique())) %>% 
  column_to_rownames("Species") %>% 
  t() %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  merge(., tb.metad)

tmp.tax.meta.table

ggscatter(tmp.tax.meta.table, x = "liver_pathology", y = "Pyruvate metabolism [PATH:ko00620]", add = "reg.line", color = "donor_ethnicity", merge = T) +
  stat_cor() +
  stat_regline_equation()

ggscatter(tmp.tax.meta.table,
          x = "colon_mucus", y = "Pyruvate metabolism [PATH:ko00620]",
          color = "donor_ethnicity", palette = "jco", # Color by groups "cyl"
          shape = "donor_ethnicity" # Change point shape by groups "cyl"
) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor(label.x = min(tmp.tax.meta.table$`colon_mucus`, na.rm = T),
           label.y = min(tmp.tax.meta.table$`Pyruvate metabolism [PATH:ko00620]`, na.rm = T), vjust = -2)






















tmp.tax.table <- 
  ls.glom.ps.filt[["Genus"]] %>% 
  subset_samples(., timepoint == "post-HFD") %>% 
  prune_taxa(taxa_sums(.) > 0, .) %>% 
  microbiome::transform(., "total")

tmp.tax.meta.table <- 
  merge(
    tax_table(tmp.tax.table), 
    otu_table(tmp.tax.table), 
    by = "row.names", all = T)[-1] %>% 
  select(Genus, 8:ncol(.)) %>% 
  filter(Genus %in% c(tmp.corr.flatted.filt$feature %>% unique())) %>% 
  column_to_rownames("Genus") %>% 
  t() %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  merge(., tb.metad)

tmp.tax.meta.table

ggscatter(tmp.tax.meta.table, x = "IL17_serum", y = "g: Akkermansia", add = "reg.line") +
  stat_cor(label.x = 100, label.y = .2) +
  stat_regline_equation(label.x = 100, label.y = .18)



ggplot(tmp.tax.meta.table, aes(x = "IL17_serum", y = "g: Akkermansia", color = "g: Akkermansia")) +
  geom_point() +
  stat_smooth(aes(fill = group, color = group), method = "lm", formula = formula) +
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = formula
  )

lm.my.xvar <- "ALT_activity"

lm.my.taxa <- tmp.corr.flatted.filt %>% filter(metad == lm.my.xvar & sig.f != "") %>% pull(feature) %>% unique()

lm.tmp.data <-
  tmp.tax.meta.table %>% 
  select(sample_id, donor_ethnicity_timepoint, any_of(lm.my.taxa), any_of(lm.my.xvar)) %>% 
  pivot_longer(cols = all_of(lm.my.taxa), names_to = "taxa", values_to = "abund") #%>% 
#mutate(abund = ifelse(abund == 0, NA_real_, abund))

lm.my.yvar = "abund"

lm.my.fml <- paste(lm.my.yvar, "~", lm.my.xvar)

lm.my.group <- rlang::sym("taxa")

lm.my.wrap.fml <- rlang::new_formula(rhs = lm.my.group, lhs = NULL)

lm.my.nest <- rlang::quo(-!!lm.my.group)

lm.tmp.data %>% 
  nest(!!lm.my.nest) %>% 
  mutate(model = map(data, ~ lm(lm.my.fml, data = .x)),
         adj.r.squared = map_dbl(model, ~ signif(summary(.x)$adj.r.squared, 5)),
         intercept = map_dbl(model, ~ signif(.x$coef[[1]],5)),
         slope = map_dbl(model, ~ signif(.x$coef[[2]], 5)),
         pvalue = map_dbl(model, ~ signif(summary(.x)$coef[2,4], 5)) 
  ) %>% 
  select(-data, -model) %>% 
  left_join(lm.tmp.data) %>% 
  filter(pvalue < 0.1) %>% 
  
  ggplot(aes_string(lm.my.xvar, lm.my.yvar)) +
  geom_point(aes(color = donor_ethnicity_timepoint)) +
  geom_smooth(method = "lm", fullrange = T) +
  facet_wrap(lm.my.wrap.fml, scales = "free_y", nrow = 1) #+
geom_text(aes(3, 40, label = paste("Adj R2 = ", adj.r.squared, "\n",
                                   "Intercept =",intercept, "\n",
                                   "Slope =", slope, "\n",
                                   "P =", pvalue)))








# facetRegression <- function(dat, xvar, yvar, group) {




# }

facetRegression(mpg, "displ", "hwy", "class")









# End ---------------------------------------------------------------------

# res <- corr.test(tmp.clinic %>% column_to_rownames("mouse_n"),
#                  method = "spearman",
#                  adjust = "fdr",
#                  use    = "pairwise")
# 
# res$r
# res$p


## CLR + Spearman correlations -------------------------------------------------


# Taxa/Count table - CLR transformation

# ls.glom.ps.filt.clr <- list(
#   
#   "Phylum" = 
#     ls.glom.ps.filt[["Phylum"]] %>% 
#     subset_samples(., timepoint == "post-HFD") %>% 
#     prune_taxa(taxa_sums(.) > 0, .) %>% 
#     microbiome::transform(., "clr"),
#   
#   "Class" = 
#     ls.glom.ps.filt[["Class"]] %>% 
#     subset_samples(., timepoint == "post-HFD") %>% 
#     prune_taxa(taxa_sums(.) > 0, .) %>% 
#     microbiome::transform(., "clr"),
#   
#   "Order" = 
#     ls.glom.ps.filt[["Order"]] %>% 
#     subset_samples(., timepoint == "post-HFD") %>% 
#     prune_taxa(taxa_sums(.) > 0, .) %>% 
#     microbiome::transform(., "clr"),
#   
#   "Family" = 
#     ls.glom.ps.filt[["Family"]] %>% 
#     subset_samples(., timepoint == "post-HFD") %>% 
#     prune_taxa(taxa_sums(.) > 0, .) %>% 
#     microbiome::transform(., "clr"),
#   
#   "Genus" = 
#     ls.glom.ps.filt[["Genus"]] %>% 
#     subset_samples(., timepoint == "post-HFD") %>% 
#     prune_taxa(taxa_sums(.) > 0, .) %>% 
#     microbiome::transform(., "clr")
#   
# )
# 
# ls.glom.tb.filt.clr <- list(
#   
#   "Phylum" = merge(tax_table(ls.glom.ps.filt.clr[["Phylum"]]), otu_table(ls.glom.ps.filt.clr[["Phylum"]]), by = "row.names", all = T)[-1] %>% select(1:2, 8:ncol(.)),
#   "Class" = merge(tax_table(ls.glom.ps.filt.clr[["Class"]]), otu_table(ls.glom.ps.filt.clr[["Class"]]), by = "row.names", all = T)[-1] %>% select(1:3, 8:ncol(.)),
#   "Order" = merge(tax_table(ls.glom.ps.filt.clr[["Order"]]), otu_table(ls.glom.ps.filt.clr[["Order"]]), by = "row.names", all = T)[-1] %>% select(1:4, 8:ncol(.)),
#   "Family" = merge(tax_table(ls.glom.ps.filt.clr[["Family"]]), otu_table(ls.glom.ps.filt.clr[["Family"]]), by = "row.names", all = T)[-1] %>% select(1:5, 8:ncol(.)),
#   "Genus" = merge(tax_table(ls.glom.ps.filt.clr[["Genus"]]), otu_table(ls.glom.ps.filt.clr[["Genus"]]), by = "row.names", all = T)[-1] %>% select(1:6, 8:ncol(.))
#   
# )
# 
# i = 5
# 
# tmp.rank <- ranks[[i]]
# 
# tmp.taxa.tb <- 
#   ls.glom.tb.filt.clr[[tmp.rank]] %>% 
#   filter(!str_detect(get(tmp.rank), " Un.\\)$")) %>% 
#   select(any_of(tmp.rank), any_of(tb.metad$sample_id)) %>% 
#   column_to_rownames(tmp.rank) %>% 
#   t() %>% 
#   data.frame(check.names = F) %>% 
#   rownames_to_column("sample_id")
# 
# tmp.corr.tb <- merge(tb.metad, tmp.taxa.tb, all.y = T)
# 
# tmp.corr.res <- 
#   corr.test(tmp.corr.tb %>% select(-c(sample_id:donor_ethnicity_timepoint)), 
#             method = "spearman", use = "pairwise", adjust = "fdr", alpha = 0.05)
# 
# tmp.corr.n <- 
#   fun.get.flatten.corr.mt(tmp.corr.res$n, tmp.corr.res$n) %>% 
#   filter(from %in% names(tb.metad)) %>% 
#   filter(!to %in% names(tb.metad)) %>% 
#   select(from, n = r) %>% 
#   distinct()
# 
# tmp.corr.flatted <-
#   fun.get.flatten.corr.mt(tmp.corr.res$r, tmp.corr.res$p) %>% 
#   filter(from %in% names(tb.metad)) %>% 
#   filter(!to %in% names(tb.metad)) %>% 
#   merge(., tmp.corr.n) %>% 
#   mutate(sig.f = case_when(
#     abs(r) >=  0.5 & p <= 0.001 ~ "***", 
#     abs(r) >=  0.5 & p <= 0.01 ~ "**", 
#     abs(r) >=  0.5 & p <= 0.05 ~ "*", 
#     .default = "")
#   ) %>% 
#   mutate(r.vjust = ifelse(p <= 0.05, .65, .5)) %>% 
#   mutate(r.tier = case_when(
#     r >=  0.75                 ~ "[0.75, 1]",
#     r >=  0.50 & r <  0.75     ~ "[0.50, 0.75)",
#     r >=  0.00 & r <  0.50     ~ "[0, 0.50)",
#     r >  -0.50 & r <  0.00     ~ "(−0.50, 0)",
#     r >  -0.75 & r <= -0.50    ~ "(−0.75, −0.50]",
#     r <= -0.75                 ~ "[−1, −0.75]",
#   )) %>% 
#   mutate(r.color = ifelse(abs(r) >= .75, "white", "black")) %>% 
#   mutate(r.tier = factor(r.tier, c("[0.75, 1]", "[0.50, 0.75)", "[0, 0.50)", "(−0.50, 0)", "(−0.75, −0.50]", "[−1, −0.75]"))) %>% 
#   mutate(from = factor(from, levels = names(tb.metad)[10:ncol(tb.metad)]))
# 
# 
# tmp.taxa.has.sig <- tmp.corr.flatted %>% filter(sig.f != "") %>% pull(to) %>% unique()
# 
# tmp.taxa.has.sig <- 
#   df.msl.filt %>% 
#   filter(rank == "Genus" & comp == "Belgian: post-HFD|Yanomami: post-HFD") %>% 
#   mutate(sig = ifelse(qval < 0.05 & abs(coef) > 1, as.character(enrich_group), "Not Sig.")) %>% 
#   arrange(sig, -abs(coef), -qval) %>% 
#   group_by(sig) %>% 
#   mutate(top = 1:n()) %>% 
#   ungroup() %>% 
#   filter(top <= 10) %>% 
#   pull(feature)
# 
# tmp.corr.flatted.filt <- 
#   tmp.corr.flatted %>% 
#   filter(to %in% tmp.taxa.has.sig)
# 
# ggplot(tmp.corr.flatted.filt, aes(x = from, y = to, fill = r.tier)) +
#   geom_tile(colour = "white", size = .25) +
#   geom_text(aes(label = sig.f, vjust = r.vjust, color = r.color), size = 3, hjust= 0.5, fontface='bold') +
#   scale_color_identity() +
#   scale_fill_manual(values = c("#67001FFF", "#B2182BFF", "#F4A582FF", "#D1E5F0FF", "#2166ACFF", "#053061FF")) + 
#   
#   #scale_fill_distiller(palette = "RdBu", na.value="#353535", breaks=c(-1, -0.5, 0, 0.5, 1), limits=c(-1,1)) +
#   labs(fill = "Coef. de correlação (r)") + xlab("") + ylab("") +
#   #guides(fill = guide_colourbar(barheight = 0.5, barwidth = 7.5, label.position = "bottom", title.position = "top", title.hjust = .5, draw.ulim = FALSE, draw.llim = FALSE)) + 
#   my.theme.1 + 
#   theme(legend.position = "right", 
#         legend.key.height = unit(3, "mm"), 
#         axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1))
