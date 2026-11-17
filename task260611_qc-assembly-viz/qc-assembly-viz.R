
# QC and Assembly Qual ---------------------------------------------------------


# Setup ------------------------------------------------------------------------


## Environment and Packages ----------------------------------------------------


rm(list = ls())
graphics.off()
gc()

set.seed("123")

library("ggtext")
library("patchwork")
library("tidyverse")
library("openxlsx")

getwd()
setwd("~/Mouse-HFD-WMS/")
getwd()

list.dirs()

path.wk.dir <- paste0(getwd(), "/task260611_qc-assembly-viz")
path.dt.dir <- paste0(getwd(), "/data")

path.fg.dir <- paste0(path.wk.dir, "/figures")
path.tb.dir <- paste0(path.wk.dir, "/tables")

for (i in c(path.fg.dir, path.tb.dir)) {
  
  if (!dir.exists(i)) { dir.create(i) }
  
}; remove(i)

## Functions -------------------------------------------------------------------

fun.get.scientific <- function(x){
  
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
  
}

fun.get.summary <- function(tb, gp_comp, gp_facet, value_var) {
  
  tb %>% 
    select(all_of(gp_comp), all_of(gp_facet), all_of(value_var)) %>% 
    rename(value_var := !!value_var) %>% 
    group_by(pick(gp_comp), pick(gp_facet)) %>% 
    summarise(temp = list(c(psych::describe(value_var), quantile(value_var, na.rm = T), IQR = EnvStats::iqr(value_var, na.rm = T)))) %>% 
    unnest_wider(temp) %>% 
    select(-c(vars, trimmed, mad)) %>% 
    mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>% 
    data.frame(check.names = F) %>% 
    arrange(pick(gp_facet), pick(gp_comp)) -> p
  
  return(p)
  
}

## Data loading ----------------------------------------------------------------

# Metadata
tb.metad.raw <- 
  read.delim(paste0(getwd(), "/misc/metadata_full.tsv")); head(tb.metad.raw)

# QC data
list.files(paste0(path.dt.dir, "/QC"))

tb.reads.qc.raw <-
  read.delim(paste0(path.dt.dir, "/QC/read_counts.tsv")); head(tb.reads.qc.raw)

tb.asm.raw <-
  read.delim(paste0(path.dt.dir, "/QC/quast.tsv")); head(tb.asm.raw)

tb.align.raw <-
  read.delim(paste0(path.dt.dir, "/QC/alignment_counts.tsv")); head(tb.align.raw)

## Graph settings --------------------------------------------------------------

# Themes

my.theme.1 <- 
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_text(size = 10, colour = 'black', face = "bold"),
        strip.background = element_rect(fill = NA, color = NA), 
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        panel.border = element_rect(color = "black", linewidth = 0, fill = NA),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_text(size = 12, face = 'bold', colour = 'black'),
        axis.title.y = element_text(size = 12, face = 'bold', colour = 'black'),
        axis.title.y.right = element_text(size = 12, face = 'bold', colour = 'black'),
        axis.text.x = element_text(size = 6, colour = 'black'),
        axis.text.y = element_text(size = 6, colour = 'black'),
        
        plot.title = element_text(size = 10, face = 'bold', colour = 'black'),
        plot.subtitle = element_text(size = 6, colour = 'black'),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

my.theme.2 <- 
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
        panel.border = element_rect(color = "black", linewidth = 0, fill = NA),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_markdown(size = 12, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 12, face = 'bold', colour = 'black'),
        axis.title.y.right = element_markdown(size = 12, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black'),
        axis.text.y = element_markdown(size = 6, colour = 'black'),
        
        plot.title = element_markdown(size = 10, face = 'bold', colour = 'black'),
        plot.subtitle = element_markdown(size = 6, colour = 'black'),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

# Colors

colors.gp <- paletteer::paletteer_d("ggthemes::Tableau_20")[1:length(unique(tb.metad.raw$assembly_id))]
names(colors.gp) <- unique(tb.metad.raw$assembly_id) %>% sort()

# Legends

gg.leg.gp <-
  ggplot(tb.metad.raw, aes(x = 1, y = assembly_id, fill = assembly_id)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "Assembly", values = colors.gp) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.gp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.gp)); gg.leg.gp

# Analysis ---------------------------------------------------------------------

## QC before/after -------------------------------------------------------------

# Data wrangling

tb.metad.raw

# Raw -> Long
tb.reads.qc.long <-
  tb.reads.qc.raw %>% 
  mutate(old_id = str_remove(ID, "_run.*")) %>% 
  group_by(old_id) %>% 
  summarise(Raw = sum(Raw, na.rm = T), 
            Trimmed = sum(Trimmed, na.rm = T), .groups = "drop") %>% 
  merge(., tb.metad.raw %>% select(old_id = raw_sample_id_1, sample_id, assembly_id)) %>% 
  arrange(sample_id) %>% 
  select(sample_id, assembly_id, Raw, Trimmed) %>% 
  pivot_longer(cols = -where(is.character), names_to = "var", values_to = "value") %>% 
  mutate(assembly_id = factor(assembly_id, levels = unique(tb.metad.raw$assembly_id) %>% sort())) %>% 
  mutate(var = factor(var, levels = c("Raw", "Trimmed"))) %>% 
  data.frame(check.names = F); head(tb.reads.qc.long)

# Summarize by assembly
tb.reads.qc.summary <-
  tb.reads.qc.long %>% 
  fun.get.summary(., gp_comp = "assembly_id", gp_facet = "var", value_var = "value") %>% 
  arrange(assembly_id); head(tb.reads.qc.summary)

tb.reads.qc.gg <-
  tb.reads.qc.long %>% 
  mutate(x_jit = as.numeric(var)*.5 - .75 + jitter(as.numeric(assembly_id), 0)) %>% 
  mutate(grouping = interaction(sample_id, assembly_id)) %>% 
  mutate(value = value/1000000)

ggplot(tb.reads.qc.gg, aes(x = assembly_id,  y = value,  group = grouping)) + 
  geom_blank() +
  geom_line(aes(x_jit), alpha = .25, linewidth = .25) +
  geom_point(aes(x_jit, fill = assembly_id, shape = var), size = 1.75, alpha = .75) + 
  scale_y_continuous(expand = expansion(mult = c(.1, .1), add = 0),
                     breaks = seq(floor(min(tb.reads.qc.gg$value) / 50) * 50, 
                                  ceiling(max(tb.reads.qc.gg$value) / 50) * 50, 25),
                     labels = round,
                     trans = scales::transform_log10()
  ) +
  scale_x_discrete(expand = expansion(mult = c(.05, .05), add = 0)) +
  scale_shape_manual(values = c(22, 21)) + 
  scale_fill_manual(values = colors.gp) + 
  xlab("") + ylab("Read Counts (million, log<sub>10</sub>)<br>") +
  my.theme.2 + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + gg.leg.gp + plot_layout(widths = c(.75, .25)) -> gg.reads.qc; gg.reads.qc

png(filename = paste0(path.fg.dir, "/QC_read-counts.png"), width = 5000 * 2, height = 5000 * .75, units = "px", res = 1200, type = "cairo")

gg.reads.qc

dev.off()

## Assemblies ------------------------------------------------------------------

tb.asm.raw %>% names()
tb.asm.raw$Assembly %>% unique()

tb.asm.long <- 
  tb.asm.raw %>% 
  pivot_longer(cols = -Assembly, names_to = "assembly_id", values_to = "value") %>% 
  filter(Assembly %in% c("# contigs (>= 0 bp)", "# contigs (>= 5000 bp)", "# contigs (>= 10000 bp)", "# contigs (>= 50000 bp)",
                         "Total length (>= 0 bp)", "Total length (>= 5000 bp)", "Total length (>= 10000 bp)", "Total length (>= 50000 bp)",
                         "Largest contig", "Total length", "N50")) %>%
  mutate(var = factor(Assembly, levels = unique(.$Assembly))) %>% 
  mutate(assembly_id = str_replace_all(assembly_id, "\\.", "-")) %>% 
  select(assembly_id, var, value) %>% 
  data.frame(check.names = F); tb.asm.long


ggplot(tb.asm.long, aes(x = value, y = assembly_id, fill = assembly_id)) +
  geom_col() + 
  geom_text(aes(label = scales::label_number(accuracy = .1, 
                                             scale_cut = scales::cut_short_scale())(value)), 
            hjust = -.15, size = 2.15) + 
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "solid") +
  facet_wrap(vars(var), scales = "free_x") +
  scale_fill_manual(values = colors.gp) + 
  scale_x_continuous(expand = expansion(mult = c(0, .5)), labels = fun.get.scientific) +
  xlab("") + ylab("") +
  my.theme.1 +
  theme(panel.spacing.x = unit(.5, "lines"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + inset_element(gg.leg.gp,align_to = "full",
                                                        left = .825, 
                                                        bottom = 0, 
                                                        right = .95, 
                                                        top = .35) -> gg.asm; gg.asm

png(filename = paste0(path.fg.dir, "/QC_assemblies.png"), width = 5000 * 2.25, height = 5000 * 1.75, units = "px", res = 1200, type = "cairo")

gg.asm

dev.off()

## Alignments ------------------------------------------------------------------

tb.align.long <- 
  tb.align.raw %>% 
  select(old_id = ID, aligned.hq, aligned.any) %>% 
  merge(., tb.metad.raw %>% select(old_id = raw_sample_id_1, sample_id, assembly_id)) %>% 
  pivot_longer(cols = -where(is.character), names_to = "var", values_to = "value") %>% 
  mutate(assembly_id = factor(assembly_id, levels = unique(tb.metad.raw$assembly_id) %>% sort())) %>% 
  mutate(var = case_when(
    var == "aligned.hq" ~ "Aligned (HQ)",
    var == "aligned.any" ~ "Aligned (Any)"
  )) %>% 
  mutate(var = factor(var, levels = c("Aligned (Any)", "Aligned (HQ)"))); head(tb.align.long)

tb.align.gg <-
  tb.align.long %>% 
  mutate(x_jit = as.numeric(var)*.25 - .5 + jitter(as.numeric(assembly_id), .5)) %>% 
  mutate(grouping = interaction(sample_id, assembly_id)) %>% 
  mutate(value = value) %>% 
  group_by(var) %>%
  mutate(mean_var = mean(value))

ggplot(tb.align.gg, aes(x = assembly_id,  y = value)) + 
  geom_blank() +
  geom_hline(aes(yintercept= mean_var), linewidth = 0.5, linetype = "dotted", color = "red") +
  geom_jitter(aes(assembly_id, fill = assembly_id), shape = 21, size = 1.5, alpha = .75, width = .15) + 
  facet_wrap(vars(var), scales = "free_y", axes = "all", ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(.1, .1), add = 0),
                     # breaks = seq(floor(min(tb.align.gg$value) / 10) * 10,
                     #              ceiling(max(tb.align.gg$value) / 10) * 10, 10),
                     # limits = c(floor(min(tb.align.gg$value) / 10) * 10,
                     #            ceiling(max(tb.align.gg$value) / 10) * 10),
                     labels = function(x) round(x, 2)
  ) +
  scale_x_discrete(expand = expansion(mult = c(.05, .05), add = 0)) +
  scale_fill_manual(values = colors.gp) + 
  xlab("") + ylab("Alignment rate (%)<br>") +
  my.theme.2 + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + gg.leg.gp + plot_layout(widths = c(.7, .3)) -> gg.align; gg.align


png(filename = paste0(path.fg.dir, "/QC_alignment.png"), width = 5000 * 1.5, height = 5000 * 1, units = "px", res = 1200, type = "cairo")

gg.align

dev.off()

# rm(list = ls(pattern = "[.]?reads[.]?"))
# rm(list = ls(pattern = "[.]?asm[.]?"))
# rm(list = ls(pattern = "[.]?align[.]?"))

## Annotations -----------------------------------------------------------------

path.dtb.dir <- paste0(getwd(), "/tables")

tb.summary <- tibble(
  
  "class" = c(rep("functional", 12), rep("taxa_contig", 8), rep("taxa_sum", 8)),
  
  "attribute" = c(
    
    "prokka: predicted gene", 
    "prokka: gene name", 
    "prokka: product", 
    "emapper: ortholog", 
    "emapper: gene name", 
    "emapper: COG", 
    "emapper: KO",
    "emapper: GO", 
    "emapper: EC", 
    "dbCAN: 1T",
    "dbCAN: 2T", 
    "dbCAN: 3T",
    
    "taxa: contigs",
    "taxa: Kingdom (c)", 
    "taxa: Phylum (c)", 
    "taxa: Class (c)", 
    "taxa: Order (c)", 
    "taxa: Family (c)", 
    "taxa: Genus (c)", 
    "taxa: Species (c)",
    
    "taxa: sum",
    "taxa: Kingdom (%)", 
    "taxa: Phylum (%)", 
    "taxa: Class (%)", 
    "taxa: Order (%)", 
    "taxa: Family (%)", 
    "taxa: Genus (%)", 
    "taxa: Species (%)"
    
  )) %>% 
  data.frame(check.names = F)

for (i in unique(tb.metad.raw$assembly_id)) {
  
  print(i) 
  
  tmp.tax <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dtb.dir}/{i}/{i}_taxonomy.xlsx"))
  tmp.annot <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dtb.dir}/{i}/{i}_annotation.xlsx"))
  tmp.cz <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dtb.dir}/{i}/{i}_cazymes.xlsx"))
  tmp.counts <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dtb.dir}/{i}/{i}_counts.xlsx"))
  
  tmp.counts.mean <-
    tmp.counts %>% 
    mutate(!!i := rowMeans(pick(where(is.numeric)), na.rm = TRUE)) %>% 
    select(contig_id, any_of(i))
  
  tmp.annot.light <- 
    tmp.annot %>% 
    select(contig_id, 
           gene_id,
           prokka_gene, 
           prokka_product, 
           emapper_seed_ortholog, 
           emapper_Preferred_name,
           emapper_COG_category, 
           emapper_GOs, emapper_EC, 
           emapper_KEGG_ko) %>%
    mutate(prokka_product = str_replace(prokka_product, "hypothetical protein", NA_character_)) %>% 
    mutate(across(starts_with("prokka_"), ~ ifelse(!is.na(.), "annotated", NA_character_))) %>% 
    mutate(across(starts_with("emapper_"), ~ ifelse(!is.na(.), "annotated", NA_character_)))
  
  tmp.cz.light <- 
    tmp.cz %>% 
    select(contig_id, gene_id, dbCAN_n_tools) %>% 
    mutate(dbCAN_1 = if_else(dbCAN_n_tools >= 1, "annotated", NA_character_)) %>% 
    mutate(dbCAN_2 = if_else(dbCAN_n_tools >= 2, "annotated", NA_character_)) %>% 
    mutate(dbCAN_3 = if_else(dbCAN_n_tools >= 3, "annotated", NA_character_))
  
  tmp.annot.cz.light <- 
    merge(tmp.annot.light, tmp.cz.light, all = T) %>% 
    merge(., tmp.counts.mean)
  
  tmp.tax.light <-
    tmp.tax %>%
    mutate(across(Kingdom:Species, ~ ifelse(str_detect(., ": Unassigned"), NA_character_, .))) %>% 
    mutate(across(Kingdom:Species, ~ ifelse(str_detect(., " Un.\\)$"), NA_character_, .))) %>%
    mutate(across(Kingdom:Species, ~ ifelse(!is.na(.), "annotated", NA_character_))) %>% 
    merge(., tmp.counts.mean)
  
  names(tmp.annot.cz.light)
  
  tmp.attr.id <- i
  
  tmp.summary <- tibble(!!i := c(
    
    tmp.annot.cz.light %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(prokka_gene)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(prokka_product)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(emapper_seed_ortholog)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(emapper_Preferred_name)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(emapper_COG_category)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(emapper_KEGG_ko)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(emapper_GOs)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(emapper_EC)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(dbCAN_1)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(dbCAN_2)) %>% nrow(),
    tmp.annot.cz.light %>% filter(!is.na(dbCAN_3)) %>% nrow(),
    
    tmp.tax.light %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Kingdom)) %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Phylum)) %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Class)) %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Order)) %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Family)) %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Genus)) %>% nrow(),
    tmp.tax.light %>% filter(!is.na(Species)) %>% nrow(),
    
    tmp.tax.light %>% pull(i) %>% sum(),
    tmp.tax.light %>% group_by(Kingdom) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Kingdom == "annotated") %>% pull(i),
    tmp.tax.light %>% group_by(Phylum) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Phylum == "annotated") %>% pull(i),
    tmp.tax.light %>% group_by(Class) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Class == "annotated") %>% pull(i),
    tmp.tax.light %>% group_by(Order) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Order == "annotated") %>% pull(i),
    tmp.tax.light %>% group_by(Family) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Family == "annotated") %>% pull(i),
    tmp.tax.light %>% group_by(Genus) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Genus == "annotated") %>% pull(i),
    tmp.tax.light %>% group_by(Species) %>% summarise(!!i := sum(pick(i), na.rm = TRUE)) %>% filter(Species == "annotated") %>% pull(i)
    
  )) %>% 
    data.frame(check.names = F)
  
  tb.summary <- bind_cols(tb.summary, tmp.summary)
  
  rm(list = ls(pattern = "^tmp."))
  
}; remove(i)

tb.summary %>% View()

tb.metad.raw %>% view()

tb.summary.func <- 
  merge(
    tb.summary %>%
      filter(class == "functional") %>% 
      select(-class) %>% 
      column_to_rownames("attribute") %>% 
      t() %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column("assembly_id") %>% 
      pivot_longer(cols = -assembly_id, names_to = "attribute", values_to = "abs"),
    tb.summary %>%
      filter(class == "functional") %>% 
      select(-class) %>% 
      column_to_rownames("attribute") %>% 
      t() %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column("assembly_id") %>% 
      mutate(across(`prokka: predicted gene`:`dbCAN: 3T`, ~ ./`prokka: predicted gene`)) %>% 
      pivot_longer(cols = -assembly_id, names_to = "attribute", values_to = "perc")
  ) %>% 
  mutate(attribute = factor(attribute, levels = unique(tb.summary$attribute)) %>% droplevels()) %>% 
  mutate(class = "functional") %>% 
  arrange(assembly_id, attribute); tb.summary.func

tb.summary.tax.contig <- 
  merge(
    tb.summary %>%
      filter(class == "taxa_contig") %>% 
      select(-class) %>% 
      column_to_rownames("attribute") %>% 
      t() %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column("assembly_id") %>% 
      pivot_longer(cols = -assembly_id, names_to = "attribute", values_to = "abs"),
    tb.summary %>%
      filter(class == "taxa_contig") %>% 
      select(-class) %>% 
      column_to_rownames("attribute") %>% 
      t() %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column("assembly_id") %>% 
      mutate(across(`taxa: contigs`:`taxa: Species (c)`, ~ ./`taxa: contigs`)) %>% 
      pivot_longer(cols = -assembly_id, names_to = "attribute", values_to = "perc")
  ) %>% 
  mutate(attribute = factor(attribute, levels = unique(tb.summary$attribute)) %>% droplevels()) %>% 
  mutate(class = "taxa_contig") %>% 
  arrange(assembly_id, attribute); tb.summary.tax.contig

tb.summary.tax.sum <- 
  merge(
    tb.summary %>%
      filter(class == "taxa_sum") %>% 
      select(-class) %>% 
      column_to_rownames("attribute") %>% 
      t() %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column("assembly_id") %>% 
      pivot_longer(cols = -assembly_id, names_to = "attribute", values_to = "abs"),
    tb.summary %>%
      filter(class == "taxa_sum") %>% 
      select(-class) %>% 
      column_to_rownames("attribute") %>% 
      t() %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column("assembly_id") %>% 
      mutate(across(`taxa: sum`:`taxa: Species (%)`, ~ ./`taxa: sum`)) %>% 
      pivot_longer(cols = -assembly_id, names_to = "attribute", values_to = "perc")
  ) %>% 
  mutate(attribute = factor(attribute, levels = unique(tb.summary$attribute)) %>% droplevels()) %>% 
  mutate(class = "taxa_sum") %>% 
  arrange(assembly_id, attribute); tb.summary.tax.sum

tb.summary.2 <- 
  bind_rows(
    tb.summary.func,
    tb.summary.tax.contig,
    tb.summary.tax.sum
  ) %>% 
  mutate(
    attribute2 = 
      case_when(
        attribute == "prokka: predicted gene" ~ "Predicted genes (w/ PROKKA)", 
        attribute == "prokka: gene name" ~ NA_character_, 
        attribute == "prokka: product" ~ NA_character_, 
        attribute == "emapper: ortholog" ~ "EggNOG orthologs", 
        attribute == "emapper: gene name" ~ NA_character_, 
        attribute == "emapper: COG" ~ "Cluster of Orthologous Groups (COG)", 
        attribute == "emapper: KO" ~ "KEGG Orthology (KO)",
        attribute == "emapper: GO" ~ "Gene Ontology (GO)", 
        attribute == "emapper: EC" ~ "Enzyme Commission (EC) numbers", 
        attribute == "dbCAN: 1T" ~ NA_character_,
        attribute == "dbCAN: 2T" ~ NA_character_, 
        attribute == "dbCAN: 3T" ~ "Carbohydrate-Active enZymes",
        
        attribute == "taxa: contigs" ~ NA_character_,
        attribute == "taxa: Kingdom (c)" ~ "Kingdom", 
        attribute == "taxa: Phylum (c)" ~ "Phylum", 
        attribute == "taxa: Class (c)" ~ "Class", 
        attribute == "taxa: Order (c)" ~ "Order", 
        attribute == "taxa: Family (c)" ~ "Family", 
        attribute == "taxa: Genus (c)" ~ "Genus", 
        attribute == "taxa: Species (c)" ~ "Species",
        
        attribute == "taxa: sum" ~ NA_character_,
        attribute == "taxa: Kingdom (%)" ~ "Kingdom", 
        attribute == "taxa: Phylum (%)" ~ "Phylum", 
        attribute == "taxa: Class (%)" ~ "Class", 
        attribute == "taxa: Order (%)" ~ "Order", 
        attribute == "taxa: Family (%)" ~ "Family", 
        attribute == "taxa: Genus (%)" ~ "Genus", 
        attribute == "taxa: Species (%)" ~ "Species"
      )
  ) %>% 
  mutate(attribute2 = factor(attribute2, levels = rev(unique(.$attribute2))) %>% droplevels()); tb.summary.2


for (i in unique(tb.summary.2$class)) {
  
  tb.summary.gg <- 
    tb.summary.2 %>% 
    filter(class == i & !is.na(attribute2)) %>% 
    droplevels()
  
  gg.summary <- 
    ggplot(tb.summary.gg, aes(x = perc, y = attribute2, fill = attribute2)) +
    geom_col() + 
    geom_text(aes(label = scales::label_number(accuracy = .1, 
                                               scale_cut = scales::cut_short_scale())(abs)), 
              hjust = -.15, size = 2.15) + 
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "solid") +
    facet_wrap(vars(assembly_id), scales = "free_x") +
    scale_fill_manual(values = rev(colorRampPalette(c(paletteer::paletteer_d("ggthemes::Jewel_Bright")))(length(levels(tb.summary.gg$attribute2))))) + 
    scale_x_continuous(
      expand = expansion(mult = c(0, .25)), 
      breaks = seq(0, 1, .25),
      labels = scales::percent_format()) +
    xlab("") + ylab("") +
    my.theme.1 +
    theme(panel.spacing.x = unit(.5, "lines"))
  
  png(filename = paste0(path.fg.dir, "/QC_annot_", i,".png"), width = 5000 * 2.25, height = 5000 * 1.4, units = "px", res = 1200, type = "cairo")
  
  print(gg.summary)
  
  dev.off()
  
}


# Saving tables -----------------------------------------------------------

tb.reads.qc.summary <-
  tb.reads.qc.long %>% 
  fun.get.summary(., gp_comp = "assembly_id", gp_facet = "var", value_var = "value") %>% 
  arrange(assembly_id); head(tb.reads.qc.summary)

tb.reads.qc.summary

tb.asm.raw

tb.align.summary <- 
  tb.align.gg %>% 
  fun.get.summary(., gp_comp = "assembly_id", gp_facet = "var", value_var = "value")
  
tb.align.summary

tb.summary.3 <- 
  tb.summary %>% 
  filter(class == "taxa_sum") %>% 
  mutate(across(where(is.numeric), ~./1000000)) %>% 
  mutate(class = "taxa_perc") %>% 
  filter(!attribute == "taxa: sum") %>% 
  bind_rows(tb.summary %>% filter(class == "functional"), .)

tb.summary.3

ls.save <- list(
  
  "QC" = tb.reads.qc.summary,
  "Assemblies" = tb.asm.raw,
  "Alignment" = tb.align.summary,
  "Annotations" = tb.summary.3
  
)

openxlsx::write.xlsx(x = ls.save, file = paste0(path.tb.dir, "/qc_asm_annotation-stats.xlsx"), rowNames = F)


# End ---------------------------------------------------------------------

