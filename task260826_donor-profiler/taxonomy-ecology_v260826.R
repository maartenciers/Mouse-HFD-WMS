
# Taxonomy/Ecology analysis ----------------------------------------------------


# Setup ------------------------------------------------------------------------


## Environment and Packages ----------------------------------------------------


rm(list = ls())
graphics.off()
gc()

set.seed("123")

library("ggrepel")
library("ggtext")
library("ggpubr")
library("ggvenn")
library("ggalluvial")
library("patchwork")
library("cowplot")
library("svglite")
library("ComplexHeatmap")
library("openxlsx")
library("xlsx")
library("DESeq2")
library("vegan")
library("phyloseq")
library("microbiome")
library("metagMisc")
library("Maaslin2")
library("tidyverse")

getwd()
setwd("~/Mouse-HFD-WMS/")
getwd()

list.dirs()

path.wk.dir <- file.path(getwd(), "task260826_donor-profiler/")
path.dt.dir <- file.path(getwd(), "tables/donor")

path.fg.dir <- file.path(path.wk.dir, "figures")
path.tb.dir <- file.path(path.wk.dir, "tables")

for (i in c(path.fg.dir, path.tb.dir)) {
  
  if (!dir.exists(i)) { dir.create(i, recursive = T) }
  
}; rm(i)

#load(paste0(path.wk.dir, "taxonomy-ecology_260826.RData"))


## Functions -------------------------------------------------------------------


fun.get.scientific <- function(x){
  
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
  
}

fun.add.small.legend <- function(myPlot, pointSize = 3, textSize = 8, spaceLegend = 0.8) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_markdown(size = textSize + 2, face = "bold"),
          legend.text  = element_markdown(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
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

fun.get.kruskal <- function(tb, gp_comp, gp_facet, value_var) {
  
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
    mutate(sig.f = case_when(adj.p.value <= 0.001 ~ "***", 
                             adj.p.value <= 0.01 ~ "**", 
                             adj.p.value <= 0.05 ~ "*", 
                             adj.p.value <= 0.1 ~ "~", 
                             .default = "ns")) %>% 
    select(-data) %>% 
    data.frame(check.names = F) -> p
  
  return(p)
  
}

fun.get.wilcox <- function(tb, gp_comp, gp_facet, value_var, 
                           p_adjust_method = "fdr", conf_level = 0.95,
                           tb.stat = NULL, alpha.global = 0.05) {
  
  gp_comp <- rlang::ensym(gp_comp); value_var <- rlang::ensym(value_var)
  gp_comp_str <- rlang::as_string(gp_comp); value_var_str <- rlang::as_string(value_var)
  
  fml <- reformulate(gp_comp_str, value_var_str)
  
  tb %>%
    nest(data = -all_of(gp_facet)) %>%
    mutate(
      results = map(data, function(d) {
        
        groups <- sort(unique(d[[gp_comp_str]]))
        n_groups <- length(groups)
        
        if (n_groups < 2) return(NULL)
        
        pairs <- if (n_groups == 2) { list(groups) } else { combn(groups, 2, simplify = FALSE) }
        
        map_dfr(pairs, function(gr) {
          sub_dat <- d[d[[gp_comp_str]] %in% gr, ]
          sub_dat[[gp_comp_str]] <- droplevels(factor(sub_dat[[gp_comp_str]]))
          
          # Wilcox test
          wt <- tryCatch(
            wilcox.test(formula = fml, data = sub_dat, 
                        conf.int = TRUE, conf.level = conf_level, correct = TRUE) %>% 
              broom::tidy(),
            error = function(e) {
              tibble(statistic = NA_real_, p.value = NA_real_,
                     conf.low = NA_real_, conf.high = NA_real_,
                     estimate = NA_real_)
            }
          )
          
          # Effect size
          eff <- tryCatch(
            rstatix::wilcox_effsize(sub_dat, formula = fml),
            error = function(e) {
              tibble(effsize = NA_real_, magnitude = NA_character_)
            }
          )
          
          tibble(
            group1 = gr[1],
            group2 = gr[2],
            statistic = wt$statistic,
            p.value = wt$p.value,
            conf.low = wt$conf.low,
            conf.high = wt$conf.high,
            estimate = wt$estimate,
            effect_effsize = eff$effsize,
            effect_magnitude = as.character(eff$magnitude)
          ) %>% 
            mutate(adj.p.value = p.adjust(p.value, method = p_adjust_method))
        })
      })
    ) %>%
    unnest(results) %>%
    mutate(
      p.full = paste0("𝘱-adj: ", ifelse(adj.p.value < 0.001, "<0.001", sprintf("%1.3f", adj.p.value))),
      sig.f = case_when(
        adj.p.value <= 0.001 ~ "***",
        adj.p.value <= 0.01  ~ "**",
        adj.p.value <= 0.05  ~ "*",
        adj.p.value <= 0.1   ~ "~",
        .default = "ns"
      )
    ) %>%
    select(all_of(gp_facet), group1, group2, statistic, p.value, adj.p.value,
           conf.low, conf.high, estimate, effect_effsize, effect_magnitude, 
           p.full, sig.f) %>%
    data.frame(check.names = FALSE) -> out
  
  # Optional: blank out post-hoc if global test is not significant
  if (!is.null(tb.stat)) {
    
    out <- out %>%
      left_join(
        tb.stat %>%
          select(all_of(gp_facet), global.p = adj.p.value),
        by = gp_facet
      ) %>%
      mutate(
        statistic = ifelse(global.p > alpha.global, NA, statistic),
        p.value = ifelse(global.p > alpha.global, NA, p.value),
        adj.p.value = ifelse(global.p > alpha.global, NA, adj.p.value),
        conf.low = ifelse(global.p > alpha.global, NA, conf.low),
        conf.high = ifelse(global.p > alpha.global, NA, conf.high),
        estimate = ifelse(global.p > alpha.global, NA, estimate),
        effect_effsize = ifelse(global.p > alpha.global, NA, effect_effsize),
        effect_magnitude = ifelse(global.p > alpha.global, NA, effect_magnitude),
        p.full = ifelse(global.p > alpha.global, NA, p.full),
        sig.f = ifelse(global.p > alpha.global, NA, sig.f))
    
  }
  
  out <- out %>% 
    mutate(across(where(is.numeric), ~ round(.x, 6))) %>%
    data.frame(check.names = FALSE)
  
  return(out)
  
}

fun.get.brackets <- function(tb, tb.compm, gp_comp, gp_facet, value_var,
                             only.sig = T, alpha = 0.1, label.type = c("stars", "p.adj"), 
                             step.mult = 0.1) {
  
  label.type <- match.arg(label.type)
  
  # base height per facet
  tb.height <- tb %>%
    group_by(.data[[gp_facet]]) %>%
    summarise(
      y.max = max(.data[[value_var]], na.rm = TRUE),
      y.min = min(.data[[value_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y.range = y.max - y.min,
      step = ifelse(y.range == 0, 0.1, y.range * step.mult)
    )
  
  tb.brackets <- tb.compm %>%
    filter(!is.na(adj.p.value)) %>%
    {if (only.sig) filter(., adj.p.value <= alpha) else .} %>%
    left_join(tb.height, by = setNames("var", gp_facet)) %>%
    group_by(.data[[gp_facet]]) %>%
    arrange(group1, group2, .by_group = TRUE) %>%
    mutate(
      y.position = y.max + row_number() * step,
      label = case_when(
        label.type == "stars" & adj.p.value <= 0.001 ~ "***",
        label.type == "stars" & adj.p.value <= 0.01  ~ "**",
        label.type == "stars" & adj.p.value <= 0.05  ~ "*",
        label.type == "stars" & adj.p.value <= 0.1  ~ "~",
        label.type == "stars" ~ "ns",
        label.type == "p.adj" & adj.p.value < 0.001 ~ "< 0.001",
        label.type == "p.adj" ~ paste0(sprintf("%1.3f", round(adj.p.value, 3)))
      )
    ) %>%
    ungroup()
  
  return(tb.brackets)
}

fun.get.brackets.2 <- function(tb, tb.compm, gp_comp, gp_facet, value_var,
                               only.sig = T, alpha = 0.1, label.type = c("stars", "p.adj"), 
                               step.mult = 0.1, sqrt_scale = FALSE) {
  
  label.type <- match.arg(label.type)
  
  # base height per facet: highest (mean + SD) across groups
  
  tb.height <- tb %>%
    group_by(.data[[gp_facet]], .data[[gp_comp]]) %>%
    summarise(
      y.bar = mean(.data[[value_var]], na.rm = TRUE) + sd(.data[[value_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(.data[[gp_facet]]) %>%
    summarise(
      y.max = max(y.bar, na.rm = TRUE),
      y.min = 0,
      .groups = "drop"
    ) %>%
    mutate(
      y.range = y.max - y.min,
      step = ifelse(y.range == 0, 0.1, y.range * step.mult)
    )
  
  tb.brackets <- tb.compm %>%
    filter(!is.na(adj.p.value)) %>%
    {if (only.sig) filter(., adj.p.value <= alpha) else .} %>%
    left_join(tb.height, by = setNames("var", gp_facet)) %>%
    group_by(.data[[gp_facet]]) %>%
    arrange(group1, group2, .by_group = TRUE) %>%
    mutate(
      y.position = if (sqrt_scale) {
        # Step evenly in sqrt space, then transform back
        (sqrt(y.max) + row_number() * sqrt(step))^2
      } else {
        y.max + row_number() * step
      },
      label = case_when(
        label.type == "stars" & adj.p.value <= 0.001 ~ "***",
        label.type == "stars" & adj.p.value <= 0.01  ~ "**",
        label.type == "stars" & adj.p.value <= 0.05  ~ "*",
        label.type == "stars" & adj.p.value <= 0.1  ~ "~",
        label.type == "stars" ~ "ns",
        label.type == "p.adj" & adj.p.value < 0.001 ~ "< 0.001",
        label.type == "p.adj" ~ paste0(sprintf("%1.3f", round(adj.p.value, 3)))
      )
    ) %>%
    ungroup()
  
  return(tb.brackets)
}

fun.upper.whisker <- function(x) {
  boxplot.stats(x)$stats[5]
}

fun.get.boxp.compm.pair <- function(tb, tb.compm, gp_comp, gp_fill = NULL, gp_col = NULL,
                                    gp_facet, value_var, only.sig = F, alpha = 0.1, 
                                    label.type = "stars", color.set = NULL) {
  
  if (is.null(gp_fill)) { gp_fill = gp_comp }
  if (is.null(gp_col)) { gp_col = gp_comp }
  
  tb.brackets <- fun.get.brackets(
    tb = tb,
    tb.compm = tb.compm,
    gp_facet = gp_facet,
    gp_comp = gp_comp,
    value_var = value_var,
    only.sig = only.sig,
    alpha = alpha,
    label.type = label.type
  )
  
  if (nrow(tb.brackets) > 0) {
    
    gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                         fill = .data[[gp_fill]], color = .data[[gp_col]])) +
      facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
      geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
      ggpubr::stat_pvalue_manual(
        tb.brackets,
        vjust = 0.1,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "label",
        tip.length = 0.01,
        bracket.size = 0.25,
        size = 3.5,
        hide.ns = only.sig,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = color.set) +  
      scale_color_manual(values = color.set) +  
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      coord_cartesian(clip = "off") +
      ylab("") + xlab("") +
      my.theme.2 + 
      theme(legend.position = "none",
            panel.spacing = unit(1, "cm"), 
            axis.text.x = element_markdown(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
    
    
    return(gg)
    
  }
  
  gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                       fill = .data[[gp_fill]], color = .data[[gp_col]])) +
    facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
    geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
    scale_fill_manual(values = color.set) +  
    scale_color_manual(values = color.set) +  
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.2 + 
    theme(legend.position = "none",
          panel.spacing = unit(1, "cm"), 
          axis.text.x = element_markdown(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
  
  return(gg)
  
}


fun.get.boxp.compm.group <- function(tb, tb.compm, tb.kruskal, gp_comp,  gp_fill = NULL, gp_col = NULL,
                                     gp_facet, value_var,
                                     only.sig = F, alpha = 0.1, label.type = "stars", color.set = NULL) {
  
  if (is.null(gp_fill)) { gp_fill = gp_comp }
  if (is.null(gp_col)) { gp_col = gp_comp }
  
  tb.brackets <- fun.get.brackets(
    tb = tb,
    tb.compm = tb.compm,
    gp_facet = gp_facet,
    gp_comp = gp_comp,
    value_var = value_var,
    only.sig = only.sig,
    alpha = alpha,
    label.type = label.type
  )
  
  tb.global <- tb.kruskal %>%
    mutate(
      p.full = paste0(
        "KW, 𝘱-adj: ",
        ifelse(adj.p.value < 0.001, "<0.001", sprintf("%1.3f", round(adj.p.value, 3)))
      )
    )
  
  if (nrow(tb.brackets) > 0) {
    
    gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                         fill = .data[[gp_fill]], color = .data[[gp_col]])) +
      facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
      geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
      geom_text(
        data = tb.global,
        aes(x = 0.5, y = Inf, label = p.full),
        inherit.aes = FALSE,
        vjust = 2,
        hjust = 0,
        size = 3.5,
        color = "black"
      ) +
      ggpubr::stat_pvalue_manual(
        tb.brackets,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "label",
        tip.length = 0.01,
        bracket.size = 0.25,
        size = 3.5,
        hide.ns = only.sig,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = color.set) +  
      scale_color_manual(values = color.set) +  
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
      coord_cartesian(clip = "off") +
      ylab("") + xlab("") +
      my.theme.2 + 
      theme(legend.position = "none",
            panel.spacing = unit(1, "cm"), 
            axis.text.x = element_markdown(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
    
    return(gg)
    
  }
  
  gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                       fill = .data[[gp_fill]], color = .data[[gp_col]])) +
    facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
    geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
    geom_text(
      data = tb.global,
      aes(x = 0.5, y = Inf, label = p.full),
      inherit.aes = FALSE,
      vjust = 2,
      hjust = 0,
      size = 2.25,
      color = "black"
    ) +
    
    scale_fill_manual(values = color.set) +  
    scale_color_manual(values = color.set) +  
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.2 + 
    theme(legend.position = "none",
          panel.spacing = unit(1, "cm"), 
          axis.text.x = element_markdown(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
  
  return(gg)
  
}


## Data loading ----------------------------------------------------------------


tb.metad <- 
  read.delim("misc/metadata_donor2.tsv") %>% 
  rename(sample_id = Sample_ID) %>% 
  arrange(Group, sample_id) %>% 
  filter(!sample_id %in% c("XAM_015", "XAM_022", "XAM_031", "XAM_074", "XAM_080")) %>% 
  mutate(across(where(is.character), .fns = ~replace_na(., "Missing data"))) %>% 
  mutate(across(.cols = -sample_id, .fns = ~str_replace_all(., "_", " "))) %>% 
  mutate(across(.cols = Group:Group_tribe, .fns = ~factor(., unique(.)))) %>% 
  mutate(Sex = factor(Sex, levels = c("Female", "Male"))) %>% 
  mutate(Age_group = factor(Age_group, levels = c("Child", "Adolescent", "Young adult", "Adult", "Elders"))) %>% 
  mutate(BMI_status = factor(BMI_status, levels = c("Underweight", "Normal", "Overweight", "Obese", "Missing data"))) %>% 
  group_by(Group_tribe) %>% 
  mutate(Group_tribe = str_glue("{Group_tribe} (n={n()})")) %>%
  group_by(Ethnicity) %>% 
  mutate(Ethnicity = str_glue("{Ethnicity} (n={n()})")) %>% 
  data.frame(check.names = F) %>% 
  mutate(across(.cols = Group:Group_tribe, .fns = ~factor(., unique(.)))) %>% 
  mutate(sample_id = str_replace(sample_id, "_", "-"))

# Annotation + Counts

ls.master <- list()

for (i in unique(tb.metad$Group)) {
  
  tmp.1 <-
    openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_counts.xlsx"), sheet = 1) %>%
    select(contig_id, any_of(tb.metad %>% filter(Group == i) %>% pull(sample_id)))
  
  tmp.2 <- 
    openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_taxonomy.xlsx"), sheet = 1)
  
  ls.master[[i]] <-
    merge(tmp.2, tmp.1, all.y = T, by = c("contig_id")) %>%
    group_by(across(tax_id:Species)) %>%
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE))) %>% 
    data.frame(check.names = F) %>%  
    mutate(across(any_of(tb.metad$sample_id), ~replace_na(., 0)))
  
}; rm(i, tmp.1, tmp.2)


tb.metad <- 
  tb.metad %>% 
  mutate(across(sample_id:Group, ~str_replace(., "BIC", "BA"))) %>% 
  mutate(across(.cols = Group:Group_tribe, .fns = ~factor(., unique(.))))

names(ls.master) <- names(ls.master) %>% str_replace(., "BIC", "BA")

names(ls.master[["BA"]]) <- names(ls.master[["BA"]]) %>% str_replace(., "BIC", "BA")

save.image(paste0(path.wk.dir, "taxonomy-ecology_260826.RData"))


# PS prep -----------------------------------------------------------------


for (i in seq_along(names(ls.master))) {
  
  if (i == 1) {
    
    tb.taxa.counts <- ls.master[[i]]
    
  } else {
    
    tb.taxa.counts <- 
      base::merge(
        tb.taxa.counts, 
        ls.master[[i]], 
        by = c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
        all = T, 
        sort = F)
    
  }
  
}; rm(i)

tb.taxa.counts %>%
  select(tax_id, any_of(tb.metad$sample_id)) %>%
  group_by(tax_id) %>%
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate("prev" = rowSums(.[tb.metad %>% pull(sample_id) %>% sort()] >= 1)) %>%
  filter(prev >= (tb.metad %>% pull(sample_id) %>% length())*.1) %>%
  select(-prev) %>%
  pull(tax_id) -> tmp.filt.id

tb.taxa.counts <- 
  tb.taxa.counts %>% 
  filter(tax_id %in% tmp.filt.id) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)))

tmp.taxa.counts <- tb.taxa.counts

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

sample_data(ps.table)

ps.table.tss <- phyloseq_standardize_otu_abundance(physeq = ps.table, method = 'total')

rm(list = ls(pattern = "^tmp."))


# Graph setup -------------------------------------------------------------


# Themes
my.theme.1 <- 
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_text(size = 10, colour = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(.1, "lines"),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(2.5, "pt"),
        
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_text(size = 10, colour = 'black'),
        axis.title.y = element_text(size = 10, colour = 'black'),
        axis.title.y.right = element_text(size = 10, colour = 'black'),
        axis.text.x = element_text(size = 8, colour = 'black'),
        axis.text.y = element_text(size = 8, colour = 'black'),
        
        plot.title = element_text(size = 10, face = 'bold', colour = 'black'),
        plot.subtitle = element_text(size = 6, colour = 'black'),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

my.theme.2 <- 
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_markdown(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_markdown(size = 10, colour = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(.1, "lines"),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(2.5, "pt"),
        
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_markdown(size = 10, colour = 'black'),
        axis.title.y = element_markdown(size = 10, colour = 'black'),
        axis.title.y.right = element_markdown(size = 10, colour = 'black'),
        axis.text.x = element_markdown(size = 8, colour = 'black'),
        axis.text.y = element_markdown(size = 8, colour = 'black'),
        
        plot.title = element_markdown(size = 10, face = 'bold', colour = 'black'),
        plot.subtitle = element_markdown(size = 6, colour = 'black'),
        plot.margin = unit(c(0.5, 0.75, 0.5, 0.75), 'cm'))

# Colors

colors.et <- c("#FFA500FF", "#0072B2FF")
names(colors.et) <- levels(tb.metad$Ethnicity); colors.et

colors.sx <- c("#7640A9FF", "#659A32FF")
names(colors.sx) <- levels(tb.metad$Sex); colors.sx

colors.ag <- c("#FDDDABFF", "#E7A974FF", "#F65620FF", "#9F2519FF",  "#7C1D6FFF")
names(colors.ag) <- levels(tb.metad$Age_group); colors.ag

colors.tax1 <- c("#B06888FF", "#5080C8FF", "#E8B800FF", "#18C898FF", "#F88048FF", "#C8C8C8FF", "#303028FF")
names(colors.tax1) <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota", "Other: Bacteria", "Other: Domains"); colors.tax1

colors.tax2 <- c("#B06888FF", "#5080C8FF", "#E8B800FF", "#18C898FF", "#F88048FF", "#C8C8C8FF", "#303028FF")
names(colors.tax2) <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota", "Other: Bacteria", "Other: Domains"); colors.tax2

shapes.gt <- c(21, 22, 23, 24)
names(shapes.gt) <- levels(tb.metad$Group_tribe)

# Legends

gg.leg.et <-
  ggplot(names(colors.et) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.et) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5), 
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

gg.leg.et <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.et)); gg.leg.et

gg.leg.sx <-
  ggplot(names(colors.sx) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.sx) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5), 
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

gg.leg.sx <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.sx)); gg.leg.sx

gg.leg.ag <-
  ggplot(names(colors.ag) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.ag) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5), 
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

gg.leg.ag <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.ag)); gg.leg.ag

gg.leg.tax1 <-
  ggplot(names(colors.tax1) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.tax1) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5), 
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

gg.leg.tax1 <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.tax1)); gg.leg.tax1

gg.leg.tax2 <-
  ggplot(names(colors.tax2) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.tax2) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5), 
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

gg.leg.tax2 <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.tax2)); gg.leg.tax2


gg.leg.shape.gt <-
  ggplot(names(shapes.gt) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., shape = .)) + 
  geom_point(size = 2, fill = "black") + 
  scale_shape_manual(values = shapes.gt) + 
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(1, 0, 1, 3)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5), 
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

gg.leg.shape.gt <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.shape.gt))

gg.leg.all <- 
  gg.leg.et + 
  gg.leg.sx + 
  gg.leg.ag + 
  gg.leg.tax1 + 
  gg.leg.tax2 + 
  gg.leg.shape.gt + 
  plot_layout(ncol = 1, heights = c(lapply(list(colors.et, 
                                                colors.sx, 
                                                colors.ag, 
                                                colors.tax1, 
                                                colors.tax2,
                                                shapes.gt), length)))

gg.leg.all

{
  png(filename = paste0(path.fg.dir, "/legends.png"), width = 5000 * .65, height = 5000 * 1.75, units = "px", res = 1200, type = "cairo")
  
  print(gg.leg.all)
  
  dev.off()
}


# Analysis ---------------------------------------------------------------------


## Taxa profile ----------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/profiles"))) { dir.create(str_glue("{path.fg.dir}/profiles"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/profiles"))) { dir.create(str_glue("{path.tb.dir}/profiles"), recursive = T) }


ls.glom.ps <- list(
  
  "Kingdom" = tax_glom(ps.table, taxrank="Kingdom"),
  "Phylum" = tax_glom(ps.table, taxrank="Phylum"),
  "Class" = tax_glom(ps.table, taxrank="Class"),
  "Order" = tax_glom(ps.table, taxrank="Order"),
  "Family" = tax_glom(ps.table, taxrank="Family"),
  "Genus" = tax_glom(ps.table, taxrank="Genus"),
  "Species" = tax_glom(ps.table, taxrank="Species")
  
)

ls.glom.tb <- list(
  
  "Kingdom" = merge(tax_table(ls.glom.ps[["Kingdom"]]), otu_table(ls.glom.ps[["Kingdom"]]), by = "row.names", all = T)[-1],
  "Phylum" = merge(tax_table(ls.glom.ps[["Phylum"]]), otu_table(ls.glom.ps[["Phylum"]]), by = "row.names", all = T)[-1],
  "Class" = merge(tax_table(ls.glom.ps[["Class"]]), otu_table(ls.glom.ps[["Class"]]), by = "row.names", all = T)[-1],
  "Order" = merge(tax_table(ls.glom.ps[["Order"]]), otu_table(ls.glom.ps[["Order"]]), by = "row.names", all = T)[-1],
  "Family" = merge(tax_table(ls.glom.ps[["Family"]]), otu_table(ls.glom.ps[["Family"]]), by = "row.names", all = T)[-1],
  "Genus" = merge(tax_table(ls.glom.ps[["Genus"]]), otu_table(ls.glom.ps[["Genus"]]), by = "row.names", all = T)[-1],
  "Species" = merge(tax_table(ls.glom.ps[["Species"]]), otu_table(ls.glom.ps[["Species"]]), by = "row.names", all = T)[-1]
  
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
    select(1:6, 8:ncol(.)),
  
  "Species" = fun.get.filt.tax.glom(my.glom.ps = subset_taxa(ls.glom.ps[["Species"]]), min.prev = 0.1, min.max.abs = 0.0001) %>% 
    select(1:7, 8:ncol(.))
  
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
  "Genus" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Genus"]]), min.prev = 0.1, min.max.abs = 0.001), !str_detect(Phylum, "p: Unassigned| Un.\\)$")),
  "Species" = subset_taxa(fun.get.filt.tax.glom.ps(my.glom.ps = subset_taxa(ls.glom.ps[["Species"]]), min.prev = 0.1, min.max.abs = 0.0001), !str_detect(Phylum, "p: Unassigned| Un.\\)$"))
  
)

ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

top.phylum <- c("Bacillota", "Bacteroidota", "Verrucomicrobiota", "Pseudomonadota", "Actinomycetota")

save.image(paste0(path.wk.dir, "taxonomy-ecology_260826.RData"))

### Profile: Samples -----------------------------------------------------------


gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.tb.glom <- ls.glom.tb.filt[[i]]
  tmp.tb.glom <- 
    tmp.tb.glom %>% 
    filter(str_detect(Phylum, "^p: ")) %>% 
    mutate(Phylum = str_remove(Phylum, "p: ")) %>% 
    mutate(taxa = .[[tmp.rank]]) %>% 
    select(Kingdom, Phylum, taxa, any_of(tb.metad$sample_id)) %>% 
    mutate(Phylum2 = ifelse(Phylum %in% top.phylum, Phylum, "Others")) %>% 
    mutate(Phylum2 = factor(Phylum2, levels = c(top.phylum, "Others"))) %>% 
    mutate(across(where(is.numeric), ~./sum(.))) %>% 
    rowwise() %>%
    mutate(mean_of_cols = mean(c_across(any_of(tb.metad$sample_id)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(taxa = make.unique(taxa) %>% 
             str_replace_all(., pattern = "__", replacement = ": ") %>% 
             str_replace_all(., pattern = "_", replacement = " ")) %>% 
    select(Kingdom, starts_with("Phylum"), taxa, any_of(tb.metad$sample_id), mean_of_cols) %>% 
    arrange(., Phylum2, -mean_of_cols); tmp.tb.glom
  
  tmp.tb.top.ind <- tmp.tb.glom %>% select(Phylum2, taxa, mean_of_cols)  %>% filter(!str_detect(taxa, "Un.\\)")) %>% arrange(-mean_of_cols) %>% pull(taxa) %>% .[1]
  
  tmp.tb.top.taxa <- 
    bind_rows(
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[1]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(10),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[2]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(5),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[3]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[4]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[5]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(3)
      
    ); tmp.tb.top.taxa
  
  tmp.tb.top.taxa <- 
    bind_rows(
      tmp.tb.top.taxa, 
      tmp.tb.top.taxa %>% 
        select(Phylum2) %>% 
        distinct() %>% 
        mutate(taxa = str_glue("Other: {Phylum2}")) %>%
        add_row(Phylum2 = "Others", taxa = "Other: Bacteria") %>% 
        add_row(Phylum2 = "Others", taxa = "Other: Domains")) %>% 
    mutate(Phylum2 = factor(factor(Phylum2, levels = c(top.phylum, "Others")))) %>% 
    arrange(Phylum2) %>% 
    mutate(taxa = factor(taxa, levels = unique(taxa)))
  
  tmp.tb.top.taxa <- 
    merge(
      tmp.tb.top.taxa, 
      data.frame(
        
        taxa  = 
          c(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[5]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == "Others") %>% pull(taxa) %>% as.character),
        
        my.colors  = c(
          
          c("#B06888FF", "#884060FF", "#D0A0A8FF", 
            "#A090C8FF", "#684870FF", "#C8C0F8FF",
            "#C24841FF", "#A40000FF", "#FF0000FF",
            "#8810B0FF", "#381858FF", "#F0D8E8FF",
            "#ED444AFF", "#C8133BFF", "#F88080FF",
            "#F8B8C0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]]))],
          
          c("#5080C8FF", "#80D8F8FF", "#182848FF", 
            "#48A8B8FF", "#507070FF",
            "#B8D8E0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]]))],
          
          c("#E8B800FF", "#A89E5EFF", "#735231FF",
            "#FAE9A0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]]))],
          
          c("#18C898FF", "#90E880FF", "#508858FF", 
            "#E0F898FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]]))],
          
          c("#F88048FF", "#A13601FF", "#FF3200FF",
            "#F8B880FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[5]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[5]]))],
          
          
          c("#C8C8C8FF", "#303028FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == "Others")) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == "Others"))]
        )
      )
    ) %>% arrange(taxa) %>% 
    
    mutate(my.colors = case_when(
      taxa == "Other: Bacillota" ~ "#F8B8C0FF",
      taxa == "Other: Bacteroidota" ~ "#B8D8E0FF",
      taxa == "Other: Verrucomicrobiota" ~ "#FAE9A0FF",
      taxa == "Other: Pseudomonadota" ~ "#E0F898FF",
      taxa == "Other: Actinomycetota" ~ "#F8B880FF",
      .default = my.colors))
  
  
  tmp.tb.glom %>% 
    mutate(taxa2 = case_when(
      Phylum2 == "Others" & str_detect(Phylum, "Unassigned") ~ "Unassigned",
      Phylum2 == "Others" & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
      Phylum2 == "Others" & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
      Phylum2 != "Others" & taxa %in% paste(levels(tmp.tb.top.taxa$taxa)) ~ taxa,
      Phylum2 != "Others" & !taxa %in% paste(levels(tmp.tb.top.taxa$taxa)) ~ str_glue("Other: {Phylum2}")
    )) %>%
    mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.top.taxa$Phylum2))) %>% 
    mutate(taxa2 = factor(taxa2, levels(tmp.tb.top.taxa$taxa))) %>% 
    select(Kingdom, starts_with("Phylum"), starts_with("taxa"), taxa, any_of(tb.metad$sample_id), mean_of_cols) %>% 
    group_by(Phylum2, taxa2) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") -> tmp.tb.glom.c
  
  tmp.tb.glom.c %>% 
    mutate(across(tb.metad$sample_id, ~as.numeric(.))) %>% 
    mutate(across(where(is.numeric), ~./sum(.))) %>% 
    filter(complete.cases(.)) %>% 
    mutate(across(tb.metad$sample_id, ~as.numeric(.))) %>% 
    merge(x = reshape2::melt(.), y = tb.metad, by.x = "variable", by.y = "sample_id") %>% 
    arrange(., taxa2, -value) -> tmp.tb.glom.top
  
  tmp.tb.top.id <- tmp.tb.glom.top %>% filter(taxa2 == tmp.tb.top.ind) %>% select(ID = variable, value) %>% pull("ID") %>% as.character()
  
  tmp.tb.glom.top <- tmp.tb.glom.top %>% mutate(ID = factor(variable, levels = tmp.tb.top.id))
  
  tmp.tb.top.taxa <- tmp.tb.top.taxa %>% filter(taxa %in% tmp.tb.glom.top$taxa2)
  
  tmp.colors <- tmp.tb.top.taxa$my.colors
  names(tmp.colors) <- tmp.tb.top.taxa$taxa
  
  if (i > 1) {
    
    tmp.ls.gg.leg.taxa <- list()
    
    for (p in levels(tmp.tb.top.taxa$Phylum2)) {
      
      tmp.tb.gg.leg.taxa <- tmp.tb.top.taxa %>% filter(Phylum2 == p)
      
      tmp.ls.gg.leg.taxa[[p]] <- ggpubr::as_ggplot(ggpubr::get_legend(
        
        ggplot(tmp.tb.gg.leg.taxa, aes(x = 1, y = taxa, fill = taxa)) + 
          geom_point(shape = 22, size = 3) + 
          scale_fill_manual(name = p, values = tmp.tb.gg.leg.taxa$my.colors) +  
          guides(shape = guide_legend(override.aes = list(size = 2)),
                 color = guide_legend(override.aes = list(size = 2)),
                 fill = guide_legend(ncol = 1, override.aes = list(size = 2))) +
          my.theme.1 + 
          theme(legend.position = "right",
                legend.title = element_text(size = 8, colour = 'black', face = "bold"),
                legend.text = element_text(size = 6, colour = 'black', margin = margin(0, 0, 0, 2)),
                legend.key.size = unit(6, 'pt'),
                legend.justification = c(0, 1))
      ))
      
    }; rm(p)
    
    gg.bp.leg <- (wrap_plots(tmp.ls.gg.leg.taxa) + plot_layout(nrow = 1))
    
  } else {
    
    gg.bp.leg <-
      ggplot(tmp.tb.top.taxa %>% group_by(taxa) %>% slice_head(n = 1), aes(x = 1, y = taxa, fill = taxa)) + 
      geom_point(shape = 22, size = 3) + 
      scale_fill_manual(name = "Phylum", values = tmp.tb.top.taxa %>% group_by(taxa) %>% slice_head(n = 1) %>% pull(my.colors)) +  
      my.theme.1 + 
      theme(legend.position = "right",
            legend.title = element_text(size = 8, colour = 'black', face = "bold"),
            legend.text = element_text(size = 7, colour = 'black', margin = margin(0, 0, 0, 2)),
            legend.key.size = unit(.5, 'lines'),
            legend.justification = c(0, .5))
    
    gg.bp.leg <- ggpubr::as_ggplot(ggpubr::get_legend(gg.bp.leg)); gg.bp.leg
    
  }
  
  gg.bp <- 
    ggplot(tmp.tb.glom.top %>% mutate(Ethnicity = str_remove(Ethnicity, " \\(.*")), aes(fill = taxa2, x = variable, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = 1) +
    facet_grid(. ~ Ethnicity, space = "free_x", scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    my.theme.2 + 
    theme(strip.background.x = element_rect(fill = NA, colour = NA), 
          strip.text.x = element_markdown(size = 8, face = "plain"),
          panel.spacing.x = unit(.25, "lines"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()); gg.bp
  
  gg.bp.plots[[i]] <- gg.bp + gg.bp.leg + plot_layout(ncol = 2, widths = c(.35, .75))
  
}; rm(i, tmp.rank)

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p
gg.bp.c <- gg.bp.plots[[2]]; gg.bp.c
gg.bp.o <- gg.bp.plots[[3]]; gg.bp.o
gg.bp.f <- gg.bp.plots[[4]]; gg.bp.f
gg.bp.g <- gg.bp.plots[[5]]; gg.bp.g


### Profile: Group -------------------------------------------------------------


gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.tb.glom <- ls.glom.tb.filt[[i]]
  tmp.tb.glom <- 
    tmp.tb.glom %>% 
    filter(str_detect(Phylum, "^p: ")) %>% 
    mutate(Phylum = str_remove(Phylum, "p: ")) %>% 
    mutate(taxa = .[[tmp.rank]]) %>% 
    select(Kingdom, Phylum, taxa, any_of(tb.metad$sample_id)) %>% 
    mutate(Phylum2 = ifelse(Phylum %in% top.phylum, Phylum, "Others")) %>% 
    mutate(Phylum2 = factor(Phylum2, levels = c(top.phylum, "Others"))) %>% 
    mutate(across(where(is.numeric), ~./sum(.))) %>% 
    rowwise() %>%
    mutate(mean_of_cols = mean(c_across(any_of(tb.metad$sample_id)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(taxa = make.unique(taxa) %>% 
             str_replace_all(., pattern = "__", replacement = ": ") %>% 
             str_replace_all(., pattern = "_", replacement = " ")) %>% 
    select(Kingdom, starts_with("Phylum"), taxa, any_of(tb.metad$sample_id), mean_of_cols) %>% 
    arrange(., Phylum2, -mean_of_cols); tmp.tb.glom
  
  tmp.tb.top.ind <- tmp.tb.glom %>% select(Phylum2, taxa, mean_of_cols)  %>% filter(!str_detect(taxa, "Un.\\)")) %>% arrange(-mean_of_cols) %>% pull(taxa) %>% .[1]
  
  tmp.tb.top.taxa <- 
    bind_rows(
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[1]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(10),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[2]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(5),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[3]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[4]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[5]]) %>% filter(!str_detect(taxa, "Un.\\)")) %>% head(3)
      
    ); tmp.tb.top.taxa
  
  tmp.tb.top.taxa <- 
    bind_rows(
      tmp.tb.top.taxa, 
      tmp.tb.top.taxa %>% 
        select(Phylum2) %>% 
        distinct() %>% 
        mutate(taxa = str_glue("Other: {Phylum2}")) %>%
        add_row(Phylum2 = "Others", taxa = "Other: Bacteria") %>% 
        add_row(Phylum2 = "Others", taxa = "Other: Domains")) %>% 
    mutate(Phylum2 = factor(factor(Phylum2, levels = c(top.phylum, "Others")))) %>% 
    arrange(Phylum2) %>% 
    mutate(taxa = factor(taxa, levels = unique(taxa)))
  
  tmp.tb.top.taxa <- 
    merge(
      tmp.tb.top.taxa, 
      data.frame(
        
        taxa  = 
          c(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[5]]) %>% pull(taxa) %>% as.character,
            tmp.tb.top.taxa %>% filter(Phylum2 == "Others") %>% pull(taxa) %>% as.character),
        
        my.colors  = c(
          
          c("#B06888FF", "#884060FF", "#D0A0A8FF", 
            "#A090C8FF", "#684870FF", "#C8C0F8FF",
            "#C24841FF", "#A40000FF", "#FF0000FF",
            "#8810B0FF", "#381858FF", "#F0D8E8FF",
            "#ED444AFF", "#C8133BFF", "#F88080FF",
            "#F8B8C0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]]))],
          
          c("#5080C8FF", "#80D8F8FF", "#182848FF", 
            "#48A8B8FF", "#507070FF",
            "#B8D8E0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]]))],
          
          c("#E8B800FF", "#A89E5EFF", "#735231FF",
            "#FAE9A0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]]))],
          
          c("#18C898FF", "#90E880FF", "#508858FF", 
            "#E0F898FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]]))],
          
          c("#F88048FF", "#A13601FF", "#FF3200FF",
            "#F8B880FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[5]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[5]]))],
          
          
          c("#C8C8C8FF", "#303028FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == "Others")) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == "Others"))]
        )
      )
    ) %>% arrange(taxa) %>% 
    
    mutate(my.colors = case_when(
      taxa == "Other: Bacillota" ~ "#F8B8C0FF",
      taxa == "Other: Bacteroidota" ~ "#B8D8E0FF",
      taxa == "Other: Verrucomicrobiota" ~ "#FAE9A0FF",
      taxa == "Other: Pseudomonadota" ~ "#E0F898FF",
      taxa == "Other: Actinomycetota" ~ "#F8B880FF",
      .default = my.colors))
  
  
  tmp.tb.glom %>% 
    mutate(taxa2 = case_when(
      Phylum2 == "Others" & str_detect(Phylum, "Unassigned") ~ "Unassigned",
      Phylum2 == "Others" & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
      Phylum2 == "Others" & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
      Phylum2 != "Others" & taxa %in% paste(levels(tmp.tb.top.taxa$taxa)) ~ taxa,
      Phylum2 != "Others" & !taxa %in% paste(levels(tmp.tb.top.taxa$taxa)) ~ str_glue("Other: {Phylum2}")
    )) %>%
    mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.top.taxa$Phylum2))) %>% 
    mutate(taxa2 = factor(taxa2, levels(tmp.tb.top.taxa$taxa))) %>% 
    select(Kingdom, starts_with("Phylum"), starts_with("taxa"), taxa, any_of(tb.metad$sample_id), mean_of_cols) %>% 
    group_by(Phylum2, taxa2) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") -> tmp.tb.glom.c
  
  tmp.tb.glom.c %>% 
    mutate(across(tb.metad$sample_id, ~as.numeric(.))) %>% 
    mutate(across(where(is.numeric), ~./sum(.))) %>% 
    filter(complete.cases(.)) %>% 
    mutate(across(tb.metad$sample_id, ~as.numeric(.))) %>% 
    merge(x = reshape2::melt(.), y = tb.metad, by.x = "variable", by.y = "sample_id") %>% 
    arrange(., taxa2, -value) -> tmp.tb.glom.top
  
  tmp.tb.top.id <- tmp.tb.glom.top %>% filter(taxa2 == tmp.tb.top.ind) %>% select(ID = variable, value) %>% pull("ID") %>% as.character()
  
  tmp.tb.glom.top <- tmp.tb.glom.top %>% mutate(ID = factor(variable, levels = tmp.tb.top.id))
  
  tmp.tb.top.taxa <- tmp.tb.top.taxa %>% filter(taxa %in% tmp.tb.glom.top$taxa2)
  
  tmp.colors <- tmp.tb.top.taxa$my.colors
  names(tmp.colors) <- tmp.tb.top.taxa$taxa
  
  if (i > 1) {
    
    tmp.ls.gg.leg.taxa <- list()
    
    for (p in levels(tmp.tb.top.taxa$Phylum2)) {
      
      tmp.tb.gg.leg.taxa <- tmp.tb.top.taxa %>% filter(Phylum2 == p)
      
      tmp.ls.gg.leg.taxa[[p]] <- ggpubr::as_ggplot(ggpubr::get_legend(
        
        ggplot(tmp.tb.gg.leg.taxa, aes(x = 1, y = taxa, fill = taxa)) + 
          geom_point(shape = 22, size = 3) + 
          scale_fill_manual(name = p, values = tmp.tb.gg.leg.taxa$my.colors) +  
          guides(shape = guide_legend(override.aes = list(size = 2)),
                 color = guide_legend(override.aes = list(size = 2)),
                 fill = guide_legend(ncol = 1, override.aes = list(size = 2))) +
          my.theme.1 + 
          theme(legend.position = "right",
                legend.title = element_text(size = 8, colour = 'black', face = "bold"),
                legend.text = element_text(size = 6, colour = 'black', margin = margin(0, 0, 0, 2)),
                legend.key.size = unit(6, 'pt'),
                legend.justification = c(0, 1))
      ))
      
    }; rm(p)
    
    gg.bp.leg <- (wrap_plots(tmp.ls.gg.leg.taxa) + plot_layout(nrow = 1))
    
  } else {
    
    gg.bp.leg <-
      ggplot(tmp.tb.top.taxa %>% group_by(taxa) %>% slice_head(n = 1), aes(x = 1, y = taxa, fill = taxa)) + 
      geom_point(shape = 22, size = 3) + 
      scale_fill_manual(name = "Phylum", values = tmp.tb.top.taxa %>% group_by(taxa) %>% slice_head(n = 1) %>% pull(my.colors)) +  
      my.theme.1 + 
      theme(legend.position = "right",
            legend.title = element_text(size = 8, colour = 'black', face = "bold"),
            legend.text = element_text(size = 7, colour = 'black', margin = margin(0, 0, 0, 2)),
            legend.key.size = unit(.5, 'lines'),
            legend.justification = c(0, .5))
    
    gg.bp.leg <- ggpubr::as_ggplot(ggpubr::get_legend(gg.bp.leg)); gg.bp.leg
    
  }
  
  gg.bp <- 
    
    tmp.tb.glom.top %>% mutate(Ethnicity = str_remove(Ethnicity, " \\(.*")) %>% 
    group_by(Ethnicity, taxa2) %>% 
    summarise(value = mean(value, na.rm = T), .groups = "drop") %>% 
    
    ggplot(., 
           aes(x = Ethnicity, stratum = taxa2, alluvium = taxa2, y = value, fill = taxa2)) + 
    geom_stratum(color = "#000000FF", alpha = 1, width = .5, linewidth = 0) +
    geom_flow(width = .5) +
    #ggplot(tmp.tb.glom.top %>% mutate(Ethnicity = str_remove(Ethnicity, " \\(.*")), aes(fill = taxa2, x = Ethnicity, y = value)) + 
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    my.theme.2 + 
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          panel.spacing.x = unit(.5, "lines")); gg.bp
  
  gg.bp.plots[[i]] <- gg.bp
  
}; rm(i, tmp.rank)

gg.bp.p2 <- gg.bp.plots[[1]]; gg.bp.p2
gg.bp.c2 <- gg.bp.plots[[2]]; gg.bp.c2
gg.bp.o2 <- gg.bp.plots[[3]]; gg.bp.o2
gg.bp.f2 <- gg.bp.plots[[4]]; gg.bp.f2
gg.bp.g2 <- gg.bp.plots[[5]]; gg.bp.g2

ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax_1p.png"), plot = gg.bp.p2 + gg.bp.p + plot_layout(widths = c(5, 95)), dpi = 600, width = 2500 * 4.5, height = 2500 * .85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax_2c.png"), plot = gg.bp.c2 + gg.bp.c + plot_layout(widths = c(5, 95)), dpi = 600, width = 2500 * 4.5, height = 2500 * .85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax_3o.png"), plot = gg.bp.o2 + gg.bp.o + plot_layout(widths = c(5, 95)), dpi = 600, width = 2500 * 4.5, height = 2500 * .85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax_4f.png"), plot = gg.bp.f2 + gg.bp.f + plot_layout(widths = c(5, 95)), dpi = 600, width = 2500 * 4.5, height = 2500 * .85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax_5g.png"), plot = gg.bp.g2 + gg.bp.g + plot_layout(widths = c(5, 95)), dpi = 600, width = 2500 * 4.5, height = 2500 * .85, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))


### Summaries ------------------------------------------------------------------


ls.glom.tb
ls.glom.tb.filt

openxlsx::write.xlsx(x = ls.glom.tb, file = paste0(path.tb.dir, "/profiles/taxa-glom-table.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = ls.glom.tb.filt, file = paste0(path.tb.dir, "/profiles/taxa-glom-filt-table.xlsx"), rowNames = F)

# Perc Class/Unclass (unfilt)

ls.glom.unclass <- list()

for (i in seq_along(ls.glom.tb)) {
  
  tmp.rank = names(ls.glom.tb[[i]])[[i]]
  
  if (i > 1) {
    
    ls.glom.unclass[[tmp.rank]] <- 
      ls.glom.tb[[i]] %>% 
      select(tmp_rank = i, any_of(tb.metad$sample_id)) %>% 
      mutate(tmp_rank = ifelse(str_detect(tmp_rank, paste0(str_to_lower(str_sub(names(ls.glom.tb[[i]])[[i]], 1, 1)), ": ")), "Classified", "Unclassified")) %>% 
      group_by(tmp_rank) %>% 
      summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
      mutate(across(where(is.numeric), ~./sum(.))) %>% 
      data.frame(check.names = F) %>% 
      rename(!!tmp.rank := tmp_rank)
    
  } else {
    
    ls.glom.unclass[[tmp.rank]] <- 
      ls.glom.tb[[i]] %>% 
      select(tmp_rank = i, any_of(tb.metad$sample_id)) %>% 
      mutate(tmp_rank = ifelse(!str_detect(tmp_rank, "d: Unassigned"), "Classified", "Unclassified")) %>% 
      group_by(tmp_rank) %>% 
      summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
      mutate(across(where(is.numeric), ~./sum(.))) %>% 
      data.frame(check.names = F) %>% 
      rename(!!tmp.rank := tmp_rank)
    
  }
  
}; rm(i)

ls.glom.unclass

# Perc Class/Unclass (groups) 

ls.glom.unclass.group <- list()

for (i in seq_along(ls.glom.unclass)) {
  
  tmp.rank = names(ls.glom.unclass[[i]])[[i]]
  
  ls.glom.unclass.group[[tmp.rank]] <- 
    ls.glom.unclass[[i]] %>% 
    pivot_longer(cols = tb.metad$sample_id, names_to = "sample_id", values_to = "perc") %>% 
    rename("status" = 1) %>% 
    mutate(rank = names(ls.glom.unclass)[[i]])
  
  
}; rm(i)

bind_rows(ls.glom.unclass.group) %>% 
  merge(tb.metad, ., by = "sample_id") %>% 
  mutate(rank = factor(rank, levels = unique(rank))) -> tb.glom.unclass.group

fun.get.summary(tb = tb.glom.unclass.group, gp_comp = "Group_tribe", gp_facet = c("status", "rank"), value_var = "perc") %>% 
  arrange(rank, Group_tribe, status) %>% 
  mutate(full_mean = str_glue("{round(mean, 4)*100} (±{round(sd, 4)*100})")) -> tb.glom.unclass.group

tb.glom.unclass.group

bind_rows(ls.glom.unclass.group) %>% 
  merge(tb.metad, ., by = "sample_id") %>% 
  mutate(rank = factor(rank, levels = unique(rank))) -> tb.glom.unclass.ethnicity

fun.get.summary(tb = tb.glom.unclass.ethnicity, gp_comp = "Ethnicity", gp_facet = c("status", "rank"), value_var = "perc") %>% 
  arrange(rank, Ethnicity, status) %>% 
  mutate(full_mean = str_glue("{round(mean, 4)*100} (±{round(sd, 4)*100})")) -> tb.glom.unclass.ethnicity

tb.glom.unclass.ethnicity

# Perc table for filt tables

ls.glom.tb.filt.perc <- list()

for (i in seq_along(ls.glom.tb.filt)) {
  
  tmp.rank = names(ls.glom.tb.filt[[i]])[[i+1]]
  
  ls.glom.tb.filt.perc[[tmp.rank]] <- 
    ls.glom.tb.filt[[i]] %>% 
    filter(str_detect(Phylum, "^p: ")) %>% 
    mutate(across(where(is.numeric), ~./sum(.))) %>% 
    data.frame(check.names = F) %>% 
    mutate(mean_of_cols = rowMeans(select(., where(is.numeric)))) %>% 
    arrange(-mean_of_cols) %>% 
    select(-mean_of_cols)
  
}

lapply(ls.glom.tb.filt.perc, head)

# Perc table for filt tables (groups)

ls.glom.tb.filt.perc.group <- list()

for (i in seq_along(ls.glom.tb.filt.perc)) {
  
  tmp.rank = names(ls.glom.tb.filt.perc)[[i]]
  
  ls.glom.tb.filt.perc.group[[tmp.rank]] <- 
    ls.glom.tb.filt.perc[[i]] %>% 
    rowwise() %>%
    mutate(mean_of_cols = mean(c_across(any_of(tb.metad$sample_id)), na.rm = TRUE)) %>%
    ungroup() %>% 
    arrange(-mean_of_cols) %>% 
    select(-mean_of_cols) %>% 
    mutate(across(where(is.character), .fns = ~factor(., unique(.)))) %>% 
    pivot_longer(cols = tb.metad$sample_id, names_to = "sample_id", values_to = "perc") %>% 
    merge(., tb.metad, by = "sample_id") %>% 
    group_by(Group_tribe, across(2:(2+i))) %>% 
    summarise(mean = mean(perc), sd = sd(perc), .groups = "drop") %>% 
    mutate(full_mean = str_glue("{round(mean, 4)*100} (±{round(sd, 4)*100})")) %>% 
    select(-c(mean,sd)) %>% 
    pivot_wider(names_from = Group_tribe, values_from = full_mean) %>% 
    select(1:(1+i), starts_with(levels(tb.metad$Group_tribe))) %>% 
    arrange(across(all_of(tmp.rank)))
  
}; rm(i)

ls.glom.tb.filt.perc.ethnicity <- list()

for (i in seq_along(ls.glom.tb.filt.perc)) {
  
  tmp.rank = names(ls.glom.tb.filt.perc)[[i]]
  
  ls.glom.tb.filt.perc.ethnicity[[tmp.rank]] <- 
    ls.glom.tb.filt.perc[[i]] %>% 
    rowwise() %>%
    mutate(mean_of_cols = mean(c_across(any_of(tb.metad$sample_id)), na.rm = TRUE)) %>%
    ungroup() %>% 
    arrange(-mean_of_cols) %>% 
    select(-mean_of_cols) %>% 
    mutate(across(where(is.character), .fns = ~factor(., unique(.)))) %>% 
    pivot_longer(cols = tb.metad$sample_id, names_to = "sample_id", values_to = "perc") %>% 
    merge(., tb.metad, by = "sample_id") %>% 
    group_by(Ethnicity, across(2:(2+i))) %>% 
    summarise(mean = mean(perc), sd = sd(perc), .groups = "drop") %>% 
    mutate(full_mean = str_glue("{round(mean, 4)*100} (±{round(sd, 4)*100})")) %>% 
    select(-c(mean,sd)) %>% 
    pivot_wider(names_from = Ethnicity, values_from = full_mean) %>% 
    select(1:(1+i), starts_with(levels(tb.metad$Ethnicity))) %>% 
    arrange(across(all_of(tmp.rank)))
  
}; rm(i)

# Saving Tables

ls.glom.unclass
tb.glom.unclass.group
tb.glom.unclass.ethnicity

ls.glom.tb.filt.perc
ls.glom.tb.filt.perc.group
ls.glom.tb.filt.perc.ethnicity

openxlsx::write.xlsx(x = ls.glom.unclass, file = paste0(path.tb.dir, "/profiles/taxa-unclass-sample.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = tb.glom.unclass.ethnicity, file = paste0(path.tb.dir, "/profiles/taxa-unclass-ethnicity.xlsx"), rowNames = F)

openxlsx::write.xlsx(x = ls.glom.tb.filt.perc, file = paste0(path.tb.dir, "/profiles/taxa-glom-perc-sample.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = ls.glom.tb.filt.perc.ethnicity, file = paste0(path.tb.dir, "/profiles/taxa-glom-perc-ethnicity.xlsx"), rowNames = F)

rm(ls.glom.unclass, ls.glom.unclass.group, tb.glom.unclass.group, tb.glom.unclass.ethnicity)
rm(ls.glom.tb.filt.perc, ls.glom.tb.filt.perc.group, ls.glom.tb.filt.perc.ethnicity)


## Alpha diversity -------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diversity"))) { dir.create(str_glue("{path.fg.dir}/diversity"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diversity"))) { dir.create(str_glue("{path.tb.dir}/diversity"), recursive = T) }

fun.get.alpha.div.ps <- function(my.ps.table.ls) {
  
  tmp.ls <- list()
  
  for (i in seq_along(my.ps.table.ls)) {
    
    tmp.ls[[names(my.ps.table.ls)[[i]]]] <- 
      microbiome::alpha(my.ps.table.ls[[i]], index = c("observed", "diversity_shannon", "diversity_gini_simpson")) %>% 
      select(observed, diversity_shannon, diversity_gini_simpson) %>% 
      rename("Richness" = 1, "Shannon index" = 2, "Gini-Simpson" = 3) %>% 
      rownames_to_column("sample_id") %>% 
      pivot_longer(cols = where(is.numeric), names_to = "var", values_to = "value") %>%
      mutate(var = factor(var, c("Richness", "Shannon index", "Gini-Simpson"))) %>%
      mutate(value = replace_na(value, 0)) %>%
      mutate(rank = names(my.ps.table.ls)[[i]]) %>%
      merge(tb.metad, .) %>%
      arrange(var) 
    
  }; rm(i)
  
  return(bind_rows(tmp.ls))
  
}

tb.adiv <- fun.get.alpha.div.ps(ls.glom.ps.filt)


### Ethnicity ------------------------------------------------------------------


ls.adiv.summ <- list()
ls.adiv.compm <- list()
ls.adiv.gg <- list()

for (i in ranks) {
  
  tmp.adiv <- tb.adiv %>% filter(rank == i) %>% mutate(Ethnicity2 = str_remove(Ethnicity, " \\(.*"))
  
  ls.adiv.summ[[i]] <- 
    fun.get.summary(tb = tmp.adiv,
                    gp_comp = "Ethnicity2", 
                    gp_facet = "var",
                    value_var = "value") %>% 
    mutate(rank = i); ls.adiv.summ[[i]]
  
  ls.adiv.compm[[i]] <- 
    fun.get.wilcox(tb = tmp.adiv, 
                   gp_comp = "Ethnicity2", 
                   gp_facet = "var", 
                   value_var = "value",
                   p_adjust_method = "fdr") %>% 
    mutate(rank = i); ls.adiv.compm[[i]]
  
  ls.adiv.gg[[i]] <- 
    fun.get.boxp.compm.pair(tb = tmp.adiv, 
                            tb.compm = ls.adiv.compm[[i]], 
                            gp_facet = "var", 
                            gp_comp = "Ethnicity2",
                            gp_fill = "Ethnicity",
                            gp_col = "Ethnicity",
                            value_var = "value",
                            only.sig = T,
                            alpha = 0.1, 
                            label.type = "star", 
                            color.set = colors.et) + 
    labs(title = paste0(i, "<br>")) + 
    theme(plot.title = element_markdown(size = 12, face = 'bold', colour = 'black'),
          axis.text.x = element_markdown(lineheight = 1.15))
  
  print(ls.adiv.gg[[i]] )
  
}; rm(i)


# Saving graphs
{
  png(filename = paste0(path.fg.dir, "/diversity/alpha-ethnicity.png"), width = 5000 * 1.35, height = 5000 * 4.55, units = "px", res = 1200, type = "cairo")
  
  print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))
  
  dev.off()
}


ls.adiv <- list(
  
  "Measures" = tb.adiv,
  
  "Summary" = ls.adiv.summ %>% bind_rows(),
  
  "Wilcoxon" = ls.adiv.compm %>% bind_rows()
  
)

openxlsx::write.xlsx(x = ls.adiv, file = paste0(path.tb.dir, "/diversity/alpha-ethnicity.xlsx"), rowNames = F)

rm(list = ls(pattern = "adiv"))


## Beta diversity ----------------------------------------------------------

my.donors <- c("BEL-051", "BEL-026", "BEL-016", "BEL-039", "BEL-005", "BEL-019", "BEL-010", 
               "CUR-009", "CUR-013", "BA-032", "BA-014", "BA-001", "BA-062", "BA-020")


# Ethnicity

fun.get.beta.div.ethnicity <- function(my.glom.ps, my.dist = "bray") {
  
  my.glom.tb <- merge(tax_table(my.glom.ps), otu_table(my.glom.ps), by = "row.names", all = T) %>% select(-c(Kingdom:Species))
  
  my.metad <- filter(tb.metad, sample_id %in% sample_names(my.glom.ps))
  
  pcoa.factors <- my.metad
  
  pcoa.values <- 
    my.glom.tb %>% 
    rename("taxa" = 1) %>% 
    column_to_rownames("taxa") %>% 
    select(my.metad$sample_id) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
    t()
  
  if (my.dist == "euclidean") {
    
    pcoa.dist <- phyloseq::distance(my.glom.ps, my.dist)
    
  } else {
    
    pcoa.dist <- phyloseq::distance(phyloseq_standardize_otu_abundance(physeq = my.glom.ps, method = 'total'), my.dist)
    
  }
  
  pcoa.cmds <- cmdscale(pcoa.dist, eig = TRUE, k = 3)
  
  pcoa.site.score <- 
    data.frame(pcoa.cmds$points) %>% 
    rename("PCoA1" = 1, "PCoA2" = 2, "PCoA3" = 3) %>% 
    rownames_to_column("sample_id") %>% 
    merge(., my.metad)
  
  pcoa.perma.comp <- adonis2(pcoa.dist ~ Ethnicity, data = my.metad, method = my.dist, permutations = 9999, seed = 123, by = "terms")
  
  pcoa.perma.df.comp <- 
    pcoa.perma.comp %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column() %>% 
    dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>% 
    mutate(Source = str_to_title(Source),
           sig = case_when(
             P <= 0.001 ~ "***",
             P <= 0.01 ~ "**",
             P <= 0.05 ~ "*",
             P <= 0.1 ~ "~",
             .default = "ns"))
  
  pcoa.perma.ph <- pairwiseAdonis::pairwise.adonis2(pcoa.dist ~ Ethnicity, data = my.metad, p.adjust.m = "fdr", perm = 9999)
  
  pcoa.perma.ph <- 
    bind_rows(pcoa.perma.ph) %>% 
    select(-parent_call) %>% 
    t() %>% 
    data.frame(check.names = T) %>% 
    rename(Model = 1, Residual = 2, Total = 3) %>% 
    rownames_to_column("comparison") %>% 
    separate_wider_delim(cols = comparison, delim = ".", names = c("comparison", "tmp")) %>% 
    pivot_longer(cols = -c(comparison, tmp), names_to = "component", values_to = "value") %>% 
    pivot_wider(names_from = "tmp", values_from = "value") %>% 
    separate_wider_delim(cols = comparison, delim = "_vs_", names = c("group1", "group2")); pcoa.perma.ph
  
  # Axes 1:2
  
  pcoa.site.score  <- 
    pcoa.site.score %>% 
    mutate(my.label = ifelse(sample_id %in% my.donors, sample_id, NA_character_)) %>% 
    mutate(PCoA1.d = ifelse(sample_id %in% my.donors, PCoA1, NA_integer_)) %>% 
    mutate(PCoA2.d = ifelse(sample_id %in% my.donors, PCoA2, NA_integer_))
  
  ggplot() +
    stat_ellipse(data = pcoa.site.score,
                 aes(x = PCoA1, y = PCoA2, color = Ethnicity, fill = Ethnicity),
                 geom = "polygon", linetype = "dashed", linewidth = 0.25, alpha = .05, show.legend = F, type = "norm") +
    geom_point(data = pcoa.site.score, 
               aes(x = PCoA1, y = PCoA2, fill = Ethnicity), 
               shape = 21, size = 2.25, alpha = .5, color = "white") +
    
    geom_text_repel(data = pcoa.site.score, 
                    aes(label = my.label, x = PCoA1, y = PCoA2), color = "black",
                    size = 2.25, segment.alpha = 1, segment.color = 'red', force = 10,
                    min.segment.length = 0.2, max.overlaps = 999, seed = 123, show.legend = F) +
    
    geom_point(data = pcoa.site.score, 
               aes(x = PCoA1.d, y = PCoA2.d, fill = Ethnicity), 
               shape = 21, size = 2.25, alpha = 1, color = "red") +
    
    annotate("richtext", x = Inf, y = -Inf, size = 2.5,
             label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp[1, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp[1, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp[1, "P"], 3)))),
             fill = NA, label.color = NA, hjust = 1, vjust = 0, color = "red") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    my.theme.2 +
    scale_color_manual(values = colors.et) +
    scale_fill_manual(values = colors.et) +
    labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
         y = paste0("PCoA2 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[2], "%)<br>")) + 
    theme(plot.subtitle = element_blank(),
          legend.position = "none") -> gg.pcoa.axes12; gg.pcoa.axes12
  
  gg.pcoa.panel <- gg.pcoa.axes12
  
  
  {
    png(filename = paste0(path.fg.dir, "/diversity/beta-div-ethnicity-panel-", str_to_lower(my.dist), ".png"), width = 5000 * 1, height = 5000 * .85, units = "px", res = 1200, type = "cairo")
    
    print(gg.pcoa.panel)
    
    dev.off()
  }
  
  pcoa.disp <- 
    betadisper(d = pcoa.dist, group = my.metad$Ethnicity, type = "centroid")$distances %>% 
    data.frame(check.names = F) %>% 
    dplyr::rename("value" = ".") %>% 
    rownames_to_column("sample_id") %>% 
    merge(my.metad, .) %>% 
    mutate(var = "Within-group<br>dispersion") %>% 
    mutate(Ethnicity2 = str_remove(Ethnicity, " \\(.*"))
  
  pcoa.disp.summary <- 
    fun.get.summary(tb = pcoa.disp, gp_comp = "Ethnicity2", gp_facet = "var", value_var = "value")
  
  pcoa.disp.posthoc <-
    fun.get.wilcox(tb = pcoa.disp, gp_comp = "Ethnicity2", gp_facet = "var", value_var = "value",
                   p_adjust_method = "fdr", conf_level = 0.95, alpha.global = 0.05)
  
  gg.pcoa.disp <- 
    fun.get.boxp.compm.pair(tb = pcoa.disp, 
                            tb.compm = pcoa.disp.posthoc,
                            gp_facet = "var", 
                            gp_comp = "Ethnicity2",
                            gp_fill = "Ethnicity",
                            gp_col = "Ethnicity",
                            value_var = "value",
                            only.sig = T,
                            alpha = 0.1, 
                            label.type = "star", color.set = colors.et) +
    theme(plot.title = element_markdown(size = 12, face = 'bold', colour = 'black'),
          axis.text.x = element_markdown(lineheight = 1.15)); gg.pcoa.disp
  
  {
    png(filename = paste0(path.fg.dir, "/diversity/beta-div-ethnicity-betadisp-", str_to_lower(my.dist), ".png"), width = 5000 * .5, height = 5000 * .7, units = "px", res = 1200, type = "cairo")
    
    print(gg.pcoa.disp)
    
    dev.off()
  }
  
  ls.pcoa <-
    
    list(
      
      "PERMANOVA" = pcoa.perma.df.comp,
      
      "PERMANOVA (post-hoc)" = pcoa.perma.ph,
      
      "Dispersion (Measures)" = pcoa.disp,
      "Dispersion (Summary)" = pcoa.disp.summary,
      "Dispersion (Wilcoxon)" = pcoa.disp.posthoc
      
    )
  
  openxlsx::write.xlsx(x = ls.pcoa, file = paste0(path.tb.dir, "/diversity/beta-div-ethnicity-", str_to_lower(my.dist), ".xlsx"), rowNames = F)
  
}

fun.get.beta.div.ethnicity(my.glom.ps = ps.table, my.dist = "bray")

rm(list = ls(pattern = "[.]?beta[.]?"))
rm(list = ls(pattern = "[.]?bdiv[.]?"))
rm(list = ls(pattern = "[.]?pcoa[.]?"))

# Confounders

fun.get.beta.div.confounders <- function(my.glom.ps, my.dist = "bray") {
  
  my.glom.tb <- merge(tax_table(my.glom.ps), otu_table(my.glom.ps), by = "row.names", all = T) %>% select(-c(Kingdom:Species))
  
  my.metad <- filter(tb.metad, sample_id %in% sample_names(my.glom.ps))
  
  pcoa.factors <- my.metad
  
  pcoa.values <- 
    my.glom.tb %>% 
    rename("taxa" = 1) %>% 
    column_to_rownames("taxa") %>% 
    select(my.metad$sample_id) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
    t()
  
  if (my.dist == "euclidean") {
    
    pcoa.dist <- phyloseq::distance(my.glom.ps, my.dist)
    
  } else {
    
    pcoa.dist <- phyloseq::distance(phyloseq_standardize_otu_abundance(physeq = my.glom.ps, method = 'total'), my.dist)
    
  }
  
  pcoa.cmds <- cmdscale(pcoa.dist, eig = TRUE, k = 3)
  
  pcoa.site.score <- 
    data.frame(pcoa.cmds$points) %>% 
    rename("PCoA1" = 1, "PCoA2" = 2, "PCoA3" = 3) %>% 
    rownames_to_column("sample_id") %>% 
    merge(., my.metad)
  
  # Sex
  
  pcoa.perma.comp.1 <- adonis2(pcoa.dist ~ Group_tribe * Sex, data = my.metad, method = my.dist, permutations = 9999, seed = 123, by = "terms")
  
  pcoa.perma.df.comp.1 <- 
    pcoa.perma.comp.1 %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column() %>% 
    dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>% 
    mutate(Source = str_to_title(Source),
           sig = case_when(
             P <= 0.001 ~ "***",
             P <= 0.01 ~ "**",
             P <= 0.05 ~ "*",
             P <= 0.1 ~ "~",
             .default = "ns"))
  
  # Axes 1:2
  
  ggplot() +
    stat_ellipse(data = pcoa.site.score,
                 aes(x = PCoA1, y = PCoA2, color = Sex, fill = Sex),
                 geom = "polygon", linetype = "dashed", linewidth = 0.25, alpha = .05, show.legend = F, type = "norm") +
    geom_point(pcoa.site.score, 
               mapping = aes(x = PCoA1, y = PCoA2, fill = Sex, shape = Ethnicity), 
               color = "white", size = 2.25, alpha = 1) +
    annotate("richtext", x = Inf, y = -Inf, size = 2.5,
             label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp.1[3, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp.1[3, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp.1[3, "P"], 3)))),
             fill = NA, label.color = NA, hjust = 1, vjust = 0, color = "red") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    my.theme.2 +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_color_manual(values = colors.sx) +
    scale_fill_manual(values = colors.sx) +
    labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
         y = paste0("PCoA2 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[2], "%)<br>")) + 
    theme(plot.subtitle = element_blank(),
          legend.position = "none") -> gg.pcoa.axes12
  
  # Age
  
  pcoa.perma.comp.2 <- adonis2(pcoa.dist ~ Group_tribe * Age_group, data = my.metad, method = my.dist, permutations = 9999, seed = 123, by = "terms")
  
  pcoa.perma.df.comp.2 <- 
    pcoa.perma.comp.2 %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column() %>% 
    dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>% 
    mutate(Source = str_to_title(Source),
           sig = case_when(
             P <= 0.001 ~ "***",
             P <= 0.01 ~ "**",
             P <= 0.05 ~ "*",
             P <= 0.1 ~ "~",
             .default = "ns"))
  
  
  # Axes 1:3
  
  ggplot() +
    stat_ellipse(data = pcoa.site.score,
                 aes(x = PCoA1, y = PCoA2, color = Age_group, fill = Age_group),
                 geom = "polygon", linetype = "dashed", linewidth = 0.25, alpha = .05, show.legend = F, type = "norm") +
    geom_point(pcoa.site.score, 
               mapping = aes(x = PCoA1, y = PCoA2, fill = Age_group, shape = Ethnicity), 
               color = "white", size = 2.25, alpha = 1) +
    annotate("richtext", x = Inf, y = -Inf, size = 2.5,
             label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp.2[3, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp.2[3, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp.2[3, "P"], 3)))),
             fill = NA, label.color = NA, hjust = 1, vjust = 0, color = "red") +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.05))) +
    my.theme.2 +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_color_manual(values = colors.ag) +
    scale_fill_manual(values = colors.ag) +
    labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
         y = paste0("PCoA2 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[2], "%)<br>")) + 
    theme(plot.subtitle = element_blank(),
          legend.position = "none") -> gg.pcoa.axes13
  
  gg.pcoa.panel <- 
    gg.pcoa.axes12 + gg.pcoa.axes13; gg.pcoa.panel
  
  
  {
    png(filename = paste0(path.fg.dir, "/diversity/beta-div-confounders-panel-", str_to_lower(my.dist), ".png"), width = 5000 * 1.65, height = 5000 * .85, units = "px", res = 1200, type = "cairo")
    
    print(gg.pcoa.panel)
    
    dev.off()
  }
  
  ls.pcoa <-
    
    list(
      
      "PERMANOVA (Sex)" = pcoa.perma.df.comp.1,
      "PERMANOVA (Age)" = pcoa.perma.df.comp.2 
    )
  
  openxlsx::write.xlsx(x = ls.pcoa, file = paste0(path.tb.dir, "/diversity/beta-div-confounders-", str_to_lower(my.dist), ".xlsx"), rowNames = F)
  
}

fun.get.beta.div.confounders(my.glom.ps = ps.table, my.dist = "bray")

rm(list = ls(pattern = "[.]?beta[.]?"))
rm(list = ls(pattern = "[.]?bdiv[.]?"))
rm(list = ls(pattern = "[.]?pcoa[.]?"))


## MaAsLin2 --------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diff_abund"))) { dir.create(str_glue("{path.fg.dir}/diff_abund"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diff_abund"))) { dir.create(str_glue("{path.tb.dir}/diff_abund"), recursive = T) }

ls.msl <- list()

# Ethnicity

comps <- c("Belgian (n=39)|Yanomami (n=78)")

for (j in seq_along(comps)) {
  
  for (i in seq_along(ranks)) {
    
    tmp.comp.factor <- "Ethnicity" 
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
               reference = paste0(tmp.comp.factor, ",", str_remove({tmp.comp}, '\\|.*')),
               min_abundance = 0,
               min_prevalence = 0.1,
               plot_heatmap = F, 
               plot_scatter = F); head(tmp.msl.comp)
    
    ls.msl[[ranks[[i]]]][[comps[[j]]]] <- 
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
  
}; rm(i)

tb.tax.ann <- 
  bind_rows(ls.tax.ann) %>% 
  mutate(Phylum2 = case_when(
    ! Phylum %in% names(colors.tax2) & str_detect(Phylum, "Unassigned") ~ "Unassigned",
    ! Phylum %in% names(colors.tax2) & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
    ! Phylum %in% names(colors.tax2) & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
    Phylum %in% names(colors.tax2) ~ Phylum)) %>%
  filter(Phylum2 != "Unassigned") %>% 
  mutate(Phylum2 = factor(Phylum2, names(colors.tax2)))

df.msl.raw <- 
  lapply(ls.msl, bind_rows) %>% bind_rows(.) %>% 
  filter(rank %in% c("Family", "Genus")) %>% 
  mutate(comp = factor(comp, levels = unique(comp))) %>% 
  mutate(enrich_group = factor(enrich_group, levels = unique(c(levels(tb.metad$Ethnicity), levels(tb.metad$Group_tribe), "Others")))) %>% 
  mutate(rank = factor(rank, levels = ranks)) %>% 
  merge(., tb.tax.ann, "feature", all.x = T) %>% 
  select(rank, Kingdom:Phylum2, feature, comp, enrich_group, coef:N.not.zero) %>% 
  arrange(rank, comp, enrich_group, -abs(coef))

df.msl.filt <-
  df.msl.raw %>% 
  filter(qval < 0.05) %>% 
  filter(enrich_group != "Others") %>% 
  filter(!str_detect(feature, " Un.\\)$"))

df.msl <- list("Raw" = df.msl.raw, "Filtered" = df.msl.filt)

openxlsx::write.xlsx(x = df.msl, file = paste0(path.tb.dir, "/diff_abund/maaslin.xlsx"), rowNames = F)


### Graph: Belgium vs. Yanomami: Genus -----------------------------------------

tmp.rank <- "Genus"

tmp.msl.features <-
  df.msl.filt %>% 
  filter(rank == tmp.rank) %>% 
  pull(feature) %>% 
  unique()

tmp.msl.features <- 
  df.msl.filt %>% 
  filter(rank == tmp.rank) %>% 
  group_by(across(rank:feature)) %>% 
  summarize(coef = mean(coef), .groups = "drop") %>% 
  arrange(-abs(coef)) %>% 
  slice_head(n = 60) %>% 
  arrange(-coef) %>% 
  pull(feature)

tmp.colors <- colors.et

tmp.tb.msl.gg <- 
  df.msl.filt %>% 
  filter(rank == tmp.rank) %>% 
  filter(feature %in% tmp.msl.features) %>% 
  mutate(feature = str_remove(feature, "^.*: ")) %>% 
  mutate(feature = factor(feature, levels = str_remove(tmp.msl.features, "^.*: "))) %>% 
  arrange(feature)

tmp.hline <- 
  tmp.tb.msl.gg %>% 
  mutate(feature.n = as.numeric(feature)) %>% 
  select(enrich_group, feature, feature.n) %>% 
  filter(feature %in% c("Campylobacter D", "Akkermansia")) %>% 
  pull(feature.n) %>% 
  mean()


# annotation bar
tmp.gg.msl.1 <- 
  ggplot(tmp.tb.msl.gg, aes(x = 1, y = feature, fill = Phylum2)) +
  geom_tile(aes(width = 0.9, height = 0.9), color = "white", size = .35) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(name = "", values = colors.tax2) +  
  my.theme.2 + 
  theme(panel.border = element_rect(colour = "#00000000", fill = NA),
        
        axis.line.y = element_line(color = "#00000000"),
        axis.line.x = element_line(color = "#00000000"),
        
        axis.text.y = element_markdown(size = 8, face = "italic"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 0.1, 0, 0), 'cm')); tmp.gg.msl.1

# heatmap
tmp.ps.tss <- 
  ls.glom.ps.filt[[tmp.rank]] %>% 
  phyloseq_standardize_otu_abundance(physeq = ., method = 'total')

tmp.dt.mean <-
  merge(tax_table(tmp.ps.tss), otu_table(tmp.ps.tss), by = "row.names", all = T) %>% 
  select(feature = any_of(tmp.rank), any_of(tb.metad$sample_id)) %>%
  column_to_rownames("feature") %>% 
  t() %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(cols = -sample_id, names_to = "feature", values_to = "perc") %>% 
  merge(tb.metad, .) %>% 
  group_by(sample_id, feature) %>% 
  summarise(perc = mean(perc, na.rm = T), .groups = "drop")

tmp.mt.mean <- 
  tmp.dt.mean %>% 
  pivot_wider(names_from = "feature", values_from = "perc") %>% 
  column_to_rownames("sample_id") %>% 
  as.matrix()

tmp.tb.clr <-
  microbiome::transform(tmp.mt.mean, transform = "clr") %>% 
  scale() %>% t() %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("feature") %>% 
  pivot_longer(cols = any_of(tb.metad$sample_id), names_to = "sample_id", values_to = "groupwise_zscore") %>% 
  data.frame(check.names = F) %>% 
  filter(feature %in% tmp.msl.features) %>% 
  mutate(feature = factor(str_remove(feature, "^.*: "), levels = levels(tmp.tb.msl.gg$feature))) %>% 
  mutate(sample_id = factor(sample_id, levels = tb.metad$sample_id))


tmp.gg.msl.2 <- 
  ggplot(data = tmp.tb.clr, mapping = aes(x = sample_id, y = feature, fill = groupwise_zscore)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("#4A7BB7FF", "#6EA6CDFF", "#FFFFFFFF", "#F67E4BFF", "#DD3D2DFF"), 
                       breaks = seq(-max(ceiling(tmp.tb.clr$groupwise_zscore)), max(ceiling(tmp.tb.clr$groupwise_zscore)), 2), 
                       limits = c(-max(abs(ceiling(tmp.tb.clr$groupwise_zscore))), 
                                  max(abs(ceiling(tmp.tb.clr$groupwise_zscore)))),
                       name = "CLR-Z") + 
  guides(fill = guide_colourbar(barwidth = 3.5, 
                                barheight = 0.25,
                                label.position = "bottom", 
                                title.position = "top", 
                                title.hjust = .5, 
                                draw.ulim = T, 
                                draw.llim = T)) +
  geom_hline(yintercept = tmp.hline, linewidth = 0.25, linetype = "dotted") +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  my.theme.1 + 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(0, 2, 0, 0), 'cm')); tmp.gg.msl.2

tmp.tb.sample_id <-
  tmp.tb.clr %>% 
  arrange(sample_id) %>% 
  select(sample_id) %>% 
  distinct() %>% 
  merge(., tb.metad) %>% 
  select(sample_id, Ethnicity) %>% 
  mutate(Ethnicity2 = str_remove(Ethnicity, " \\(.*")) %>% 
  mutate(my.y = "Population")

tmp.gg.msl.3 <-
  ggplot(tmp.tb.sample_id, aes(y = my.y, x = sample_id, fill = Ethnicity)) +
  geom_tile() +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(name = "", values = colors.et) +  
  labs(title = "Top 60 most variable Genera<br>Standardized CLR-transformed data (CLR-Z)") +
  my.theme.2 + 
  theme(panel.border = element_rect(colour = "#00000000", fill = NA),
        
        axis.line.y = element_line(color = "#00000000"),
        axis.line.x = element_line(color = "#00000000"),
        
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        
        plot.title = element_markdown(size = 10, face = 'plain', colour = 'black', hjust = .5),
        plot.margin = unit(c(0, 2, .1, 0), 'cm')); tmp.gg.msl.3

# main plot
tmp.gg.msl.4 <-
  ggplot(tmp.tb.msl.gg, aes(x = coef, y = feature, color = enrich_group, fill = enrich_group)) + 
  ggstats::geom_stripped_rows(color = NA) +
  geom_tile(aes(width = 0, height = 0.9), color = "#00000000", fill = "#00000000", size = .35) +
  geom_errorbar(aes(xmin = coef - stderr, xmax = coef + stderr), 
                color = "black", width = .25, linewidth = .25, position = position_dodge(.9)) + 
  #geom_point(size = 1.25, alpha = 1, shape = 16) + 
  geom_bar(position = position_dodge(), stat = "identity", width = .75, color = "black", linewidth = .15) +
  geom_hline(yintercept = tmp.hline, linewidth = 0.25, linetype = "dotted") +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "solid") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)), 
                     limits = ~ c(-1, 1) * max(abs(tmp.tb.msl.gg$coef)+abs(tmp.tb.msl.gg$stderr), na.rm = T),
                     breaks = c(-max(ceiling(tmp.tb.clr$groupwise_zscore)), 0, max(ceiling(tmp.tb.clr$groupwise_zscore)))) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_fill_manual(values = c(tmp.colors)) +  
  scale_color_manual(values = c(tmp.colors)) +  
  xlab("<br>MaAsLin (β)<br>coefficient") + 
  my.theme.2 +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, .5, 0, 0), 'cm')); tmp.gg.msl.4

plot_spacer() + tmp.gg.msl.3 + plot_spacer() + 
  tmp.gg.msl.1 + tmp.gg.msl.2 + tmp.gg.msl.4 + 
  plot_layout(nrow = 2, heights = c(.01, .99), widths = c(.025, .825, .15))

{
  png(filename = paste0(path.fg.dir, "/diff_abund/maaslin-panel-a.png"), width = 5000 * 1.15, height = 5000 * 2.25, units = "px", res = 1200, type = "cairo")
  
  print(
    (plot_spacer() + theme_void()) + tmp.gg.msl.3 + (plot_spacer() + theme_void()) + 
      tmp.gg.msl.1 + tmp.gg.msl.2 + tmp.gg.msl.4 + 
      plot_layout(nrow = 2, heights = c(.01, .99), widths = c(.025, .8, .175)))
  
  dev.off()
  }


rm(list = ls(pattern = "[.]?tmp[.]?"))


### Graph: Belgium vs. Yanomami: Genus -----------------------------------------


tmp.rank <- "Family"

df.msl.raw %>% 
  filter(rank == tmp.rank) %>% 
  mutate(sig = ifelse(qval < 0.05, as.character(enrich_group), "Not Sig.")) %>% 
  mutate(sig = factor(sig, c(names(colors.et), "Not Sig."))) %>% 
  droplevels() %>% 
  arrange(sig, -abs(coef), -qval) %>% 
  group_by(sig) %>% 
  mutate(top = 1:n()) %>% 
  ungroup() %>% 
  mutate(my.label = ifelse(top <= 10 & sig != "Not Sig.", str_remove(feature, "^.*: "), NA_character_)) %>% 
  mutate(my.nudge.x = case_when( !is.na(my.label) & coef > 0 ~ Inf,
                                 !is.na(my.label) & coef < 0 ~ -Inf,
                                 .default = NA_integer_ )) %>% 
  mutate(my.hjust = case_when( !is.na(my.label) & coef > 0 ~ "left",
                               !is.na(my.label) & coef < 0 ~ "right",
                               .default = NA_character_ )) -> tmp.tb.msl.gg

tmp.gg.vp <-
  tmp.tb.msl.gg %>%
  filter(!is.na(coef)) %>% 
  ggplot(., aes(x = coef, y = -log10(qval), fill = sig)) + 
  #geom_text(aes(x = Inf, y = -log10(0.05), label = "qval=0.05"), size = 1.75, color = "red", hjust = 1, vjust = 3, check_overlap = T) + 
  geom_text_repel(aes(label = my.label, nudge_x = my.nudge.x), color = "black",
                  direction = "y", size = 2, segment.alpha = .5, segment.color = 'black', segment.size = 0.25,
                  min.segment.length = 0.2, max.overlaps = 999, force = 20, force_pull = 0.5, seed = 123) +
  geom_point(shape = 21, alpha = 1, size = 2, color = "white") + 
  geom_hline(yintercept = -log10(0.05), color = c("red"), linetype = "dotted", linewidth = .25) + 
  ylab("-log<sub>10</sub> (<i>q</i>-value)<br>") +
  xlab("<br>MaAsLin (β) coefficient") +
  scale_x_continuous(limits = ~ c(-1.5, 1.5) * max(abs(.x))) + 
  scale_fill_manual(name = "Group", values = c(colors.et, "Not Sig." = "gray75")) +
  coord_cartesian(clip = "off") +
  my.theme.2 +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()); tmp.gg.vp


{
  png(filename = paste0(path.fg.dir, "/diff_abund/maaslin-panel-b.png"), width = 5000 * 1, height = 5000 * .75, units = "px", res = 1200, type = "cairo")
  
  print(tmp.gg.vp)
  
  dev.off()
}

rm(list = ls(pattern = "[.]?tmp[.]?"))


save.image(paste0(path.wk.dir, "taxonomy-ecology_260826.RData"))


# End --------------------------------------------------------------------------

