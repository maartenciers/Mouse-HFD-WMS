

# Contribution Plots ------------------------------------------------------


# Setup ------------------------------------------------------------------------


## Environment and Packages ----------------------------------------------------


rm(list = ls())
graphics.off()
gc()

set.seed("123")

library("ggrepel")
library("ggtext")
library("ggpubr")
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

path.wk.dir <- paste0(getwd(), "/task260708_pw-contrib")
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
          )
        })
      })
    ) %>%
    unnest(results) %>%
    mutate(
      adj.p.value = p.adjust(p.value, method = p_adjust_method),
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
                             step.mult = 0.08) {
  
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
        label.type == "p.adj" & adj.p.value < 0.001 ~ "𝘱-adj < 0.001",
        label.type == "p.adj" ~ paste0("𝘱-adj = ", sprintf("%1.3f", round(adj.p.value, 3)))
      )
    ) %>%
    ungroup()
  
  return(tb.brackets)
}

fun.get.brackets.2 <- function(tb, tb.compm, gp_comp, gp_facet, value_var,
                               only.sig = T, alpha = 0.1, label.type = c("stars", "p.adj"), 
                               step.mult = 0.08, sqrt_scale = FALSE) {
  
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
        label.type == "p.adj" & adj.p.value < 0.001 ~ "𝘱-adj < 0.001",
        label.type == "p.adj" ~ paste0("𝘱-adj = ", sprintf("%1.3f", round(adj.p.value, 3)))
      )
    ) %>%
    ungroup()
  
  return(tb.brackets)
}

fun.upper.whisker <- function(x) {
  boxplot.stats(x)$stats[5]
}

fun.get.boxp.compm1 <- function(tb, tb.compm, gp_comp, gp_fill = NULL, gp_col = NULL,
                                gp_facet, value_var, only.sig = F, alpha = 0.1, 
                                label.type = "stars") {
  
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
      geom_boxplot(width = .65, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
      ggpubr::stat_pvalue_manual(
        tb.brackets,
        vjust = -0.5,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "label",
        tip.length = 0.01,
        bracket.size = 0.25,
        size = 2.75,
        hide.ns = only.sig,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
      scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      coord_cartesian(clip = "off") +
      ylab("") + xlab("") +
      my.theme.1 + 
      theme(legend.position = "none",
            panel.spacing = unit(1, "cm"), 
            axis.text.x = element_text(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
    
    return(gg)
    
  }
  
  gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                       fill = .data[[gp_fill]], color = .data[[gp_col]])) +
    facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
    geom_boxplot(width = .65, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
    scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.1 + 
    theme(legend.position = "none",
          panel.spacing = unit(1, "cm"), 
          axis.text.x = element_text(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
  
  return(gg)
  
}


fun.get.boxp.compm3 <- function(tb, tb.compm, tb.kruskal, gp_comp,  gp_fill = NULL, gp_col = NULL,
                                gp_facet, value_var,
                                only.sig = F, alpha = 0.1, label.type = "stars") {
  
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
    label.type = label.type,
    step.mult = 0.15
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
      geom_boxplot(width = .65, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
      geom_text(
        data = tb.global,
        aes(x = 0.5, y = Inf, label = p.full),
        inherit.aes = FALSE,
        vjust = 2,
        hjust = 0,
        size = 2.75,
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
        size = 2.75,
        hide.ns = only.sig,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
      scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      coord_cartesian(clip = "off") +
      ylab("") + xlab("") +
      my.theme.1 + 
      theme(legend.position = "none",
            panel.spacing = unit(1, "cm"), 
            axis.text.x = element_text(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
    
    return(gg)
    
  }
  
  
  gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                       fill = .data[[gp_fill]], color = .data[[gp_col]])) +
    facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
    geom_boxplot(width = .65, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
    geom_text(
      data = tb.global,
      aes(x = 0.5, y = Inf, label = p.full),
      inherit.aes = FALSE,
      vjust = 2,
      hjust = 0,
      size = 2.75,
      color = "black"
    ) +
    
    scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.1 + 
    theme(legend.position = "none",
          panel.spacing = unit(1, "cm"), 
          axis.text.x = element_text(size = 8, colour = 'black', angle = 45, hjust = 1, vjust = 1))
  
  return(gg)
  
}


## Data loading ----------------------------------------------------------------


### HFD data -------------------------------------------------------------------


# Metadata
tb.metad <- 
  read.delim(paste0(getwd(), "/misc/metadata_full.tsv")) %>% 
  filter(experiment == my.exp) %>% 
  mutate(timepoint = case_when(
    experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
    experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD",
  )) %>% 
  mutate(donor_ethnicity = ifelse(donor_ethnicity == "Western", "Belgian", donor_ethnicity)) %>% 
  select(sample_id, assembly_id, donor_id:cage_id, mouse_n, timepoint) %>% 
  mutate(donor_ethnicity = factor(donor_ethnicity, c("Belgian", "Yanomami"))) %>% 
  mutate(timepoint = factor(timepoint, c("pre-HFD", "post-HFD"))) %>% 
  arrange(donor_ethnicity, timepoint, mouse_n) %>% 
  mutate(donor_ethnicity_timepoint = str_glue("{donor_ethnicity}: {timepoint}")) %>% 
  mutate(donor_ethnicity_timepoint = factor(donor_ethnicity_timepoint, unique(.$donor_ethnicity_timepoint))) %>% 
  mutate(donor_mouse_n = str_remove_all(sample_id, "_T[12]$")); head(tb.metad)

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


### Donor data -----------------------------------------------------------------


## Graph settings --------------------------------------------------------------


# Themes
my.theme.1 <- 
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_text(size = 10, colour = 'black'),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.clip = "off",
        
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
        strip.clip = "off",
        
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

colors.gp <- c("#FFA500FF", "#0072B2FF")
names(colors.gp) <- levels(tb.metad$donor_ethnicity)

colors.gp.tp <- c("#FFA500FF", "#FF5500FF", "#0072B2FF", "#2F00B2FF")
names(colors.gp.tp) <- levels(tb.metad$donor_ethnicity_timepoint)

colors.tax <- c("#B06888FF", "#5080C8FF", "#E8B800FF", "#18C898FF", "#F88048FF", "#C8C8C8FF", "#303028FF")
names(colors.tax) <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota", "Other: Bacteria", "Other: Domains")

# Legends

gg.leg.gp <-
  ggplot(tb.metad, aes(x = 1, y = donor_ethnicity, fill = donor_ethnicity)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.gp) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.gp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.gp)); gg.leg.gp

gg.leg.gp.tp <-
  ggplot(tb.metad, aes(x = 1, y = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.gp.tp) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.gp.tp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.gp.tp)); gg.leg.gp.tp

gg.leg.tax <-
  ggplot(names(colors.tax) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.tax) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.tax <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.tax)); gg.leg.tax

gg.leg.all <- gg.leg.gp + gg.leg.gp.tp + gg.leg.tax + plot_layout(ncol = 1, heights = c(3, 5, 8))
gg.leg.all

{
  png(filename = paste0(path.fg.dir, "/legends.png"), width = 5000 * .65, height = 5000 * 1, units = "px", res = 1200, type = "cairo")
  
  print(gg.leg.all)
  
  dev.off()
}


# Analysis ---------------------------------------------------------------------


## New taxa profile ------------------------------------------------------------


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


### Samples --------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/profiles"))) { dir.create(str_glue("{path.fg.dir}/profiles"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/profiles"))) { dir.create(str_glue("{path.tb.dir}/profiles"), recursive = T) }

gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.tb.glom <- ls.glom.tb.filt[[i]]
  tmp.tb.glom <- 
    tmp.tb.glom %>% 
    filter(str_detect(Phylum, "^p: ")) %>% 
    mutate(Phylum = str_remove(Phylum, "p: ")) %>% 
    mutate(taxa = .[[tmp.rank]]) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
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
  
  tmp.tb.top.ind <- tmp.tb.glom %>% select(Phylum2, taxa, mean_of_cols) %>% arrange(-mean_of_cols) %>% pull(taxa) %>% .[1]
  
  tmp.tb.top.taxa <- 
    bind_rows(
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[1]]) %>% head(15),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[2]]) %>% head(5),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[3]]) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[4]]) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[5]]) %>% head(3)
      
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
    ) %>% arrange(taxa)
  
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
    ggplot(tmp.tb.glom.top, aes(fill = taxa2, x = mouse_n, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    facet_grid(. ~ donor_ethnicity + timepoint, space = "free_x", scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    guides(fill=guide_legend(ncol = ifelse(length(unique(tmp.tb.glom.top$taxa)) > 10, 1, 1))) +
    my.theme.2 + 
    theme(panel.spacing = unit(.1, "lines"), 
          strip.background = element_rect(fill = "grey85", colour = NA),
          strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
          axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)); gg.bp
  
  gg.bp.plots[[i]] <- gg.bp + gg.bp.leg + plot_layout(ncol = 2, widths = c(.3, .7))
  
}; rm(i)

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p
gg.bp.c <- gg.bp.plots[[2]]; gg.bp.c
gg.bp.o <- gg.bp.plots[[3]]; gg.bp.o
gg.bp.f <- gg.bp.plots[[4]]; gg.bp.f
gg.bp.g <- gg.bp.plots[[5]]; gg.bp.g

ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_1p.png"), plot = gg.bp.p, dpi = 1200, width = 5000 * 3.5, height = 5000 * .9, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_2c.png"), plot = gg.bp.c, dpi = 1200, width = 5000 * 3.5, height = 5000 * .9, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_3o.png"), plot = gg.bp.o, dpi = 1200, width = 5000 * 3.5, height = 5000 * .9, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_4f.png"), plot = gg.bp.f, dpi = 1200, width = 5000 * 3.5, height = 5000 * .9, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_5g.png"), plot = gg.bp.g, dpi = 1200, width = 5000 * 3.5, height = 5000 * .9, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))


### Group ----------------------------------------------------------------------


gg.bp.plots <- list()
ls.top.taxa <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.tb.glom <- ls.glom.tb.filt[[i]]
  tmp.tb.glom <- 
    tmp.tb.glom %>% 
    filter(str_detect(Phylum, "^p: ")) %>% 
    mutate(Phylum = str_remove(Phylum, "p: ")) %>% 
    mutate(taxa = .[[tmp.rank]]) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
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
  
  tmp.tb.top.ind <- tmp.tb.glom %>% select(Phylum2, taxa, mean_of_cols) %>% arrange(-mean_of_cols) %>% pull(taxa) %>% .[1]
  
  tmp.tb.top.taxa <- 
    bind_rows(
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[1]]) %>% head(15),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[2]]) %>% head(5),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[3]]) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[4]]) %>% head(3),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[5]]) %>% head(3)
      
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
    ) %>% arrange(taxa)
  
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
  
  tmp.tb.top.taxa <- tmp.tb.top.taxa %>% filter(taxa %in% tmp.tb.glom.top$taxa2)
  
  tmp.colors <- tmp.tb.top.taxa$my.colors
  names(tmp.colors) <- tmp.tb.top.taxa$taxa
  
  ls.top.taxa[[tmp.rank]] <- tmp.tb.top.taxa
  
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
    ggplot(tmp.tb.glom.top %>% 
             mutate(group = paste0(str_sub(donor_ethnicity, 1, 1), "_", str_remove_all(variable, ".*_"))) %>% 
             mutate(d_ethnicity = paste0(str_sub(donor_ethnicity, 1, 1))) %>% 
             mutate(t_point = paste0(str_remove_all(variable, ".*_"))), 
           aes(fill = taxa2, x = group, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    facet_grid(. ~ d_ethnicity + t_point, scales = "free", space = "free") +
    scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_fill_manual(values = tmp.colors) +  
    labs(fill = tmp.rank) + xlab("") + ylab("Relative abundance (%)<br>") +
    guides(fill=guide_legend(ncol = ifelse(length(unique(tmp.tb.glom.top$taxa2)) > 10, 1, 1))) +
    my.theme.2 + 
    theme(panel.spacing = unit(.1, "lines"), 
          strip.background = element_rect(fill = "grey85", colour = NA),
          strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  gg.bp.plots[[i]] <- gg.bp + gg.bp.leg + plot_layout(ncol = 2, widths = c(.1, .9))
  
}; rm(i)

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p
gg.bp.c <- gg.bp.plots[[2]]; gg.bp.c
gg.bp.o <- gg.bp.plots[[3]]; gg.bp.o
gg.bp.f <- gg.bp.plots[[4]]; gg.bp.f
gg.bp.g <- gg.bp.plots[[5]]; gg.bp.g

ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_1p.png"), plot = gg.bp.p, dpi = 1200, width = 5000 * 3, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_2c.png"), plot = gg.bp.c, dpi = 1200, width = 5000 * 3, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_3o.png"), plot = gg.bp.o, dpi = 1200, width = 5000 * 3, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_4f.png"), plot = gg.bp.f, dpi = 1200, width = 5000 * 3, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_5g.png"), plot = gg.bp.g, dpi = 1200, width = 5000 * 3, height = 5000 * 0.85, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))

ls.top.taxa


## Taxa contribution -----------------------------------------------------------


### HFD data -------------------------------------------------------------------


# Annotation + Counts

ls.counts <- list()

for (i in unique(tb.metad$assembly_id)) {
  
  ls.counts[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_counts.xlsx"), sheet = 1)
  
}; rm(i)

ls.taxa <- list()

for (i in unique(tb.metad$assembly_id)) {
  
  ls.taxa[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_taxonomy.xlsx"), sheet = 1)
  
}; rm(i)

ls.ann <- list()

for (i in unique(tb.metad$assembly_id)) {
  
  ls.ann[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_annotation.xlsx"), sheet = 2)
  
}; rm(i)

tb.tax <-
  lapply(ls.taxa, . %>% select(tax_id:Species)) %>%
  bind_rows() %>%
  distinct()

tb.ann <-
  lapply(ls.ann, . %>% select(ortholog_id = emapper_seed_ortholog, ko_id = emapper_KEGG_ko)) %>%
  bind_rows() %>%
  distinct() %>%
  mutate(ko_id = str_remove_all(ko_id, "ko:")) %>%
  filter(!is.na(ko_id))

ls.master <- list()

for (i in 1:length(ls.counts)) {
  
  ls.master[[i]] <-
    reduce(
      list(
        ls.ann[[i]] %>% select(contig_id, gene_id, ortholog_id = emapper_seed_ortholog),
        ls.taxa[[i]] %>% select(contig_id, tax_id),
        ls.counts[[i]]
      ), full_join, by = c("contig_id")) %>%
    filter(ortholog_id %in% tb.ann$ortholog_id) %>%
    distinct() %>%
    select(-c(contig_id, gene_id)) %>%
    group_by(ortholog_id, tax_id) %>%
    summarise(across(where(is.numeric), ~sum(.x, na.rm = T)), .groups = "drop") %>%
    mutate("prev" = rowSums(.[tb.metad %>% filter(assembly_id == unique(tb.metad$assembly_id)[[i]]) %>% pull(sample_id) %>% sort()] >= 1)) %>%
    filter(prev >= (tb.metad %>% filter(assembly_id == unique(tb.metad$assembly_id)[[i]]) %>% pull(sample_id) %>% length())*.1) %>%
    select(-prev)
  
}; rm(i)

tb.work <-
  reduce(ls.master, full_join, by = c("ortholog_id", "tax_id")) %>%
  mutate(across(any_of(tb.metad$sample_id), ~replace_na(., 0))) %>%
  select(ends_with("_id"), any_of(tb.metad$sample_id))

rm(ls.counts, ls.ann, ls.taxa, ls.master)
gc()


### Info data ------------------------------------------------------------------


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


### Loading LEfSe DAA tables ---------------------------------------------------

tmp.path <- path.tb.dir %>% str_replace(., "task260708_pw-contrib", "task260619_functional-profiler")

df.lefse.raw <-
  openxlsx::read.xlsx(xlsxFile = paste0(tmp.path, "/diff_abund/lefse.xlsx"), rowNames = F)

df.lefse.filt <- df.lefse.raw %>% filter(ef_lda > 2, padj < 0.05)

df.ancom.raw <-
  openxlsx::read.xlsx(xlsxFile = paste0(tmp.path, "/diff_abund/ancom.xlsx"), rowNames = F)

df.ancom.filt <- 
  df.ancom.raw %>% filter(padj < 0.05) %>% 
  filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_class == "Metabolism")

tmp.to.filter <-
  c(
    "Biosynthesis of various plant secondary metabolites [PATH:ko00999]",
    "Drug metabolism - cytochrome P450 [PATH:ko00982]",
    "Flavone and flavonol biosynthesis [PATH:ko00944]",
    "Glycosphingolipid biosynthesis - ganglio series [PATH:ko00604]",
    "Glycosphingolipid biosynthesis - globo and isoglobo series [PATH:ko00603]",
    "Insect hormone biosynthesis [PATH:ko00981]",
    "Isoquinoline alkaloid biosynthesis [PATH:ko00950]",
    "Mannose type O-glycan biosynthesis [PATH:ko00515]",
    "Metabolism of xenobiotics by cytochrome P450 [PATH:ko00980]",
    "Photosynthesis [PATH:ko00195]",
    "Steroid hormone biosynthesis [PATH:ko00140]",
    "Tropane, piperidine and pyridine alkaloid biosynthesis [PATH:ko00960]",
    "Various types of N-glycan biosynthesis [PATH:ko00513]",
    "Primary bile acid biosynthesis [PATH:ko00120]"
  )


df.diffabund.pws <- 
  merge(df.lefse.raw %>% 
        filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_class == "Metabolism") %>%
        mutate(ef_lda = ifelse(enrich_group == str_remove(comp, "\\|.*"), -ef_lda, ef_lda)) %>% 
        select(ko_met_class = ko_subclass, ko_pathway = feature, enrich_group, 
               lefse_ef_lda = ef_lda, lefse_pval = pvalue, lefse_padj = padj),
      
      df.ancom.raw %>% 
        filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_class == "Metabolism") %>% 
        select(ko_met_class = ko_subclass, ko_pathway = feature, ancomBC_ef_W = ef_W, ancomBC_pval = pvalue, ancomBC_padj = padj)) %>%
  
  arrange(ko_met_class, enrich_group, desc(abs(lefse_ef_lda))) %>% 
  mutate(enrich_group = ifelse(abs(lefse_ef_lda) >= 2 & lefse_padj <= 0.05, enrich_group, "Not Sig"))
  
openxlsx::write.xlsx(x = df.diffabund.pws, file = paste0(path.tb.dir, "/contrib/diffabund_pw.xlsx"), rowNames = F)


### Graph prep -----------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/contrib"))) { dir.create(str_glue("{path.fg.dir}/contrib"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/contrib"))) { dir.create(str_glue("{path.tb.dir}/contrib"), recursive = T) }

tmp.tb.funtax.top.genus <- 
  ls.top.taxa[["Genus"]] %>% 
  filter(!str_detect(taxa, "^Other:")) %>% 
  select(-my.colors) %>% 
  select(Phylum = Phylum2, Genus = taxa) %>% 
  mutate(across(where(is.factor), ~ as.character(.)))

# Getting the most functionally active genus

tmp.tb.funtax.top.genus <- 
  bind_rows(
    tmp.tb.funtax.top.genus, 
    tmp.tb.funtax.top.genus %>% 
      select(Phylum) %>% 
      distinct() %>% 
      mutate(Genus = str_glue("Other: {str_remove(Phylum, 'p: ')}")) %>%
      add_row(Phylum = "Others", Genus = "Other: Bacteria") %>% 
      add_row(Phylum = "Others", Genus = "Other: Domains")
  ) %>% 
  mutate(Phylum = str_remove(Phylum, 'p: ')) %>% 
  mutate(Phylum = factor(Phylum, 
                         levels = c("Bacillota", 
                                    "Bacteroidota", 
                                    "Verrucomicrobiota", 
                                    "Pseudomonadota", 
                                    "Actinomycetota",
                                    "Others"))) %>% 
  arrange(Phylum) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus)))

tmp.tb.funtax.top.genus %>% pull(Phylum) %>% table()

tmp.tb.funtax.top.genus <- 
  bind_cols(
    tmp.tb.funtax.top.genus, 
    data.frame(my.colors  = c(
      
      colorRampPalette(c("#B06888FF", "#884060FF", "#D0A0A8FF", 
                         "#A090C8FF", "#684870FF", "#C8C0F8FF",
                         "#C24841FF", "#A40000FF", "#FF0000FF",
                         "#8810B0FF", "#381858FF", "#F0D8E8FF",
                         "#ED444AFF", "#C8133BFF", "#F88080FF",
                         "#F8B8C0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Bacillota"))),
      
      colorRampPalette(c("#5080C8FF", "#80D8F8FF", "#182848FF", 
                         "#48A8B8FF", "#507070FF",
                         "#B8D8E0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Bacteroidota"))),
      
      colorRampPalette(c("#E8B800FF", "#A89E5EFF", "#735231FF", 
                         "#FAE9A0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Verrucomicrobiota"))),
      
      colorRampPalette(c("#18C898FF", "#90E880FF", "#508858FF", 
                         "#E0F898FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Pseudomonadota"))),
      
      colorRampPalette(c("#F88048FF", "#A13601FF", "#FF3200FF",
                         "#F8B880FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Actinomycetota"))),
      
      colorRampPalette(c("#C8C8C8FF", "#303028FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Others"))))
      
    )
  )

tmp.colors.funtax <- tmp.tb.funtax.top.genus$my.colors
names(tmp.colors.funtax) <- tmp.tb.funtax.top.genus$Genus

# Creating funtax top gen legend

tmp.ls.gg.leg.funtax.gen <- list()

for (i in levels(tmp.tb.funtax.top.genus$Phylum)) {
  
  tmp.tb.gg.leg.funtax <- tmp.tb.funtax.top.genus %>% filter(Phylum == i)
  
  tmp.ls.gg.leg.funtax.gen[[i]] <- ggpubr::as_ggplot(ggpubr::get_legend(
    
    ggplot(tmp.tb.gg.leg.funtax, aes(x = 1, y = Genus, fill = Genus)) + 
      geom_point(shape = 22, size = 3) + 
      scale_fill_manual(name = i, values = tmp.tb.gg.leg.funtax$my.colors) +  
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
  
}

tmp.gg.leg.funtax <- (wrap_plots(tmp.ls.gg.leg.funtax.gen) + plot_layout(nrow = 1))
tmp.gg.leg.funtax

{
  png(filename = paste0(path.fg.dir, "/legends2.png"), width = 5000 * 2.5, height = 5000 * .5, units = "px", res = 1200, type = "cairo")
  
  print(tmp.gg.leg.funtax)
  
  dev.off()
}

tmp.gg.leg.funtax <- (wrap_plots(tmp.ls.gg.leg.funtax.gen) + plot_layout(ncol = 1, heights = c(10, 5, 3, 4, 4, 3)))
tmp.gg.leg.funtax

{
  png(filename = paste0(path.fg.dir, "/legends3.png"), width = 5000 * .5, height = 5000 * 1.15, units = "px", res = 1200, type = "cairo")

  print(tmp.gg.leg.funtax)
  
dev.off()
}

tb.ko.tax.counts <-
  tb.work %>%
  merge(tb.ann, ., all.y = T) %>%
  separate_longer_delim(., ko_id, delim = ",") %>%
  merge(., tb.aux.ann.ko.def, all.x = T) %>%
  select(ortholog_id, starts_with("ko_"), tax_id, any_of(tb.metad %>% filter(timepoint == "post-HFD") %>% pull(sample_id))) %>%
  filter(!is.na(ko_pathway)) %>%
  select(ortholog_id, ko_class, ko_subclass, ko_pathway, tax_id, any_of(tb.metad$sample_id)) %>%
  distinct() %>%
  group_by(ko_class, ko_subclass, ko_pathway, tax_id) %>%
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  data.frame(check.names = F) %>% 
  mutate(across(where(is.numeric), ~./sum(.)))

tmp.scs <- 
  df.lefse.filt %>% 
  filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_class == "Metabolism") %>% 
  pull(ko_subclass) %>% 
  unique()

ls.funtax.pw <- list()

for (sc in tmp.scs) {
  
  # Selecting set:
  tmp.funtax.pws <- 
    df.lefse.filt %>% 
    filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_subclass == sc) %>%
    pull(feature) %>% 
    unique() %>% 
    sort()
  
  tmp.funtax.pws <- tmp.funtax.pws[!tmp.funtax.pws %in% tmp.to.filter]
  
  tmp.lefse <- 
    df.lefse.filt %>% 
    filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_subclass == sc) %>%
    mutate(sig.f = case_when(padj <= 0.001 ~ "***",
                             padj <= 0.01 ~ "**",
                             padj <= 0.05 ~ "*",
                             padj <= 0.1 ~ "~",
                             .default = NA_character_)) %>% 
    mutate(my.label = ifelse(enrich_group == str_remove(comp, "\\|.*"), str_glue("−{round(ef_lda, 2)} {sig.f}"), str_glue("{round(ef_lda, 2)} {sig.f}")))
  
  # Genus agglom
  tmp.funtax.sc.genus <- 
    tb.ko.tax.counts %>% 
    filter(ko_pathway %in% tmp.funtax.pws) %>% 
    merge(., tb.tax, by = "tax_id", all.x = T) %>% 
    select(ko_class, ko_subclass, ko_pathway, Kingdom, Phylum, Genus, any_of(tb.metad$sample_id)) %>% 
    group_by(across(-any_of(tb.metad$sample_id))) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  
  # Top genus
  
  tmp.funtax.sc.top.genus <- 
    tmp.funtax.sc.genus %>% 
    mutate(Phylum2 = ifelse(Phylum %in% paste("p:", levels(tmp.tb.funtax.top.genus$Phylum)), str_remove(Phylum, "p: "), "Others")) %>%
    filter(!str_detect(Phylum, "Unassigned")) %>% 
    mutate(Genus2 = case_when(
      Phylum2 == "Others" & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
      Phylum2 == "Others" & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
      Phylum2 != "Others" & Genus %in% paste(levels(tmp.tb.funtax.top.genus$Genus) %>% str_remove(., " \\([0-9]\\)$")) ~ Genus,
      Phylum2 != "Others" & !Genus %in% paste(levels(tmp.tb.funtax.top.genus$Genus) %>% str_remove(., " \\([0-9]\\)$")) ~ str_glue("Other: {Phylum2}")
    )) %>%
    group_by(Phylum2, Genus2, ko_class, ko_subclass, ko_pathway) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
    merge(., 
          data.frame(Genus2 = levels(tmp.tb.funtax.top.genus$Genus) %>% str_remove(., " \\([0-9]\\)$"),
                     Genus3 = levels(tmp.tb.funtax.top.genus$Genus)))
  
  # Sample
  
  tmp.ancom.filt.pws <- df.ancom.filt$feature %>% unique()
  
  tmp.tb.funtax.sc.prof.ra <- 
    tmp.funtax.sc.top.genus %>% 
    mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
    mutate(Genus3 = factor(Genus3, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "ra") %>% 
    merge(tb.metad, .) %>% 
    mutate(donor_ethnicity2 = str_sub(donor_ethnicity, 1, 1)) %>% 
    group_by(ko_pathway, donor_ethnicity2, Genus3) %>%
    summarize(ra = mean(ra, na.rm = TRUE), .groups = 'drop') %>%
    mutate(ko_pathway2 = 
             ifelse(ko_pathway %in% c(tmp.ancom.filt.pws), 
                    str_glue("<sup>†</sup> {ko_pathway}"), 
                    ko_pathway) %>% 
             str_replace(., " \\[PATH:", "<br>["))
  
  tmp.tb.funtax.sc.summary <- 
    tmp.tb.funtax.sc.prof.ra %>% 
    group_by(donor_ethnicity2, ko_pathway, ko_pathway2) %>% 
    summarise(ra = sum(ra, na.rm = T), .groups = "drop") %>% 
    merge(., tmp.lefse %>% select(ko_pathway = feature, enrich_group, ef_lda, my.label) %>% distinct()) %>% 
    group_by(ko_pathway, ko_pathway2) %>% 
    mutate(my.label.pos = max(ra, na.rm = T)) %>% 
    mutate(my.mean = mean(ra, na.rm = T)) %>% 
    arrange(-my.mean) %>% 
    arrange(enrich_group) %>%
    mutate(ko_pathway2 = factor(ko_pathway2, levels = unique(.$ko_pathway2)))
  
  ggplot(tmp.tb.funtax.sc.prof.ra %>% 
           mutate(ko_pathway2 = factor(ko_pathway2, levels = levels(tmp.tb.funtax.sc.summary$ko_pathway2))), 
         aes(y = donor_ethnicity2)) + 
    geom_bar(aes(x = ra, fill = Genus3), stat = "identity", position = "stack", width = .75) +
    geom_text(data = tmp.tb.funtax.sc.summary, aes(x = my.label.pos, y = 1.5, label = my.label), vjust = .5, hjust = -.25, size = 2.5, color = "black", check_overlap = T) +
    facet_grid(ko_pathway2 ~ ., scales = "free", space = "free", switch = "y") +
    scale_x_continuous(expand = expansion(mult = c(0, .2), add = 0), labels = scales::percent_format(accuracy = .01, suffix = "")) +
    scale_y_discrete(expand = expansion(mult = 0.025, add = .75)) +
    scale_fill_manual(values = tmp.colors.funtax) +  
    xlab("<br>Relative abundance (%)") + ylab("") +
    labs(title = paste0(sc, "<br>")) +
    coord_transform(x = "sqrt", clip = "off") +
    my.theme.2 + 
    theme(strip.text.y.left = element_markdown(size = 8, colour = 'black', face = "plain", angle = 0, hjust = 1, lineheight = 1.15),
          strip.placement = "outside",
          panel.spacing = unit(.1, "lines"),
          plot.title = element_markdown(size = 10, face = 'bold', hjust = 1, colour = 'black'), 
          plot.margin = unit(c(0.5, 1, 0.5, 0.75), 'cm'),
          axis.text.x = element_markdown(size = 7, colour = 'black', angle = 45, hjust = 1, vjust = 1)) -> tmp.gg.funtax; tmp.gg.funtax
  
  ls.funtax.pw[[sc]] <- tmp.gg.funtax
  
  tmp.funtax.sc.genus %>% 
    rowwise() %>% 
    mutate(mean_of_cols = (mean(c_across(any_of(tb.metad$sample_id)), na.rm = T))) %>%
    ungroup() %>% 
    arrange(ko_pathway, -mean_of_cols) %>% 
    select(-mean_of_cols) %>% 
    data.frame(check.names = F) -> tmp.1
  
  tmp.funtax.sc.genus %>% 
    pivot_longer(cols = any_of(tb.metad$sample_id), names_to = "sample_id", values_to = "ra") %>% 
    merge(., tb.metad %>% select(sample_id, donor_ethnicity), all.x = T) %>% 
    group_by(across(-c(sample_id, ra))) %>% 
    summarise(ra = sum(ra, na.rm = T)) %>% 
    pivot_wider(names_from = "donor_ethnicity", values_from = "ra") %>% 
    rowwise() %>% 
    mutate(mean_of_cols = (mean(c_across(any_of(tb.metad$donor_ethnicity)), na.rm = T))) %>%
    ungroup() %>% 
    arrange(ko_pathway, -mean_of_cols) %>% 
    select(-mean_of_cols) %>% 
    data.frame(check.names = F) -> tmp.2
  
  tmp.2 %>% 
    group_by(across(starts_with("ko_"))) %>% 
    slice_head(n = 10) -> tmp.3
  
  openxlsx::write.xlsx(x = list("sample" = tmp.1, "group" = tmp.2, "top5" = tmp.3), file = paste0(path.tb.dir, "/contrib/contrib_pw-", str_replace_all(str_to_lower(sc), " ", "_"), ".xlsx"), rowNames = F)
  
}

ls.funtax.pw

my.hei <- 
  df.lefse.filt %>% 
  filter(comp == "Belgian: post-HFD|Yanomami: post-HFD" & ko_class == "Metabolism") %>% 
  filter(!feature %in% tmp.to.filter) %>% 
  mutate(ko_subclass = factor(ko_subclass, tmp.scs)) %>% 
  arrange(ko_subclass) %>% 
  select(ko_subclass, feature) %>%
  distinct() %>% 
  pull(ko_subclass) %>% 
  table()

unname(my.hei)

tmp.gg.funtax <- (wrap_plots(ls.funtax.pw) + plot_layout(ncol = 1, heights = unname(my.hei)))
tmp.gg.funtax

{
  png(filename = paste0(path.fg.dir, "/contrib/contrib_pws.png"), width = 2500 * 1.35, height = 2500 * 10, units = "px", res = 600, type = "cairo")
  
  print(tmp.gg.funtax)
  
  dev.off()
}


# End ---------------------------------------------------------------------

