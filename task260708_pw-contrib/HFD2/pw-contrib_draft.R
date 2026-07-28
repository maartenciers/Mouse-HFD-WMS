
# Functional -------------------------------------------------------------------


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
library("microbiomeMarker")
library("Maaslin2")
library("tidyverse")

getwd()
setwd("~/YanomamiMouseHFD/")
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
        label.type == "p.adj" & adj.p.value < 0.001 ~ "< 0.001",
        label.type == "p.adj" ~ paste0(sprintf("%1.3f", round(adj.p.value, 3)))
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
      geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
      ggpubr::stat_pvalue_manual(
        tb.brackets,
        vjust = -1,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "label",
        tip.length = 0.01,
        bracket.size = 0.25,
        size = 2.5,
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
            legend.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
    
    return(gg)
    
  }
  
  gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                       fill = .data[[gp_fill]], color = .data[[gp_col]])) +
    facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
    geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
    scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.1 + 
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
  
  return(gg)
  
}


fun.get.boxp.compm2 <- function(tb, tb.compm, gp_comp, gp_fill = NULL, gp_col = NULL,
                                gp_facet, group_var, value_var, only.sig = F, alpha = 0.1, 
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
      geom_line(aes(group = .data[[group_var]]), color = "black", alpha = .25, linewidth = .25, 
                position = position_jitter(w = 0.1, h = 0, seed = 123)) +
      geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
      ggpubr::stat_pvalue_manual(
        tb.brackets,
        vjust = -1,
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        label = "label",
        tip.length = 0.01,
        bracket.size = 0.25,
        size = 2.5,
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
            legend.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
    
    return(gg)
    
  }
  
  gg <- ggplot(tb, aes(x = .data[[gp_comp]], y = .data[[value_var]],
                       fill = .data[[gp_fill]], color = .data[[gp_col]])) +
    facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
    geom_line(aes(group = .data[[group_var]]), color = "black", alpha = .25, linewidth = .25, 
              position = position_jitter(w = 0.1, h = 0, seed = 123)) +
    geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
    scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.1 + 
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
  
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
        size = 2,
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
        size = 2.5,
        hide.ns = only.sig,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
      scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
      coord_cartesian(clip = "off") +
      ylab("") + xlab("") +
      my.theme.1 + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
    
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
      size = 2,
      color = "black"
    ) +
    
    scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
    coord_cartesian(clip = "off") +
    ylab("") + xlab("") +
    my.theme.1 + 
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
  
  return(gg)
  
}


## Data loading ----------------------------------------------------------------


### HFD data -------------------------------------------------------------------


# Metadata
tb.metad <-
  read.delim(paste0(getwd(), "/misc/metadata_full.tsv")) %>%
  filter(experiment == my.exp) %>%
  mutate(
    timepoint = case_when(
      experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
      experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD")
  ) %>%
  mutate(donor_ethnicity = ifelse(donor_ethnicity == "Western", "Belgian", donor_ethnicity)) %>%
  select(sample_id, assembly_id, donor_id:cage_id, mouse_n, timepoint) %>%
  mutate(donor_ethnicity = factor(donor_ethnicity, c("Belgian", "Yanomami"))) %>%
  mutate(timepoint = factor(timepoint, c("pre-HFD", "post-HFD"))) %>%
  arrange(donor_ethnicity, timepoint, mouse_n) %>%
  mutate(donor_ethnicity_timepoint = str_glue("{donor_ethnicity}: {timepoint}")) %>%
  mutate(donor_ethnicity_timepoint = factor(donor_ethnicity_timepoint, unique(.$donor_ethnicity_timepoint))) %>% 
  mutate(donor_mouse_n = str_remove_all(sample_id, "_T[12]$")); head(tb.metad)

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


### Info -----------------------------------------------------------------------


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


# Loading DAA tables ------------------------------------------------------


tmp.path <- path.tb.dir %>% str_replace(., "task260708_pw-contrib", "task260619_functional-profiler")

df.lefse.raw <-
  openxlsx::read.xlsx(xlsxFile = paste0(tmp.path, "/diff_abund/lefse.xlsx"), rowNames = F)

df.lefse.filt <- df.lefse.raw %>% filter(ef_lda > 2, padj < 0.05)



## Graph settings --------------------------------------------------------------


# Themes
my.theme.1 <- 
  theme_bw(base_family = "Liberation Sans") + 
  theme(legend.position = "none", 
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0.5),
        
        strip.text = element_text(size = 10, colour = 'black', face = "bold"),
        strip.background = element_rect(fill = NA, colour = NA),
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0, fill = NA),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_text(size = 10, colour = 'black'),
        axis.title.y = element_text(size = 10, colour = 'black'),
        axis.title.y.right = element_text(size = 10, colour = 'black'),
        axis.text.x = element_text(size = 7, colour = 'black'),
        axis.text.y = element_text(size = 7, colour = 'black'),
        
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
        strip.background = element_rect(fill = NA, colour = NA),
        
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0, fill = NA),
        
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_markdown(size = 10, colour = 'black'),
        axis.title.y = element_markdown(size = 10, colour = 'black'),
        axis.title.y.right = element_markdown(size = 10, colour = 'black'),
        axis.text.x = element_markdown(size = 7, colour = 'black'),
        axis.text.y = element_markdown(size = 7, colour = 'black'),
        
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

colors.ko <- c("#F4C40FFF", "#E6352FFF", "#3D79F3FF", "#34A74BFF")
names(colors.ko) <- levels(tb.aux.ann.ko$ko_class)

colors.ko.met <- c("#F4C40FFF", "#FE9B00FF", "#D8443CFF", "#633372FF", "#9F5691FF", "#DE597CFF", "#E6A2A6FF",  "#92C051FF", "#2B9B81FF", "#78A5C3FF", "#1F6E9CFF")
names(colors.ko.met) <- levels(tb.aux.ann.ko.met$ko_subclass)

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

gg.leg.ko <-
  ggplot(tb.aux.ann.ko, aes(x = 1, y = ko_class, fill = ko_class)) +
  geom_point(shape = 22, size = 3) +
  scale_fill_manual(name = "KEGG annotation", values = colors.ko) +
  my.theme.1 +
  theme(legend.position = "right",
        legend.title = element_text(size = 8, colour = 'black', face = "bold"),
        legend.text = element_text(size = 7, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))
gg.leg.ko <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.ko)); gg.leg.ko

gg.leg.ko.met <-
  ggplot(tb.aux.ann.ko.met, aes(x = 1, y = ko_subclass, fill = ko_subclass)) +
  geom_point(shape = 22, size = 3) +
  scale_fill_manual(name = "KEGG Metabolism annotation", values = colors.ko.met) +
  my.theme.1 +
  theme(legend.position = "right",
        legend.title = element_text(size = 8, colour = 'black', face = "bold"),
        legend.text = element_text(size = 7, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))
gg.leg.ko.met <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.ko.met)); gg.leg.ko.met

gg.leg.all <- gg.leg.gp + gg.leg.gp.tp + gg.leg.tax + gg.leg.ko + gg.leg.ko.met + plot_layout(ncol = 1, heights = c(3, 5, 8, 5, 12))
gg.leg.all

png(filename = paste0(path.fg.dir, "/legends.png"), width = 5000 * .65, height = 5000 * 1.5, units = "px", res = 1200, type = "cairo")

gg.leg.all

dev.off()


# -------------------------------------------------------------------------


names(colors.tax) <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota", "Other: Bacteria", "Other: Domains")

tb.ko.tax.counts <-
  tb.work %>%
  merge(tb.ann, ., all.y = T) %>% 
  separate_longer_delim(., ko_id, delim = ",") %>% 
  distinct() %>% 
  merge(., tb.aux.ann.ko.def, by = "ko_id", all.x = T) %>% 
  select(ortholog_id, ko_subclass, tax_id, any_of(tb.metad %>% filter(timepoint == "post-HFD") %>% pull(sample_id))) %>% 
  distinct() %>% 
  group_by(ko_subclass, tax_id) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  data.frame(check.names = F) %>% 
  mutate(across(where(is.numeric), ~./sum(.)))

tb.ko.tax.counts %>% filter(ko_subclass == "Carbohydrate metabolism") %>% select(any_of(tb.metad$sample_id)) %>% colSums()


#for (tmp.funtax.sc in rel.funtax.fun %>% filter(ko_class == "Metabolism") %>% pull(ko_subclass) %>% unique() %>% sort()) {

tmp.funtax.sc = "Lipid metabolism"

print(tmp.funtax.sc)

tmp.funtax.fun <- tb.aux.ann.ko.def %>% filter(ko_subclass == tmp.funtax.sc)

# Getting the most functionally active genus

tmp.tb.funtax.genus <- 
  tb.ko.tax.counts %>% 
  filter(ko_subclass == tmp.funtax.sc) %>% 
  merge(., tb.tax, by = "tax_id", all.x = T) %>% 
  mutate(across(any_of(tb.metad$sample_id), ~replace_na(., 0))) %>% 
  select(Kingdom:Genus, any_of(tb.metad$sample_id)) %>% 
  group_by(across(-any_of(tb.metad$sample_id))) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  rowwise() %>%
  mutate(mean_of_cols = mean(c_across(any_of(tb.metad$sample_id)), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(-mean_of_cols) %>% 
  filter(mean_of_cols > 0) %>%
  select(mean_of_cols, everything()); head(tmp.tb.funtax.genus)

tmp.tb.funtax.genus.perc <- 
  tmp.tb.funtax.genus %>% 
  mutate(across(where(is.numeric), ~./sum(.))) %>% 
  rowwise() %>%
  mutate(mean_of_cols = (mean(c_across(any_of(tb.metad$sample_id))))*100, na.rm = TRUE) %>%
  ungroup() %>%
  arrange(-mean_of_cols) %>% 
  select(mean_of_cols, everything())

tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Bacillota") # 10
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Bacteroidota") # 5
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Verrucomicrobiota") # 2
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Pseudomonadota") # 2
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Actinomycetota") # 2

tmp.tb.funtax.top.genus <- 
  bind_rows(
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Bacillota") %>% 
      head(15),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Bacteroidota") %>% 
      head(5),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Verrucomicrobiota") %>% 
      head(3),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Pseudomonadota") %>% 
      head(3),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Actinomycetota") %>% 
      head(3)
    
  ); tmp.tb.funtax.top.genus

tmp.tb.funtax.genus.perc %>%
  filter(Genus %in% tmp.tb.funtax.top.genus$Genus) %>%
  select(any_of(tb.metad$sample_id)) %>%
  colSums() * 100

tmp.tb.funtax.top.genus <- 
  bind_rows(
    tmp.tb.funtax.top.genus, 
    tmp.tb.funtax.top.genus %>% 
      select(Phylum) %>% 
      distinct() %>% 
      mutate(Genus = str_glue("Other: {str_remove(Phylum, 'p: ')}")) %>%
      add_row(Phylum = "Others", Genus = "Other: Bacteria") %>% 
      add_row(Phylum = "Others", Genus = "Other: Domains") %>% 
      add_row(Phylum = "Others", Genus = "Unassigned")
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
                         "#F87890FF", "#980000FF", "#F04030FF",
                         "#B860B8FF", "#802880FF", "#E098D0FF",
                         "#F05870FF", "#782840FF", "#F87890FF",
                         "#F8B8C0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Bacillota"))),
      
      colorRampPalette(c("#5080C8FF", "#4858A0FF", "#A8D0F8FF", "#48A8B8FF", "#507070FF",
                         "#B8D8E0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Bacteroidota"))),
      
      colorRampPalette(c("#E8B800FF", 
                         "#FAE9A0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Verrucomicrobiota"))),
      
      colorRampPalette(c("#18C898FF", "#90E880FF", "#508858FF", 
                         "#E0F898FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Pseudomonadota"))),
      
      colorRampPalette(c("#F88048FF", "#A13601FF", 
                         "#F8B880FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Actinomycetota"))),
      
      colorRampPalette(c("#C8C8C8FF", "#888890FF", "#303028FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Others"))))
      
    )
  )

tmp.colors.funtax <- tmp.tb.funtax.top.genus$my.colors
names(tmp.colors.funtax) <- tmp.tb.funtax.top.genus$Genus

# Creating funtax top gen legend

tmp.ls.gg.leg.funtax.gen <- list()
tmp.ls.hei.leg.funtax.gen <- list()

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
  
  tmp.ls.hei.leg.funtax.gen[[i]] <- nrow(tmp.tb.gg.leg.funtax) + 3
  
}

tmp.gg.leg.funtax <- (wrap_plots(tmp.ls.gg.leg.funtax.gen) + plot_layout(nrow = 1))
tmp.gg.leg.funtax

# Genus agglom

tmp.funtax.sc.genus <- 
  tb.ko.tax.counts %>% 
  filter(ko_subclass == tmp.funtax.sc) %>% 
  merge(., tb.tax, by = "tax_id", all.x = T) %>% 
  select(ko_subclass, Kingdom, Phylum, Genus, any_of(tb.metad$sample_id)) %>% 
  group_by(across(-any_of(tb.metad$sample_id))) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# Top genus

tmp.funtax.sc.top.genus <- 
  tmp.funtax.sc.genus %>% 
    mutate(Phylum2 = ifelse(Phylum %in% paste("p:", levels(tmp.tb.funtax.top.genus$Phylum)), str_remove(Phylum, "p: "), "Others")) %>%
    mutate(Genus2 = case_when(
      Phylum2 == "Others" & str_detect(Phylum, "Unassigned") ~ "Unassigned",
      Phylum2 == "Others" & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
      Phylum2 == "Others" & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
      Phylum2 != "Others" & Genus %in% paste(levels(tmp.tb.funtax.top.genus$Genus)) ~ Genus,
      Phylum2 != "Others" & !Genus %in% paste(levels(tmp.tb.funtax.top.genus$Genus)) ~ str_glue("Other: {Phylum2}")
    )) %>%
    group_by(Phylum2, Genus2, ko_subclass) %>% 
    summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# Sample

tmp.tb.funtax.sc.prof.cpm <- 
  tmp.funtax.sc.top.genus %>% 
  mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
  mutate(Genus2 = factor(Genus2, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "cpm") %>% 
  merge(tb.metad, .)

tmp.tb.funtax.sc.prof.perc <- 
  tmp.funtax.sc.top.genus %>% 
  mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
  mutate(Genus2 = factor(Genus2, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
  mutate(across(.cols = any_of(tb.metad$sample_id), .fns = ~./sum(.))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "perc") %>% 
  merge(tb.metad, .)

tmp.tb.funtax.sc.prof <- merge(tmp.tb.funtax.sc.prof.cpm, tmp.tb.funtax.sc.prof.perc)

ggplot(tmp.tb.funtax.sc.prof, aes(x = sample_id, y = cpm, fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0), labels = fun.get.scientific) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Counts per million (CPM) mapped reads<br>") +
  labs(title = paste0(tmp.funtax.sc, "<br>")) +
  my.theme.1 + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'),
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 3.25, colour = 'black', angle = 75, hjust = 1, vjust = 1)) #-> tmp.gg.funtax.sc.barprof.abs

ggplot(tmp.tb.funtax.sc.prof, aes(x = sample_id, y = perc, fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme.1 + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'), 
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 3.25, colour = 'black', angle = 75, hjust = 1, vjust = 1)) #-> tmp.gg.funtax.sc.barprof.rel





ggplot(tb.metad, aes(x = sample_id, y = 1, fill = Group_tribe)) + 
  geom_tile() + 
  geom_text(data = tb.metad %>% group_by(Group_tribe) %>% mutate(x_pos = n()/2), 
            aes(label = Group_tribe, x = x_pos), hjust = 0.5, size = 2.5, fontface = "bold", check_overlap = TRUE) + 
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0)) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  facet_grid(. ~ Group_tribe, scales = "free", space = "free") +
  labs(title = paste0(tmp.funtax.sc, "<br>")) +
  scale_fill_manual(values = colors.tr) + 
  theme_void(base_family = "Liberation Sans") + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black'),
        panel.border = element_rect(color = "black", linewidth = .25, fill = NA),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0.5, 0.1, 0.05, 0.75), 'cm')) -> tmp.gg.funtax.sc.barprof.bar

tmp.gg.funtax.sc.barprof.bar + (plot_spacer() + theme_void()) + tmp.gg.funtax.sc.barprof.abs + tmp.gg.leg.funtax + plot_layout(ncol = 2, heights = c(0.075, 1), widths = c(.4, .6)) -> tmp.gg.funtax.sc.barprof.a
tmp.gg.funtax.sc.barprof.bar + (plot_spacer() + theme_void()) + tmp.gg.funtax.sc.barprof.rel + tmp.gg.leg.funtax + plot_layout(ncol = 2, heights = c(0.075, 1), widths = c(.4, .6)) -> tmp.gg.funtax.sc.barprof.b

png(filename = paste0("task260106_funcProfile/figures/taxa_contrib/ko_subclass/bp-ko-sc-taxa_", str_remove_all(str_replace_all(str_remove_all(str_to_lower(tmp.funtax.sc), "[[:punct:]]"), " ", "-"), "-[[:alnum:]]$"), ".png"), width = 5000 * 3.5, height = 5000 * 1.5, units = "px", res = 1200, type = "cairo")

print((tmp.gg.funtax.sc.barprof.a / tmp.gg.funtax.sc.barprof.b))

dev.off()

#ggsave(filename = paste0("task260106_funcProfile/figures/taxa_contrib/ko_subclass/bp-ko-sc-taxa_", str_remove_all(str_replace_all(str_remove_all(str_to_lower(tmp.funtax.sc), "[[:punct:]]"), " ", "-"), "-[[:alnum:]]$"), ".svg"), plot = (tmp.gg.funtax.sc.barprof.a / tmp.gg.funtax.sc.barprof.b), dpi = 1200, width = 5000 * 3.5, height = 5000 * 1.5, units = "px")
#ggsave(filename = paste0("task260106_funcProfile/figures/taxa_contrib/ko_subclass/bp-ko-sc-taxa_", str_remove_all(str_replace_all(str_remove_all(str_to_lower(tmp.funtax.sc), "[[:punct:]]"), " ", "-"), "-[[:alnum:]]$"), ".png"), plot = (tmp.gg.funtax.sc.barprof.a / tmp.gg.funtax.sc.barprof.b), dpi = 1200, width = 5000 * 3.5, height = 5000 * 1.5, units = "px")

# Group

tmp.tb.funtax.sc.prof.cpm <- 
  tmp.tb.funtax.sc.top.genus %>% 
  mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
  mutate(Genus2 = factor(Genus2, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "cpm") %>% 
  merge(tb.metad, .) %>% 
  group_by(Group_tribe, Phylum2, Genus2) %>% 
  summarise(cpm = mean(cpm, na.rm = T), .groups = "drop")

tmp.tb.funtax.sc.prof.perc <- 
  tmp.tb.funtax.sc.top.genus %>% 
  mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
  mutate(Genus2 = factor(Genus2, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
  mutate(across(.cols = tb.metad$sample_id, .fns = ~./sum(.))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "perc") %>% 
  merge(tb.metad, .) %>% 
  group_by(Group_tribe, Phylum2, Genus2) %>% 
  summarise(perc = mean(perc, na.rm = T), .groups = "drop") %>% 
  group_by(Group_tribe) %>% 
  mutate(perc = perc/sum(perc, na.rm = T))

tmp.tb.funtax.sc.prof <- merge(tmp.tb.funtax.sc.prof.cpm, tmp.tb.funtax.sc.prof.perc)

ggplot(tmp.tb.funtax.sc.prof, aes(x = Group_tribe, y = cpm, fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ Group_tribe, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0), labels = scientific) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Counts per million (CPM) mapped reads<br>") +
  my.theme + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'), 
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> tmp.gg.funtax.sc.barprof.abs; tmp.gg.funtax.sc.barprof.abs

ggplot(tmp.tb.funtax.sc.prof, aes(x = Group_tribe, y = perc, fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ Group_tribe, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'), 
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> tmp.gg.funtax.sc.barprof.rel; tmp.gg.funtax.sc.barprof.rel

ggplot(tb.metad, aes(x = Group_tribe, y = 1, fill = Group_tribe)) + 
  geom_tile() + 
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0)) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  facet_grid(. ~ Group_tribe, scales = "free", space = "free") +
  scale_fill_manual(values = colors.tr) + 
  coord_cartesian(clip = "off") +
  labs(title = paste0(tmp.funtax.sc, "<br>")) +
  theme_void(base_family = "Liberation Sans") + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black', hjust = 0),
        panel.border = element_rect(color = "black", linewidth = .25, fill = NA),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0.5, 0.1, 0.05, 0.75), 'cm')) -> tmp.gg.funtax.sc.barprof.bar

tmp.gg.funtax.sc.barprof.bar + (plot_spacer() + theme_void()) + tmp.gg.funtax.sc.barprof.abs + tmp.gg.leg.funtax + plot_layout(ncol = 2, heights = c(0.075, 1), widths = c(.15, .85)) -> tmp.gg.funtax.sc.barprof.a
tmp.gg.funtax.sc.barprof.bar + (plot_spacer() + theme_void()) + tmp.gg.funtax.sc.barprof.rel + tmp.gg.leg.funtax + plot_layout(ncol = 2, heights = c(0.075, 1), widths = c(.15, .85)) -> tmp.gg.funtax.sc.barprof.b

png(filename = paste0("task260106_funcProfile/figures/taxa_contrib/ko_subclass/bp-ko-sc-taxa-group_", str_remove_all(str_replace_all(str_remove_all(str_to_lower(tmp.funtax.sc), "[[:punct:]]"), " ", "-"), "-[[:alnum:]]$"), ".png"), width = 5000 * 2.5, height = 5000 * 1.5, units = "px", res = 1200, type = "cairo")

print((tmp.gg.funtax.sc.barprof.a / tmp.gg.funtax.sc.barprof.b))

dev.off()

#ggsave(filename = paste0("task260106_funcProfile/figures/taxa_contrib/ko_subclass/bp-ko-sc-taxa-group_", str_remove_all(str_replace_all(str_remove_all(str_to_lower(tmp.funtax.sc), "[[:punct:]]"), " ", "-"), "-[[:alnum:]]$"), ".svg"), plot = (tmp.gg.funtax.sc.barprof.a / tmp.gg.funtax.sc.barprof.b), dpi = 1200, width = 5000 * 2.5, height = 5000 * 1.5, units = "px")
#ggsave(filename = paste0("task260106_funcProfile/figures/taxa_contrib/ko_subclass/bp-ko-sc-taxa-group_", str_remove_all(str_replace_all(str_remove_all(str_to_lower(tmp.funtax.sc), "[[:punct:]]"), " ", "-"), "-[[:alnum:]]$"), ".png"), plot = (tmp.gg.funtax.sc.barprof.a / tmp.gg.funtax.sc.barprof.b), dpi = 1200, width = 5000 * 2.5, height = 5000 * 1.5, units = "px")

}




tb.ko.tax.counts <-
  tb.work %>%
  merge(tb.ann, ., all.y = T) %>% 
  separate_longer_delim(., ko_id, delim = ",") %>% 
  distinct() %>% 
  merge(., tb.aux.ann.ko.def, by = "ko_id", all.x = T) %>% 
  select(ortholog_id, ko_class, ko_subclass, ko_pathway, tax_id, any_of(tb.metad %>% filter(timepoint == "post-HFD") %>% pull(sample_id))) %>% 
  distinct() %>% 
  group_by(ko_class, ko_subclass, ko_pathway, tax_id) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  data.frame(check.names = F) %>% 
  mutate(across(where(is.numeric), ~./sum(.)))

tb.ko.tax.counts %>% filter(ko_pathway == "Propanoate metabolism [PATH:ko00640]") %>% select(any_of(tb.metad$sample_id)) %>% colSums()+1 %>% lo



df.lefse.filt %>% filter(ko_subclass == "Carbohydrate metabolism") %>% pull(feature) %>% unique() %>% sort()

tmp.funtax.pw = "Propanoate metabolism [PATH:ko00640]" 

print(tmp.funtax.pw)

tmp.funtax.fun <- tb.aux.ann.ko.def %>% filter(ko_pathway == tmp.funtax.pw)

# Getting the most functionally active genus

tmp.tb.funtax.genus <- 
  tb.ko.tax.counts %>% 
  filter(ko_pathway == tmp.funtax.pw) %>% 
  merge(., tb.tax, by = "tax_id", all.x = T) %>% 
  mutate(across(any_of(tb.metad$sample_id), ~replace_na(., 0))) %>% 
  select(Kingdom:Genus, any_of(tb.metad$sample_id)) %>% 
  group_by(across(-any_of(tb.metad$sample_id))) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  rowwise() %>%
  mutate(mean_of_cols = mean(c_across(any_of(tb.metad$sample_id)), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(-mean_of_cols) %>% 
  filter(mean_of_cols > 0) %>%
  select(mean_of_cols, everything()); head(tmp.tb.funtax.genus)

tmp.tb.funtax.genus.perc <- 
  tmp.tb.funtax.genus %>% 
  mutate(across(where(is.numeric), ~./sum(.))) %>% 
  rowwise() %>%
  mutate(mean_of_cols = (mean(c_across(any_of(tb.metad$sample_id))))*100, na.rm = TRUE) %>%
  ungroup() %>%
  arrange(-mean_of_cols) %>% 
  select(mean_of_cols, everything())

tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Bacillota") # 10
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Bacteroidota") # 5
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Verrucomicrobiota") # 2
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Pseudomonadota") # 2
tmp.tb.funtax.genus.perc %>% select(Kingdom, Phylum, Genus) %>% filter(Phylum == "p: Actinomycetota") # 2

tmp.tb.funtax.top.genus <- 
  bind_rows(
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Bacillota") %>% 
      head(15),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Bacteroidota") %>% 
      head(5),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Verrucomicrobiota") %>% 
      head(3),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Pseudomonadota") %>% 
      head(3),
    
    tmp.tb.funtax.genus.perc %>% 
      select(Phylum, Genus) %>% 
      filter(Phylum == "p: Actinomycetota") %>% 
      head(3)
    
  ); tmp.tb.funtax.top.genus

tmp.tb.funtax.genus.perc %>%
  filter(Genus %in% tmp.tb.funtax.top.genus$Genus) %>%
  select(any_of(tb.metad$sample_id)) %>%
  colSums() * 100

tmp.tb.funtax.top.genus <- 
  bind_rows(
    tmp.tb.funtax.top.genus, 
    tmp.tb.funtax.top.genus %>% 
      select(Phylum) %>% 
      distinct() %>% 
      mutate(Genus = str_glue("Other: {str_remove(Phylum, 'p: ')}")) %>%
      add_row(Phylum = "Others", Genus = "Other: Bacteria") %>% 
      add_row(Phylum = "Others", Genus = "Other: Domains") %>% 
      add_row(Phylum = "Others", Genus = "Unassigned")
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
                         "#F87890FF", "#980000FF", "#F04030FF",
                         "#B860B8FF", "#802880FF", "#E098D0FF",
                         "#F05870FF", "#782840FF", "#F87890FF",
                         "#F8B8C0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Bacillota"))),
      
      colorRampPalette(c("#5080C8FF", "#4858A0FF", "#A8D0F8FF", "#48A8B8FF", "#507070FF",
                         "#B8D8E0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Bacteroidota"))),
      
      colorRampPalette(c("#E8B800FF", 
                         "#FAE9A0FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Verrucomicrobiota"))),
      
      colorRampPalette(c("#18C898FF", "#90E880FF", "#508858FF", 
                         "#E0F898FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Pseudomonadota"))),
      
      colorRampPalette(c("#F88048FF", "#A13601FF", 
                         "#F8B880FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Actinomycetota"))),
      
      colorRampPalette(c("#C8C8C8FF", "#888890FF", "#303028FF"))(nrow(tmp.tb.funtax.top.genus %>% filter(Phylum == "Others"))))
      
    )
  )

tmp.colors.funtax <- tmp.tb.funtax.top.genus$my.colors
names(tmp.colors.funtax) <- tmp.tb.funtax.top.genus$Genus

# Creating funtax top gen legend

tmp.ls.gg.leg.funtax.gen <- list()
tmp.ls.hei.leg.funtax.gen <- list()

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
  
  tmp.ls.hei.leg.funtax.gen[[i]] <- nrow(tmp.tb.gg.leg.funtax) + 3
  
}

tmp.gg.leg.funtax <- (wrap_plots(tmp.ls.gg.leg.funtax.gen) + plot_layout(nrow = 1))
tmp.gg.leg.funtax

# Genus agglom
tmp.funtax.sc.genus <- 
  tb.ko.tax.counts %>% 
  filter(ko_pathway == tmp.funtax.pw) %>% 
  merge(., tb.tax, by = "tax_id", all.x = T) %>% 
  select(ko_subclass = ko_pathway, Kingdom, Phylum, Genus, any_of(tb.metad$sample_id)) %>% 
  group_by(across(-any_of(tb.metad$sample_id))) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# Top genus

tmp.funtax.sc.top.genus <- 
  tmp.funtax.sc.genus %>% 
  mutate(Phylum2 = ifelse(Phylum %in% paste("p:", levels(tmp.tb.funtax.top.genus$Phylum)), str_remove(Phylum, "p: "), "Others")) %>%
  mutate(Genus2 = case_when(
    Phylum2 == "Others" & str_detect(Phylum, "Unassigned") ~ "Unassigned",
    Phylum2 == "Others" & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
    Phylum2 == "Others" & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
    Phylum2 != "Others" & Genus %in% paste(levels(tmp.tb.funtax.top.genus$Genus)) ~ Genus,
    Phylum2 != "Others" & !Genus %in% paste(levels(tmp.tb.funtax.top.genus$Genus)) ~ str_glue("Other: {Phylum2}")
  )) %>%
  group_by(Phylum2, Genus2, ko_subclass) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# Sample

tmp.tb.funtax.sc.prof.cpm <- 
  tmp.funtax.sc.top.genus %>% 
  mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
  mutate(Genus2 = factor(Genus2, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "cpm") %>% 
  merge(tb.metad, .)

tmp.tb.funtax.sc.prof.perc <- 
  tmp.funtax.sc.top.genus %>% 
  mutate(Phylum2 = factor(Phylum2, levels(tmp.tb.funtax.top.genus$Phylum))) %>% 
  mutate(Genus2 = factor(Genus2, levels(tmp.tb.funtax.top.genus$Genus))) %>% 
  mutate(across(.cols = any_of(tb.metad$sample_id), .fns = ~./sum(.))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "perc") %>% 
  merge(tb.metad, .)

tmp.tb.funtax.sc.prof <- merge(tmp.tb.funtax.sc.prof.cpm, tmp.tb.funtax.sc.prof.perc)

ggplot(tmp.tb.funtax.sc.prof, aes(x = donor_ethnicity, y = cpm, fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
  facet_grid(donor_ethnicity ~ ., scales = "free", space = "free") +
  #scale_y_continuous(transform = "log2") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Counts per million (CPM) mapped reads<br>") +
  labs(title = paste0(tmp.funtax.pw, "<br>")) +
  coord_flip() + 
  my.theme.1 + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'),
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 3.25, colour = 'black', angle = 75, hjust = 1, vjust = 1)) #-> tmp.gg.funtax.sc.barprof.abs

ggplot(tmp.tb.funtax.sc.prof, aes(x = donor_ethnicity, y = log2(cpm+1), fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(donor_ethnicity ~ ., scales = "free", space = "free") +
  #scale_y_continuous(transform = "log2") +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Counts per million (CPM) mapped reads<br>") +
  labs(title = paste0(tmp.funtax.pw, "<br>")) +
  coord_flip() + 
  my.theme.1 + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'),
        plot.title = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 3.25, colour = 'black', angle = 75, hjust = 1, vjust = 1)) #-> tmp.gg.funtax.sc.barprof.abs



ggplot(tmp.tb.funtax.sc.prof, aes(x = sample_id, y = perc, fill = Genus2)) + 
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = tmp.colors.funtax) +  
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme.1 + 
  theme(strip.text = element_blank(),
        panel.spacing = unit(.1, "lines"),
        plot.margin = unit(c(0, 0.1, 0.5, 0.75), 'cm'), 
        axis.title.y = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 8, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 3.25, colour = 'black', angle = 75, hjust = 1, vjust = 1)) #-> tmp.gg.funtax.sc.barprof.rel




# THIS IS THE RIGHT ONE

tb.work.fun.ko <-
  tb.work %>%
  select(ortholog_id, tax_id, any_of(tb.metad$sample_id)) %>%
  group_by(ortholog_id, tax_id) %>%
  summarise(across(where(is.numeric), ~sum(.x, na.rm = T)), .groups = "drop") %>%
  merge(tb.ann, ., all.y = T) %>%
  separate_longer_delim(., ko_id, delim = ",") %>%
  merge(., tb.aux.ann.ko.def, all.x = T) %>%
  select(ortholog_id, tax_id, starts_with("ko_"), any_of(tb.metad$sample_id)) %>%
  filter(!is.na(ko_class))

tb.ko.pw.counts <-
  tb.work.fun.ko %>%
  filter(!is.na(ko_pathway)) %>%
  select(ortholog_id, ko_pathway, tax_id, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_pathway, tax_id) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_pathway) %>%
  data.frame(check.names = F)

tb.ko.pw.counts %>% filter(ko_pathway == "ABC transporters [PATH:ko02010]") %>% select(any_of(tb.metad$sample_id)) %>% colSums()

####






tb.ko.tax.counts <-
  tb.work %>%
  merge(tb.ann, ., all.y = T) %>% 
  separate_longer_delim(., ko_id, delim = ",") %>% 
  distinct() %>% 
  merge(., tb.aux.ann.ko.def, by = "ko_id", all.x = T) %>% 
  #select(ortholog_id, ko_pathway, tax_id, any_of(tb.metad %>% filter(timepoint == "post-HFD") %>% pull(sample_id))) %>% 
  select(ortholog_id, ko_pathway, tax_id, any_of(tb.metad %>% pull(sample_id))) %>% 
  distinct() %>% 
  group_by(ko_pathway, tax_id) %>% 
  summarise(across(.cols = any_of(tb.metad$sample_id), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  data.frame(check.names = F) -> tmp.0

tmp.0 %>% filter(ko_pathway == "ABC transporters [PATH:ko02010]") %>% select(any_of(tb.metad$sample_id)) %>% colSums()





tb.ko.pw.counts %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "cpm") %>% 
  merge(., tb.metad) -> tmp.1

tb.ko.pw.counts <-
  tb.work.fun.ko %>%
  filter(!is.na(ko_pathway)) %>%
  select(ortholog_id, ko_pathway, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_pathway) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_pathway) %>%
  data.frame(check.names = F) %>% 
  mutate(across(.cols = tb.metad$sample_id, .fns = ~./sum(.)))

tb.ko.pw.counts %>% 
  pivot_longer(cols = where(is.numeric), names_to = "sample_id", values_to = "cpm") %>% 
  merge(., tb.metad) -> tmp.1


fun.get.summary(tmp.1, gp_comp = "donor_ethnicity_timepoint", gp_facet = "ko_pathway", value_var = "cpm") %>% View()
