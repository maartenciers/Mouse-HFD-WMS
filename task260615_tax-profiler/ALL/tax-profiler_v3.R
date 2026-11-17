
# Taxonomy and Ecology ---------------------------------------------------------


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

path.wk.dir <- paste0(getwd(), "/task260615_tax-profiler")
path.dt.dir <- paste0(getwd(), "/tables")

for (i in c(path.wk.dir)) {
  
  if (!dir.exists(i)) { dir.create(i) }
  
}; rm(i)

my.exp <- "ALL"

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

fun.get.boxp.compm1 <- function(tb, tb.compm, gp_comp, gp_fill = NULL, gp_col = NULL, gp_shape = NULL,
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
      {if(is.null(gp_shape)) geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
      {if(!is.null(gp_shape)) geom_point(aes(shape = .data[[gp_shape]]), size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
      {if(!is.null(gp_shape)) scale_shape_manual(values = shapes.exp)} +
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
    {if(is.null(gp_shape)) geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
    {if(!is.null(gp_shape)) geom_point(aes(shape = .data[[gp_shape]]), size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
    {if(!is.null(gp_shape)) scale_shape_manual(values = shapes.exp)} +    
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

fun.get.boxp.compm3 <- function(tb, tb.compm, tb.kruskal, gp_comp,  gp_fill = NULL, gp_col = NULL, gp_shape = NULL,
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
      {if(is.null(gp_shape)) geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
      {if(!is.null(gp_shape)) geom_point(aes(shape = .data[[gp_shape]]), size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
      {if(!is.null(gp_shape)) scale_shape_manual(values = shapes.exp)} +
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
    {if(is.null(gp_shape)) geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
    {if(!is.null(gp_shape)) geom_point(aes(shape = .data[[gp_shape]]), size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
    {if(!is.null(gp_shape)) scale_shape_manual(values = shapes.exp)} +
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
  filter(experiment_type == "high-fat diet") %>% 
  mutate(timepoint = case_when(
    experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
    experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD",
  )) %>% 
  mutate(mouse_n = str_glue("{donor_id}")) %>% 
  mutate(donor_ethnicity = ifelse(donor_ethnicity == "Western", "Belgian", donor_ethnicity)) %>% 
  select(experiment, sample_id, assembly_id, donor_id:cage_id, mouse_n, timepoint) %>% 
  mutate(donor_ethnicity = factor(donor_ethnicity, c("Belgian", "Yanomami"))) %>% 
  mutate(timepoint = factor(timepoint, c("pre-HFD", "post-HFD"))) %>% 
  mutate(experiment = factor(experiment, levels = unique(.$experiment) %>% sort())) %>% 
  arrange(experiment, donor_ethnicity, timepoint, mouse_n) %>% 
  mutate(donor_ethnicity_timepoint = str_glue("{donor_ethnicity}: {timepoint}")) %>% 
  mutate(donor_ethnicity_timepoint = factor(donor_ethnicity_timepoint, unique(.$donor_ethnicity_timepoint))) %>% 
  mutate(mouse_n = str_replace(mouse_n, "Donor-2-MCBL1", "BEL-051")) %>%
  mutate(mouse_n = factor(mouse_n, unique(.$mouse_n))); head(tb.metad)

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
  pivot_longer(cols = any_of(tb.metad$sample_id), names_to = "sample_id", values_to = "cpm") %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  merge(tb.metad %>% select(-c(cage_id)), ., by = "sample_id") %>% 
  group_by(across(-c(sample_id, cpm))) %>% 
  summarise(cpm = mean(cpm, na.rm = T), .groups = "drop") %>% 
  mutate(sample_id = ifelse(timepoint == "pre-HFD", str_glue("{mouse_n}_T1"), str_glue("{mouse_n}_T2"))) %>% 
  select(tax_id:Species, sample_id, cpm) %>% 
  pivot_wider(names_from = "sample_id", values_from = "cpm") %>% 
  data.frame(check.names = F)

tb.metad <- 
  tb.metad %>% 
  select(-c(sample_id, cage_id)) %>% 
  mutate(sample_id = ifelse(timepoint == "pre-HFD", str_glue("{mouse_n}_T1"), str_glue("{mouse_n}_T2"))) %>% 
  select(experiment, sample_id, everything()) %>% 
  distinct() %>% 
  mutate(donor_id = str_replace(donor_id, "Donor-2-MCBL1", "BEL-051")) %>%
  group_by(donor_ethnicity) %>%
  mutate(pair = paste0("pair_", match(donor_id, unique(donor_id)))) %>%
  ungroup() %>% 
  data.frame(check.names = F)

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

colors.exp <- c("#D9D9D9FF", "#969696FF", "#525252FF", "#252525FF", "#000000FF")
names(colors.exp) <- paste0("HFD", 1:5)

shapes.exp <- c(21, 22, 23, 24, 25, 3, 4)
names(shapes.exp) <- unique(tb.metad$pair)

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

gg.leg.exp.shape <-
  ggplot(tb.metad %>% filter(!is.na(pair)), aes(x = 1, y = pair, shape = pair)) + 
  geom_point(size = 2, fill = "black") + 
  scale_shape_manual(values = shapes.exp) + 
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(2, 0, 2, 5)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.exp.shape <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.exp.shape))


gg.leg.all <- gg.leg.gp + gg.leg.gp.tp + gg.leg.tax + gg.leg.exp.shape + plot_layout(ncol = 1, heights = c(3, 5, 8, 8))
gg.leg.all

{
  png(filename = paste0(path.fg.dir, "/legends.png"), width = 5000 * .65, height = 5000 * 1, units = "px", res = 1200, type = "cairo")
  
  print(gg.leg.all)
  
  dev.off()
}


# Analysis ---------------------------------------------------------------------


## Taxa profile ----------------------------------------------------------------


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

save.image(paste0(path.wk.dir, "/", my.exp,  "/tax-profiler_v3.RData"))


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
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[1]]) %>% head(9),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[2]]) %>% head(5),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[3]]) %>% head(2),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[4]]) %>% head(2),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[5]]) %>% head(2)
      
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
          
          c("#B06888FF", "#884060FF", "#D0A0A8FF", "#A090C8FF", "#684870FF", 
            "#C8C0F8FF", "#F87890FF", "#980000FF", "#F04030FF", 
            "#F8B8C0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]]))],
          
          c("#5080C8FF", "#4858A0FF", "#A8D0F8FF", "#48A8B8FF", "#507070FF",
            "#B8D8E0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]]))],
          
          c("#E8B800FF", "#A89E5EFF",
            "#FAE9A0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]]))],
          
          c("#18C898FF", "#90E880FF", 
            "#E0F898FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]]))],
          
          c("#F88048FF", "#A13601FF", 
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
    mutate(across(tb.metad %>% pull(sample_id), ~as.numeric(.))) %>% 
    mutate(across(where(is.numeric), ~./sum(.))) %>% 
    filter(complete.cases(.)) %>% 
    mutate(across(tb.metad %>% pull(sample_id), ~as.numeric(.))) %>% 
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
    ggplot(tmp.tb.glom.top %>% 
             mutate(group = paste0(str_sub(donor_ethnicity, 1, 1), "_", str_remove_all(variable, ".*_"))) %>% 
             mutate(d_ethnicity = paste0(str_sub(donor_ethnicity, 1, 1))) %>% 
             mutate(t_point = paste0(str_remove_all(variable, ".*_"))), aes(fill = taxa2, x = mouse_n, y = value)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    facet_grid(. ~ d_ethnicity + t_point, scales = "free", space = "free") +
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
gg.bp.s <- gg.bp.plots[[6]]; gg.bp.s

ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_1p.png"), plot = gg.bp.p, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.8, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_2c.png"), plot = gg.bp.c, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.8, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_3o.png"), plot = gg.bp.o, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.8, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_4f.png"), plot = gg.bp.f, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.8, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_5g.png"), plot = gg.bp.g, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.8, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_6s.png"), plot = gg.bp.s, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.8, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))


### Group ----------------------------------------------------------------------


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
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[1]]) %>% head(9),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[2]]) %>% head(5),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[3]]) %>% head(2),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[4]]) %>% head(2),
      
      tmp.tb.glom %>% select(Phylum2, taxa) %>% filter(Phylum2 == top.phylum[[5]]) %>% head(2)
      
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
          
          c("#B06888FF", "#884060FF", "#D0A0A8FF", "#A090C8FF", "#684870FF", 
            "#C8C0F8FF", "#F87890FF", "#980000FF", "#F04030FF", 
            "#F8B8C0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[1]]))],
          
          c("#5080C8FF", "#4858A0FF", "#A8D0F8FF", "#48A8B8FF", "#507070FF",
            "#B8D8E0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[2]]))],
          
          c("#E8B800FF", "#A89E5EFF",
            "#FAE9A0FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[3]]))],
          
          c("#18C898FF", "#90E880FF", 
            "#E0F898FF")[if(nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]])) == 0) 0 else 1:nrow(tmp.tb.top.taxa %>% filter(Phylum2 == top.phylum[[4]]))],
          
          c("#F88048FF", "#A13601FF", 
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
gg.bp.s <- gg.bp.plots[[6]]; gg.bp.s

ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_1p.png"), plot = gg.bp.p, dpi = 1200, width = 5000 * 3, height = 5000 * 0.65, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_2c.png"), plot = gg.bp.c, dpi = 1200, width = 5000 * 3, height = 5000 * 0.65, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_3o.png"), plot = gg.bp.o, dpi = 1200, width = 5000 * 3, height = 5000 * 0.65, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_4f.png"), plot = gg.bp.f, dpi = 1200, width = 5000 * 3, height = 5000 * 0.65, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_5g.png"), plot = gg.bp.g, dpi = 1200, width = 5000 * 3, height = 5000 * 0.65, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_6s.png"), plot = gg.bp.s, dpi = 1200, width = 5000 * 3, height = 5000 * 0.65, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))

openxlsx::write.xlsx(x = ls.glom.tb.filt, file = paste0(path.tb.dir, "/profiles/taxa-glom-filt-table.xlsx"), rowNames = F)


### Summaries ------------------------------------------------------------------


ls.glom.tb
ls.glom.tb.filt

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
  
}

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
  
  
}

bind_rows(ls.glom.unclass.group) %>% 
  merge(tb.metad, ., by = "sample_id") %>% 
  mutate(rank = factor(rank, levels = unique(rank))) -> tb.glom.unclass.group

fun.get.summary(tb = tb.glom.unclass.group, gp_comp = c("experiment", "donor_ethnicity_timepoint"), gp_facet = c("status", "rank"), value_var = "perc") %>% 
  arrange(rank, donor_ethnicity_timepoint, status) %>% 
  mutate(full_mean = str_glue("{round(mean, 4)*100} (±{round(sd, 4)*100})")) -> tb.glom.unclass.group

tb.glom.unclass.group

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
    group_by(experiment, donor_ethnicity_timepoint, across(2:(2+i))) %>% 
    summarise(mean = mean(perc), sd = sd(perc), .groups = "drop") %>% 
    mutate(full_mean = str_glue("{round(mean, 4)*100} (±{round(sd, 4)*100})")) %>% 
    select(-c(mean,sd)) %>% 
    pivot_wider(names_from = c(experiment, donor_ethnicity_timepoint), values_from = full_mean) %>% 
    select(1:(1+i), starts_with(levels(tb.metad$experiment))) %>% 
    arrange(across(all_of(tmp.rank)))
  
}

# Saving Tables

ls.glom.unclass
tb.glom.unclass.group

ls.glom.tb.filt.perc
ls.glom.tb.filt.perc.group

openxlsx::write.xlsx(x = ls.glom.unclass, file = paste0(path.tb.dir, "/profiles/taxa-unclass-sample.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = tb.glom.unclass.group, file = paste0(path.tb.dir, "/profiles/taxa-unclass-group.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = ls.glom.tb.filt.perc, file = paste0(path.tb.dir, "/profiles/taxa-glom-perc-sample.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = ls.glom.tb.filt.perc.group, file = paste0(path.tb.dir, "/profiles/taxa-glom-perc-group.xlsx"), rowNames = F)

rm(ls.glom.unclass, tb.glom.unclass.group)
rm(ls.glom.tb.filt.perc, ls.glom.tb.filt.perc.group)


## Alpha diversity -------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diversity"))) { dir.create(str_glue("{path.fg.dir}/diversity"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diversity"))) { dir.create(str_glue("{path.tb.dir}/diversity"), recursive = T) }

fun.get.alpha.div.ps <- function(my.ps.table.ls) {
  
  tmp.ls <- list()
  
  for (i in seq_along(my.ps.table.ls)) {
    
    tmp.ls[[names(my.ps.table.ls)[[i]]]] <- 
      microbiome::alpha(my.ps.table.ls[[i]], index = c("observed", "diversity_shannon", "diversity_gini_simpson")) %>% 
      select(observed, diversity_shannon, diversity_gini_simpson) %>% 
      rename("Observed" = 1, "Shannon index" = 2, "Gini-Simpson" = 3) %>% 
      rownames_to_column("sample_id") %>% 
      pivot_longer(cols = where(is.numeric), names_to = "var", values_to = "value") %>%
      mutate(var = factor(var, c("Observed", "Shannon index", "Gini-Simpson"))) %>%
      mutate(value = replace_na(value, 0)) %>%
      mutate(rank = names(my.ps.table.ls)[[i]]) %>%
      merge(tb.metad, .) %>%
      arrange(var) 
    
  }; rm(i)
  
  return(bind_rows(tmp.ls))
  
}

tb.adiv <- fun.get.alpha.div.ps(ls.glom.ps.filt)

# T1

my.cond <- "pre-HFD"

ls.adiv.summ <- list()
ls.adiv.compm <- list()
ls.adiv.gg <- list()

for (i in ranks) {
  
  tmp.adiv <- tb.adiv %>% filter(timepoint == my.cond & rank == i)
  
  ls.adiv.summ[[i]] <- 
    fun.get.summary(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value") %>% 
    mutate(rank = i); ls.adiv.summ[[i]]
  
  ls.adiv.compm[[i]] <- 
    fun.get.wilcox(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                   p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05) %>% 
    mutate(rank = i) %>% 
    suppressWarnings(); ls.adiv.compm[[i]]
  
  ls.adiv.gg[[i]] <- 
    fun.get.boxp.compm1(tb = tmp.adiv, 
                        tb.compm = ls.adiv.compm[[i]], 
                        gp_facet = "var", 
                        gp_comp = "donor_ethnicity_timepoint",
                        gp_shape = "pair",
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") + 
    labs(title = paste0(i, "\n")) + ylab(paste0(my.cond, "\n")) +  
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]] )
  
}; rm(i)


ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(timepoint == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

{
  png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 4.75, units = "px", res = 1200, type = "cairo")
  
  print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))
  
  dev.off()
}

# T2

my.cond <- "post-HFD"

ls.adiv.summ <- list()
ls.adiv.compm <- list()
ls.adiv.gg <- list()

for (i in ranks) {
  
  tmp.adiv <- tb.adiv %>% filter(timepoint == my.cond & rank == i) %>% droplevels()
  
  ls.adiv.summ[[i]] <- 
    fun.get.summary(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value") %>% 
    mutate(rank = i); ls.adiv.summ[[i]]
  
  ls.adiv.compm[[i]] <- 
    fun.get.wilcox(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                   p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05) %>% 
    mutate(rank = i) %>% 
    suppressWarnings(); ls.adiv.compm[[i]]
  
  ls.adiv.gg[[i]] <- 
    fun.get.boxp.compm1(tb = tmp.adiv, 
                        tb.compm = ls.adiv.compm[[i]], 
                        gp_facet = "var", 
                        gp_comp = "donor_ethnicity_timepoint",
                        gp_shape = "pair",
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") + 
    labs(title = paste0(i, "\n")) + ylab(paste0(my.cond, "\n")) +  
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]] )
  
}; rm(i)


ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(timepoint == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

{
  png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 4.75, units = "px", res = 1200, type = "cairo")
  
  print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))
  
  dev.off()
}


# BEL

my.cond <- "Belgian"

ls.adiv.summ <- list()
ls.adiv.compm <- list()
ls.adiv.gg <- list()

for (i in ranks) {
  
  tmp.adiv <- tb.adiv %>% filter(donor_ethnicity == my.cond & rank == i) %>% droplevels()
  
  ls.adiv.summ[[i]] <- 
    fun.get.summary(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value") %>% 
    mutate(rank = i); ls.adiv.summ[[i]]
  
  ls.adiv.compm[[i]] <- 
    fun.get.wilcox(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                   p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05) %>% 
    mutate(rank = i) %>% 
    suppressWarnings(); ls.adiv.compm[[i]]
  
  # BUG: https://github.com/tidyverse/ggplot2/issues/3535
  
  ls.adiv.gg[[i]] <- 
    fun.get.boxp.compm1(tb = tmp.adiv, 
                        tb.compm = ls.adiv.compm[[i]], 
                        gp_facet = "var", 
                        gp_comp = "donor_ethnicity_timepoint",
                        gp_shape = "pair",
                        #group_var = "mouse_n",
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") + 
    labs(title = paste0(i, "\n")) + ylab(paste0(my.cond, "\n")) +  
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]] )
  
}; rm(i)


ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(donor_ethnicity == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

{
  png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 4.75, units = "px", res = 1200, type = "cairo")
  
  print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))
  
  dev.off()
}


# YANOM

my.cond <- "Yanomami"

ls.adiv.summ <- list()
ls.adiv.compm <- list()
ls.adiv.gg <- list()

for (i in ranks) {
  
  tmp.adiv <- tb.adiv %>% filter(donor_ethnicity == my.cond & rank == i) %>% droplevels()
  
  ls.adiv.summ[[i]] <- 
    fun.get.summary(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value") %>% 
    mutate(rank = i); ls.adiv.summ[[i]]
  
  ls.adiv.compm[[i]] <- 
    fun.get.wilcox(tb = tmp.adiv, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                   p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05) %>% 
    mutate(rank = i) %>% 
    suppressWarnings(); ls.adiv.compm[[i]]
  
  ls.adiv.gg[[i]] <- 
    fun.get.boxp.compm1(tb = tmp.adiv, 
                        tb.compm = ls.adiv.compm[[i]], 
                        gp_facet = "var", 
                        gp_comp = "donor_ethnicity_timepoint",
                        gp_shape = "pair",
                        #group_var = "mouse_n",
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") + 
    labs(title = paste0(i, "\n")) + ylab(paste0(my.cond, "\n")) +  
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]] )
  
}; rm(i)


ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(donor_ethnicity == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

{
  png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 4.75, units = "px", res = 1200, type = "cairo")
  
  print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))
  
  dev.off()
}


rm(list = ls(pattern = "[.]?alpha[.]?"))
rm(list = ls(pattern = "[.]?adiv[.]?"))
rm(my.cond)


## Beta diversity --------------------------------------------------------------

fun.get.beta.div <- function(my.glom.ps) {
  
  my.glom.tb <- merge(tax_table(my.glom.ps), otu_table(my.glom.ps), by = "row.names", all = T) %>% select(-c(Kingdom:Species))
  
  my.metad <- filter(tb.metad, sample_id %in% sample_names(my.glom.ps))
  
  for (i in seq_along(bdiv.indexes)) {
    
    print(names(bdiv.indexes)[i])
    
    pcoa.factors <- my.metad
    
    pcoa.values <- 
      my.glom.tb %>% 
      rename("taxa" = 1) %>% 
      column_to_rownames("taxa") %>% 
      select(my.metad$sample_id) %>% 
      mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
      t()
    
    if (bdiv.indexes[[i]] == "euclidean") {
      
      pcoa.dist <- phyloseq::distance(my.glom.ps, bdiv.indexes[[i]])
      
    } else {
      
      pcoa.dist <- phyloseq::distance(phyloseq_standardize_otu_abundance(physeq = my.glom.ps, method = 'total'), bdiv.indexes[[i]])
      
    }
    
    pcoa.cmds <- cmdscale(pcoa.dist, eig = TRUE, k = 3)
    
    pcoa.site.score <- 
      data.frame(pcoa.cmds$points) %>% 
      rename("PCoA1" = 1, "PCoA2" = 2, "PCoA3" = 3) %>% 
      rownames_to_column("sample_id") %>% 
      merge(., my.metad)
    
    #pcoa.permutations <- with(my.metad, how(nperm = 999, blocks = cage_id))
    
    pcoa.perma.comp1 <- adonis2(pcoa.dist ~ donor_ethnicity + timepoint, data = my.metad, method = bdiv.indexes[[i]], permutations = 999, seed = 123)
    
    pcoa.perma.df.comp1 <- 
      pcoa.perma.comp1 %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column() %>% 
      dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>% 
      mutate(Source = str_to_title(Source),
             sig = case_when(P <= 0.001 ~ "***", P <= 0.01 ~ "**", P <= 0.05 ~ "*", P <= 0.1 ~ ".", .default = " "))
    
    
    pcoa.perma.ph.comp1.1 <- pairwiseAdonis::pairwise.adonis2(pcoa.dist ~ donor_ethnicity_timepoint, strata = "experiment", data = my.metad, p.adjust.m = "fdr", perm = 999)
    
    pcoa.perma.ph.comp1.1 <- bind_rows(pcoa.perma.ph.comp1.1) %>% 
      select(-parent_call) %>% 
      t() %>% 
      data.frame(check.names = T) %>% 
      rename(Model = 1, Residual = 2, Total = 3) %>% 
      rownames_to_column("comparison") %>% 
      separate_wider_delim(cols = comparison, delim = ".", names = c("comparison", "tmp")) %>% 
      pivot_longer(cols = -c(comparison, tmp), names_to = "component", values_to = "value") %>% 
      pivot_wider(names_from = "tmp", values_from = "value") %>% 
      separate_wider_delim(cols = comparison, delim = "_vs_", names = c("group1", "group2")) %>% 
      mutate(strata = "experiment"); pcoa.perma.ph.comp1.1
    
    pcoa.perma.ph.comp1.2 <- pairwiseAdonis::pairwise.adonis2(pcoa.dist ~ donor_ethnicity_timepoint, data = my.metad, p.adjust.m = "fdr", perm = 999)
    
    pcoa.perma.ph.comp1.2 <- bind_rows(pcoa.perma.ph.comp1.2) %>% 
      select(-parent_call) %>% 
      t() %>% 
      data.frame(check.names = T) %>% 
      rename(Model = 1, Residual = 2, Total = 3) %>% 
      rownames_to_column("comparison") %>% 
      separate_wider_delim(cols = comparison, delim = ".", names = c("comparison", "tmp")) %>% 
      pivot_longer(cols = -c(comparison, tmp), names_to = "component", values_to = "value") %>% 
      pivot_wider(names_from = "tmp", values_from = "value") %>% 
      separate_wider_delim(cols = comparison, delim = "_vs_", names = c("group1", "group2")) %>% 
      mutate(strata = NA_character_); pcoa.perma.ph.comp1.2
    
    pcoa.perma.ph <- bind_rows(pcoa.perma.ph.comp1.1, pcoa.perma.ph.comp1.2)
    
    # Axes 1:2
    
    ggplot() +
      stat_ellipse(data = pcoa.site.score,
                   aes(x = PCoA1, y = PCoA2, color = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint),
                   geom = "polygon", linetype = "blank", alpha = .15, show.legend = F, type = "norm") +
      geom_point(pcoa.site.score, 
                 mapping = aes(x = PCoA1, y = PCoA2, fill = donor_ethnicity_timepoint, 
                               color = donor_ethnicity_timepoint, shape = pair), 
                 size = 2, alpha = 1) +
      annotate("richtext", x = Inf, y = -Inf, size = 2.5,
               label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp1[1, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp1[1, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp1[1, "P"], 3)))),
               fill = NA, label.color = NA, hjust = 1, vjust = 0, color = "red") +
      geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
      my.theme.2 +
      scale_shape_manual(values = shapes.exp) +
      scale_color_manual(values = colors.gp.tp) +
      scale_fill_manual(values = colors.gp.tp) +
      labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
           y = paste0("PCoA2 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[2], "%)<br>")) + 
      theme(plot.subtitle = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0.1, 0.1, 0.1), 'cm')) -> gg.pcoa.axes12; gg.pcoa.axes12
    
    # Axes 1:3
    
    ggplot() +
      stat_ellipse(data = pcoa.site.score,
                   aes(x = PCoA1, y = PCoA3, color = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint),
                   geom = "polygon", linetype = "blank", alpha = .15, show.legend = F, type = "norm") +
      geom_point(pcoa.site.score, 
                 mapping = aes(x = PCoA1, y = PCoA3, fill = donor_ethnicity_timepoint, 
                               color = donor_ethnicity_timepoint, shape = pair), 
                 size = 2, alpha = 1) +
      geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
      my.theme.2 +
      scale_shape_manual(values = shapes.exp) +
      scale_color_manual(values = colors.gp.tp) +
      scale_fill_manual(values = colors.gp.tp) +
      labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
           y = paste0("PCoA3 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[3], "%)<br>")) + 
      theme(plot.subtitle = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0.1, 0.1, 0.1), 'cm')) -> gg.pcoa.axes13; gg.pcoa.axes13
    
    
    gg.pcoa.panel <- gg.pcoa.axes12 + gg.pcoa.axes13
    
    {
      png(filename = paste0(path.fg.dir, "/diversity/beta-div-pcoa-", str_to_lower(bdiv.indexes[[i]]), ".png"), width = 5000 * 1.5, height = 5000 * .75, units = "px", res = 1200, type = "cairo")
      
      print(gg.pcoa.panel)
      
      dev.off()
    }
    
    
    pcoa.disp.comp2 <- 
      betadisper(d = pcoa.dist, group = my.metad$donor_ethnicity_timepoint, type = "centroid")$distances %>% 
      data.frame(check.names = F) %>% 
      dplyr::rename("value" = ".") %>% 
      rownames_to_column("sample_id") %>% 
      merge(my.metad, .) %>% 
      mutate(var = "Within-group dispersion")
    
    pcoa.disp.summary.comp2 <- 
      fun.get.summary(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value")
    
    pcoa.disp.kruskal.comp2 <-
      fun.get.kruskal(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value")
    
    pcoa.disp.posthoc.comp2 <-
      fun.get.wilcox(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                     p_adjust_method = "fdr", conf_level = 0.95, tb.stat = pcoa.disp.kruskal.comp2, alpha.global = 0.05)
    
    gg.pcoa.disp.compm.comp2 <- 
      fun.get.boxp.compm3(tb = pcoa.disp.comp2, 
                          #tb.compm = pcoa.disp.compm.comp2,
                          tb.compm = pcoa.disp.posthoc.comp2,
                          tb.kruskal = pcoa.disp.kruskal.comp2,
                          gp_facet = "var", 
                          gp_comp = "donor_ethnicity_timepoint",
                          gp_fill = "donor_ethnicity_timepoint",
                          gp_shape = "pair",
                          value_var = "value",
                          only.sig = T, 
                          alpha = 0.1, 
                          label.type = "stars") + 
      ylab(paste0("Within-group dispersion\n")) +  
      theme(strip.text = element_blank())
    
    pcoa.shift.mt <- as.matrix(pcoa.dist)
    
    pcoa.shift.vc <- sapply(unique(tb.metad$donor_id), function(x) {
      
      pre  <- tb.metad$sample_id[tb.metad$donor_id == x & tb.metad$timepoint == "pre-HFD"]
      post <- tb.metad$sample_id[tb.metad$donor_id == x & tb.metad$timepoint == "post-HFD"]
      
      pcoa.shift.mt[pre, post]
      
    })
    
    pcoa.shift.dt <- 
      data.frame("donor_id" = names(pcoa.shift.vc), "value" = unname(pcoa.shift.vc)) %>% 
      merge(tb.metad %>% select(experiment, donor_id, donor_ethnicity, pair) %>% distinct(), .) %>% 
      mutate(var = "pre- and post-HFD\ndistance")
    
    pcoa.shift.summary.comp2 <- 
      fun.get.summary(tb = pcoa.shift.dt, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value")
    
    pcoa.shift.compm.comp2 <-
      fun.get.wilcox(tb = pcoa.shift.dt, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value",
                     p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05)
    
    gg.pcoa.shift.compm.comp2 <- 
      fun.get.boxp.compm1(tb = pcoa.shift.dt, 
                          tb.compm = pcoa.shift.compm.comp2,
                          #tb.compm = pcoa.disp.posthoc.comp2,
                          #tb.kruskal = pcoa.disp.kruskal.comp2,
                          gp_facet = "var", 
                          gp_comp = "donor_ethnicity",
                          #gp_fill = "donor_ethnicity_timepoint",
                          gp_shape = "pair",
                          value_var = "value",
                          only.sig = T,
                          alpha = 1, 
                          label.type = "p.adj") +
      ylab(paste0("pre- and post-HFD distance\n")) +  
      theme(strip.text = element_blank())
    
    gg.pcoa.panel <- 
      gg.pcoa.disp.compm.comp2 + gg.pcoa.shift.compm.comp2 + 
      plot_layout(nrow = 1, widths = c(65, 35))
    
    {
      png(filename = paste0(path.fg.dir, "/diversity/beta-div-disp-", str_to_lower(bdiv.indexes[[i]]), ".png"), width = 5000 * 1.05, height = 5000 * .8, units = "px", res = 1200, type = "cairo")
      
      print(gg.pcoa.panel)
      
      dev.off()
      
      }
    
    ls.pcoa <-
      
      list(
        "PERMANOVA" = pcoa.perma.df.comp1,
        
        "PERMANOVA (post-hoc)" = pcoa.perma.ph,
        
        "Dispersion (Measures)" = pcoa.disp.comp2,
        "Dispersion (Summary)" = pcoa.disp.summary.comp2,
        "Dispersion (Kruskal)" = pcoa.disp.kruskal.comp2,
        "Dispersion (post-hoc)" = pcoa.disp.posthoc.comp2,
        
        "pre-post dist (Measures)" = pcoa.shift.dt,
        "pre-post dist (Summary)" = pcoa.shift.summary.comp2,
        "pre-post dist (Wilcox)" = pcoa.shift.compm.comp2
        
      )
    
    openxlsx::write.xlsx(x = ls.pcoa, file = paste0(path.tb.dir, "/diversity/beta-div-", str_to_lower(bdiv.indexes[[i]]), ".xlsx"), rowNames = F)
    
  }; rm(i)
  
}

bdiv.indexes <- c("bray")
names(bdiv.indexes) <- c("Bray—Curtis")

fun.get.beta.div(ps.table)

rm(list = ls(pattern = "[.]?beta[.]?"))
rm(list = ls(pattern = "[.]?bdiv[.]?"))
rm(list = ls(pattern = "[.]?pcoa[.]?"))


## MaAsLin2 ---------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diff_abund"))) { dir.create(str_glue("{path.fg.dir}/diff_abund"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diff_abund"))) { dir.create(str_glue("{path.tb.dir}/diff_abund"), recursive = T) }

# Comps:			       
# 1. T1: Bel vs. Yanom	 
# 2. T2: Bel vs. Bic		
# 3. Bel: T1 vs. T2		
# 4. Yanom: T1 vs. T2		

ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

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
  
}

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

openxlsx::write.xlsx(x = df.msl, file = paste0(path.tb.dir, "/diff_abund/maaslin.xlsx"), rowNames = F)

# Graph: Volcano Plot

for (i in c(1:4)) {
  
  for (r in ranks[4:5]) {
    
    df.msl.raw %>% 
      filter(comp == unique(.$comp)[[i]]) %>% 
      filter(rank == r) %>% 
      mutate(sig = ifelse(qval < 0.05 & abs(coef) > 1, as.character(enrich_group), "Not Sig.")) %>% 
      mutate(sig = factor(sig, c(names(colors.gp.tp), "Not Sig."))) %>% 
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
      geom_text_repel(aes(label = my.label, nudge_x = my.nudge.x, hjust = my.hjust, color = sig),
                      direction = "y", size = 2, segment.alpha = .5, segment.color = 'black', segment.size = 0.25,
                      min.segment.length = 0.2, max.overlaps = 999, force = 20, force_pull = 0.5, seed = 123) +
      geom_point(shape = 21, alpha = 1, size = 2, color = "white") + 
      geom_hline(yintercept = -log10(0.05), color = c("black"), linetype = "dotted", linewidth = .25) + 
      geom_vline(xintercept = c(-1, 1), color = c("black"), linetype = "dotted", linewidth = .25) +
      ylab("-log<sub>10</sub> (<i>q</i>-value)<br>") +
      xlab("<br>MaAsLin (β) coefficient") +
      scale_x_continuous(limits = ~ c(-1.5, 1.5) * max(abs(.x))) + 
      scale_fill_manual(name = "Group", values = c(colors.gp.tp, "Not Sig." = "gray75")) +
      scale_color_manual(values = c(colors.gp.tp)) +
      coord_cartesian(clip = "off") +
      my.theme.2 +
      theme(legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank())
    
    png(filename = paste0(path.fg.dir, "/diff_abund/maaslin-vp-comp", i, "_", str_to_lower(r), ".png"), width = 5000 * .85, height = 5000 * .75, units = "px", res = 1200, type = "cairo")
    
    print(tmp.gg.vp)
    
    dev.off()
    
  }
}

rm(i, r)
rm(list = ls(pattern = "[.]?tmp[.]?"))

# Graph: Heatmap

fun.get.heatmap.clr.z <- function(my.msl.table, my.rank, my.top) {
  
  my.metad <- tb.metad %>% arrange(donor_ethnicity_timepoint)
  
  top.taxa <- 
    my.msl.table %>% 
    filter(rank == my.rank) %>% 
    arrange(-abs(coef)) %>% 
    pull(feature) %>% 
    unique() %>% 
    head(my.top)
  
  tmp.ps.clr <-
    ls.glom.ps.filt[[my.rank]] %>% 
    microbiome::transform(., "clr")
  
  tmp.ps.clr <- prune_taxa(tax_table(tmp.ps.clr)[, my.rank] %in% top.taxa, tmp.ps.clr)
  
  tmp.tb.clr <-
    merge(tax_table(tmp.ps.clr), otu_table(tmp.ps.clr), by = "row.names", all = T) %>% 
    select(feature = my.rank, any_of(my.metad$sample_id)) %>% 
    merge(., my.msl.table %>% select(Phylum, Phylum2, feature), all.x = T) %>% 
    distinct()
  
  hp.anno.col <-
    my.metad %>% 
    column_to_rownames("sample_id") %>% 
    select(donor_ethnicity_timepoint, mouse_n)
  
  hp.anno.row <-
    tmp.tb.clr %>% 
    select(Phylum, Phylum2, feature) %>% 
    column_to_rownames("feature")
  
  hp.values <- 
    tmp.tb.clr %>% 
    select(feature, any_of(my.metad$sample_id)) %>% 
    column_to_rownames("feature") %>%
    t() %>% 
    scale() %>% 
    t()
  
  hp.breaks <- max(abs(min(hp.values, na.rm = T)), max(hp.values, na.rm = T))
  
  hp.col.fun.1 <- 
    colorRamp2::colorRamp2(c(seq(-hp.breaks, hp.breaks, .25)),
                           colorRampPalette(c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#FFFFFFFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"))(length(seq(-hp.breaks, hp.breaks, .25))))
  
  hp.a <- 
    ComplexHeatmap::rowAnnotation(
      `Phylum` = 
        hp.anno.row$Phylum2, 
      col = 
        list(`Phylum` = colors.tax), 
      annotation_legend_param = 
        list(`Phylum` = list(labels = names(colors.tax), at = names(colors.tax))), 
      show_annotation_name = FALSE,
      simple_anno_size = unit(.25, "cm"),
      show_legend = c(TRUE))
  
  hp.b <- 
    ComplexHeatmap::columnAnnotation(
      
      `Condition` = 
        hp.anno.col$donor_ethnicity_timepoint,
      
      col = 
        list(
          `Experiment` = colors.exp,
          `Condition` = colors.gp.tp 
        ),
      show_annotation_name = FALSE,
      simple_anno_size = unit(.25, "cm"),
      show_legend = TRUE)
  
  # Plot Row
  ComplexHeatmap::Heatmap(hp.values, 
                          name = "CLR-Z",  
                          col = hp.col.fun.1,
                          cluster_rows = T, 
                          cluster_columns = F,
                          row_names_gp = gpar(fontsize = 9),
                          column_names_gp = gpar(fontsize = 9),
                          column_title = str_glue("Top {nrow(hp.values)} most variable taxa at {my.rank} level\nStandardized CLR-transformed data (CLR-Z)"), 
                          column_title_gp = gpar(fontsize = 12),
                          column_labels = my.metad$mouse_n,
                          show_column_names = T,
                          heatmap_legend_param = 
                            list(color_bar = "continuous",
                                 ncol = 1,
                                 legend_direction = "vertical", 
                                 legend_height = unit(3, "cm"),
                                 grid_width = unit(.25, "cm"),
                                 title_position = "leftcenter-rot",
                                 title_gp = gpar(fontsize = 10),
                                 labels_gp = gpar(fontsize = 8)), 
                          left_annotation = hp.a,
                          top_annotation = hp.b) -> hp.plot
  
  png(filename = paste0(path.fg.dir, "/diff_abund/heatmap-clr-z_", str_to_lower(my.rank),".png"), width = 5000 * 2.25, height = 5000 * (1 + (nrow(hp.values)/60)), units = "px", res = 1200, type = "cairo")
  
  draw(hp.plot, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 10, 1, 50), "mm"))
  
  dev.off()
  
  # Plot Col
  ComplexHeatmap::Heatmap(hp.values, 
                          name = "CLR-Z",  
                          col = hp.col.fun.1,
                          cluster_rows = T, 
                          cluster_columns = T,
                          row_names_gp = gpar(fontsize = 9),
                          column_names_gp = gpar(fontsize = 9),
                          column_title = str_glue("Top {nrow(hp.values)} most variable taxa at {my.rank} level\nStandardized CLR-transformed data (CLR-Z)"), 
                          column_title_gp = gpar(fontsize = 12),
                          column_labels = my.metad$mouse_n,
                          show_column_names = T,
                          heatmap_legend_param = 
                            list(color_bar = "continuous",
                                 ncol = 1,
                                 legend_direction = "vertical", 
                                 legend_height = unit(3, "cm"),
                                 grid_width = unit(.25, "cm"),
                                 title_position = "leftcenter-rot",
                                 title_gp = gpar(fontsize = 10),
                                 labels_gp = gpar(fontsize = 8)), 
                          left_annotation = hp.a,
                          top_annotation = hp.b) -> hp.plot
  
  png(filename = paste0(path.fg.dir, "/diff_abund/heatmap2-clr-z_", str_to_lower(my.rank),".png"), width = 5000 * 2.25, height = 5000 * (1.2 + (nrow(hp.values)/60)), units = "px", res = 1200, type = "cairo")
  
  draw(hp.plot, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 10, 1, 50), "mm"))
  
  dev.off()
  
}

fun.get.heatmap.clr.z(df.msl.raw, "Genus", 50)
fun.get.heatmap.clr.z(df.msl.raw, "Family", 30)

graphics.off()
gc()

save.image(paste0(path.wk.dir, "/", my.exp,  "/tax-profiler_v3.RData"))


# End ---------------------------------------------------------------------

