
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

path.wk.dir <- paste0(getwd(), "/task260708_DiD")
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


# Metadata
tb.metad <- 
  read.delim(paste0(getwd(), "/misc/metadata_full.tsv")) %>% 
  filter(experiment_type == "high-fat diet") %>% 
  mutate(timepoint = case_when(
    experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
    experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD",
  )) %>% 
  mutate(mouse_n = str_glue("{experiment}_{donor_id}")) %>% 
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

# MaAsLin2 tables

tmp.path <- path.tb.dir %>% str_replace(., "task260708_DiD", "task260615_tax-profiler")

ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

comps <- c("Belgian: pre-HFD|Yanomami: pre-HFD",
           "Belgian: post-HFD|Yanomami: post-HFD",
           "Belgian: pre-HFD|Belgian: post-HFD",
           "Yanomami: pre-HFD|Yanomami: post-HFD")

top.phylum <- c("p: Bacillota", "p: Bacteroidota", "p: Verrucomicrobiota", "p: Pseudomonadota", "p: Actinomycetota", "Other: Bacteria", "Other: Domains")

df.msl.raw <- 
  openxlsx::read.xlsx(xlsxFile = paste0(tmp.path, "/diff_abund/maaslin.xlsx"), rowNames = F, sheet = 1) %>% 
  mutate(rank = factor(rank, ranks)) %>% 
  mutate(comp = factor(comp, comps)) %>% 
  mutate(Phylum2 = factor(Phylum2, top.phylum)) %>% 
  mutate(enrich_group = factor(enrich_group, levels(tb.metad$donor_ethnicity_timepoint)))


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
names(colors.tax) <- top.phylum

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


# -------------------------------------------------------------------------

tmp.rank = "Genus"

res_g1 <- df.msl.raw %>% filter(rank == tmp.rank & comp == "Belgian: pre-HFD|Belgian: post-HFD")
res_g2 <- df.msl.raw %>% filter(rank == tmp.rank & comp == "Yanomami: pre-HFD|Yanomami: post-HFD")

signif_g1 <- res_g1 %>% filter(qval < 0.05 & abs(coef) > 1) %>% pull(feature)
signif_g2 <- res_g2 %>% filter(qval < 0.05 & abs(coef) > 1) %>% pull(feature)

shared_features <- intersect(res_g1$feature, res_g2$feature)

concordance <-
  merge(
    res_g1 %>% filter(feature %in% shared_features) %>% select(feature, coef_g1 = coef),
    res_g2 %>% filter(feature %in% shared_features) %>% select(feature, coef_g2 = coef)
  )

# Classify regulatory pattern in one pass
# Boolean membership vectors computed once for clarity
in_g1 <- concordance$feature %in% signif_g1
in_g2 <- concordance$feature %in% signif_g2
in_both <- in_g1 & in_g2


concordance <-
  concordance %>% 
  mutate(regulation = 
           case_when(
             in_g1 & !in_both & concordance$coef_g1 < 0 ~ "Belgian only, pre-HFD higher",
             in_g1 & !in_both & concordance$coef_g1 > 0 ~ "Belgian only, post-HFD higher",
             in_g2 & !in_both & concordance$coef_g2 < 0 ~ "Yanomami only, pre-HFD higher",
             in_g2 & !in_both & concordance$coef_g2 > 0 ~ "Yanomami only, post-HFD higher",
             in_both & sign(coef_g1) == sign(coef_g2) & concordance$coef_g1 > 0 ~ "Concordant, post-HFD higher",
             in_both & sign(coef_g1) == sign(coef_g2) & concordance$coef_g1 < 0 ~ "Concordant, pre-HFD higher",
             in_both & sign(coef_g1) != sign(coef_g2) ~ "Discordant",
             .default = "NS"
           )) %>% 
  mutate(regulation = factor(regulation, c("Concordant, post-HFD higher", 
                                           "Concordant, pre-HFD higher",
                                           "Discordant",
                                           "Belgian only, pre-HFD higher",
                                           "Belgian only, post-HFD higher",
                                           "Yanomami only, pre-HFD higher",
                                           "Yanomami only, post-HFD higher",
                                           "NS")))

table(concordance$regulation)


# Top genes to label
# Subset to concordant genes only
concordant_genes <- 
  concordance %>% 
  filter(!str_detect(regulation, "NS")) %>% 
  rowwise() %>% 
  mutate(mean_coef = mean(abs(coef_g1) + abs(coef_g2)) / 2) %>% ungroup() %>% 
  arrange(-mean_coef) %>%
  mutate(top = 1:n()) %>% 
  mutate(feature = str_remove(feature, "^.*: ")) %>% 
  mutate(my.label = ifelse(top < 30, feature, NA_character_)) %>% 
  mutate(my.nudge.x = case_when( !is.na(my.label) & coef_g1 > 0 & coef_g2 > 0 ~ Inf,
                                 !is.na(my.label) & coef_g1 < 0 & coef_g2 < 0 ~ -Inf,
                                 !is.na(my.label) & coef_g1 > 0 & coef_g2 < 0 ~ Inf,
                                 !is.na(my.label) & coef_g1 < 0 & coef_g2 > 0 ~ -Inf,
                                 .default = NA_integer_ )) %>% 
  mutate(my.nudge.y = case_when( !is.na(my.label) & coef_g1 > 0 & coef_g2 > 0 ~ Inf,
                                 !is.na(my.label) & coef_g1 < 0 & coef_g2 < 0 ~ -Inf,
                                 !is.na(my.label) & coef_g1 > 0 & coef_g2 < 0 ~ -Inf,
                                 !is.na(my.label) & coef_g1 < 0 & coef_g2 > 0 ~ Inf,
                                 .default = NA_integer_ ))



levels(concordance$regulation)

tmp.colors <- c("#343D1FFF", "#8BA749FF", "#8E2322FF", "#FFA500FF", "#FF5500FF", "#0072B2FF", "#2F00B2FF", "grey70")
names(tmp.colors) <- levels(concordance$regulation)

gg.did <-
  ggplot(concordance, aes(x = coef_g1, y = coef_g2, fill = regulation)) +
  
  geom_text_repel(data = concordant_genes, 
                  aes(label = my.label, nudge_x = my.nudge.x, color = regulation),
                  direction = "y", size = 2.25, segment.alpha = .5, segment.color = 'black', segment.size = 0.25,
                  min.segment.length = 0.2, max.overlaps = 999, force = 20, force_pull = 0.5, seed = 123) +
  
  geom_point(data = . %>% filter(regulation == "NS"), shape = 21, color = "white", size = 2) +
  geom_point(data = . %>% filter(regulation != "NS"), shape = 21, color = "white", size = 2.5) +
  
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +

  scale_x_continuous(limits = ~ c(-1.1, 1.1) * max(abs(.x))) +
  scale_y_continuous(limits = ~ c(-1.1, 1.1) * max(abs(.x))) +
  
  scale_fill_manual(values = tmp.colors) +
  scale_color_manual(values = tmp.colors) +
  
  labs(x = "<br>MaAsLin (β) coefficient<br>Belgian: pre- vs post-HFD", 
       y = "MaAsLin (β) coefficient<br>Yanomami: pre- vs post-HFD<br>") +
  
  my.theme.2 +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_markdown(size = 10, colour = 'black', lineheight = 1.25),
        axis.title.y = element_markdown(size = 10, colour = 'black', lineheight = 1.25)); gg.did

gg.leg.tmp <-
  ggplot(names(tmp.colors) %>% data.frame() %>% filter(. != "NS") %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = tmp.colors) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.tmp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.tmp)); gg.leg.tmp


gg.did + gg.leg.tmp + plot_layout(ncol = 2, widths = c(.35, .65)) -> gg.did.leg; gg.did.leg

{
  png(filename = paste0(path.fg.dir, "/DiD_genus.png"), width = 5000 * 2.5, height = 5000 * .85, units = "px", res = 1200, type = "cairo")
  
  print((gg.did.leg))
  
  dev.off()
}

openxlsx::write.xlsx(x = concordance, file = paste0(path.tb.dir, "/maaslin-concordance-table.xlsx"), rowNames = F)


# End ---------------------------------------------------------------------

