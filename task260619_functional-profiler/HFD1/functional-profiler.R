
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

path.wk.dir <- paste0(getwd(), "/task260619_functional-profiler")
path.dt.dir <- paste0(getwd(), "/tables")

for (i in c(path.wk.dir)) {
  
  if (!dir.exists(i)) { dir.create(i) }
  
}; rm(i)

my.exp <- "HFD1"

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
      geom_point(shape = 21, color = "white", size = 2.5, position = position_jitter(w = 0.2, h = 0, seed = 123)) +  
      #eom_jitter(shape = 16, alpha = 1, size = 1, width = .25, height = 0) + 
      #geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      stat_summary(aes(group = .data[[gp_comp]]), geom = "crossbar", fun = "mean", color = "black", width = 0.5, linewidth = .25) +
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
    geom_point(shape = 21, color = "white", size = 2.5, position = position_jitter(w = 0.2, h = 0, seed = 123)) +  
    #eom_jitter(shape = 16, alpha = 1, size = 1, width = .25, height = 0) + 
    #geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    stat_summary(aes(group = .data[[gp_comp]]), geom = "crossbar", fun = "mean", color = "black", width = 0.5, linewidth = .25) +
    
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
                position = position_jitter(w = 0.2, h = 0, seed = 123)) +
      geom_point(shape = 21, color = "white", size = 2.5, position = position_jitter(w = 0.2, h = 0, seed = 123)) +  
      #eom_jitter(shape = 16, alpha = 1, size = 1, width = .25, height = 0) + 
      #geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      stat_summary(aes(group = .data[[gp_comp]]), geom = "crossbar", fun = "mean", color = "black", width = 0.5, linewidth = .25) +
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
              position = position_jitter(w = 0.2, h = 0, seed = 123)) +
    geom_point(shape = 21, color = "white", size = 2.5, position = position_jitter(w = 0.2, h = 0, seed = 123)) +  
    #eom_jitter(shape = 16, alpha = 1, size = 1, width = .25, height = 0) + 
    #geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    stat_summary(aes(group = .data[[gp_comp]]), geom = "crossbar", fun = "mean", color = "black", width = 0.5, linewidth = .25) +
    
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
      geom_point(shape = 21, color = "white", size = 2.5, position = position_jitter(w = 0.2, h = 0, seed = 123)) +  
      #eom_jitter(shape = 16, alpha = 1, size = 1, width = .25, height = 0) + 
      #geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
      stat_summary(aes(group = .data[[gp_comp]]), geom = "crossbar", fun = "mean", color = "black", width = 0.5, linewidth = .25) +
      #stat_summary(fun = "mean", color = "black", fill = "white", shape = 23, size = 0.25, stroke = 0.25) +
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
    geom_point(shape = 21, color = "white", size = 2.5, position = position_jitter(w = 0.2, h = 0, seed = 123)) +  
    #eom_jitter(shape = 16, alpha = 1, size = 1, width = .25, height = 0) + 
    #geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, linewidth = 0.25) +
    stat_summary(aes(group = .data[[gp_comp]]), geom = "crossbar", fun = "mean", color = "black", width = 0.5, linewidth = .25) +
    #stat_summary(fun = "mean", color = "black", fill = "white", shape = 23, size = 0.25, stroke = 0.25) +
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



### Donors ---------------------------------------------------------------------



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

# KO phyloseq table

tb.ko.counts <-
  tb.work %>%
  select(ortholog_id, any_of(tb.metad$sample_id)) %>%
  merge(tb.ann, ., all.y = T) %>%
  select(-ortholog_id) %>%
  separate_longer_delim(., ko_id, delim = ",") %>%
  group_by(ko_id) %>%
  summarise(across(where(is.numeric), ~sum(.x, na.rm = T)), .groups = "drop") %>%
  data.frame(check.names = F)

tmp.tb.metad <- tb.metad %>% rename(ID = sample_id)
tmp.tb.count <- tb.ko.counts %>% select(ID = ko_id, tmp.tb.metad$ID)

tmp.tb.taxon <-
  tb.ko.counts %>%
  mutate(Kingdom = ko_id,
         Phylum = ko_id,
         Class = ko_id,
         Order = ko_id,
         Family = ko_id,
         Genus = ko_id,
         Species = ko_id) %>%
  select(ID = ko_id, Kingdom:Species)

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
ps.table.tss <- phyloseq_standardize_otu_abundance(physeq = ps.table, method = 'total')

rm(list = ls(pattern = "^tmp."))


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
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        panel.border = element_rect(color = "black", linewidth = 0, fill = NA),
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_text(size = 10, face = 'bold', colour = 'black'),
        axis.title.y = element_text(size = 10, face = 'bold', colour = 'black'),
        axis.title.y.right = element_text(size = 10, face = 'bold', colour = 'black'),
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
        strip.background = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray95", linewidth = 0.5, linetype = 1),
        panel.border = element_rect(color = "black", linewidth = 0, fill = NA),
        axis.ticks = element_line(color = "black", linewidth = .25),
        axis.ticks.length = unit(1, "pt"),
        axis.line.x = element_line(colour = 'black', linewidth = 0.25),
        axis.line.y = element_line(colour = 'black', linewidth = 0.25),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.y.right = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black'),
        axis.text.y = element_markdown(size = 6, colour = 'black'),
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


# Diversity --------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diversity"))) { dir.create(str_glue("{path.fg.dir}/diversity"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diversity"))) { dir.create(str_glue("{path.tb.dir}/diversity"), recursive = T) }


### Venn Diagram ----------------------------------------------------------------


fun.get.venn <- function(my.tb.glom, my.column, my.title = "") {
  
  tb.venn <-
    my.tb.glom %>%
    select(any_of(my.column), where(is.numeric)) %>%
    column_to_rownames(my.column) %>%
    t() %>%
    data.frame(check.names = F) %>%
    rownames_to_column("sample_id") %>%
    merge(tb.metad, .) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    group_by(donor_ethnicity_timepoint) %>%
    summarise(across(.cols = where(is.numeric), ~ sum(.x, na.rm = T)), .groups = "drop") %>%
    column_to_rownames("donor_ethnicity_timepoint") %>%
    t() %>%
    data.frame(check.names = F) %>%
    mutate(across(levels(tb.metad$donor_ethnicity_timepoint), ~ifelse(. > 0, 1, 0))) %>%
    rownames_to_column(var = "ID")
  
  ls.venn <- list(
    `Belgian: pre-HFD` = tb.venn %>% filter(`Belgian: pre-HFD` != 0) %>% pull(ID),
    `Belgian: post-HFD` = tb.venn %>% filter(`Belgian: post-HFD` != 0) %>% pull(ID),
    `Yanomami: pre-HFD` = tb.venn %>% filter(`Yanomami: pre-HFD` != 0) %>% pull(ID),
    `Yanomami: post-HFD` = tb.venn %>% filter(`Yanomami: post-HFD` != 0) %>% pull(ID)
  )
  
  gg.venn <-
    ggvenn(ls.venn,
           fill_color = c(colors.gp.tp[["Belgian: pre-HFD"]],
                          colors.gp.tp[["Belgian: post-HFD"]],
                          colors.gp.tp[["Yanomami: pre-HFD"]],
                          colors.gp.tp[["Yanomami: post-HFD"]]),
           stroke_color = "#00000000",
           stroke_size = 0.5,
           set_name_size = 0,
           fill_alpha = .35,
           text_size = 3.25,
           show_percentage = F) +
    labs(title = my.title) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)); gg.venn
  
  ls.venn.exc <-
    list(
      `Belgian__pre-HFD` = setdiff(ls.venn[[1]], Reduce(union, list(ls.venn[[2]], ls.venn[[3]], ls.venn[[4]]))),
      `Belgian__post-HFD` = setdiff(ls.venn[[2]], Reduce(union, list(ls.venn[[1]], ls.venn[[3]], ls.venn[[4]]))),
      `Yanomami__pre-HFD` = setdiff(ls.venn[[3]], Reduce(union, list(ls.venn[[2]], ls.venn[[1]], ls.venn[[4]]))),
      `Yanomami__post-HFD` = setdiff(ls.venn[[4]], Reduce(union, list(ls.venn[[2]], ls.venn[[3]], ls.venn[[1]])))
    )
  
  return(list(gg.venn, ls.venn.exc))
  
}

gg.venn.1 <- fun.get.venn(tb.ko.counts, "ko_id", "KEGG Orthologs (prev > 10%)")

gg.venn.2 <- fun.get.venn(tb.ko.counts %>%
                            mutate("prev" = rowSums(.[tb.metad %>% pull(sample_id) %>% sort()] >= 1)) %>%
                            filter(prev >= (tb.metad %>% pull(sample_id) %>% length())*.25) %>%
                            select(-prev), "ko_id", "KEGG Orthologs (prev > 25%)")

gg.venn <- gg.venn.1[[1]] + gg.venn.2[[1]] + plot_layout(nrow = 1); gg.venn

png(filename = paste0(path.fg.dir, "/diversity/venn-ko.png"), width = 5000 * 2, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

gg.venn

dev.off()

openxlsx::write.xlsx(x = gg.venn.1[[2]], file = paste0(path.tb.dir, "/diversity/venn-10perc.xlsx"), rowNames = F)
openxlsx::write.xlsx(x = gg.venn.2[[2]], file = paste0(path.tb.dir, "/diversity/venn-25perc.xlsx"), rowNames = F)

rm(list = ls(pattern = "[.]?venn[.]?"))


## Alpha diversity -------------------------------------------------------------


ls.glom.ps <- list("KEGG Orthologs" = ps.table)

ls.glom.tb <- list(
  
  "KEGG Orthologs" = merge(tax_table(ls.glom.ps[["KEGG Orthologs"]]), otu_table(ls.glom.ps[["KEGG Orthologs"]]), by = "row.names", all = T)[-1]
  
)

ranks <- c("KEGG Orthologs")

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


tb.adiv <- fun.get.alpha.div.ps(ls.glom.ps)

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
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") +
    labs(title = paste0(i, "\n")) +
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]])
  
}; rm(i)

ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(timepoint == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))

dev.off()

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
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") +
    labs(title = paste0(i, "\n")) +
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]])
  
}; rm(i)

ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(timepoint == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))

dev.off()

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
  
  ls.adiv.gg[[i]] <-
    fun.get.boxp.compm2(tb = tmp.adiv,
                        tb.compm = ls.adiv.compm[[i]],
                        gp_facet = "var",
                        gp_comp = "donor_ethnicity_timepoint",
                        group_var = "mouse_n",
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") +
    labs(title = paste0(i, "\n")) +
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  print(ls.adiv.gg[[i]])
  
}; rm(i)

ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(donor_ethnicity == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))

dev.off()

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
    fun.get.boxp.compm2(tb = tmp.adiv,
                        tb.compm = ls.adiv.compm[[i]],
                        gp_facet = "var",
                        gp_comp = "donor_ethnicity_timepoint",
                        group_var = "mouse_n",
                        value_var = "value",
                        only.sig = T,
                        alpha = 1,
                        label.type = "p.adj") +
    labs(title = paste0(i, "\n")) +
    theme(plot.title = element_text(size = 12, face = 'bold', colour = 'black'))
  
  print(ls.adiv.gg[[i]])
  
}; rm(i)

ls.adiv.res <- list(
  
  "Raw" = tb.adiv %>% filter(donor_ethnicity == my.cond) %>% droplevels(),
  "Summary" = bind_rows(ls.adiv.summ),
  "Wilcox" = bind_rows(ls.adiv.compm)
  
)

openxlsx::write.xlsx(x = ls.adiv.res, file = paste0(path.tb.dir, "/diversity/alpha-div-", my.cond, ".xlsx"), rowNames = F)

png(filename = paste0(path.fg.dir, "/diversity/alpha-div-", my.cond, ".png"), width = 5000 * 1.1, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

print(wrap_plots(ls.adiv.gg) + plot_layout(ncol = 1))

dev.off()

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
    
    pcoa.permutations <- with(my.metad, how(nperm = 999, blocks = cage_id))
    
    pcoa.perma.comp1 <- adonis2(pcoa.dist ~ donor_ethnicity + timepoint, data = my.metad, method = bdiv.indexes[[i]], permutations = pcoa.permutations, seed = 123)
    
    pcoa.perma.df.comp1 <-
      pcoa.perma.comp1 %>%
      data.frame(check.names = F) %>%
      rownames_to_column() %>%
      dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>%
      mutate(Source = str_to_title(Source),
             sig = case_when(P <= 0.001 ~ "***", P <= 0.01 ~ "**", P <= 0.05 ~ "*", P <= 0.1 ~ ".", .default = " "))
    
    pcoa.perma.ph.comp1.1 <- pairwiseAdonis::pairwise.adonis2(pcoa.dist ~ donor_ethnicity_timepoint, strata = "cage_id", data = my.metad, p.adjust.m = "fdr", perm = 999)
    
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
      mutate(strata = "cage_id"); pcoa.perma.ph.comp1.1
    
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
      # stat_ellipse(data = pcoa.site.score,
      #              aes(x = PCoA1, y = PCoA2, fill = donor_ethnicity_timepoint),
      #              geom = "polygon", linetype = 0, alpha = .25, show.legend = F) +
      geom_line(pcoa.site.score, mapping = aes(group = mouse_n, x = PCoA1, y = PCoA2), color = "black", alpha = .25, linewidth = .25) +
      geom_point(pcoa.site.score, mapping = aes(x = PCoA1, y = PCoA2, fill = donor_ethnicity_timepoint), shape = 21, color = "white", size = 3) +
      annotate("richtext", x = Inf, y = -Inf, size = 2.5,
               label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp1[1, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp1[1, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp1[1, "P"], 3)))),
               fill = NA, label.color = NA, hjust = 1, vjust = 0) +
      geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
      my.theme.2 +
      scale_fill_manual(values = colors.gp.tp) +
      labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
           y = paste0("PCoA2 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[2], "%)<br>")) +
      theme(plot.subtitle = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0.1, 0.1, 0.1), 'cm')) -> gg.pcoa.axes12
    
    gg.pcoa.axes12.param <- ggplot_build(gg.pcoa.axes12)
    
    gg.pcoa.axes12.param.ylimits <- gg.pcoa.axes12.param$layout$panel_params[[1]]$y$continuous_range
    
    gg.pcoa.axes12.param.xlimits <- gg.pcoa.axes12.param$layout$panel_params[[1]]$x$continuous_range
    
    pcoa.site.score.longer.axes12 <-
      pcoa.site.score %>%
      pivot_longer(cols = c("PCoA1", "PCoA2"), names_to = "PC", values_to = "value") %>%
      distinct()
    
    pcoa.site.score.summary.axes12 <-
      pcoa.site.score.longer.axes12 %>%
      group_by(donor_ethnicity_timepoint, PC) %>%
      summarise(temp = list(c(psych::describe(value), quantile(value, na.rm = T), IQR = EnvStats::iqr(value, na.rm = T)))) %>%
      unnest_wider(temp) %>%
      select(-c(vars, trimmed, mad)) %>%
      mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>%
      data.frame(check.names = F) %>%
      arrange(PC, donor_ethnicity_timepoint)
    
    pcoa.site.score.kruskal.axes12 <-
      fun.get.kruskal(tb = pcoa.site.score.longer.axes12, gp_comp = "donor_ethnicity_timepoint", gp_facet = "PC", value_var = "value") %>%
      arrange(PC)
    
    pcoa.site.score.posthoc.axes12 <-
      fun.get.wilcox(tb = pcoa.site.score.longer.axes12,
                     tb.stat = pcoa.site.score.kruskal.axes12,
                     gp_comp = "donor_ethnicity_timepoint", gp_facet = "PC", value_var = "value") %>%
      arrange(PC)
    
    pcoa.site.score.longer.axes12 %>%
      filter(PC == "PCoA1") %>%
      ggplot(., aes(x = donor_ethnicity_timepoint, y = value, fill = donor_ethnicity_timepoint, color = donor_ethnicity_timepoint)) +
      geom_boxplot(width = 0.66, color = "black", outlier.size = 1, outlier.shape = 21, linewidth = 0.25) +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(limits = c(gg.pcoa.axes12.param.xlimits[1], gg.pcoa.axes12.param.xlimits[2]), expand = c(0, 0)) +
      scale_fill_manual(values = colors.gp.tp) +
      coord_flip() +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), 'cm')) -> gg.pcoa.pc1.compm.axes12
    
    pcoa.site.score.longer.axes12 %>%
      filter(PC == "PCoA2") %>%
      ggplot(., aes(x = donor_ethnicity_timepoint, y = value, fill = donor_ethnicity_timepoint, color = donor_ethnicity_timepoint)) +
      geom_boxplot(width = 0.66, color = "black", outlier.size = 1, outlier.shape = 21, linewidth = 0.25) +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(limits = c(gg.pcoa.axes12.param.ylimits[1], gg.pcoa.axes12.param.ylimits[2]), expand = c(0, 0)) +
      scale_color_manual(values = colors.gp.tp) +
      scale_fill_manual(values = colors.gp.tp) +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), 'cm')) -> gg.pcoa.pc2.compm.axes12
    
    gg.pcoa.axes12.panel <- gg.pcoa.pc1.compm.axes12 + plot_spacer() + gg.pcoa.axes12 + gg.pcoa.pc2.compm.axes12 + plot_layout(ncol = 2, heights = c(.1, 1), widths = c(1, .1))
    
    # Axes 1:3
    
    ggplot() +
      # stat_ellipse(data = pcoa.site.score,
      #              aes(x = PCoA1, y = PCoA2, fill = donor_ethnicity_timepoint),
      #              geom = "polygon", linetype = 0, alpha = .25, show.legend = F) +
      geom_line(pcoa.site.score, mapping = aes(group = mouse_n, x = PCoA1, y = PCoA3), color = "black", alpha = .25, linewidth = .25) +
      geom_point(pcoa.site.score, mapping = aes(x = PCoA1, y = PCoA3, fill = donor_ethnicity_timepoint), shape = 21, color = "white", size = 3) +
      # annotate("richtext", x = Inf, y = -Inf, size = 2.5,
      #          label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp1[1, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp1[1, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp1[1, "P"], 3)))),
      #          fill = NA, label.color = NA, hjust = 1, vjust = 0) +
      geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted") +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
      my.theme.2 +
      scale_fill_manual(values = colors.gp.tp) +
      labs(x = paste0("<br>PCoA1 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[1], "%)"),
           y = paste0("PCoA3 (", round(pcoa.cmds$eig / sum(pcoa.cmds$eig) * 100, 2)[3], "%)<br>")) +
      theme(plot.subtitle = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0.1, 0.1, 0.1), 'cm')) -> gg.pcoa.axes13
    
    gg.pcoa.axes13.param <- ggplot_build(gg.pcoa.axes13)
    
    gg.pcoa.axes13.param.ylimits <- gg.pcoa.axes13.param$layout$panel_params[[1]]$y$continuous_range
    
    gg.pcoa.axes13.param.xlimits <- gg.pcoa.axes13.param$layout$panel_params[[1]]$x$continuous_range
    
    pcoa.site.score.longer.axes13 <-
      pcoa.site.score %>%
      pivot_longer(cols = c("PCoA1", "PCoA3"), names_to = "PC", values_to = "value") %>%
      distinct()
    
    pcoa.site.score.summary.axes13 <-
      pcoa.site.score.longer.axes13 %>%
      group_by(donor_ethnicity_timepoint, PC) %>%
      summarise(temp = list(c(psych::describe(value), quantile(value, na.rm = T), IQR = EnvStats::iqr(value, na.rm = T)))) %>%
      unnest_wider(temp) %>%
      select(-c(vars, trimmed, mad)) %>%
      mutate(across(where(is.numeric), ~case_when(is.nan(.) ~ NA_integer_, is.infinite(.) ~ NA_integer_,.default = .))) %>%
      data.frame(check.names = F) %>%
      arrange(PC, donor_ethnicity_timepoint)
    
    pcoa.site.score.kruskal.axes13 <-
      fun.get.kruskal(tb = pcoa.site.score.longer.axes13, gp_comp = "donor_ethnicity_timepoint", gp_facet = "PC", value_var = "value") %>%
      arrange(PC)
    
    pcoa.site.score.posthoc.axes13 <-
      fun.get.wilcox(tb = pcoa.site.score.longer.axes13,
                     tb.stat = pcoa.site.score.kruskal.axes13,
                     gp_comp = "donor_ethnicity_timepoint", gp_facet = "PC", value_var = "value") %>%
      arrange(PC)
    
    pcoa.site.score.longer.axes13 %>%
      filter(PC == "PCoA1") %>%
      ggplot(., aes(x = donor_ethnicity_timepoint, y = value, fill = donor_ethnicity_timepoint, color = donor_ethnicity_timepoint)) +
      geom_boxplot(width = 0.66, color = "black", outlier.size = 1, outlier.shape = 21, linewidth = 0.25) +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(limits = c(gg.pcoa.axes13.param.xlimits[1], gg.pcoa.axes13.param.xlimits[2]), expand = c(0, 0)) +
      scale_fill_manual(values = colors.gp.tp) +
      coord_flip() +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), 'cm')) -> gg.pcoa.pc1.compm.axes13
    
    pcoa.site.score.longer.axes13 %>%
      filter(PC == "PCoA3") %>%
      ggplot(., aes(x = donor_ethnicity_timepoint, y = value, fill = donor_ethnicity_timepoint, color = donor_ethnicity_timepoint)) +
      geom_boxplot(width = 0.66, color = "black", outlier.size = 1, outlier.shape = 21, linewidth = 0.25) +
      geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dotted") +
      scale_y_continuous(limits = c(gg.pcoa.axes13.param.ylimits[1], gg.pcoa.axes13.param.ylimits[2]), expand = c(0, 0)) +
      scale_color_manual(values = colors.gp.tp) +
      scale_fill_manual(values = colors.gp.tp) +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), 'cm')) -> gg.pcoa.pc2.compm.axes13
    
    gg.pcoa.axes13.panel <- gg.pcoa.pc1.compm.axes13 + plot_spacer() + gg.pcoa.axes13 + gg.pcoa.pc2.compm.axes13 + plot_layout(ncol = 2, heights = c(.1, 1), widths = c(1, .1))
    
    pcoa.disp.comp2 <- 
      betadisper(d = pcoa.dist, group = my.metad$donor_ethnicity, type = "centroid")$distances %>% 
      data.frame(check.names = F) %>% 
      dplyr::rename("value" = ".") %>% 
      rownames_to_column("sample_id") %>% 
      merge(my.metad, .) %>% 
      mutate(var = "Within-group\ndispersion")
    
    pcoa.disp.summary.comp2 <- 
      fun.get.summary(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value")
    
    # pcoa.disp.kruskal.comp2 <- 
    #   fun.get.kruskal(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value")
    # 
    # pcoa.disp.posthoc.comp2 <- 
    #   fun.get.wilcox(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value",
    #                  p_adjust_method = "fdr", conf_level = 0.95, tb.stat = pcoa.disp.kruskal.comp2, alpha.global = 0.05)
    
    pcoa.disp.compm.comp2 <-
      fun.get.wilcox(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value",
                     p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05)
    
    gg.pcoa.disp.compm.comp2 <- 
      fun.get.boxp.compm1(tb = pcoa.disp.comp2, 
                          tb.compm = pcoa.disp.compm.comp2,
                          #tb.compm = pcoa.disp.posthoc.comp2,
                          #tb.kruskal = pcoa.disp.kruskal.comp2,
                          gp_facet = "var", 
                          gp_comp = "donor_ethnicity",
                          gp_fill = "donor_ethnicity_timepoint",
                          value_var = "value",
                          only.sig = T,
                          alpha = 1, 
                          label.type = "p.adj")
    
    gg.pcoa.panel <-
      gg.pcoa.pc1.compm.axes12 + plot_spacer() + plot_spacer() + gg.pcoa.pc1.compm.axes13 + plot_spacer() + plot_spacer() +
      gg.pcoa.axes12 + gg.pcoa.pc2.compm.axes12 + plot_spacer() + gg.pcoa.axes13 + gg.pcoa.pc2.compm.axes13 + gg.pcoa.disp.compm.comp2 + plot_layout(ncol = 6, heights = c(.1, 1), widths = c(1, .1, .05, 1, .1, .5))
    
    print(gg.pcoa.panel)
    
    pcoa.shift.mt <- as.matrix(pcoa.dist)
    
    pcoa.shift.vc <- sapply(unique(tb.metad$donor_mouse_n), function(x) {
      
      pre  <- tb.metad$sample_id[tb.metad$donor_mouse_n == x & tb.metad$timepoint == "pre-HFD"]
      post <- tb.metad$sample_id[tb.metad$donor_mouse_n == x & tb.metad$timepoint == "post-HFD"]
      
      pcoa.shift.mt[pre, post]
      
    }) %>% 
      unlist()
    
    pcoa.shift.dt <- 
      data.frame("donor_mouse_n" = names(pcoa.shift.vc), "value" = pcoa.shift.vc) %>% 
      merge(tb.metad %>% select(donor_mouse_n, donor_ethnicity) %>% distinct(), .) %>% 
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
                          value_var = "value",
                          only.sig = T,
                          alpha = 1, 
                          label.type = "p.adj")
    
    gg.pcoa.panel <- 
      gg.pcoa.pc1.compm.axes12 + plot_spacer() + plot_spacer() + gg.pcoa.pc1.compm.axes13 + plot_spacer() + plot_spacer() + plot_spacer() +
      gg.pcoa.axes12 + gg.pcoa.pc2.compm.axes12 + plot_spacer() + gg.pcoa.axes13 + gg.pcoa.pc2.compm.axes13 + gg.pcoa.disp.compm.comp2 + gg.pcoa.shift.compm.comp2 + plot_layout(ncol = 7, heights = c(.1, 1), widths = c(1, .1, .05, 1, .1, .5, .5))
    
    print(gg.pcoa.panel)
    
    #Saving graphs
    png(filename = paste0(path.fg.dir, "/diversity/beta-div-", str_to_lower(bdiv.indexes[[i]]), ".png"), width = 5000 * 2.75, height = 5000 * 1, units = "px", res = 1200, type = "cairo")
    
    print(gg.pcoa.panel)
    
    dev.off()
    
    ls.pcoa <-
      
      list(
        "PERMANOVA" = pcoa.perma.df.comp1,
        
        "PERMANOVA (post-hoc)" = pcoa.perma.ph,
        
        "PC1 (Summary)" = pcoa.site.score.posthoc.axes12 %>% filter(PC == "PCoA1"),
        "PC2 (Summary)" = pcoa.site.score.posthoc.axes12 %>% filter(PC == "PCoA2"),
        "PC3 (Summary)" = pcoa.site.score.posthoc.axes13 %>% filter(PC == "PCoA3"),
        
        "Dispersion (Measures)" = pcoa.disp.comp2,
        "Dispersion (Summary)" = pcoa.disp.summary.comp2,
        "Dispersion (Wilcox)" = pcoa.disp.compm.comp2,
        
        "pre-post dist (Measures)" = pcoa.shift.dt,
        "pre-post dist (Summary)" = pcoa.shift.summary.comp2,
        "pre-post dist (Wilcox)" = pcoa.shift.compm.comp2
        
      )
    
    openxlsx::write.xlsx(x = ls.pcoa, file = paste0(path.tb.dir, "/diversity/beta-div-", str_to_lower(bdiv.indexes[[i]]), ".xlsx"), rowNames = F)
    
  }; rm(i)
  
}

bdiv.indexes <- c("bray", "jaccard")
names(bdiv.indexes) <- c("Bray—Curtis", "Jaccard")
fun.get.beta.div(ps.table)

bdiv.indexes <- c("euclidean")
names(bdiv.indexes) <- c("Euclidean")
fun.get.beta.div(microbiome::transform(ps.table, "clr"))

rm(list = ls(pattern = "[.]?beta[.]?"))
rm(list = ls(pattern = "[.]?bdiv[.]?"))
rm(list = ls(pattern = "[.]?pcoa[.]?"))


# Profiles ---------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/profiles"))) { dir.create(str_glue("{path.fg.dir}/profiles"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/profiles"))) { dir.create(str_glue("{path.tb.dir}/profiles"), recursive = T) }

tb.work.fun.ko <-
  tb.work %>%
  select(ortholog_id, any_of(tb.metad$sample_id)) %>%
  group_by(ortholog_id) %>%
  summarise(across(where(is.numeric), ~sum(.x, na.rm = T)), .groups = "drop") %>%
  merge(tb.ann, ., all.y = T) %>%
  separate_longer_delim(., ko_id, delim = ",") %>%
  merge(., tb.aux.ann.ko.def, all.x = T) %>%
  select(ortholog_id, starts_with("ko_"), any_of(tb.metad$sample_id)) %>%
  filter(!is.na(ko_class))

openxlsx::write.xlsx(x = tb.work.fun.ko, file = paste0(path.tb.dir, "/profiles/tb-work-ko.xlsx"), rowNames = F)

names(tb.work.fun.ko)


## KEGG class -------------------------------------------------------------------


# Sample

tb.ann.ko.prof.1 <-
  tb.work.fun.ko %>%
  select(ortholog_id, ko_class, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_class) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_class) %>%
  data.frame(check.names = F)

tb.ann.ko.prof.2 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class %in% levels(tb.aux.ann.ko$ko_class)) %>%
  mutate(ko_class = factor(ko_class, levels(tb.aux.ann.ko$ko_class))) %>%
  pivot_longer(cols = -ko_class, names_to = "sample_id", values_to = "cpm") %>%
  merge(tb.metad, .)

tb.ann.ko.prof.3 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class %in% levels(tb.aux.ann.ko$ko_class)) %>%
  mutate(ko_class = factor(ko_class, levels(tb.aux.ann.ko$ko_class))) %>%
  mutate(across(.cols = tb.metad$sample_id, .fns = ~./sum(.))) %>%
  pivot_longer(cols = -ko_class, names_to = "sample_id", values_to = "perc") %>%
  merge(tb.metad, .)

tb.ann.ko.prof.4 <- merge(tb.ann.ko.prof.2, tb.ann.ko.prof.3)

ggplot(tb.ann.ko.prof.4, aes(x = mouse_n, y = cpm, fill = ko_class)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity + timepoint, space = "free_x", scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0), labels = fun.get.scientific) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko) +
  xlab("") + ylab("CPM mapped reads<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)) -> gg.ko.barprof.abs; gg.ko.barprof.abs

ggplot(tb.ann.ko.prof.4, aes(x = mouse_n, y = perc, fill = ko_class)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity + timepoint, space = "free_x", scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko) +
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)) -> gg.ko.barprof.rel; gg.ko.barprof.rel

gg.ko.barprof.abs + gg.leg.ko + plot_layout(ncol = 2, widths = c(.3, .7)) -> gg.ko.barprof.a; gg.ko.barprof.a
gg.ko.barprof.rel + gg.leg.ko + plot_layout(ncol = 2, widths = c(.3, .7)) -> gg.ko.barprof.b; gg.ko.barprof.b

(gg.ko.barprof.a / gg.ko.barprof.b)

png(filename = paste0(path.fg.dir, "/profiles/bp-ko-class_sample.png"), width = 5000 * 3.5, height = 5000 * 1.5, units = "px", res = 1200, type = "cairo")

print((gg.ko.barprof.a / gg.ko.barprof.b))

dev.off()

rm(list = ls(pattern = "[.]?prof[.]?"))
rm(list = ls(pattern = "[.]?barprof[.]?"))

# Group

tb.ann.ko.prof.1 <-
  tb.work.fun.ko %>%
  select(ortholog_id, ko_class, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_class) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_class) %>%
  data.frame(check.names = F)

tb.ann.ko.prof.2 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class %in% levels(tb.aux.ann.ko$ko_class)) %>%
  mutate(ko_class = factor(ko_class, levels(tb.aux.ann.ko$ko_class))) %>%
  pivot_longer(cols = -ko_class, names_to = "sample_id", values_to = "cpm") %>%
  merge(tb.metad, .) %>%
  group_by(donor_ethnicity_timepoint, ko_class) %>%
  summarise(cpm = mean(cpm, na.rm = T), .groups = "drop")

tb.ann.ko.prof.3 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class %in% levels(tb.aux.ann.ko$ko_class)) %>%
  mutate(ko_class = factor(ko_class, levels(tb.aux.ann.ko$ko_class))) %>%
  pivot_longer(cols = -ko_class, names_to = "sample_id", values_to = "perc") %>%
  merge(tb.metad, .) %>%
  group_by(sample_id, donor_ethnicity_timepoint, ko_class) %>%
  summarise(perc = sum(perc, na.rm = T), .groups = "drop") %>%
  group_by(donor_ethnicity_timepoint, ko_class) %>%
  summarise(perc = mean(perc, na.rm = T), .groups = "drop") %>%
  group_by(donor_ethnicity_timepoint) %>%
  mutate(perc = perc/sum(perc, na.rm = T))

tb.ann.ko.prof.4 <- merge(tb.ann.ko.prof.2, tb.ann.ko.prof.3)

ggplot(tb.ann.ko.prof.4 %>%
         mutate(d_ethnicity = paste0(str_sub(donor_ethnicity_timepoint, 1, 1))) %>%
         mutate(t_point = ifelse(str_detect(donor_ethnicity_timepoint, "pre-HFD"), "T1", "T2")),
       aes(fill = ko_class, x = donor_ethnicity_timepoint, y = cpm)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
  facet_grid(. ~ d_ethnicity + t_point, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0), labels = fun.get.scientific) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko) +
  xlab("") + ylab("CPM mapped reads<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> gg.ko.barprof.abs; gg.ko.barprof.abs

ggplot(tb.ann.ko.prof.4 %>%
         mutate(d_ethnicity = paste0(str_sub(donor_ethnicity_timepoint, 1, 1))) %>%
         mutate(t_point = ifelse(str_detect(donor_ethnicity_timepoint, "pre-HFD"), "T1", "T2")),
       aes(fill = ko_class, x = donor_ethnicity_timepoint, y = perc)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
  facet_grid(. ~ d_ethnicity + t_point, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko) +
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> gg.ko.barprof.rel; gg.ko.barprof.rel

gg.ko.barprof.abs + gg.leg.ko + plot_layout(ncol = 2, widths = c(.1, .9)) -> gg.ko.barprof.a; gg.ko.barprof.a
gg.ko.barprof.rel + gg.leg.ko + plot_layout(ncol = 2, widths = c(.1, .9)) -> gg.ko.barprof.b; gg.ko.barprof.b

(gg.ko.barprof.a / gg.ko.barprof.b)

png(filename = paste0(path.fg.dir, "/profiles/bp-ko-class_group.png"), width = 5000 * 3, height = 5000 * 1.3, units = "px", res = 1200, type = "cairo")

print((gg.ko.barprof.a / gg.ko.barprof.b))

dev.off()

rm(list = ls(pattern = "[.]?prof[.]?"))
rm(list = ls(pattern = "[.]?barprof[.]?"))


## KEGG metab -------------------------------------------------------------------


# Sample

tb.ann.ko.prof.1 <-
  tb.work.fun.ko %>%
  select(ortholog_id, ko_class, ko_subclass, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_class, ko_subclass) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_class, ko_subclass) %>%
  data.frame(check.names = F)

tb.ann.ko.prof.2 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class == "Metabolism") %>%
  mutate(ko_subclass = factor(ko_subclass, unique(.$ko_subclass))) %>%
  pivot_longer(cols = -c(ko_class, ko_subclass), names_to = "sample_id", values_to = "cpm") %>%
  merge(tb.metad, .)

tb.ann.ko.prof.3 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class == "Metabolism") %>%
  mutate(ko_subclass = factor(ko_subclass, unique(.$ko_subclass))) %>%
  mutate(across(.cols = tb.metad$sample_id, .fns = ~./sum(.))) %>%
  pivot_longer(cols = -c(ko_class, ko_subclass), names_to = "sample_id", values_to = "perc") %>%
  merge(tb.metad, .)

tb.ann.ko.prof.4 <- merge(tb.ann.ko.prof.2, tb.ann.ko.prof.3)

ggplot(tb.ann.ko.prof.4, aes(x = mouse_n, y = cpm, fill = ko_subclass)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity + timepoint, space = "free_x", scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0), labels = fun.get.scientific) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko.met) +
  xlab("") + ylab("CPM mapped reads<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)) -> gg.ko.barprof.abs

ggplot(tb.ann.ko.prof.4, aes(x = mouse_n, y = perc, fill = ko_subclass)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'sum', width = .95) +
  facet_grid(. ~ donor_ethnicity + timepoint, space = "free_x", scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko.met) +
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)) -> gg.ko.barprof.rel

gg.ko.barprof.abs + gg.leg.ko.met + plot_layout(ncol = 2, widths = c(.3, .7)) -> gg.ko.barprof.a; gg.ko.barprof.a
gg.ko.barprof.rel + gg.leg.ko.met + plot_layout(ncol = 2, widths = c(.3, .7)) -> gg.ko.barprof.b; gg.ko.barprof.b

(gg.ko.barprof.a / gg.ko.barprof.b)

png(filename = paste0(path.fg.dir, "/profiles/bp-ko-class-met_sample.png"), width = 5000 * 3.5, height = 5000 * 1.5, units = "px", res = 1200, type = "cairo")

print((gg.ko.barprof.a / gg.ko.barprof.b))

dev.off()

rm(list = ls(pattern = "[.]?prof[.]?"))
rm(list = ls(pattern = "[.]?barprof[.]?"))

# Group

tb.ann.ko.prof.1 <-
  tb.work.fun.ko %>%
  select(ortholog_id, ko_class, ko_subclass, tb.metad$sample_id) %>%
  distinct() %>%
  group_by(ko_class, ko_subclass) %>%
  summarise(across(.cols = tb.metad$sample_id, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(ko_class, ko_subclass) %>%
  data.frame(check.names = F)

tb.ann.ko.prof.2 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class == "Metabolism") %>%
  mutate(ko_subclass = factor(ko_subclass, unique(.$ko_subclass))) %>%
  pivot_longer(cols = -c(ko_class, ko_subclass), names_to = "sample_id", values_to = "cpm") %>%
  merge(tb.metad, .) %>%
  group_by(donor_ethnicity_timepoint, ko_class, ko_subclass) %>%
  summarise(cpm = mean(cpm, na.rm = T), .groups = "drop")

tb.ann.ko.prof.3 <-
  tb.ann.ko.prof.1 %>%
  filter(ko_class == "Metabolism") %>%
  mutate(ko_subclass = factor(ko_subclass, unique(.$ko_subclass))) %>%
  pivot_longer(cols = -c(ko_class, ko_subclass), names_to = "sample_id", values_to = "perc") %>%
  merge(tb.metad, .) %>%
  group_by(sample_id, donor_ethnicity_timepoint, ko_class, ko_subclass) %>%
  summarise(perc = sum(perc, na.rm = T), .groups = "drop") %>%
  group_by(donor_ethnicity_timepoint, ko_class, ko_subclass) %>%
  summarise(perc = mean(perc, na.rm = T), .groups = "drop") %>%
  group_by(donor_ethnicity_timepoint) %>%
  mutate(perc = perc/sum(perc, na.rm = T))

tb.ann.ko.prof.4 <- merge(tb.ann.ko.prof.2, tb.ann.ko.prof.3)

ggplot(tb.ann.ko.prof.4 %>%
         mutate(d_ethnicity = paste0(str_sub(donor_ethnicity_timepoint, 1, 1))) %>%
         mutate(t_point = ifelse(str_detect(donor_ethnicity_timepoint, "pre-HFD"), "T1", "T2")),
       aes(fill = ko_subclass, x = donor_ethnicity_timepoint, y = cpm)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
  facet_grid(. ~ d_ethnicity + t_point, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0), labels = fun.get.scientific) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko.met) +
  xlab("") + ylab("CPM mapped reads<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> gg.ko.barprof.abs

ggplot(tb.ann.ko.prof.4 %>%
         mutate(d_ethnicity = paste0(str_sub(donor_ethnicity_timepoint, 1, 1))) %>%
         mutate(t_point = ifelse(str_detect(donor_ethnicity_timepoint, "pre-HFD"), "T1", "T2")),
       aes(fill = ko_subclass, x = donor_ethnicity_timepoint, y = perc)) +
  geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
  facet_grid(. ~ d_ethnicity + t_point, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = colors.ko.met) +
  xlab("") + ylab("Relative abundance (%)<br>") +
  my.theme.1 +
  theme(panel.spacing = unit(.1, "lines"),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> gg.ko.barprof.rel

gg.ko.barprof.abs + gg.leg.ko.met + plot_layout(ncol = 2, widths = c(.1, .9)) -> gg.ko.barprof.a; gg.ko.barprof.a
gg.ko.barprof.rel + gg.leg.ko.met + plot_layout(ncol = 2, widths = c(.1, .9)) -> gg.ko.barprof.b; gg.ko.barprof.b

(gg.ko.barprof.a / gg.ko.barprof.b)


png(filename = paste0(path.fg.dir, "/profiles/bp-ko-class-met_group.png"), width = 5000 * 3, height = 5000 * 1.3, units = "px", res = 1200, type = "cairo")

print((gg.ko.barprof.a / gg.ko.barprof.b))

dev.off()

rm(list = ls(pattern = "[.]?prof[.]?"))
rm(list = ls(pattern = "[.]?barprof[.]?"))


# Diff Abund -------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diff_abund"))) { dir.create(str_glue("{path.fg.dir}/diff_abund"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diff_abund"))) { dir.create(str_glue("{path.tb.dir}/diff_abund"), recursive = T) }


## KOs -------------------------------------------------------------------------


ranks <- c("ko_id")

comps <- c("Belgian: pre-HFD|Yanomami: pre-HFD",
           "Belgian: post-HFD|Yanomami: post-HFD",
           "Belgian: pre-HFD|Belgian: post-HFD",
           "Yanomami: pre-HFD|Yanomami: post-HFD")

ls.tb <- list("ko_id" = tb.ko.counts)
ls.msl <- list()

for (j in seq_along(comps)) {
  
  for (i in seq_along(ranks)) {
    
    tmp.comp.factor <- "donor_ethnicity_timepoint"
    
    tmp.comp <- comps[[j]]
    
    tmp.meta <- tb.metad %>% filter(!!as.symbol(tmp.comp.factor) %in% unlist(str_split(tmp.comp, "\\|"))) %>% column_to_rownames("sample_id")
    
    tmp.feature.table <-
      ls.tb[[i]] %>%
      mutate("prev" = rowSums(.[tmp.meta %>% rownames()] != 0)) %>%
      filter(prev >= (tmp.meta %>% rownames() %>% length())*.1) %>%
      select(-prev) %>%
      mutate(tmpID = paste0("tmpID.", 1:n()))
    
    tmp.rel.feature.id <- tmp.feature.table %>% select(tmpID, any_of(ranks[[i]]))
    
    tmp.feature.table <- tmp.feature.table %>% column_to_rownames("tmpID") %>% select(any_of(rownames(tmp.meta)))
    
    tmp.msl.comp <-
      Maaslin2(input_data = tmp.feature.table,
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
      merge(tmp.rel.feature.id, tmp.msl.comp$results, by.x = "tmpID", by.y = "feature") %>%
      mutate(enrich_group = ifelse(coef < 0, str_remove({tmp.comp}, '\\|.*'), str_remove({tmp.comp}, '.*\\|'))) %>%
      mutate(comp = tmp.comp,
             rank = ranks[[i]]) %>%
      select(rank, feature = ranks[[i]], enrich_group, coef, stderr, pval, qval, N, N.not.zero, comp)
    
  }
  
}

rm(i, j)
rm(comps)
rm(list = ls(pattern = "[.]?tmp[.]?"))

df.msl.raw <-
  bind_rows(ls.msl) %>%
  mutate(comp = factor(comp, levels = unique(comp))) %>%
  mutate(enrich_group = factor(enrich_group, levels = unique(c(levels(tb.metad$donor_ethnicity_timepoint))))) %>%
  mutate(rank = factor(rank, levels = ranks)) %>%
  select(rank, feature, comp, enrich_group, coef:N.not.zero) %>%
  arrange(rank, comp, enrich_group, -abs(coef))

df.msl.filt <-
  df.msl.raw %>%
  filter(qval < 0.05 & abs(coef) > 1)

df.msl <- list("Raw" = df.msl.raw, "Filtered" = df.msl.filt)

openxlsx::write.xlsx(x = df.msl, file = paste0(path.tb.dir, "/diff_abund/maaslin.xlsx"), rowNames = F)


# Volcano Plot

my.daa.tb <-
  df.msl.raw %>%
  mutate(sig = case_when(qval < 0.05 & abs(coef) > 1 ~ enrich_group, .default = "Not Sig.")) %>%
  mutate(sig = factor(sig, c(names(colors.gp.tp), "Not Sig.")))

tb.daa.labels <-
  my.daa.tb %>%
  filter(sig != "Not Sig.") %>%
  select(comp, sig) %>%
  table() %>%
  data.frame() %>%
  filter(Freq != 0) %>%
  mutate(
    my.x = case_when(
      str_remove(sig, "\\|.*") == str_remove_all(comp, "\\|.*") ~ -3.162278,
      str_remove(sig, "\\|.*") == str_remove_all(comp, ".*\\|") ~ 3.162278),
    my.hjust = case_when(
      str_remove(sig, "\\|.*") == str_remove_all(comp, "\\|.*") ~ 0.5,
      str_remove(sig, "\\|.*") == str_remove_all(comp, ".*\\|") ~ 0.5)
  ) %>%
  mutate(sig = factor(sig, levels = c(names(colors.gp.tp))),
         comp = factor(comp, levels = levels(my.daa.tb$comp))) %>%
  arrange(comp, sig); head(tb.daa.labels)

my.daa.tb %>%
  filter(!is.na(coef)) %>%
  ggplot(., aes(x = ifelse(coef > 0, sqrt(coef), -sqrt(abs(coef))), y = -log10(qval), color = sig)) +
  geom_point(shape = 16, alpha = 1, size = 1, fill = "white") +
  geom_hline(yintercept = -log10(0.05), color = c("black"), linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(-1, 1), color = c("black"), linetype = "dashed", linewidth = .25) +
  geom_richtext(mapping = aes(label = str_glue("**{Freq}**"), x = my.x, y = Inf, hjust = my.hjust),
                data = tb.daa.labels, size = 3.5, vjust = 2.5, color = "black",
                inherit.aes = T, label.padding = unit(2, "pt"), label.r = unit(0, "pt"), fill = "white", label.color = "gray80", show.legend = F) +
  facet_grid(. ~ comp, axes = "all_y") +
  ylab("-log<sub>10</sub> (adjusted <i>p</i>-value)<br>") +
  xlab("<br>|MaAsLin (β) coefficient|") +
  scale_x_continuous(limits = ~ c(-1.1, 1.1) * max(abs(.x)), labels = function(x) round(x^2, 0), breaks = c(-5.477226, -3.162278, -1, 0, 1, 3.162278, 5.477226)) +
  scale_y_continuous(expand = expansion(mult = c(0.001, 0.075)), trans = "sqrt", breaks = c(0, 2, 6, 10, 14, 18)^2) +
  scale_color_manual(name = "Group", values = c(colors.gp.tp, "Not Sig." = "gray75")) +
  my.theme.2 +
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        axis.title.x = element_markdown(size = 10, face = 'bold', colour = 'black'),
        axis.title.y = element_markdown(size = 10, face = 'bold', colour = 'black'),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), 'cm')) -> gg.daa.vp.plot.1; gg.daa.vp.plot.1

ggplot(tb.daa.labels, aes(x = sig, y = 1)) +
  geom_tile(fill = "white") +
  geom_text(aes(label = sig), size = 3, fontface = "bold", check_overlap = F) +
  facet_grid(. ~ comp, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  scale_x_discrete(expand = expansion(mult = 0.1, add = 0)) +
  theme_void(base_family = "Liberation Sans") +
  theme(strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.01, 0.01, 0.1, 0.01), 'cm')) -> gg.daa.vp.labels.1; gg.daa.vp.labels.1


gg.daa.vp.labels.1 / gg.daa.vp.plot.1 + plot_layout(heights = c(0.15, 1)) -> gg.daa.vp.1; gg.daa.vp.1

png(filename = paste0(path.fg.dir, "/diff_abund/daa-volcano.png"), width = 5000 * 3, height = 5000 * 0.75, units = "px", res = 1200, type = "cairo")

print(gg.daa.vp.1)

dev.off()

rm(list = ls(pattern = "[.]?daa[.]?"))
rm(list = ls(pattern = "[.]?vp[.]?"))

# Divergent plot: class

tb.daa.dp <-
  df.msl.raw %>%
  mutate(sig = case_when(qval < 0.05 & abs(coef) > 1 ~ enrich_group, .default = "Not Sig.")) %>%
  mutate(sig = factor(sig, c(names(colors.gp.tp), "Not Sig."))) %>%
  select(-rank) %>%
  rename(fun_id = feature)

tb.daa.ann <- tb.aux.ann.ko.def

tb.daa.dp.ann <-
  tb.daa.dp %>%
  merge(.,
        tb.daa.ann %>%
          select(fun_id = 1, fun_class = 2, fun_subclass = 3) %>%
          filter(!is.na(fun_subclass)) %>%
          distinct(), by = "fun_id", all.x = T)

tb.daa.dp.ann %>%
  select(sig, comp, fun_class, fun_subclass) %>%
  filter(!is.na(fun_subclass)) %>%
  group_by(sig, comp, fun_class, fun_subclass) %>%
  summarise(freq = n(), .groups = "drop") %>%
  data.frame(.) %>%
  arrange(comp, fun_class, fun_subclass) -> tb.daa.dp.freq; tb.daa.dp.freq

for (i in seq_along(levels(tb.daa.dp.freq$comp))) {
  
  tb.daa.dp.freq %>%
    filter(as.numeric(comp) == i) %>%
    pivot_wider(., names_from = sig, values_from = freq) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    pivot_longer(cols = where(is.numeric), names_to = "sig", values_to = "freq") %>%
    mutate(my.x = case_when(
      sig == "Not Sig." ~ 0,
      sig == str_remove(comp, "\\|.*") ~ -freq,
      sig == str_remove(comp, ".*\\|") ~ freq)) %>%
    mutate(sig2 = str_replace(sig, ": ", "<br>")) -> tb.daa.dp.freq.2; tb.daa.dp.freq.2
  
  tmp.threshold <- max(tb.daa.dp.freq.2 %>% filter(sig != "Not Sig.") %>% pull(freq)) * .33
  
  tb.daa.dp.freq.2 %>%
    filter(sig == str_remove(comp, "\\|.*")) %>%
    mutate(my.label.pos = ifelse(my.x > -tmp.threshold, 2.25, -0.75)) %>%
    mutate(my.label.col = ifelse(my.x > -tmp.threshold, "black", "white")) %>%
    ggplot(., aes(x = my.x, y = fun_subclass, fill = sig)) +
    geom_col() +
    facet_grid(fun_class ~ sig2,  scales = "free", space = "free_y") +
    geom_text(aes(hjust = my.label.pos, label = freq, color = my.label.col), size = 2) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(limits = c(-max(tmp.threshold, (tb.daa.dp.freq.2 %>%
                                                         filter(sig != "Not Sig.") %>%
                                                         pull(freq) %>%
                                                         max(.))), 0),
                       expand = expansion(mult = c(0.035, 0)),
                       labels = function(x) abs(x)) +
    scale_fill_manual(name = "", values = c(colors.gp, colors.gp.tp)) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "gray85", fill = NA),
          panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "gray75"),
          panel.grid.major.y = element_line(size = 0.25, linetype = 'dotted', colour = "gray75"),
          strip.text.x = element_markdown(size = 9, colour = 'black', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.ticks.length = unit(2.5, "pt"),
          axis.text.y = element_markdown(size = 6.5, colour = 'black'),
          plot.margin = unit(c(0, -.5, 0, 0), 'cm')) -> gg.daa.dp.a
  
  tb.daa.dp.freq.2 %>%
    filter(str_detect(sig, "Not Sig.")) %>%
    ggplot(., aes(x = 1, y = fun_subclass)) +
    geom_col(fill = "gray") +
    facet_grid(fun_class ~ sig2,  scales = "free", space = "free_y") +
    geom_text(aes(label = freq, x = .5), size = 1.75) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "gray85", fill = NA),
          axis.line.y = element_blank(),
          strip.text = element_markdown(size = 9, colour = '#00000000', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_markdown(size = 6, colour = '#00000000'),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(color = "#00000000"),
          axis.ticks.length = unit(2.5, "pt"),
          plot.margin = unit(c(0, -.5, 0, -.5), 'cm')) -> gg.daa.dp.b
  
  tb.daa.dp.freq.2 %>%
    filter(sig == str_remove(comp, ".*\\|")) %>%
    mutate(my.label.pos = ifelse(my.x < tmp.threshold, -1.25, 1.5)) %>%
    mutate(my.label.col = ifelse(my.x < tmp.threshold, "black", "white")) %>%
    ggplot(., aes(x = my.x, y = fun_subclass, fill = sig)) +
    geom_col() +
    facet_grid(fun_class ~ sig2,  scales = "free", space = "free_y") +
    geom_text(aes(hjust = my.label.pos, label = freq, color = my.label.col), size = 2) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(limits = c(0, max(tmp.threshold, (tb.daa.dp.freq.2 %>%
                                                           filter(sig != "Not Sig.") %>%
                                                           pull(freq) %>%
                                                           max(.)))), expand = expansion(mult = c(0, 0.035)), labels = function(x) abs(x)) +
    scale_fill_manual(name = "", values = c(colors.gp, colors.gp.tp)) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "gray85", fill = NA),
          panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "gray75"),
          panel.grid.major.y = element_line(size = 0.25, linetype = 'dotted', colour = "gray75"),
          axis.line.y = element_blank(),
          strip.text = element_markdown(size = 9, colour = 'black', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(2.5, "pt"),
          plot.margin = unit(c(0, 0, 0, -.5), 'cm')) -> gg.daa.dp.c
  
  tb.daa.dp.freq.2 %>%
    filter(str_detect(sig, "Not Sig.")) %>%
    ggplot(., aes(x = 1, y = fun_subclass, fill = fun_class)) +
    geom_col() +
    facet_grid(fun_class ~ sig2,  scales = "free", space = "free_y") +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(name = "", values = c(colors.ko)) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "#00000000", fill = NA),
          axis.line.y = element_blank(),
          axis.line.x = element_line(color = "#00000000"),
          strip.text = element_markdown(size = 9, colour = '#00000000', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_markdown(size = 6, colour = '#00000000'),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(color = "#00000000"),
          axis.ticks.length = unit(2.5, "pt"),
          plot.margin = unit(c(0, .5, 0, -.5), 'cm')) -> gg.daa.dp.d
  
  #print(gg.daa.dp.a + gg.daa.dp.b + gg.daa.dp.c + gg.daa.dp.d + gg.leg.ko + plot_layout(nrow = 1, widths = c(.75, .25, .75, .05, 2.5), axis_titles = "collect_x"))
  
  png(filename = paste0(path.fg.dir, "/diff_abund/daa-divplot-ko-class_comp", i, ".png"), width = 5000 * 2, height = 5000 * 1, units = "px", res = 1200, type = "cairo")
  
  print(gg.daa.dp.a + gg.daa.dp.b + gg.daa.dp.c + gg.daa.dp.d + gg.leg.ko + plot_layout(nrow = 1, widths = c(.75, .25, .75, .05, 2.5), axis_titles = "collect_x"))
  
  dev.off()
  
}; rm(i)

rm(list = ls(pattern = "[.]?daa[.]?"))
rm(list = ls(pattern = "[.]?dp[.]?"))
rm(list = ls(pattern = "[.]?tmp[.]?"))

# Divergent plot: metabolism

tb.daa.dp <-
  df.msl.raw %>%
  mutate(sig = case_when(qval < 0.05 & abs(coef) > 1 ~ enrich_group, .default = "Not Sig.")) %>%
  mutate(sig = factor(sig, c(names(colors.gp.tp), "Not Sig."))) %>%
  select(-rank) %>%
  rename(fun_id = feature)

tb.daa.ann <- tb.aux.ann.ko.def

tb.daa.dp.ann <-
  tb.daa.dp %>%
  merge(.,
        tb.daa.ann %>%
          select(fun_id = 1, fun_class = 2, fun_subclass = 3, fun_pathway = 4) %>%
          filter(!is.na(fun_pathway)) %>%
          distinct(), by = "fun_id", all.x = T)

tb.daa.dp.ann %>%
  select(sig, comp, fun_class, fun_subclass, fun_pathway) %>%
  filter(fun_class == "Metabolism") %>%
  group_by(sig, comp, fun_class, fun_subclass, fun_pathway) %>%
  summarise(freq = n(), .groups = "drop") %>%
  data.frame(.) %>%
  arrange(comp, fun_class, fun_subclass, fun_pathway) -> tb.daa.dp.freq

tb.daa.dp.freq.top <-
  tb.daa.dp.freq %>%
  select(fun_class:freq) %>%
  group_by(across(-freq)) %>%
  summarise(across(.cols = freq, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  group_by(fun_subclass) %>%
  mutate(threshold = quantile(freq, .75)) %>%
  ungroup() %>%
  filter(freq >= threshold) %>% pull(fun_pathway) %>% unique()

for (i in seq_along(levels(tb.daa.dp.freq$comp))) {
  
  tb.daa.dp.freq %>%
    filter(as.numeric(comp) == i) %>%
    mutate(fun_pathway = ifelse(fun_pathway %in% tb.daa.dp.freq.top, fun_pathway, str_glue("Other: {fun_subclass}"))) %>%
    group_by(across(sig:fun_pathway)) %>%
    summarise(across(.cols = freq, ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(., names_from = sig, values_from = freq) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    pivot_longer(cols = where(is.numeric), names_to = "sig", values_to = "freq") %>%
    mutate(my.x = case_when(
      sig == "Not Sig." ~ 0,
      sig == str_remove(comp, "\\|.*") ~ -freq,
      sig == str_remove(comp, ".*\\|") ~ freq)) %>%
    mutate(sig2 = str_replace(sig, ": ", "<br>")) -> tb.daa.dp.freq.2; tb.daa.dp.freq.2
  
  tmp.threshold <- max(tb.daa.dp.freq.2 %>% filter(sig != "Not Sig.") %>% pull(freq)) * .33
  
  tb.daa.dp.freq.2 %>%
    filter(sig == str_remove(comp, "\\|.*")) %>%
    mutate(my.label.pos = ifelse(my.x > -tmp.threshold, 2.25, -0.75)) %>%
    mutate(my.label.col = ifelse(my.x > -tmp.threshold, "black", "white")) %>%
    ggplot(., aes(x = my.x, y = fun_pathway, fill = sig)) +
    geom_col() +
    facet_grid(fun_subclass ~ sig2,  scales = "free", space = "free_y") +
    geom_text(aes(hjust = my.label.pos, label = freq, color = my.label.col), size = 2) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(limits = c(-max(tmp.threshold, (tb.daa.dp.freq.2 %>%
                                                         filter(sig != "Not Sig.") %>%
                                                         pull(freq) %>%
                                                         max(.))), 0),
                       expand = expansion(mult = c(0.035, 0)),
                       labels = function(x) abs(x)) +
    scale_fill_manual(name = "", values = c(colors.gp, colors.gp.tp)) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "gray85", fill = NA),
          panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "gray75"),
          panel.grid.major.y = element_line(size = 0.25, linetype = 'dotted', colour = "gray75"),
          strip.text.x = element_markdown(size = 9, colour = 'black', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.ticks.length = unit(2.5, "pt"),
          axis.text.y = element_markdown(size = 6.5, colour = 'black'),
          plot.margin = unit(c(0, -.5, 0, 0), 'cm')) -> gg.daa.dp.a; gg.daa.dp.a
  
  tb.daa.dp.freq.2 %>%
    filter(str_detect(sig, "Not Sig.")) %>%
    ggplot(., aes(x = 1, y = fun_pathway)) +
    geom_col(fill = "gray") +
    facet_grid(fun_subclass ~ sig2,  scales = "free", space = "free_y") +
    geom_text(aes(label = freq, x = .5), size = 1.75) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "gray85", fill = NA),
          axis.line.y = element_blank(),
          strip.text = element_markdown(size = 9, colour = '#00000000', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_markdown(size = 6, colour = '#00000000'),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(color = "#00000000"),
          axis.ticks.length = unit(2.5, "pt"),
          plot.margin = unit(c(0, -.5, 0, -.5), 'cm')) -> gg.daa.dp.b; gg.daa.dp.b
  
  tb.daa.dp.freq.2 %>%
    filter(sig == str_remove(comp, ".*\\|")) %>%
    mutate(my.label.pos = ifelse(my.x < tmp.threshold, -1.25, 1.5)) %>%
    mutate(my.label.col = ifelse(my.x < tmp.threshold, "black", "white")) %>%
    ggplot(., aes(x = my.x, y = fun_pathway, fill = sig)) +
    geom_col() +
    facet_grid(fun_subclass ~ sig2,  scales = "free", space = "free_y") +
    geom_text(aes(hjust = my.label.pos, label = freq, color = my.label.col), size = 2) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(limits = c(0, max(tmp.threshold, (tb.daa.dp.freq.2 %>%
                                                           filter(sig != "Not Sig.") %>%
                                                           pull(freq) %>%
                                                           max(.)))), expand = expansion(mult = c(0, 0.035)), labels = function(x) abs(x)) +
    scale_fill_manual(name = "", values = c(colors.gp, colors.gp.tp)) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "gray85", fill = NA),
          panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "gray75"),
          panel.grid.major.y = element_line(size = 0.25, linetype = 'dotted', colour = "gray75"),
          axis.line.y = element_blank(),
          strip.text = element_markdown(size = 9, colour = 'black', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(2.5, "pt"),
          plot.margin = unit(c(0, 0, 0, -.5), 'cm')) -> gg.daa.dp.c; gg.daa.dp.c
  
  tb.daa.dp.freq.2 %>%
    filter(str_detect(sig, "Not Sig.")) %>%
    ggplot(., aes(x = 1, y = fun_pathway, fill = fun_subclass)) +
    geom_col() +
    facet_grid(fun_subclass ~ sig2,  scales = "free", space = "free_y") +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(name = "", values = c(colors.ko.met)) +
    scale_color_identity(name = "") +
    xlab("") + ylab("") +
    my.theme.2 +
    theme(legend.position = "none",
          panel.spacing.x = unit(0, "mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.border = element_rect(colour = "#00000000", fill = NA),
          axis.line.y = element_blank(),
          axis.line.x = element_line(color = "#00000000"),
          strip.text = element_markdown(size = 9, colour = '#00000000', face = "bold", lineheight = 1.25),
          strip.text.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_markdown(size = 6, colour = '#00000000'),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(color = "#00000000"),
          axis.ticks.length = unit(2.5, "pt"),
          plot.margin = unit(c(0, .5, 0, -.5), 'cm')) -> gg.daa.dp.d; gg.daa.dp.d
  
  #print(gg.daa.dp.a + gg.daa.dp.b + gg.daa.dp.c + gg.daa.dp.d + gg.leg.ko.met + plot_layout(nrow = 1, widths = c(.75, .25, .75, .05, 2.5), axis_titles = "collect_x"))
  
  png(filename = paste0(path.fg.dir, "/diff_abund/daa-divplot-ko-class-met_comp", i, ".png"), width = 5000 * 2, height = 5000 * 2, units = "px", res = 1200, type = "cairo")
  
  print(gg.daa.dp.a + gg.daa.dp.b + gg.daa.dp.c + gg.daa.dp.d + gg.leg.ko.met + plot_layout(nrow = 1, widths = c(.75, .25, .75, .05, 2.5), axis_titles = "collect_x"))
  
  dev.off()
  
}; rm(i)

rm(list = ls(pattern = "[.]?daa[.]?"))
rm(list = ls(pattern = "[.]?dp[.]?"))
rm(list = ls(pattern = "[.]?tmp[.]?"))


## KO pathway -------------------------------------------------------------------


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

#

ls.ps <-
  list(
    "ko_pathway" = ps.table.ko.pw
  )

ls.lf <- list()

ranks <- c("ko_pathway")
comps <- c("Belgian: pre-HFD|Yanomami: pre-HFD",
           "Belgian: post-HFD|Yanomami: post-HFD",
           "Belgian: pre-HFD|Belgian: post-HFD",
           "Yanomami: pre-HFD|Yanomami: post-HFD")

for (j in seq_along(comps)) {
  
  for (i in seq_along(ranks)) {
    
    tmp.comp.factor <- "donor_ethnicity_timepoint"
    
    tmp.comp <- comps[[j]]
    
    tmp.ps.table <-
      merge_phyloseq(subset_samples(ls.ps[[i]], donor_ethnicity_timepoint == unlist(str_split(tmp.comp, "\\|"))[[1]]),
                     subset_samples(ls.ps[[i]], donor_ethnicity_timepoint == unlist(str_split(tmp.comp, "\\|"))[[2]])) %>%
      prune_samples(sample_sums(.) > 0, .)
    
    tmp.lf.comp <-
      tmp.ps.table %>%
      run_lefse(taxa_rank = "Kingdom",
                wilcoxon_cutoff = Inf,
                group = tmp.comp.factor,
                kw_cutoff = Inf,
                multigrp_strat = TRUE,
                lda_cutoff = 0) %>%
      marker_table() %>%
      data.frame() %>%
      remove_rownames() %>%
      mutate(comp = tmp.comp, padj = p.adjust(pvalue, method = "fdr")) %>%
      mutate(rank = ranks[[i]])
    
    ls.lf[[ranks[[i]]]][[j]] <- tmp.lf.comp
    
  }
  
}; rm(i, j)

df.lefse.raw <-
  bind_rows(ls.lf) %>%
  mutate(comp = factor(comp, levels = unique(comp)))

df.lefse.raw <-
  merge(tb.aux.ann.ko.def %>% select(ko_class:ko_pathway) %>% distinct(), df.lefse.raw, by.x = "ko_pathway", by.y = "feature", all.y = T) %>%
  select(ko_class, ko_subclass, feature = ko_pathway, everything()) %>%
  arrange(comp, enrich_group, -ef_lda)

df.lefse.filt <- df.lefse.raw %>% filter(ef_lda > 2, padj < 0.05)


df.lefse <- list("Raw" = df.lefse.raw, "Filtered LDA>2" = df.lefse.filt)

openxlsx::write.xlsx(x = df.lefse, file = paste0(path.tb.dir, "/diff_abund/lefse.xlsx"), rowNames = F)

for (i in c(1:4)) {
  
  ls.gg.lefse <- list()
  
  ls.hei.lefse <- list()
  
  df.lefse.filt %>%
    filter(comp == unique(.$comp)[[i]]) %>%
    filter(ko_class == "Metabolism") %>%
    arrange(desc(enrich_group), ef_lda) %>%
    mutate(tmp = str_remove(comp, "\\|.*")) %>%
    mutate(ef_lda = ifelse(enrich_group == tmp, -ef_lda, ef_lda)) %>%
    mutate(my.xpos = ifelse(enrich_group == tmp, 0.3, -0.3)) %>%
    mutate(my.hjust = ifelse(enrich_group == tmp, 0, 1)) %>%
    select(-tmp) %>%
    #group_by(enrich_group) %>%
    #slice_tail(n = 20) %>%
    #ungroup() %>%
    data.frame(check.names = F) %>%
    arrange(desc(enrich_group), abs(ef_lda)) %>%
    mutate(feature = str_replace(feature, "\\[PATH:", "[")) %>%
    mutate(
      tmp1 = str_trunc(feature, width = 42, side = "right"),
      tmp2 = str_replace_all(feature, ".* \\[", "["),
      feature = ifelse(str_detect(tmp1, "\\["), feature, str_glue("{tmp1} {tmp2}")),
      feature = str_replace(feature, "\\.\\.\\.", "[...]")) %>%
    select(-starts_with("tmp")) %>%
    mutate(feature = make.unique(feature)) %>%
    mutate(feature = factor(feature, levels = unique(feature))) -> tmp.tb.lefse.gg
  
  if (nrow(tmp.tb.lefse.gg) > 0) {
    
    ls.hei.lefse[[i]] <- tmp.tb.lefse.gg %>% nrow()
    
    tmp.gg.dual.1 <-
      ggplot(tmp.tb.lefse.gg, aes(x = ef_lda, y = feature, fill = enrich_group)) +
      facet_grid(ko_subclass ~ .,  scales = "free", space = "free_y") +
      geom_col(color = "black") +
      geom_vline(xintercept = 0, color = "gray10", linetype = "solid", linewidth = .25) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = ~ c(-1.05, 1.05) * max(abs(tmp.tb.lefse.gg$ef_lda), na.rm = T)) +
      scale_y_discrete(expand = expansion(mult = 0, add = .5), labels = function(x) str_replace(x, "\\]\\..*$", "]")) +
      scale_fill_manual(name = "", values = c(colors.gp.tp)) +
      coord_cartesian(clip = "off") +
      xlab("<br>LDA score (log<sub>10</sub>)") + ylab("<br>") +
      my.theme.2 +
      theme(panel.grid.major.x = element_line(size = 0.25, linetype = 'solid', colour = "gray75"),
            panel.grid.major.y = element_line(size = 0.25, linetype = 'dotted', colour = "gray75"),
            axis.text.y = element_markdown(size = 7, colour = 'black'),
            strip.text.y = element_blank(),
            plot.margin = unit(c(0, 0, 0, .5), 'cm'))
    
    tmp.gg.dual.2 <-
      ggplot(tmp.tb.lefse.gg, aes(x = 1, y = feature, fill = ko_subclass)) +
      facet_grid(ko_subclass ~ .,  scales = "free", space = "free_y") +
      geom_col() +
      scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_fill_manual(name = "", values = colors.ko.met) +
      scale_color_identity(name = "") +
      xlab("") + ylab("") +
      my.theme.2 +
      theme(panel.border = element_rect(colour = "#00000000", fill = NA),
            strip.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_line(color = "#00000000"),
            axis.text.x = element_markdown(size = 6, colour = '#00000000'),
            axis.ticks.x = element_line(color = "#00000000"),
            plot.margin = unit(c(0, .5, 0, -.5), 'cm'))
    
    tmp.gg.dual.1 + tmp.gg.dual.2 + gg.leg.ko.met + plot_layout(nrow = 1, widths = c(.75, .04, 2.5))
    
    if (nrow(tmp.tb.lefse.gg) > 15) {
      
      ls.gg.lefse[[i]] <- tmp.gg.dual.1 + tmp.gg.dual.2 + gg.leg.ko.met + plot_layout(nrow = 1, widths = c(.75, .04, 2.5))
      
    } else {
      
      ls.gg.lefse[[i]] <- tmp.gg.dual.1 + tmp.gg.dual.2 + (plot_spacer() + theme_void()) + plot_layout(nrow = 1, widths = c(.75, .04, 2.5))
      
    }
    
  }
  
  if (length(ls.hei.lefse) > 0 & sum(unlist(ls.hei.lefse)) < 50) {
    
    png(filename = paste0(path.fg.dir, "/diff_abund/daa-lefse-comp", i, ".png"), width = 5000 * 1.75, height = 170 * sum(unlist(ls.hei.lefse)), units = "px", res = 1200, type = "cairo")
    
    print(wrap_plots(unlist(ls.gg.lefse)) + plot_layout(ncol = 1, heights = c(unlist(ls.hei.lefse))))
    
    dev.off()
    
  }
  
  if (length(ls.hei.lefse) > 0 & sum(unlist(ls.hei.lefse)) > 50) {
    
    png(filename = paste0(path.fg.dir, "/diff_abund/daa-lefse-comp", i, ".png"), width = 5000 * 1.75, height = 150 * sum(unlist(ls.hei.lefse)), units = "px", res = 1200, type = "cairo")
    
    print(wrap_plots(unlist(ls.gg.lefse)) + plot_layout(ncol = 1, heights = c(unlist(ls.hei.lefse))))
    
    dev.off()
    
  }
  
}; rm(i)

rm(list = ls(pattern = "^tmp."))


## Heatmap ---------------------------------------------------------------------

# 1

my.lefse.table = df.lefse.raw
my.ps.table = ps.table.ko.pw
my.top = 50

top.taxa <-
  my.lefse.table %>%
  filter(ko_class == "Metabolism") %>%
  arrange(-abs(ef_lda)) %>%
  pull(feature) %>%
  unique() %>%
  head(my.top)

tmp.ps.clr <-
  prune_taxa(tax_table(my.ps.table)[, "Species"] %in% top.taxa, my.ps.table) %>%
  microbiome::transform(., "clr")

tmp.tb.clr <-
  merge(tax_table(tmp.ps.clr), otu_table(tmp.ps.clr), by = "row.names", all = T) %>%
  select(feature = Row.names, any_of(tb.metad$sample_id)) %>%
  merge(., my.lefse.table %>% select(ko_class, ko_subclass, feature), all.x = T) %>%
  distinct()

hp.anno.col <-
  tb.metad %>%
  column_to_rownames("sample_id") %>%
  select(donor_ethnicity_timepoint, mouse_n)

hp.anno.row <-
  tmp.tb.clr %>%
  select(ko_class, ko_subclass, feature) %>%
  column_to_rownames("feature")

hp.values <-
  tmp.tb.clr %>%
  select(feature, any_of(tb.metad$sample_id)) %>%
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
    `Metabolism class` =
      hp.anno.row$ko_subclass,
    col =
      list(`Metabolism class` = colors.ko.met),
    annotation_legend_param =
      list(`Metabolism class` = list(labels = names(colors.ko.met), at = names(colors.ko.met))),
    show_annotation_name = FALSE,
    simple_anno_size = unit(.25, "cm"),
    show_legend = c(TRUE))

hp.b <-
  ComplexHeatmap::columnAnnotation(
    `Condition` =
      hp.anno.col$donor_ethnicity_timepoint,
    col =
      list(`Condition` = colors.gp.tp),
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
                        column_title = str_glue("Top {nrow(hp.values)} most variable pathways\nStandardized CLR-transformed data (CLR-Z)"),
                        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                        show_column_names = F,
                        heatmap_legend_param =
                          list(color_bar = "continuous",
                               ncol = 1,
                               legend_direction = "vertical",
                               legend_height = unit(3, "cm"),
                               grid_width = unit(.25, "cm"),
                               title_position = "leftcenter-rot",
                               title_gp = gpar(fontsize = 10, fontface="bold"),
                               labels_gp = gpar(fontsize = 8)),
                        left_annotation = hp.a,
                        top_annotation = hp.b) -> hp.plot

draw(hp.plot, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 10, 1, 50), "mm"))

png(filename = paste0(path.fg.dir, "/diff_abund/heatmap-clr-z_top50.png"), width = 5000 * 2.5, height = 5000 * (.9 + (nrow(hp.values)/60)), units = "px", res = 1200, type = "cairo")

draw(hp.plot, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 10, 1, 50), "mm"))

dev.off()

# 2
my.lefse.table = df.lefse.raw
my.ps.table = ps.table.ko.pw
my.top = 10

my.classes = c("Lipid metabolism",
               "Carbohydrate metabolism", 
               "Metabolism of cofactors and vitamins",
               "Amino acid metabolism",
               "Energy metabolism",
               "Glycan biosynthesis and metabolism")

for (my.class in my.classes) {
  
  top.taxa <-
    my.lefse.table %>%
    filter(ko_subclass == my.class) %>%
    arrange(-abs(ef_lda)) %>%
    pull(feature) %>%
    unique() %>%
    head(my.top)
  
  tmp.ps.clr <-
    prune_taxa(tax_table(my.ps.table)[, "Species"] %in% top.taxa, my.ps.table) %>%
    microbiome::transform(., "clr")
  
  tmp.tb.clr <-
    merge(tax_table(tmp.ps.clr), otu_table(tmp.ps.clr), by = "row.names", all = T) %>%
    select(feature = Row.names, any_of(tb.metad$sample_id)) %>%
    merge(., my.lefse.table %>% select(ko_class, ko_subclass, feature), all.x = T) %>%
    distinct()
  
  hp.anno.col <-
    tb.metad %>%
    column_to_rownames("sample_id") %>%
    select(donor_ethnicity_timepoint, mouse_n)
  
  hp.anno.row <-
    tmp.tb.clr %>%
    select(ko_class, ko_subclass, feature) %>%
    column_to_rownames("feature")
  
  hp.values <-
    tmp.tb.clr %>%
    select(feature, any_of(tb.metad$sample_id)) %>%
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
      `Metabolism class` =
        hp.anno.row$ko_subclass,
      col =
        list(`Metabolism class` = colors.ko.met),
      annotation_legend_param =
        list(`Metabolism class` = list(labels = names(colors.ko.met), at = names(colors.ko.met))),
      show_annotation_name = F,
      simple_anno_size = unit(.25, "cm"),
      show_legend = c(F))
  
  hp.b <-
    ComplexHeatmap::columnAnnotation(
      `Condition` =
        hp.anno.col$donor_ethnicity_timepoint,
      col =
        list(`Condition` = colors.gp.tp),
      show_annotation_name = F,
      simple_anno_size = unit(.25, "cm"),
      show_legend = F)
  
  # Plot Row
  ComplexHeatmap::Heatmap(hp.values,
                          name = "CLR-Z",
                          col = hp.col.fun.1,
                          cluster_rows = T,
                          cluster_columns = F,
                          row_names_gp = gpar(fontsize = 9),
                          column_title = str_glue("Top {nrow(hp.values)} most variable {my.class} pathways\nStandardized CLR-transformed data (CLR-Z)"),
                          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                          show_column_names = F,
                          heatmap_legend_param =
                            list(color_bar = "continuous",
                                 ncol = 1,
                                 legend_direction = "vertical",
                                 legend_height = unit(3, "cm"),
                                 grid_width = unit(.25, "cm"),
                                 title_position = "leftcenter-rot",
                                 title_gp = gpar(fontsize = 10, fontface="bold"),
                                 labels_gp = gpar(fontsize = 8)),
                          left_annotation = hp.a,
                          top_annotation = hp.b) -> hp.plot
  
  draw(hp.plot, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 10, 1, 50), "mm"))
  
  png(filename = paste0(path.fg.dir, "/diff_abund/heatmap-clr-z_", str_replace_all(str_to_lower(my.class), " ", "-"),".png"), width = 5000 * 2, height = 5000 * .5, units = "px", res = 1200, type = "cairo")
  
  draw(hp.plot, merge_legend = TRUE, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(1, 10, 1, 50), "mm"))
  
  dev.off()
  
}


# End ---------------------------------------------------------------------

