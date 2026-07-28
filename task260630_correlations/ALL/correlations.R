
# Correlations -----------------------------------------------------------------


# Setup ------------------------------------------------------------------------


## Environment and Packages ----------------------------------------------------


rm(list = ls())
graphics.off()
gc()

set.seed("123")

library("ggrepel")
library("ggtext")
library("ggpubr")
library("ggcorrplot")
library("patchwork")
library("cowplot")
library("svglite")
library("ComplexHeatmap")
library("openxlsx")
library("xlsx")
library("DESeq2")
library("vegan")
library("psych")
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

path.wk.dir <- paste0(getwd(), "/task260630_correlations/")
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


fun.get.flatten.corr.mt <- function(cormat, pmat) {
  
  # cormat : matrix of the correlation coefficients
  # pmat : matrix of the correlation p-values
  
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[row(cormat)[ut]],
    to = rownames(cormat)[col(cormat)[ut]],
    r  = cormat[ut],
    p = pmat[ut]
  )
}

## Data loading ----------------------------------------------------------------


### Tax data -------------------------------------------------------------------


# Metadata
tb.metad <- 
  read.delim(paste0(getwd(), "/misc/metadata_full.tsv")) %>% 
  filter(experiment_type == "high-fat diet") %>% 
  mutate(timepoint = case_when(
    experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
    experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD",
  )) %>% 
  mutate(mouse_n = ifelse(is_pooled == "N", str_glue("{experiment}_{donor_id}_{mouse_n}"), str_glue("{experiment}_{donor_id}"))) %>% 
  mutate(donor_ethnicity = ifelse(donor_ethnicity == "Western", "Belgian", donor_ethnicity)) %>% 
  select(experiment, sample_id, assembly_id, donor_id:cage_id, mouse_n, timepoint) %>% 
  mutate(donor_ethnicity = factor(donor_ethnicity, c("Belgian", "Yanomami"))) %>% 
  mutate(timepoint = factor(timepoint, c("pre-HFD", "post-HFD"))) %>% 
  mutate(experiment = factor(experiment, levels = unique(.$experiment) %>% sort())) %>% 
  arrange(experiment, donor_ethnicity, timepoint, mouse_n) %>% 
  mutate(donor_ethnicity_timepoint = str_glue("{donor_ethnicity}: {timepoint}")) %>% 
  mutate(donor_ethnicity_timepoint = factor(donor_ethnicity_timepoint, unique(.$donor_ethnicity_timepoint))) %>% 
  mutate(mouse_n = str_replace(mouse_n, "Donor-2-MCBL1", "D2M-001")) %>%
  mutate(mouse_n = factor(mouse_n, unique(.$mouse_n))) %>% 
  mutate(group_id = str_remove(mouse_n, "_M[0-9]+$")) %>% 
  mutate(group_id = factor(group_id, unique(.$group_id))); head(tb.metad)

tmp.clinic <-
  openxlsx::read.xlsx(xlsxFile = paste0(getwd(), "/misc/clinical-data-HFD.xlsx"), sheet = 1) %>% 
  select(!contains("_log10")) %>% 
  mutate(FFA = ifelse(is.na(FFA_1_2), FFA_3_4_5, FFA_1_2)) %>% 
  select(-c(donor_ethnicity:mouse_id), -c(weight_pre, weight_post, liver_weight, eWAT_weight, FFA_1_2, FFA_3_4_5)) %>%
  mutate(mouse_n = str_replace(mouse_n, "D2M", "D2M-001")) %>% 
  mutate(group_id = str_remove(mouse_n, "_M[0-9]+$")) %>% 
  mutate(group_id = factor(group_id, unique(.$group_id)))

tmp.order <- c(
  # Whole-animal phenotype
  "weight_delta",
  
  # Adiposity
  "eWAT_relative_weight",
  "adipocyte_area",
  
  # Glucose homeostasis
  "fasting_blood_glucose",
  "insulin",
  
  # Energy homeostasis hormones
  "leptin",
  "ghrelin",
  
  # Adipokines / metabolic mediators
  "resistin",
  "PAI_1",
  "FFA",
  
  # Systemic inflammation
  "MCP1_serum",
  "TNF_serum",
  "IL6_serum",
  "IL17_serum",
  "IFNg_serum",
  
  # Liver phenotype
  "liver_relative_weight",
  "ALT_activity",
  "ORO_liver",
  "liver_pathology",
  
  # Gut phenotype
  "colon_mucus"
)

tmp.clinic <- 
  tmp.clinic %>% 
  select(mouse_n, group_id, all_of(tmp.order))
  

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
  select(-c(sample_id, mouse_n)) %>% 
  filter(timepoint == "post-HFD") %>% 
  group_by(across(-c(cpm))) %>% 
  summarise(cpm = mean(cpm, na.rm = T), .groups = "drop") %>% 
  select(tax_id:Species, group_id, cpm) %>% 
  pivot_wider(names_from = "group_id", values_from = "cpm") %>% 
  data.frame(check.names = F)

tb.metad <- 
  merge(tb.metad, tmp.clinic, all.x = T) %>% 
  select(sample_id, everything()) %>% 
  pivot_longer(cols = any_of(tmp.order), names_to = "tmp", values_to = "value") %>% 
  select(-c(sample_id, cage_id, mouse_n)) %>% 
  filter(timepoint == "post-HFD") %>%
  group_by(across(-c(value))) %>% 
  summarise(value = mean(value, na.rm = T), .groups = "drop") %>% 
  mutate(tmp = factor(tmp, levels = tmp.order)) %>% 
  arrange(tmp) %>% 
  pivot_wider(names_from = "tmp", values_from = "value") %>% 
  data.frame(check.names = F)

tmp.tb.metad <- tb.metad %>% rename(ID = group_id)

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

colors.msl <- c("#67001FFF", "#B2182BFF", "#F4A582FF", "#D1E5F0FF", "#2166ACFF", "#053061FF")
names(colors.msl) <- c("[1, Inf)", "[0.5, 1)", "[0, 0.5)", "(−0.5, 0)", "(−1, −0.5]", "(-Inf, −1]")

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

gg.leg.msl <-
  ggplot(names(colors.msl) %>% data.frame() %>% mutate(across(everything(), .fns = ~factor(., unique(.)))), aes(x = 1, y = ., fill = .)) +
  geom_point(shape = 22, size = 3) +
  scale_fill_manual(name = "MaAsLin (β) coefficient", values = colors.msl) +
  my.theme.1 +
  theme(legend.position = "right",
        legend.title = element_text(size = 8, colour = 'black', face = "bold"),
        legend.text = element_text(size = 7, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.msl <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.msl)); gg.leg.msl


gg.leg.all <- gg.leg.gp + gg.leg.gp.tp + gg.leg.tax +  gg.leg.msl + plot_layout(ncol = 1, heights = c(3, 3, 8, 7))
gg.leg.all

png(filename = paste0(path.fg.dir, "/legends.png"), width = 5000 * .65, height = 5000 * 1.75, units = "px", res = 1200, type = "cairo")

gg.leg.all

dev.off()


# Correlations -----------------------------------------------------------------


## Taxa - Meta -----------------------------------------------------------------


# Glom tables

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

# MaAsLin2

tmp.comps <- c("Belgian: post-HFD|Yanomami: post-HFD")

tmp.ls.msl.meta <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.meta.var <- names(tb.metad)[10:ncol(tb.metad)]
  
  for (j in seq_along(tmp.meta.var)) {
    
    tmp.var <- tmp.meta.var[[j]]
    
    tmp.comp.factor <- "donor_ethnicity_timepoint" 
    tmp.comp <- tmp.comps[[1]]
    
    tmp.meta <- 
      tb.metad %>% 
      filter(!!as.symbol(tmp.comp.factor) %in% unlist(str_split(tmp.comp, "\\|"))) %>% 
      column_to_rownames("group_id") %>% 
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
               min_abundance = 0.01,
               #min_prevalence = 0.5,
               plot_heatmap = F, 
               plot_scatter = F); head(tmp.msl.comp)
    
    tmp.ls.msl.meta[[ranks[[i]]]][[j]] <- 
      merge(tmp.rel.tax.id, tmp.msl.comp$results, by.x = "tmpID", by.y = "feature") %>% 
      merge(., tmp.prev.table) %>% 
      mutate(comp = tmp.comp, rank = ranks[[i]]) %>% 
      select(rank, feature = ranks[[i]], metadata, coef, stderr, pval, qval, 
             N, N.not.zero, any_of(unlist(str_split(tmp.comp, pattern = "\\|"))), comp) 
    
  }
  
}; rm(i, j)

# Annotating and Filtering

tmp.ls.tax.ann <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]] 
  
  tmp.ls.tax.ann[[i]] <- 
    tax_table(ps.table) %>% 
    data.frame(check.names = F) %>% 
    remove_rownames() %>% 
    mutate(feature = .[[tmp.rank]]) %>% 
    select(Kingdom, Phylum, feature) %>% 
    distinct()
  
}; rm(i, tmp.rank)

tmp.tb.tax.ann <- 
  bind_rows(tmp.ls.tax.ann) %>% 
  mutate(Phylum2 = case_when(
    ! Phylum %in% names(colors.tax) & str_detect(Phylum, "Unassigned") ~ "Unassigned",
    ! Phylum %in% names(colors.tax) & !str_detect(Kingdom, "d: Bacteria") ~ "Other: Domains",
    ! Phylum %in% names(colors.tax) & str_detect(Kingdom, "d: Bacteria") ~ "Other: Bacteria",
    Phylum %in% names(colors.tax) ~ Phylum)) %>%
  filter(Phylum2 != "Unassigned") %>% 
  mutate(Phylum2 = factor(Phylum2, names(colors.tax)))

df.msl.taxa.meta.raw <- 
  bind_rows(tmp.ls.msl.meta) %>% 
  mutate(comp = factor(comp, levels = unique(comp))) %>% 
  mutate(rank = factor(rank, levels = ranks)) %>% 
  mutate(metadata = factor(metadata, levels = names(tb.metad)[10:ncol(tb.metad)])) %>% 
  merge(., tmp.tb.tax.ann, "feature", all.x = T) %>% 
  select(rank, Kingdom:Phylum2, feature, metadata, coef:comp) %>% 
  arrange(rank, metadata, -abs(coef), qval) %>% 
  filter(!str_detect(feature, "Un.\\)$"))

df.msl.taxa.meta.filt <- df.msl.taxa.meta.raw %>% filter(qval < 0.1)

df.msl.taxa.meta <- list("Raw" = df.msl.taxa.meta.raw, "Filtered" = df.msl.taxa.meta.filt)

openxlsx::write.xlsx(x = df.msl.taxa.meta, file = paste0(path.tb.dir, "/maaslin_taxa-meta.xlsx"), rowNames = F)

# Graph

tmp.corr.flatted <-
  df.msl.taxa.meta.raw %>% 
  filter(rank == "Genus") %>% 
  mutate(sig.f = case_when(
    qval <= 0.1 ~ "*", 
    #qval > 0.1 & pval < 0.05 ~ "+",
    .default = "")
  ) %>% 
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

tmp.taxa.has.sig <- tmp.corr.flatted %>% filter(sig.f != "") %>% pull(feature) %>% unique()

tmp.corr.flatted.filt <- 
  tmp.corr.flatted %>% 
  #filter(feature %in% tmp.taxa.has.sig) %>% 
  rename(metad = metadata) %>% 
  mutate(feature2 = str_glue("{feature} (B={`Belgian: post-HFD`};Y={`Yanomami: post-HFD`})"))

tmp.gg.msl.1 <- 
  ggplot(tmp.corr.flatted.filt, aes(x = metad, y = feature2, fill = r.tier)) +
  geom_tile(colour = "white", size = .25) +
  geom_text(aes(label = sig.f, color = r.color), size = 2.5, hjust= .5, vjust = .75, , fontface='bold') +
  scale_color_identity() +
  scale_fill_manual(values = colors.msl) + 
  xlab("") + ylab("") +
  my.theme.1 + 
  theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(0, 0, 0, .5), 'cm')); tmp.gg.msl.1

tmp.gg.msl.2 <-
  ggplot(tmp.corr.flatted.filt, aes(x = 1, y = feature2, fill = Phylum2)) +
  geom_tile(colour = "white", size = .25) +
  scale_fill_manual(name = "", values = colors.tax) +
  scale_color_identity(name = "") +
  xlab("") + ylab("") +
  my.theme.1 +
  theme(panel.border = element_rect(colour = "#00000000", fill = NA),
        strip.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = "#00000000"),
        axis.text.x = element_text(colour = '#00000000'),
        axis.ticks.x = element_line(color = "#00000000"),
        plot.margin = unit(c(0, .5, 0, -1), 'cm'))

print(tmp.gg.msl.1 + plot_spacer() + tmp.gg.msl.2 + (gg.leg.msl/gg.leg.tax) + plot_layout(nrow = 1, widths = c(.75, -.1, .04, 1.25)))

png(filename = paste0(path.fg.dir, "/maaslin_taxa-meta_full.png"), width = 5000 * 2, height = 160 * (length(unique(tmp.corr.flatted.filt$feature)) + 7), units = "px", res = 1200, type = "cairo")

print(tmp.gg.msl.1 + plot_spacer() + tmp.gg.msl.2 + (gg.leg.msl/gg.leg.tax) + plot_layout(nrow = 1, widths = c(.75, -.1, .04, 1.25)))

dev.off()

tmp.corr.flatted.filt <- 
  tmp.corr.flatted %>% 
  filter(feature %in% tmp.taxa.has.sig) %>% 
  rename(metad = metadata) %>% 
  mutate(feature2 = str_glue("{feature} (B={`Belgian: post-HFD`};Y={`Yanomami: post-HFD`})"))

tmp.gg.msl.1 <- 
  ggplot(tmp.corr.flatted.filt, aes(x = metad, y = feature2, fill = r.tier)) +
  geom_tile(colour = "white", size = .25) +
  geom_text(aes(label = sig.f, color = r.color), size = 2.5, hjust= .5, vjust = .75, , fontface='bold') +
  scale_color_identity() +
  scale_fill_manual(values = colors.msl) + 
  xlab("") + ylab("") +
  my.theme.1 + 
  theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1),
        plot.margin = unit(c(0, 0, 0, .5), 'cm')); tmp.gg.msl.1

tmp.gg.msl.2 <-
  ggplot(tmp.corr.flatted.filt, aes(x = 1, y = feature2, fill = Phylum2)) +
  geom_tile(colour = "white", size = .25) +
  scale_fill_manual(name = "", values = colors.tax) +
  scale_color_identity(name = "") +
  xlab("") + ylab("") +
  my.theme.1 +
  theme(panel.border = element_rect(colour = "#00000000", fill = NA),
        strip.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = "#00000000"),
        axis.text.x = element_text(colour = '#00000000'),
        axis.ticks.x = element_line(color = "#00000000"),
        plot.margin = unit(c(0, .5, 0, -1), 'cm'))

print(tmp.gg.msl.1 + plot_spacer() + tmp.gg.msl.2 + (gg.leg.msl+gg.leg.tax) + plot_layout(nrow = 1, widths = c(.75, -.1, .04, 1.25)))

png(filename = paste0(path.fg.dir, "/maaslin_taxa-meta_sig.png"), width = 5000 * 2, height = 180 * (length(unique(tmp.corr.flatted.filt$feature)) + 7), units = "px", res = 1200, type = "cairo")

print(tmp.gg.msl.1 + plot_spacer() + tmp.gg.msl.2 + (gg.leg.msl+gg.leg.tax) + plot_layout(nrow = 1, widths = c(.75, -.1, .04, 1.25)))

dev.off()

rm(i, j)
rm(list = ls(pattern = "[.]?tmp[.]?"))


# End ---------------------------------------------------------------------

