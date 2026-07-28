
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
library("Maaslin2")
library("tidyverse")

getwd()
setwd("~/YanomamiMouseHFD/")
getwd()

list.dirs()

path.wk.dir <- paste0(getwd(), "/task260624_tax-profiler_w-donors")
path.dt.dir <- paste0(getwd(), "/tables")

for (i in c(path.wk.dir)) {
  
  if (!dir.exists(i)) { dir.create(i) }
  
}; remove(i)

my.exp <- "HFD1"

path.fg.dir <- paste0(path.wk.dir, "/", my.exp ,"/figures")
path.tb.dir <- paste0(path.wk.dir, "/", my.exp,  "/tables")

for (i in c(path.fg.dir, path.tb.dir)) {
  
  if (!dir.exists(i)) { dir.create(i, recursive = T) }
  
}; remove(i)


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
                                gp_shape = NULL, gp_facet, value_var,
                                only.sig = F, alpha = 0.1, label.type = "stars", step.mult = 0.08) {
  
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
    step.mult = step.mult
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
      {if(is.null(gp_shape)) geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
      {if(!is.null(gp_shape)) geom_point(aes(shape = .data[[gp_shape]]), size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
      {if(!is.null(gp_shape)) scale_shape_manual(values = shapes.exp)} +
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
    {if(is.null(gp_shape)) geom_point(shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
    {if(!is.null(gp_shape)) geom_point(aes(shape = .data[[gp_shape]]), size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123))} +
    {if(!is.null(gp_shape)) scale_shape_manual(values = shapes.exp)} +
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


# fun.get.barp.compm3 <- function(tb, tb.compm, tb.kruskal, tb.summ, gp_comp,  gp_fill = NULL, gp_col = NULL,
#                                 gp_facet, value_var,
#                                 only.sig = F, alpha = 0.1, label.type = "stars") {
#   
#   if (is.null(gp_fill)) { gp_fill = gp_comp }
#   if (is.null(gp_col)) { gp_col = gp_comp }
#   
#   tb.brackets <- fun.get.brackets(
#     tb = tb,
#     tb.compm = tb.compm,
#     gp_facet = gp_facet,
#     gp_comp = gp_comp,
#     value_var = value_var,
#     only.sig = only.sig,
#     alpha = alpha,
#     label.type = label.type,
#     step.mult = .5
#   )
#   
#   tb.global <- tb.kruskal %>%
#     mutate(
#       p.full = paste0(
#         "KW, 𝘱-adj: ",
#         ifelse(adj.p.value < 0.001, "<0.001", sprintf("%1.3f", round(adj.p.value, 3)))
#       )
#     )
#   
#   if (nrow(tb.brackets) > 0) {
#     
#     gg <- ggplot() +
#       facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
#       
#       geom_bar(data = tb.summ, mapping = aes(x = .data[[gp_comp]], y = .data[["mean"]],
#                                    fill = .data[[gp_fill]], color = .data[[gp_col]]),
#                stat="identity", alpha = .35, position=position_dodge()) +
#       # geom_errorbar(data = tb.summ, 
#       #               mapping = aes(x = .data[[gp_comp]],
#       #                             ymin=.data[["mean"]]-.data[["sd"]], ymax=.data[["mean"]]+.data[["sd"]]), 
#       #               width=.2, position=position_dodge(.9)) + 
#     
#       geom_point(tb, mapping = aes(x = .data[[gp_comp]], y = .data[[value_var]],
#                          fill = .data[[gp_fill]], color = .data[[gp_col]]),
#                  shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
#       geom_text(
#         data = tb.global,
#         aes(x = 0.5, y = Inf, label = p.full),
#         inherit.aes = FALSE,
#         vjust = 2,
#         hjust = 0,
#         size = 2,
#         color = "black"
#       ) +
#       ggpubr::stat_pvalue_manual(
#         tb.brackets,
#         xmin = "group1",
#         xmax = "group2",
#         y.position = "y.position",
#         label = "label",
#         tip.length = 0.01,
#         bracket.size = 0.25,
#         size = 2.5,
#         hide.ns = only.sig,
#         inherit.aes = FALSE
#       ) +
#       scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
#       scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
#       scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
#       coord_cartesian(clip = "off") +
#       ylab("") + xlab("") +
#       my.theme.1 + 
#       theme(legend.position = "none",
#             legend.title = element_blank(),
#             axis.text.x = element_blank(),
#             axis.ticks.x = element_blank(),
#             plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
#     
#     return(gg)
#     
#   }
#   
#   
#   gg <- ggplot() +
#     facet_wrap(stats::as.formula(paste("~", gp_facet)), scales = "free_y", nrow = 1) +
#     
#     geom_bar(data = tb.summ, mapping = aes(x = .data[[gp_comp]], y = .data[["mean"]],
#                                  fill = .data[[gp_fill]], color = .data[[gp_col]]),
#              stat="identity", alpha = .35, position=position_dodge()) +
#     # geom_errorbar(data = tb.summ, 
#     #               mapping = aes(x = .data[[gp_comp]], 
#     #                             ymin=.data[["mean"]]-.data[["sd"]], ymax=.data[["mean"]]+.data[["sd"]]), 
#     #               width=.2, position=position_dodge(.9)) + 
#     
#     geom_point(data = tb, mapping = aes(x = .data[[gp_comp]], y = .data[[value_var]],
#                        fill = .data[[gp_fill]], color = .data[[gp_col]]),
#                shape = 16, size = 1, position = position_jitter(w = 0.1, h = 0, seed = 123)) +  
#     geom_text(
#       data = tb.global,
#       aes(x = 0.5, y = Inf, label = p.full),
#       inherit.aes = FALSE,
#       vjust = 2,
#       hjust = 0,
#       size = 2,
#       color = "black"
#     ) +
#     
#     scale_fill_manual(values = c(colors.gp, colors.gp.tp)) +  
#     scale_color_manual(values = c(colors.gp, colors.gp.tp)) +  
#     scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
#     coord_cartesian(clip = "off") +
#     ylab("") + xlab("") +
#     my.theme.1 + 
#     theme(legend.position = "none",
#           legend.title = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
#   
#   return(gg)
#   
# }



## Data loading ----------------------------------------------------------------


# Metadata
tb.metad <- 
  bind_rows(
    read.delim(paste0(getwd(), "/misc/metadata_full.tsv")),
    read.delim(paste0(getwd(), "/misc/metadata_gav.tsv")),
    read.delim(paste0(getwd(), "/misc/metadata_donor.tsv"))
  ) %>% 
  filter(experiment == my.exp) %>% 
  mutate(timepoint = case_when(
    experiment_type == "high-fat diet" & timepoint == "T1" ~ "pre-HFD",
    experiment_type == "high-fat diet" & timepoint == "T2" ~ "post-HFD",
    .default = experiment_type)
  ) %>% 
  mutate(donor_ethnicity = ifelse(donor_ethnicity == "Western", "Belgian", donor_ethnicity)) %>% 
  select(sample_id, assembly_id, experiment_type, donor_id:cage_id, mouse_n, timepoint) %>% 
  mutate(donor_ethnicity = factor(donor_ethnicity, c("Belgian", "Yanomami"))) %>% 
  mutate(timepoint = factor(timepoint, c("donor", "gavage mix", "pre-HFD", "post-HFD"))) %>% 
  mutate(experiment_type = factor(experiment_type, c("donor", "gavage mix", "high-fat diet"))) %>% 
  arrange(donor_ethnicity, experiment_type, timepoint, mouse_n) %>% 
  mutate(donor_ethnicity_timepoint = 
           ifelse(experiment_type == "high-fat diet", str_glue("{donor_ethnicity}: {timepoint}"), str_glue("{donor_ethnicity}: pre-HFD"))) %>% 
  mutate(donor_ethnicity_timepoint = factor(donor_ethnicity_timepoint, unique(.$donor_ethnicity_timepoint))); head(tb.metad)


### Mice ------------------------------------------------------------------------


# Taxonomy + Counts

mc.ls.counts <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "high-fat diet") %>% pull(assembly_id))) {
  
  mc.ls.counts[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_counts.xlsx"), sheet = 1)
  
}; rm(i)

mc.ls.taxa <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "high-fat diet") %>% pull(assembly_id))) {
  
  mc.ls.taxa[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/{i}/{i}_taxonomy.xlsx"), sheet = 1)
  
}; rm(i)

mc.ls.taxa.counts <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "high-fat diet") %>% pull(assembly_id))) {
  
  mc.ls.taxa.counts[[i]] <- 
    merge(mc.ls.taxa[[i]], mc.ls.counts[[i]]) %>% 
    group_by(across(tax_id:Species)) %>%
    summarise(across(.cols = any_of(tb.metad %>% filter(experiment_type == "high-fat diet") %>% pull(sample_id)), ~ sum(.x, na.rm = TRUE))) %>% 
    data.frame(check.names = F) %>% 
    mutate("prev" = rowSums(.[tb.metad %>% filter(assembly_id == i) %>% pull(sample_id)] != 0)) %>%
    filter(prev >= (tb.metad %>% filter(assembly_id == i) %>% pull(sample_id) %>% length())*.1) %>%
    select(-prev)
  
}; rm(i)

mc.tmp.taxa.counts <-
  reduce(mc.ls.taxa.counts, full_join) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))


### Gavages -----------------------------------------------------------------


# Taxonomy + Counts

gv.ls.counts <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "gavage mix") %>% pull(assembly_id))) {
  
  gv.ls.counts[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/gav/{i}/{i}_counts.xlsx"), sheet = 1)
  
}; rm(i)

gv.ls.taxa <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "gavage mix") %>% pull(assembly_id))) {
  
  gv.ls.taxa[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/gav/{i}/{i}_taxonomy.xlsx"), sheet = 1)
  
}; rm(i)

gv.ls.taxa.counts <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "gavage mix") %>% pull(assembly_id))) {
  
  gv.ls.taxa.counts[[i]] <- 
    merge(gv.ls.taxa[[i]], gv.ls.counts[[i]]) %>% 
    group_by(across(tax_id:Species)) %>%
    summarise(across(.cols = any_of(tb.metad %>% filter(experiment_type == "gavage mix") %>% pull(sample_id)), ~ sum(.x, na.rm = TRUE))) %>% 
    data.frame(check.names = F)
  
}; rm(i)

gv.tmp.taxa.counts <-
  reduce(gv.ls.taxa.counts, full_join) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))


### Donors ---------------------------------------------------------------------


# Taxonomy + Counts

dn.ls.counts <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "donor") %>% pull(assembly_id))) {
  
  dn.ls.counts[[i]] <- 
    openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/donor/{i}/{i}_counts.xlsx"), sheet = 1) %>% 
    select(contig_id, any_of(tb.metad %>% filter(experiment_type == "donor" & assembly_id == i) %>% pull(sample_id))) %>% 
    mutate("prev" = rowSums(.[tb.metad %>% filter(experiment_type == "donor" & assembly_id == i) %>% pull(sample_id)] != 0)) %>%
    filter(prev > 0) %>%
    select(-prev)
  
}; rm(i)

dn.ls.taxa <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "donor") %>% pull(assembly_id))) {
  
  dn.ls.taxa[[i]] <- openxlsx::read.xlsx(xlsxFile = str_glue("{path.dt.dir}/donor/{i}/{i}_taxonomy.xlsx"), sheet = 1)
  
}; rm(i)

dn.ls.taxa.counts <- list()

for (i in unique(tb.metad %>% filter(experiment_type == "donor") %>% pull(assembly_id))) {
  
  dn.ls.taxa.counts[[i]] <- 
    merge(dn.ls.taxa[[i]], dn.ls.counts[[i]]) %>% 
    group_by(across(tax_id:Species)) %>%
    summarise(across(.cols = any_of(tb.metad %>% filter(experiment_type == "donor") %>% pull(sample_id)), ~ sum(.x, na.rm = TRUE))) %>% 
    data.frame(check.names = F)
  
}; rm(i)

dn.tmp.taxa.counts <-
  reduce(dn.ls.taxa.counts, full_join) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))


### PSfying --------------------------------------------------------------------


tmp.taxa.counts <- 
  reduce(list(mc.tmp.taxa.counts, gv.tmp.taxa.counts, dn.tmp.taxa.counts), full_join) %>% 
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

rm(list = ls(pattern = "[.]?ls."))
rm(list = ls(pattern = "[.]?tmp."))


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

colors.exp <- c("#CA925BFF", "#CFC99FFF", "#89AAC6FF")
names(colors.exp) <- levels(tb.metad$experiment_type)

shapes.exp <- c(21, 22, 24)
names(shapes.exp) <- levels(tb.metad$experiment_type)

# Legends

gg.leg.gp <-
  ggplot(tb.metad %>% filter(!is.na(donor_ethnicity)), aes(x = 1, y = donor_ethnicity, fill = donor_ethnicity)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.gp) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.gp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.gp))

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

gg.leg.gp.tp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.gp.tp))

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

gg.leg.tax <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.tax))

gg.leg.exp <-
  ggplot(tb.metad %>% filter(!is.na(experiment_type)), aes(x = 1, y = experiment_type, fill = experiment_type)) + 
  geom_point(shape = 22, size = 3) + 
  scale_fill_manual(name = "", values = colors.exp, labels = str_to_sentence(names(colors.exp))) +  
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10, colour = 'black', face = "bold"),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(0, 0, 0, 2)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.exp <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.exp))

gg.leg.exp.shape <-
  ggplot(tb.metad %>% filter(!is.na(experiment_type)), aes(x = 1, y = experiment_type, shape = experiment_type)) + 
  geom_point(size = 2, fill = "black") + 
  scale_shape_manual(values = shapes.exp) + 
  my.theme.1 + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = 'black', margin = margin(2, 0, 2, 5)),
        legend.key.size = unit(.5, 'lines'),
        legend.justification = c(0, .5))

gg.leg.exp.shape <- ggpubr::as_ggplot(ggpubr::get_legend(gg.leg.exp.shape))

gg.leg.all <- gg.leg.gp + gg.leg.gp.tp + gg.leg.tax + gg.leg.exp + gg.leg.exp.shape + plot_layout(ncol = 1, heights = c(3, 5, 8, 4, 4))
gg.leg.all

png(filename = paste0(path.fg.dir, "/legends.png"), width = 5000 * .65, height = 5000 * 1, units = "px", res = 1200, type = "cairo")

gg.leg.all

dev.off()


# Analysis ---------------------------------------------------------------------


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

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

top.phylum <- c("Bacillota", "Bacteroidota", "Verrucomicrobiota", "Pseudomonadota", "Actinomycetota")


## FMT engraftment -------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/engraftment"))) { dir.create(str_glue("{path.fg.dir}/engraftment"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/engraftment"))) { dir.create(str_glue("{path.tb.dir}/engraftment"), recursive = T) }

fun.get.core <- function(my.glom.tb.ls, my.metad, my.rank, my.group) {
  
  my.metad = my.metad %>% filter(donor_ethnicity == my.group)
  my.glom = my.glom.tb.ls[[my.rank]]
  
  if (my.metad %>% filter(experiment_type == "donor") %>% pull("sample_id") %>% length() > 0) {
    
    tmp.1 <- 
      my.glom %>% 
      select(taxa = my.rank, any_of(my.metad %>% filter(experiment_type == "donor") %>% pull("sample_id"))) %>% 
      filter(!str_detect(taxa, "Un.\\)")) %>% 
      rowwise() %>%
      mutate(mean_of_cols = mean(c_across(any_of(my.metad %>% filter(experiment_type == "donor") %>% pull("sample_id"))), na.rm = TRUE)) %>%
      ungroup() %>% 
      arrange(-mean_of_cols) %>% 
      mutate(cumsum_of_means = cumsum(mean_of_cols))
    
    tmp.sum <- sum(tmp.1$mean_of_cols) * .99
    
    return(tmp.1 %>% filter(cumsum_of_means < tmp.sum) %>% pull(taxa))
    
  }
  
  tmp.1 <- 
    my.glom %>% 
    select(taxa = my.rank, any_of(my.metad %>% filter(experiment_type == "gavage mix") %>% pull("sample_id"))) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
    rowwise() %>%
    mutate(mean_of_cols = mean(c_across(any_of(my.metad %>% filter(experiment_type == "gavage mix") %>% pull("sample_id"))), na.rm = TRUE)) %>%
    ungroup() %>% 
    arrange(-mean_of_cols) %>% 
    mutate(cumsum_of_means = cumsum(mean_of_cols))
  
  tmp.sum <- sum(tmp.1$mean_of_cols) * .99
  
  return(tmp.1 %>% filter(cumsum_of_means < tmp.sum) %>% pull(taxa))
  
}

# fun.get.venn <- function(my.glom.tb.ls, my.metad, my.rank, my.group) {
#   
#   my.metad = my.metad %>% filter(donor_ethnicity == my.group)
#   
#   fun.get.core(my.glom.tb.ls, my.metad, my.rank, my.group)
#   
#   tb.venn <- 
#     my.glom.tb.ls[[my.rank]] %>% 
#     mutate(taxa = .[[my.rank]]) %>% 
#     filter(taxa %in% fun.get.core(my.glom.tb.ls, my.metad, my.rank, my.group)) %>% 
#     select(taxa, where(is.numeric)) %>% 
#     column_to_rownames("taxa") %>% 
#     t() %>% 
#     data.frame(check.names = F) %>% 
#     rownames_to_column("sample_id") %>% 
#     merge(my.metad, .) %>% 
#     mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
#     group_by(experiment_type) %>% 
#     summarise(across(.cols = where(is.numeric), ~ sum(.x, na.rm = T)), .groups = "drop") %>% 
#     column_to_rownames("experiment_type") %>% 
#     t() %>% 
#     data.frame(check.names = F) %>% 
#     mutate(across(levels(my.metad$experiment_type), ~ifelse(. > 0, 1, 0))) %>% 
#     rownames_to_column("ID")
#   
#   ls.venn <- lst(
#     !!levels(my.metad$experiment_type)[[1]] := tb.venn %>% filter(!!sym(levels(my.metad$experiment_type)[1]) != 0) %>% pull(ID),
#     !!levels(my.metad$experiment_type)[[2]] := tb.venn %>% filter(!!sym(levels(my.metad$experiment_type)[2]) != 0) %>% pull(ID),
#     !!levels(my.metad$experiment_type)[[3]] := tb.venn %>% filter(!!sym(levels(my.metad$experiment_type)[3]) != 0) %>% pull(ID)
#   )
#   
#   gg.venn <- 
#     ggvenn(ls.venn, 
#            fill_color = c(colors.exp[[1]],
#                           colors.exp[[2]],
#                           colors.exp[[3]]), 
#            stroke_color = "white", 
#            stroke_size = 0.5, 
#            set_name_size = 0, 
#            fill_alpha = .5, 
#            text_size = 3.5,  
#            show_percentage = F) +  
#     labs(title = str_glue("{my.group}: {my.rank}")) + 
#     theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)); gg.venn
#   
#   return(gg.venn)
#   
# }
# 
# 
# for (r in ranks) {
#   
#   print(fun.get.venn(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r, my.group = "Yanomami"))
#   print(fun.get.venn(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r, my.group = "Belgian"))
#   
# }; rm(r)


fun.get.rare.bp <- function(my.glom.tb.ls, my.metad, my.rank, my.group) {
  
  my.metad = my.metad %>% filter(donor_ethnicity == my.group)
  
  my.top <- fun.get.core(my.glom.tb.ls, my.metad, my.rank, my.group)
  
  my.glom.tb.ls[[my.rank]] %>% 
    mutate(taxa = .[[my.rank]]) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
    mutate(class = ifelse(taxa %in% my.top, "core", "rare")) %>% 
    mutate(class = factor(class, c("rare", "core"))) %>% 
    select(taxa, class, any_of(my.metad$sample_id)) %>% 
    pivot_longer(cols = -c(taxa, class), names_to = "sample_id", values_to = "tpm") %>% 
    mutate(prev = ifelse(tpm > 0, 1, 0)) %>% 
    group_by(class, sample_id) %>% 
    summarise(tpm = sum(tpm, na.rm = T),
              prev = sum(prev, na.rm = T), .groups = "drop") %>%
    group_by(sample_id) %>% 
    mutate(label_pos = ifelse(class == "rare", sum(tpm, na.rm = T), sum(tpm, na.rm = T)/2)) %>% 
    ungroup() %>% 
    mutate(label_col = ifelse(class == "rare", "#762A83FF", "#033800FF")) %>% 
    mutate(label_size = ifelse(class == "rare", 2, 2)) %>% 
    mutate(label_vjust = ifelse(class == "rare", -1, .5)) %>% 
    merge(my.metad, .) %>% 
    mutate(my.facet.1 = case_when(
      experiment_type == "donor" ~ "",
      experiment_type == "gavage mix" ~ "",
      experiment_type == "high-fat diet" ~ "HFD"
    )) %>% 
    mutate(my.facet.2 = case_when(
      experiment_type == "donor" ~ "D",
      experiment_type == "gavage mix" ~ "G",
      experiment_type == "high-fat diet" ~ timepoint
    )) %>% 
    mutate(my.facet.2 = factor(my.facet.2, c("D", "G", "pre-HFD", "post-HFD"))) %>% 
    mutate(my.id = case_when(
      experiment_type == "donor" ~ "Stool",
      experiment_type == "gavage mix" ~ "Gavage",
      experiment_type == "high-fat diet" ~ mouse_n
    )) %>% 
    
    ggplot(., aes(fill = class, x = my.id, y = tpm)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    facet_grid(donor_id ~ my.facet.1 + my.facet.2, space = "free_x", scales = "free_x") +
    geom_text(aes(label = prev, y = label_pos, color = label_col, 
                  size = label_size, vjust = label_vjust)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1), add = 0), 
                       limits = c(0, 1e6),
                       breaks = seq(0, 1e6, 2e5),
                       label = scales::label_number(accuracy = 1, 
                                                    scale_cut = scales::cut_short_scale())) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_size_identity() +
    scale_color_identity() +
    scale_fill_manual(values = c("#9970ABFF", "#5AAE61FF")) +
    
    labs(title = str_glue("{my.group}: {my.rank}")) + xlab("") +  ylab("CPM mapped reads<br>") +
    my.theme.2 + 
    theme(panel.spacing = unit(.2, "lines"), 
          strip.background = element_rect(fill = "grey85", colour = NA),
          strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
          axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)) -> gg.bp
  
  return(gg.bp)
  
}


for (r in ranks) {
  
  png(filename = paste0(path.fg.dir, "/engraftment/core-rare-", str_to_lower(r), ".png"), width = 5000 * 2.5, height = 5000 * .75, units = "px", res = 1200, type = "cairo")
  
  print(
    fun.get.rare.bp(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r, my.group = "Belgian") + 
      fun.get.rare.bp(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r, my.group = "Yanomami") +
      plot_layout(ncol = 2, widths = c(tb.metad %>% filter(donor_ethnicity == "Belgian") %>% nrow(),
                                       tb.metad %>% filter(donor_ethnicity == "Yanomami") %>% nrow()))
  )
  
  dev.off()
  
}; rm(r)

fun.get.rare.bp.core <- function(my.glom.tb.ls, my.metad, my.rank, my.group) {
  
  my.metad = my.metad %>% filter(donor_ethnicity == my.group)
  
  my.top <- fun.get.core(my.glom.tb.ls, my.metad, my.rank, my.group)
  
  my.glom.tb.ls[[my.rank]] %>% 
    mutate(taxa = .[[my.rank]]) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
    mutate(class = ifelse(taxa %in% my.top, "core", "rare")) %>% 
    filter(class == "core") %>% 
    select(taxa, class, any_of(my.metad$sample_id)) %>% 
    pivot_longer(cols = -c(taxa, class), names_to = "sample_id", values_to = "tpm") %>% 
    mutate(prev = ifelse(tpm > 0, 1, 0)) %>% 
    group_by(class, sample_id) %>% 
    summarise(tpm = sum(tpm, na.rm = T),
              prev = sum(prev, na.rm = T), .groups = "drop") %>%
    group_by(sample_id) %>% 
    mutate(perc = prev/length(my.top)) %>% 
    mutate(label_pos = perc/2) %>% 
    ungroup() %>% 
    mutate(label_col = "#033800FF") %>% 
    mutate(label_size = 2) %>% 
    mutate(label_vjust = .5) %>% 
    merge(my.metad, .) %>% 
    mutate(my.facet.1 = case_when(
      experiment_type == "donor" ~ "",
      experiment_type == "gavage mix" ~ "",
      experiment_type == "high-fat diet" ~ "HFD"
    )) %>% 
    mutate(my.facet.2 = case_when(
      experiment_type == "donor" ~ "D",
      experiment_type == "gavage mix" ~ "G",
      experiment_type == "high-fat diet" ~ timepoint
    )) %>% 
    mutate(my.facet.2 = factor(my.facet.2, c("D", "G", "pre-HFD", "post-HFD"))) %>% 
    mutate(my.id = case_when(
      experiment_type == "donor" ~ "Stool",
      experiment_type == "gavage mix" ~ "Gavage",
      experiment_type == "high-fat diet" ~ mouse_n
    )) %>% 
    group_by(my.facet.2) %>% 
    mutate(avg_line_pos = ifelse(experiment_type == "high-fat diet", mean(perc), NA_integer_)) %>%
    mutate(avg_label_pos = ifelse(experiment_type == "high-fat diet", max(perc), NA_integer_)) %>%
    mutate(avg_line_label = ifelse(experiment_type == "high-fat diet", str_glue("Avg.: {round(mean(perc)*100, 2)}%"), NA_character_)) %>%
    
    
    ggplot(., aes(fill = class, x = my.id, y = perc)) + 
    geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
    geom_hline(aes(yintercept = avg_line_pos), linewidth = .5, linetype = "dashed",
               color = "#BE3428FF") +
    facet_grid(donor_id ~ my.facet.1 + my.facet.2, space = "free_x", scales = "free_x") +
    
    geom_text(aes(label = avg_line_label, y = avg_label_pos, x = 1, size = 2.5),
              hjust = .1, vjust = -1, fontface = "bold", color = "#BE3428FF") +
    
    geom_text(aes(label = prev, y = label_pos, color = label_col, 
                  size = label_size, vjust = label_vjust), check_overlap = T) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.2), add = 0), 
                       limits = c(0, 1),
                       breaks = seq(0, 1, .20),
                       labels = scales::percent_format(accuracy = 1, suffix = "")) +
    scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
    scale_size_identity() +
    scale_color_identity() +
    scale_fill_manual(values = c("#5AAE61FF")) +
    
    labs(title = str_glue("{my.group}: {my.rank}")) + xlab("") +  ylab("Retained core taxa (%)<br>") +
    my.theme.2 + 
    theme(panel.spacing = unit(.2, "lines"), 
          strip.background = element_rect(fill = "grey85", colour = NA),
          strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
          axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1)) -> gg.bp; gg.bp
  
  return(gg.bp)
  
}

for (r in ranks) {
  
  png(filename = paste0(path.fg.dir, "/engraftment/core-retained-", str_to_lower(r), ".png"), width = 5000 * 2.5, height = 5000 * .75, units = "px", res = 1200, type = "cairo")
  
  print(
    fun.get.rare.bp.core(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r, my.group = "Belgian") + 
      fun.get.rare.bp.core(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r, my.group = "Yanomami") +
      plot_layout(ncol = 2, widths = c(tb.metad %>% filter(donor_ethnicity == "Belgian") %>% nrow(),
                                       tb.metad %>% filter(donor_ethnicity == "Yanomami") %>% nrow()))
  )
  
  dev.off()
  
}; rm(r)


fun.get.rare.bp.core.compm <- function(my.glom.tb.ls, my.metad, my.rank) {
  
  my.group = "Belgian"
  
  tmp.my.metad = my.metad %>% filter(donor_ethnicity == my.group)
  
  my.top <- fun.get.core(my.glom.tb.ls, tmp.my.metad, my.rank, my.group)
  
  my.glom.tb.ls[[my.rank]] %>% 
    mutate(taxa = .[[my.rank]]) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
    mutate(class = ifelse(taxa %in% my.top, "core", "rare")) %>% 
    filter(class == "core") %>% 
    select(taxa, class, any_of(tmp.my.metad$sample_id)) %>% 
    pivot_longer(cols = -c(taxa, class), names_to = "sample_id", values_to = "tpm") %>% 
    mutate(prev = ifelse(tpm > 0, 1, 0)) %>% 
    group_by(class, sample_id) %>% 
    summarise(tpm = sum(tpm, na.rm = T),
              prev = sum(prev, na.rm = T), .groups = "drop") %>%
    group_by(sample_id) %>% 
    mutate(value = prev/length(my.top)) %>% 
    ungroup() %>% 
    merge(tmp.my.metad, .) %>% 
    data.frame(check.names = F) %>% 
    filter(experiment_type == "high-fat diet") %>% 
    mutate(var = my.rank) -> tmp.tb.engraf.1; tmp.tb.engraf.1
  
  my.group = "Yanomami"
  
  tmp.my.metad = my.metad %>% filter(donor_ethnicity == my.group)
  
  my.top <- fun.get.core(my.glom.tb.ls, tmp.my.metad, my.rank, my.group)
  
  my.glom.tb.ls[[my.rank]] %>% 
    mutate(taxa = .[[my.rank]]) %>% 
    filter(!str_detect(taxa, "Un.\\)")) %>% 
    mutate(class = ifelse(taxa %in% my.top, "core", "rare")) %>% 
    filter(class == "core") %>% 
    select(taxa, class, any_of(tmp.my.metad$sample_id)) %>% 
    pivot_longer(cols = -c(taxa, class), names_to = "sample_id", values_to = "tpm") %>% 
    mutate(prev = ifelse(tpm > 0, 1, 0)) %>% 
    group_by(class, sample_id) %>% 
    summarise(tpm = sum(tpm, na.rm = T),
              prev = sum(prev, na.rm = T), .groups = "drop") %>%
    group_by(sample_id) %>% 
    mutate(value = prev/length(my.top)) %>% 
    ungroup() %>% 
    merge(tmp.my.metad, .) %>% 
    data.frame(check.names = F) %>% 
    filter(experiment_type == "high-fat diet") %>% 
    mutate(var = my.rank) -> tmp.tb.engraf.2; tmp.tb.engraf.2
  
  
  tmp.tb.engraf <- bind_rows(tmp.tb.engraf.1, tmp.tb.engraf.2)
  
  tmp.engraf.summary <- 
    fun.get.summary(tb = tmp.tb.engraf, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value") %>% 
    mutate(sd = ifelse(is.na(sd), 0, sd))
  
  tmp.engraf.kruskal <-
    fun.get.kruskal(tb = tmp.tb.engraf, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value")
  
  tmp.engraf.posthoc <-
    fun.get.wilcox(tb = tmp.tb.engraf, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                   p_adjust_method = "fdr", conf_level = 0.95, tb.stat = tmp.engraf.kruskal, alpha.global = 0.05)
  
  gg.engraf.compm <- 
    fun.get.boxp.compm3(tb = tmp.tb.engraf, 
                        tb.compm = tmp.engraf.posthoc,
                        tb.kruskal = tmp.engraf.kruskal,
                        #tb.summ = tmp.engraf.summary,
                        gp_facet = "var", 
                        gp_comp = "donor_ethnicity_timepoint",
                        gp_fill = "donor_ethnicity_timepoint",
                        value_var = "value",
                        only.sig = T,
                        alpha = 0.1, 
                        label.type = "p.adj", 
                        step.mult = 0.2) + 
    xlab("") + ylab("FMT engraftment rate (%)<br>") +   
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.35)), 
                       labels = scales::percent_format(accuracy = 1, suffix = "")) +
    my.theme.2 +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()); gg.engraf.compm
  
  ls.pcoa <-
    
    list(
      
      "Engraft (Measures)" = tmp.tb.engraf,
      "Engraft (Summary)" = tmp.engraf.summary,
      "Engraft (Kruskal)" = tmp.engraf.kruskal,
      "Engraft (post-hoc)" = tmp.engraf.posthoc
      
    )
  
  openxlsx::write.xlsx(x = ls.pcoa, file = paste0(path.tb.dir, "/engraftment/engraftment-rate_", str_to_lower(my.rank), ".xlsx"), rowNames = F)
  
  return(gg.engraf.compm)
  
}

ls.engraf.gg <- list()

for (r in ranks) {
  
  ls.engraf.gg[[r]] <- fun.get.rare.bp.core.compm(my.glom.tb.ls = ls.glom.tb, my.metad = tb.metad, my.rank = r)
  
}; rm(r)

png(filename = paste0(path.fg.dir, "/engraftment/core-retained-compm.png"), width = 5000 * 2.5, height = 5000 * .75, units = "px", res = 1200, type = "cairo")

print(wrap_plots(unlist(ls.engraf.gg)) +plot_layout(nrow = 1))

dev.off()


## Profiles --------------------------------------------------------------------


top.phylum <- c("Bacillota", "Bacteroidota", "Verrucomicrobiota", "Pseudomonadota", "Actinomycetota")


### Samples --------------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/profiles"))) { dir.create(str_glue("{path.fg.dir}/profiles"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/profiles"))) { dir.create(str_glue("{path.tb.dir}/profiles"), recursive = T) }

gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.tb.glom <- ls.glom.tb[[tmp.rank]]
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
    merge(x = tb.metad, y = reshape2::melt(.), by.x = "sample_id", by.y = "variable") %>% 
    arrange(., taxa2, -value) %>% 
    
    mutate(my.facet.1 = case_when(
      experiment_type == "donor" ~ "",
      experiment_type == "gavage mix" ~ "",
      experiment_type == "high-fat diet" ~ "HFD"
    )) %>% 
    mutate(my.facet.2 = case_when(
      experiment_type == "donor" ~ "D",
      experiment_type == "gavage mix" ~ "G",
      experiment_type == "high-fat diet" ~ timepoint
    )) %>% 
    mutate(my.facet.2 = factor(my.facet.2, c("D", "G", "pre-HFD", "post-HFD"))) %>% 
    mutate(my.id = case_when(
      experiment_type == "donor" ~ "Stool",
      experiment_type == "gavage mix" ~ "Gavage",
      experiment_type == "high-fat diet" ~ mouse_n
    )) %>% 
    arrange(donor_ethnicity, my.facet.2) -> tmp.tb.glom.top
  
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
  
  ls.gg.bp <- list()
  ls.gg.wd <- list()
  
  for (d in unique(tmp.tb.glom.top$donor_id)) {
    
    ls.gg.bp[[d]] <- 
      ggplot(tmp.tb.glom.top %>% filter(donor_id == d), aes(fill = taxa2, x = my.id, y = value)) + 
      geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
      facet_grid(donor_id ~ my.facet.1 + my.facet.2, space = "free_x", scales = "free_x") +
      scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
      scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
      scale_fill_manual(values = tmp.colors) +  
      labs(fill = tmp.rank, title = str_glue("{tmp.tb.glom.top %>% filter(donor_id == d) %>% pull(donor_ethnicity) %>% unique()}")) + xlab("") + ylab("Relative abundance (%)<br>") +
      guides(fill=guide_legend(ncol = ifelse(length(unique(tmp.tb.glom.top$taxa)) > 10, 1, 1))) +
      my.theme.2 + 
      theme(panel.spacing = unit(.1, "lines"), 
            strip.background = element_rect(fill = "grey85", colour = NA),
            strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
            axis.text.x = element_markdown(size = 6, colour = 'black', angle = 75, hjust = 1, vjust = 1))
    
    ls.gg.wd[[d]] <- 
      tmp.tb.glom.top %>% 
      filter(donor_id == d) %>% 
      mutate(tmp = paste(my.facet.2, my.id)) %>% 
      pull(tmp) %>% 
      unique() %>% 
      length()
    
  }
  
  gg.bp <- wrap_plots(unlist(ls.gg.bp)) + plot_layout(nrow = 1); gg.bp
  
  gg.bp.plots[[i]] <- gg.bp + gg.bp.leg + plot_layout(nrow = 1, widths = c(c(unlist(ls.gg.wd)), (sum(unlist(ls.gg.wd)) * 70) / 30))
  
}; rm(i, d)

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p
gg.bp.c <- gg.bp.plots[[2]]; gg.bp.c
gg.bp.o <- gg.bp.plots[[3]]; gg.bp.o
gg.bp.f <- gg.bp.plots[[4]]; gg.bp.f
gg.bp.g <- gg.bp.plots[[5]]; gg.bp.g

ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_1p.png"), plot = gg.bp.p, dpi = 1200, width = 5000 * 4, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_2c.png"), plot = gg.bp.c, dpi = 1200, width = 5000 * 4, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_3o.png"), plot = gg.bp.o, dpi = 1200, width = 5000 * 4, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_4f.png"), plot = gg.bp.f, dpi = 1200, width = 5000 * 4, height = 5000 * 0.85, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-sample_5g.png"), plot = gg.bp.g, dpi = 1200, width = 5000 * 4, height = 5000 * 0.85, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))


### Group ----------------------------------------------------------------------


gg.bp.plots <- list()

for (i in seq_along(ranks)) {
  
  tmp.rank <- ranks[[i]]
  
  tmp.tb.glom <- ls.glom.tb[[tmp.rank]]
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
    merge(x = tb.metad, y = reshape2::melt(.), by.x = "sample_id", by.y = "variable") %>% 
    arrange(., taxa2, -value) %>% 
    
    mutate(my.facet.1 = case_when(
      experiment_type == "donor" ~ "",
      experiment_type == "gavage mix" ~ "",
      experiment_type == "high-fat diet" ~ "HFD"
    )) %>% 
    mutate(my.facet.2 = case_when(
      experiment_type == "donor" ~ "D",
      experiment_type == "gavage mix" ~ "G",
      experiment_type == "high-fat diet" ~ paste0(str_remove_all(sample_id, ".*_"))
    )) %>% 
    mutate(my.facet.2 = factor(my.facet.2, c("D", "G", "T1", "T2"))) %>% 
    mutate(my.id = case_when(
      experiment_type == "donor" ~ "Stool",
      experiment_type == "gavage mix" ~ "Gavage",
      experiment_type == "high-fat diet" ~ mouse_n
    )) %>% 
    arrange(donor_ethnicity, my.facet.2) -> tmp.tb.glom.top
  
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
  
  ls.gg.bp <- list()
  ls.gg.wd <- list()
  
  for (d in unique(tmp.tb.glom.top$donor_id)) {
    
    ls.gg.bp[[d]] <-
      ggplot(tmp.tb.glom.top %>% filter(donor_id == d), 
             aes(fill = taxa2, x = my.facet.2, y = value)) + 
      geom_bar(position = position_stack(reverse = F), stat = 'summary', fun = 'mean', width = .95) +
      facet_grid(donor_id ~ my.facet.2, space = "free", scales = "free") +
      scale_y_continuous(expand = expansion(mult = 0.0025, add = 0), labels = scales::percent_format(accuracy = 1, suffix = "")) +
      scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
      scale_fill_manual(values = tmp.colors) +  
      labs(fill = tmp.rank, title = str_glue("{tmp.tb.glom.top %>% filter(donor_id == d) %>% pull(donor_ethnicity) %>% unique()}")) + xlab("") + ylab("Relative abundance (%)<br>") +
      guides(fill=guide_legend(ncol = ifelse(length(unique(tmp.tb.glom.top$taxa2)) > 10, 1, 1))) +
      my.theme.2 + 
      theme(panel.spacing = unit(.1, "lines"), 
            strip.background = element_rect(fill = "grey85", colour = NA),
            strip.text = element_markdown(size = 8, colour = 'black', face = "bold"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    
    ls.gg.wd[[d]] <- 
      tmp.tb.glom.top %>% 
      filter(donor_id == d) %>% 
      pull(my.facet.2) %>% 
      unique() %>% 
      length()
    
  }
  
  gg.bp <- wrap_plots(unlist(ls.gg.bp)) + plot_layout(nrow = 1); gg.bp
  
  gg.bp.plots[[i]] <- gg.bp + gg.bp.leg + plot_layout(nrow = 1, widths = c(c(unlist(ls.gg.wd)), (sum(unlist(ls.gg.wd)) * 80) / 20))
  
}; rm(i, d)

gg.bp.p <- gg.bp.plots[[1]]; gg.bp.p
gg.bp.c <- gg.bp.plots[[2]]; gg.bp.c
gg.bp.o <- gg.bp.plots[[3]]; gg.bp.o
gg.bp.f <- gg.bp.plots[[4]]; gg.bp.f
gg.bp.g <- gg.bp.plots[[5]]; gg.bp.g


ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_1p.png"), plot = gg.bp.p, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.75, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_2c.png"), plot = gg.bp.c, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.75, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_3o.png"), plot = gg.bp.o, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.75, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_4f.png"), plot = gg.bp.f, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.75, units = "px")
ggsave(filename = paste0(path.fg.dir, "/profiles/bp-tax-top-group_5g.png"), plot = gg.bp.g, dpi = 1200, width = 5000 * 3.5, height = 5000 * 0.75, units = "px")

rm(list = ls(pattern = "^tmp."))
rm(list = ls(pattern = "^gg.bp[.]?"))

#openxlsx::write.xlsx(x = ls.glom.tb.filt, file = paste0(path.tb.dir, "/profiles/taxa-glom-filt-table.xlsx"), rowNames = F)


## Beta diversity --------------------------------------------------------------


if (!dir.exists(str_glue("{path.fg.dir}/diversity"))) { dir.create(str_glue("{path.fg.dir}/diversity"), recursive = T) }
if (!dir.exists(str_glue("{path.tb.dir}/diversity"))) { dir.create(str_glue("{path.tb.dir}/diversity"), recursive = T) }

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
    
    pcoa.perma.comp1 <- adonis2(pcoa.dist ~ donor_ethnicity + timepoint, data = my.metad, method = bdiv.indexes[[i]], permutations = 999, seed = 123)
    
    pcoa.perma.df.comp1 <- 
      pcoa.perma.comp1 %>% 
      data.frame(check.names = F) %>% 
      rownames_to_column() %>% 
      dplyr::rename("Source" = 1, "DF" = 2, "SS" = 3, "F" = 5, "P" = 6) %>% 
      mutate(Source = str_to_title(Source),
             sig = case_when(P <= 0.001 ~ "***", P <= 0.01 ~ "**", P <= 0.05 ~ "*", P <= 0.1 ~ ".", .default = " "))
    
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
    
    pcoa.perma.ph <- bind_rows(pcoa.perma.ph.comp1.2)
    
    # Axes 1:2
    
    ggplot() +
      # stat_ellipse(data = pcoa.site.score,
      #              aes(x = PCoA1, y = PCoA2, color = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint),
      #              geom = "polygon", linetype = "blank", alpha = .15, show.legend = F, type = "norm") +
      geom_point(pcoa.site.score, 
                 mapping = aes(x = PCoA1, y = PCoA2, color = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint, shape = experiment_type), 
                 size = 2, alpha = .75) +
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
      #              aes(x = PCoA1, y = PCoA3, color = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint),
      #              geom = "polygon", linetype = "blank", alpha = .15, show.legend = F, type = "norm") +
      geom_point(pcoa.site.score, 
                 mapping = aes(x = PCoA1, y = PCoA3, color = donor_ethnicity_timepoint, fill = donor_ethnicity_timepoint, shape = experiment_type), 
                 size = 2, alpha = .75) +
      # annotate("richtext", x = Inf, y = -Inf, size = 2.5,
      #          label = paste0("PERMANOVA: *R*<sup>2</sup> = ", round(pcoa.perma.df.comp1[1, "R2"], 2), "; *p* = ", ifelse(pcoa.perma.df.comp1[1, "P"] < 0.001, "<0.001", sprintf("%1.3f", round(pcoa.perma.df.comp1[1, "P"], 3)))),
      #          fill = NA, label.color = NA, hjust = 1, vjust = 0) +
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
      betadisper(d = pcoa.dist, group = my.metad$donor_ethnicity_timepoint, type = "centroid")$distances %>% 
      data.frame(check.names = F) %>% 
      dplyr::rename("value" = ".") %>% 
      rownames_to_column("sample_id") %>% 
      merge(my.metad, .) %>% 
      mutate(var = "Within-group\ndispersion")
    
    pcoa.disp.summary.comp2 <- 
      fun.get.summary(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value")
    
    pcoa.disp.kruskal.comp2 <-
      fun.get.kruskal(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value")
    
    pcoa.disp.posthoc.comp2 <-
      fun.get.wilcox(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity_timepoint", gp_facet = "var", value_var = "value",
                     p_adjust_method = "fdr", conf_level = 0.95, tb.stat = pcoa.disp.kruskal.comp2, alpha.global = 0.05)
    
    # pcoa.disp.compm.comp2 <-
    #   fun.get.wilcox(tb = pcoa.disp.comp2, gp_comp = "donor_ethnicity", gp_facet = "var", value_var = "value",
    #                  p_adjust_method = "fdr", conf_level = 0.95, tb.stat = NULL, alpha.global = 0.05)
    
    gg.pcoa.disp.compm.comp2 <- 
      fun.get.boxp.compm3(tb = pcoa.disp.comp2, 
                          #tb.compm = pcoa.disp.compm.comp2,
                          tb.compm = pcoa.disp.posthoc.comp2,
                          tb.kruskal = pcoa.disp.kruskal.comp2,
                          gp_facet = "var", 
                          gp_comp = "donor_ethnicity_timepoint",
                          gp_fill = "donor_ethnicity_timepoint",
                          gp_shape = "experiment_type",
                          value_var = "value",
                          only.sig = T,
                          alpha = 0.1, 
                          label.type = "p.adj", 
                          step.mult = 0.1)
    
    gg.pcoa.panel <- 
      gg.pcoa.pc1.compm.axes12 + plot_spacer() + plot_spacer() + gg.pcoa.pc1.compm.axes13 + plot_spacer() + plot_spacer() +
      gg.pcoa.axes12 + gg.pcoa.pc2.compm.axes12 + plot_spacer() + gg.pcoa.axes13 + gg.pcoa.pc2.compm.axes13 + gg.pcoa.disp.compm.comp2 + plot_layout(ncol = 6, heights = c(.1, 1), widths = c(1, .1, .05, 1, .1, .5))
    
    print(gg.pcoa.panel)
    
    
    #Saving graphs
    png(filename = paste0(path.fg.dir, "/diversity/beta-div-", str_to_lower(bdiv.indexes[[i]]), ".png"), width = 5000 * 2.25, height = 5000 * 1, units = "px", res = 1200, type = "cairo")
    
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
        "Dispersion (Kruskal)" = pcoa.disp.kruskal.comp2,
        "Dispersion (post-hoc)" = pcoa.disp.posthoc.comp2
        
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

# End ---------------------------------------------------------------------

