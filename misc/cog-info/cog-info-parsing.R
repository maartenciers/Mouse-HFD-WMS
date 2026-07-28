
# COG info parsing --------------------------------------------------------

# Setup ------------------------------------------------------------------------

rm(list = ls())

library("tidyverse")

getwd()
setwd("~/YanomamiMicrobiome/misc/cog-info/")
getwd()

set.seed("123")

list.files()

# COG letters
cog.fun.tab <- read.delim(file = "cog-24.fun.tab", header = F, col.names = c("cog.letter", "cog.class", "cog.color", "cog.cat"))

cog.fun.tab.1 <- 
  cog.fun.tab %>% 
  filter(str_detect(cog.letter, pattern = "[:digit:]")) %>% 
  select("cog.class.n" = 1, "cog.class" = 2)

cog.fun.tab.2 <-
  cog.fun.tab %>% 
  filter(str_detect(cog.letter, pattern = "[:LETTER:]")) %>% 
  mutate(cog.class = case_when(
    cog.class == 1 ~ cog.fun.tab.1 %>% filter(cog.class.n == 1) %>% pull(cog.class),
    cog.class == 2 ~ cog.fun.tab.1 %>% filter(cog.class.n == 2) %>% pull(cog.class),
    cog.class == 3 ~ cog.fun.tab.1 %>% filter(cog.class.n == 3) %>% pull(cog.class),
    cog.class == 4 ~ cog.fun.tab.1 %>% filter(cog.class.n == 4) %>% pull(cog.class))) %>% 
  mutate(cog.color = str_replace(cog.color, "Mobilome: prophages, transposons", "9CFC8C")) %>% 
  mutate(cog.color = str_glue("#{cog.color}")) %>% 
  mutate(cog.cat = str_replace(cog.cat, "^$", "Mobilome: prophages, transposons"))


# COG families
cog.def.tab <- read.delim(file = "cog-24.def.tab", header = F, col.names = c("cog.id", "cog.letter", "cog.name", "cog.gene.name", "cog.path", "cog.pubmed.id", "cog.pdb.id"))
cog.map.tab <- read.delim(file = "cog-24.mapping.tab", header = T, col.names = c("cog.id", "cog.letter", "cog.gene.name", "cog.name.short", "cog.uniprot.id", "cog.uniprot.name", "cog.aa.len"))

cog.def.map.tab <- merge(cog.def.tab, cog.map.tab, all.x = T)

# COG pathway
cog.path.tab <- read.delim(file = "cog-24.pathways.tab", header = F, col.names = c("cog.path", "cog.id", "cog.letter", "cog.gene.name", "cog.name", "cog.ec"))

cog.path.tab.1 <- cog.path.tab %>% filter(cog.path != "")

# Saving

write.table(x = cog.fun.tab.2, file = "cog-info-class.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = cog.def.map.tab, file = "cog-info-def.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(x = cog.path.tab.1, file = "cog-info-pw.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
