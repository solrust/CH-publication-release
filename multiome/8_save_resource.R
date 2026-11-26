library(dplyr)
library(tidyr)

hspc.labels <- c(`1`="HSC", `2`="LMPP-1", `3`="MLP", `4`="Pro-B", `5`="pDC",
                 `6`="LMPP-2", `7`="early GMP", `8`="late GMP", `9`="EMPP/MkP", `10`="EryP")

meta <- read.csv("data/processed/metainfo.txt", sep = '\t')
head(meta)

meta2 <- meta %>%
  mutate(batch=sample,
         sample=paste(batch, genotype, sep="_"),
         population=hspc.labels[as.character(projection)])

cell.summary.multiome <- meta2 %>%
  filter(!is.na(projection)) %>%
  count(batch, sample, genotype, projection, population, name="n_cells") %>%
  complete(nesting(batch, sample, genotype),
           projection=as.numeric(names(hspc.labels)),
           fill=list(n_cells=0)) %>%
  mutate(population=hspc.labels[as.character(projection)]) %>%
  group_by(batch, sample, genotype) %>%
  mutate(total_cells_sample=sum(n_cells),
         fraction_of_sample=n_cells/total_cells_sample) %>%
  ungroup() %>%
  arrange(batch, genotype, projection)

cell.summary.multiome
cell.summary.multiome <- cell.summary.multiome %>%
  select(batch, sample, genotype, cluster=projection, population,
         n_cells, total_cells_sample, fraction_of_sample)

write.csv(cell.summary.multiome, "../resource/multiome_cell_counts.csv", row.names = FALSE)
