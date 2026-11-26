plotEnrichment_new <- function(pathway, stats, gseaParam = 1, ticksSize = 0.2, spreadES1 = 1, lsize = 1){
  pd <- plotEnrichmentData(pathway = pathway, stats = stats,
                           gseaParam = gseaParam)
  with(pd, ggplot(data = curve) + geom_line(aes(x = rank, y = ES), size = lsize,
                                            color = "green") + geom_segment(data = ticks, mapping = aes(x = rank,
                                                                                                        y = -spreadES1/16, xend = rank, yend = spreadES1/16), linewidth = ticksSize) +
         geom_hline(yintercept = posES, colour = "red", linetype = "dashed") +
         geom_hline(yintercept = negES, colour = "red", linetype = "dashed") +
         geom_hline(yintercept = 0, colour = "black") + theme(panel.background = element_blank(),
                                                              panel.grid.major = element_line(color = "grey92")) +
         labs(x = "rank", y = "enrichment score"))
}
