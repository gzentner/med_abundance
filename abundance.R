library(readxl)
library(tidyverse)

med <- read_tsv("core_mediator_complex_annotations.txt")
abundance <- read_xlsx("1-s2.0-S240547121730546X-mmc5.xlsx", skip = 1)

abundance_trimmed <- select(abundance, -c(1,3:6))

pivoted <- pivot_longer(abundance_trimmed, -c(1), names_to = "study", values_to = "quant") %>%
    filter(`Standard Name` %in% med$gene) %>%
    left_join(., med, by = c("Standard Name" = "gene")) %>%
    mutate(module = factor(module, levels = c("scaffold",
                                              "head",
                                              "middle",
                                              "tail"))) 
gene_order <- c("RGR1", 
                "MED6", "MED8", "MED11", "SRB4", "SRB5", "SRB2", "SRB6",
                "MED1", "MED4", "NUT1", "MED7", "CSE2", "NUT2", "ROX3", "SRB7", "SOH1",
                "MED2", "PGD1", "GAL11", "SIN4")

aka_order <- c("Med14",
               "Med6", "Med8", "Med11", "Med17", "Med18", "Med20", "Med22",
               "Med1", "Med4", "Med5", "Med7", "Med9", "Med10", "Med19", "Med21", "Med31",
               "Med2", "Med3", "Med15", "Med16")

pivoted$aka <- factor(pivoted$aka, levels = aka_order)

p <- ggplot(pivoted, aes(x = aka, y = quant, fill = module)) +
    stat_boxplot(geom = "errorbar", width = 0.4, na.rm = T, color = "black") +
    geom_jitter(position = position_jitter(height = 0, 
                                           width = 0.3), 
                size = 1, 
                na.rm = T) +
    geom_boxplot(width = 0.75, 
                 alpha = 0.5, 
                 outlier.size = 1, 
                 outlier.alpha = 1,
                 outlier.color = "black",
                 color = "black",
                 na.rm = T) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(color = "black", size = 16),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18)) +
    xlab("") +
    ylab("Molecules per cell") +
    scale_fill_viridis_d() 

ggsave(p, filename = "abundance_boxplot.png", width = 9, height = 6, units = "in")
