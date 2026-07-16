## ============================================================
## fig_urbanisation_firstinvaded.R
## Supplementary Figure S7: Urbanisation composition of
## first-invaded regions vs the at-risk background, with a test
## of association.
## ============================================================

library(dplyr)
library(ggplot2)
library(ggpubr)     # theme_classic2()

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!exists("analysis_cases"))
  analysis_cases <- readRDS("analysis_cases_with_cats.rds")

## Eurostat DEGURBA coding: 1 = cities (urban), 2 = towns/suburbs
## (intermediate), 3 = rural.
urbn_labels <- c("1" = "Urban", "2" = "Intermediate", "3" = "Rural")

## One row per region: its urbanisation type and whether it was ever
## first-invaded during the study window. (Each region has at most one
## event == 1, its first outbreak.)
region_level <- analysis_cases %>%
  group_by(NUTS_ID) %>%
  summarise(URBN_TYPE = dplyr::first(URBN_TYPE),
            invaded   = as.integer(any(event == 1)),
            .groups   = "drop") %>%
  filter(!is.na(URBN_TYPE)) %>%
  mutate(URBN = factor(urbn_labels[as.character(URBN_TYPE)],
                       levels = c("Urban", "Intermediate", "Rural")))

n_dropped <- sum(is.na(analysis_cases %>%
                         distinct(NUTS_ID, URBN_TYPE) %>% pull(URBN_TYPE)))
if (n_dropped > 0)
  message(n_dropped, " region(s) dropped for missing URBN_TYPE.")

## Association test: does urbanisation differ between invaded and
## non-invaded regions? (2 x 3 contingency table.)
tab  <- table(region_level$URBN, region_level$invaded)
gof  <- suppressWarnings(chisq.test(tab))
p_lab <- ifelse(gof$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", gof$p.value))

## Composition: first-invaded vs full at-risk background.
comp <- bind_rows(
  region_level %>% filter(invaded == 1) %>% count(URBN) %>%
    mutate(Group = "First-invaded regions"),
  region_level %>% count(URBN) %>%
    mutate(Group = "All at-risk regions (background)")
) %>%
  group_by(Group) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(Group = factor(Group, levels = c("First-invaded regions",
                                           "All at-risk regions (background)")))

Fig_Urbanisation <-
  ggplot(comp, aes(x = URBN, y = 100 * prop, fill = Group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = paste0("n=", n)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("First-invaded regions" = "#CB181D",
                               "All at-risk regions (background)" = "#BDBDBD")) +
  annotate("text", x = 0.6, y = max(100 * comp$prop) * 1.12,
           label = paste0("Chi-square test: ", p_lab), hjust = 0, size = 3.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Urbanisation category", y = "Percentage of regions (%)",
       fill = NULL) +
  theme_classic2(base_size = 13) +
  theme(legend.position = "top")

Fig_Urbanisation
ggsave(file.path(output_dir, "FigureS7_Urbanisation_firstinvaded.pdf"),
       Fig_Urbanisation, width = 7, height = 5)

## Console summary
cat("\nUrbanisation of first-invaded regions vs background:\n")
print(comp %>% mutate(prop = round(100 * prop, 1)) %>%
        select(Group, URBN, n, `%` = prop))
print(gof)
