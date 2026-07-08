## ============================================================
## fig_urbanisation_by_phase.R
## Phase-stratified urbanisation composition of first-invaded
## regions (green urbanisation palette; year-range x-axis).
## ============================================================

library(dplyr)
library(ggplot2)
library(ggpubr)     # theme_classic2()

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!exists("analysis_cases"))
  analysis_cases <- readRDS("analysis_cases_with_cats.rds")

urbn_labels <- c("1" = "Urban", "2" = "Intermediate", "3" = "Rural")
phase_levels <- c("Early", "Middle", "Pre-COVID", "Post-COVID")

## Green urbanisation palette, consistent with the main-text encoding.
green3 <- c("Urban" = "#006D2C", "Intermediate" = "#41AB5D", "Rural" = "#A1D99B")

## First-invaded regions with urbanisation and phase of first outbreak.
first_invaded <- analysis_cases %>%
  filter(event == 1, !is.na(URBN_TYPE)) %>%
  mutate(URBN  = factor(urbn_labels[as.character(URBN_TYPE)],
                        levels = c("Urban", "Intermediate", "Rural")),
         phase = factor(invasion_phase, levels = phase_levels))

## Phase -> "Name (year0-year1)" labels, computed from the data.
phase_years <- analysis_cases %>%
  filter(!is.na(invasion_phase)) %>%
  transmute(phase = factor(invasion_phase, levels = phase_levels), year) %>%
  group_by(phase) %>%
  summarise(y0 = min(year), y1 = max(year), .groups = "drop") %>%
  # Early is labelled by its definitional lower bound (study start, 2006);
  # the modelling data begin in 2007, so no 2006 events enter this analysis.
  mutate(y0 = ifelse(as.character(phase) == "Early", 2006, y0),
         lab = paste0(as.character(phase), "\n(", y0, "\u2013", y1, ")"))
phase_labels <- setNames(phase_years$lab, as.character(phase_years$phase))

## Phase x urbanisation table and association test.
tab <- table(first_invaded$phase, first_invaded$URBN)
cat("Phase x urbanisation counts:\n"); print(tab)
cs      <- suppressWarnings(chisq.test(tab))
min_exp <- min(cs$expected)
if (min_exp < 5) {
  test <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 1e4))
  test_name <- "Monte Carlo chi-square test"
} else {
  test <- cs; test_name <- sprintf("Chi-square test (df = %d)", cs$parameter)
}
p_lab <- ifelse(test$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", test$p.value))

## Proportional composition within each phase.
comp <- first_invaded %>%
  count(phase, URBN) %>%
  group_by(phase) %>%
  mutate(prop = n / sum(n), n_phase = sum(n)) %>%
  ungroup()
phase_totals <- comp %>% distinct(phase, n_phase)

Fig_Urban_Phase <-
  ggplot(comp, aes(x = phase, y = 100 * prop, fill = URBN)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5), size = 3, colour = "white") +
  geom_text(data = phase_totals, aes(x = phase, y = 103,
            label = paste0("n=", n_phase)), inherit.aes = FALSE, size = 3.2) +
  scale_fill_manual(values = green3, name = "Urbanisation") +
  scale_x_discrete(labels = phase_labels) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  annotate("text", x = 0.55, y = 108,
           label = paste0(sub(" \\(.*", "", test_name), ": ", p_lab),
           hjust = 0, size = 3.6) +
  labs(x = "Invasion phase", y = "Composition of first-invaded regions (%)") +
  theme_classic2(base_size = 13) +
  theme(legend.position = "top")

Fig_Urban_Phase
ggsave(file.path(output_dir, "Fig_Urbanisation_by_phase.pdf"),
       Fig_Urban_Phase, width = 8, height = 5)
