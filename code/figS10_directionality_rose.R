###############################################################################
# figS10_directionality_rose.R
#
# Supplementary Figure S10. Direction of first invasion relative to the nearest
# previously invaded region (response to Reviewer 3.2).
#
#   A  16-sector rose of the initial great-circle bearing from each event's
#      nearest source to the newly invaded region, every first invasion from
#      2007 to 2026, with the mean resultant vector.
#   B  The same rose within each display period: 2006-2010, 2011-2015,
#      2016-2020 and 2021-2026.
#
# This is descriptive, so the 2025 and 2026 seasons belong here: no model is
# fitted. The three regions invaded in 2006 have no previously invaded source and
# so no bearing; the 2006-2010 panel therefore holds the 2007-2010 events.
#
# Bearings are computed exactly as in 08_revision_analyses.R, section E
# (centroids in EPSG:3035, transformed to WGS84, initial great-circle bearing),
# and the fitted-period bearings are checked against rev_event_bearings.csv.
# 2025 and 2026 sources come from the prospective scoring tables, where the
# nearest source was taken from the invasion history before each season.
#
# Writes rev_directionality_by_period.csv and rev_directionality_rose_2007_2026.csv
# for the appendix tables.
#
#   bash run_fig.sh s10
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("figS10")
on.exit(lg$close(), add = TRUE)

SECTORS <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
             "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
PERIODS <- c("2006-2010", "2011-2015", "2016-2020", "2021-2026")
period4 <- function(y) as.character(cut(y, c(2005, 2010, 2015, 2020, 2026), labels = PERIODS))

## ---- Centroids and bearings, as in 08 section E -------------------------------
cent_ll <- suppressWarnings(nuts %>% st_transform(3035) %>% st_make_valid() %>%
                              st_centroid() %>% st_transform(4326))
cent_df <- data.frame(NUTS_ID = cent_ll$NUTS_ID, st_coordinates(cent_ll)) %>%
  rename(lon = X, lat = Y)

bearing_deg <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  dl <- (lon2 - lon1) * r
  y <- sin(dl) * cos(lat2 * r)
  x <- cos(lat1 * r) * sin(lat2 * r) - sin(lat1 * r) * cos(lat2 * r) * cos(dl)
  (atan2(y, x) / r + 360) %% 360
}

events <- bind_rows(
  ac %>% filter(event == 1) %>% transmute(NUTS_ID, year, nearest_source),
  inv25 %>% transmute(NUTS_ID, year = 2025L, nearest_source),
  inv26 %>% transmute(NUTS_ID, year = 2026L, nearest_source)
) %>%
  left_join(cent_df, by = "NUTS_ID") %>% rename(lon_to = lon, lat_to = lat) %>%
  left_join(cent_df %>% rename(nearest_source = NUTS_ID, lon_from = lon, lat_from = lat),
            by = "nearest_source") %>%
  mutate(bearing = bearing_deg(lon_from, lat_from, lon_to, lat_to),
         sector = SECTORS[((floor((bearing + 11.25) / 22.5)) %% 16) + 1],
         period = period4(year))

n_missing <- sum(is.na(events$bearing))
message("first invasions 2007-2026: ", nrow(events), "; without a computable bearing: ", n_missing)
stopifnot(n_missing == 0, !anyNA(events$period))

## Fitted-period bearings must reproduce 08.
ref <- read.csv(file.path(REV, "rev_event_bearings.csv"), stringsAsFactors = FALSE)
chk <- inner_join(events %>% filter(year <= 2024) %>% select(NUTS_ID, bearing),
                  ref %>% select(NUTS_ID, bearing_ref = bearing), by = "NUTS_ID")
stopifnot(nrow(chk) == nrow(ref), max(abs(chk$bearing - chk$bearing_ref)) < 0.1)
message("fitted-period bearings match 08 for all ", nrow(chk), " events")

rose_counts <- function(d) {
  d %>%
    count(sector, name = "n") %>%
    right_join(data.frame(sector = SECTORS), by = "sector") %>%
    mutate(n = ifelse(is.na(n), 0L, n),
           angle = (match(sector, SECTORS) - 1) * 22.5,
           pct = 100 * n / sum(n))
}

## Rayleigh test, the same function as 08.
rayleigh <- function(theta_deg) {
  th <- theta_deg * pi / 180
  n  <- length(th)
  C  <- sum(cos(th)); S <- sum(sin(th))
  R  <- sqrt(C^2 + S^2)
  data.frame(n = n,
             mean_bearing_deg = round((atan2(S, C) * 180 / pi + 360) %% 360, 1),
             resultant_length = round(R / n, 3),
             rayleigh_z = round(R^2 / n, 3),
             p_value = signif(exp(sqrt(1 + 4 * n + 4 * (n^2 - R^2)) - (1 + 2 * n)), 3))
}

dirtab <- bind_rows(
  rayleigh(events$bearing) %>% mutate(stratum = "All first invasions 2007-2026", .before = 1),
  bind_rows(lapply(PERIODS, function(p)
    rayleigh(events$bearing[events$period == p]) %>% mutate(stratum = p, .before = 1))))
print(dirtab, row.names = FALSE)
write.csv(dirtab, file.path(REV, "rev_directionality_by_period.csv"), row.names = FALSE)
write.csv(rose_counts(events) %>% select(sector, n, pct) %>% mutate(pct = round(pct, 1)),
          file.path(REV, "rev_directionality_rose_2007_2026.csv"), row.names = FALSE)

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

## Rose builder. Every aesthetic is a column of the layer's own data.
rose_plot <- function(counts, stats, value = c("n", "pct")) {
  value <- match.arg(value)
  counts$y <- counts[[value]]
  ymax <- max(counts$y) * 1.08
  arrow_df <- data.frame(x = stats$mean, y0 = 0, y1 = stats$R * ymax)
  rings <- data.frame(y = pretty(c(0, ymax), n = 3)) %>% filter(y > 0, y <= ymax)
  rings$x <- 90

  ggplot(counts) +
    geom_hline(data = rings, aes(yintercept = y), colour = RULE, linewidth = 0.25) +
    geom_col(aes(x = angle, y = y), width = 22.5, fill = TEAL, colour = "white",
             linewidth = 0.3) +
    geom_segment(data = arrow_df, aes(x = x, xend = x, y = y0, yend = y1),
                 colour = INK, linewidth = 0.55,
                 arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")) +
    geom_text(data = rings, aes(x = x, y = y, label = y), colour = MID,
              size = 1.9, vjust = -0.3, hjust = -0.1, family = FONT) +
    scale_x_continuous(limits = c(-11.25, 348.75),
                       breaks = c(0, 90, 180, 270),
                       labels = c("N", "E", "S", "W")) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    coord_polar(start = -11.25 * pi / 180, direction = 1) +
    theme_wnv_plot() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(colour = INK, size = 7.5),
          panel.grid.major.x = element_line(colour = RULE, linewidth = 0.2),
          plot.title = element_text(size = 7.5, colour = INK, hjust = 0.5),
          plot.subtitle = element_text(size = 6.2, colour = MID, hjust = 0.5))
}

## Panel A, every first invasion 2007-2026, counts.
row_all <- dirtab[1, ]
pA <- rose_plot(rose_counts(events),
                list(mean = row_all$mean_bearing_deg, R = row_all$resultant_length), "n") +
  labs(title = sprintf("All first invasions, 2007-2026 (n = %d)", row_all$n),
       subtitle = sprintf("Mean bearing %.1f°, resultant length %.3f, Rayleigh %s",
                          row_all$mean_bearing_deg, row_all$resultant_length,
                          fmt_p(row_all$p_value)))

## Panel B, by display period, percentages within period.
period_plots <- lapply(PERIODS, function(p) {
  row <- dirtab[dirtab$stratum == p, ]
  rose_plot(rose_counts(events[events$period == p, ]),
            list(mean = row$mean_bearing_deg, R = row$resultant_length), "pct") +
    labs(title = sprintf("%s (n = %d)", p, row$n),
         subtitle = sprintf("%.0f°, R = %.3f, %s", row$mean_bearing_deg,
                            row$resultant_length, fmt_p(row$p_value)))
})
pB <- wrap_plots(period_plots, ncol = 2)

figS10 <- wrap_elements(pA) + wrap_elements(pB) +
  plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = list(c("A", "B")),
                  theme = theme(plot.margin = margin(4, 4, 4, 4))) &
  theme(plot.tag = element_text(size = 10, colour = INK, face = "plain"))

save_fig(figS10, "FigureS10_directionality_rose", WIDTH_2COL, 105)

## Quadrant counts for the caption, same rule as 08.
events <- events %>% mutate(quadrant = case_when(
  bearing >= 315 | bearing < 45  ~ "Northward",
  bearing >= 45  & bearing < 135 ~ "Eastward",
  bearing >= 135 & bearing < 225 ~ "Southward",
  TRUE ~ "Westward"))
quad <- events %>% count(quadrant) %>% mutate(pct = 100 * n / sum(n))
qtxt <- function(q) sprintf("%d (%.1f%%)", quad$n[quad$quadrant == q], quad$pct[quad$quadrant == q])
per <- function(p) { r <- dirtab[dirtab$stratum == p, ]
  sprintf("%s, n = %d, mean bearing %.0f degrees, R = %.3f, %s", p, r$n, r$mean_bearing_deg,
          r$resultant_length, fmt_p(r$p_value)) }

write_caption("FigureS10_directionality_rose", paste(
  "Supplementary Figure S10. Direction of first West Nile virus invasion relative",
  "to the nearest previously invaded region, 2007 to 2026. For each first invasion",
  "the initial great-circle bearing was computed from the centroid of the nearest",
  "region already invaded at the start of that season to the centroid of the newly",
  "invaded region, and bearings were grouped into 16 sectors of 22.5 degrees, north",
  "at the top and bearings increasing clockwise. The analysis is descriptive and",
  "includes the 2025 and 2026 seasons; the three regions invaded in 2006 had no",
  "previously invaded region and so have no bearing.",
  sprintf("(A) All %d first invasions; bar length is the number of events in each", row_all$n),
  "sector. The arrow is the mean resultant vector: it points along the mean bearing,",
  "and its length is the mean resultant length on a scale on which the outer edge of",
  "the rose equals 1, so a short arrow means the bearings are close to uniformly",
  "spread.",
  sprintf("Mean bearing %.1f degrees, resultant length %.3f, Rayleigh test %s.",
          row_all$mean_bearing_deg, row_all$resultant_length, fmt_p(row_all$p_value)),
  sprintf("By quadrant, westward %s, southward %s, northward %s and eastward %s.",
          qtxt("Westward"), qtxt("Southward"), qtxt("Northward"), qtxt("Eastward")),
  "(B) The same rose within each period, with bar length the percentage of that",
  "period's first invasions so that periods of different size can be compared:",
  paste0(per(PERIODS[1]), "; ", per(PERIODS[2]), "; ", per(PERIODS[3]), "; ",
         per(PERIODS[4]), "."),
  "Data for 2026 are to 10 September 2026."))

message("figS10 done")
