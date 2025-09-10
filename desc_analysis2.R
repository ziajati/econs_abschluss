# ---------------------------------------------
# Deskriptive Analyse – EZB & Taylor-Regel
# ---------------------------------------------
# Voraussetzungen:

library(tidyverse)
library(lubridate)
library(gt)
library(scales)
library(patchwork)
library(zoo)

# ---- 1) Daten laden & vorbereiten ----
df_raw <- read_csv("taylor_rule_interpolated.csv",
                   show_col_types = FALSE)

df <- df_raw %>%
  mutate(
    datum = parse_date_time(datum, orders = "b Y"),     # "Dec 2001" -> Datum
    jahr  = year(datum),
    monat = month(datum, label = TRUE, abbr = TRUE)
  ) %>%
  arrange(datum)

# Numerische Variablen (für Statistiken/Plots)
num_vars <- c("ezb_leitzins", "hicp_inflation", "inflationsabweichung", "output_gap")

# ---- 2) Überblick (Struktur & Zeitraum) ----
cat("Zeitraum: ", format(min(df$datum), "%b %Y"), "bis", format(max(df$datum), "%b %Y"), "\n")
glimpse(df)

# ---- 3) Deskriptive Statistiken (gt) ----
# Hilfsfunktion: Kennzahlen je Variable
descr_tbl <- function(data, vars) {
  data %>%
    summarise(
      across(all_of(vars),
             list(
               n      = ~sum(!is.na(.)),
               mean   = ~mean(., na.rm = TRUE),
               sd     = ~sd(., na.rm = TRUE),
               min    = ~min(., na.rm = TRUE),
               q25    = ~quantile(., 0.25, na.rm = TRUE),
               median = ~median(., na.rm = TRUE),
               q75    = ~quantile(., 0.75, na.rm = TRUE),
               max    = ~max(., na.rm = TRUE)
             ),
             .names = "{.col}__{.fn}")
    ) %>%
    pivot_longer(everything(),
                 names_to = c("variable", ".value"),
                 names_sep = "__") %>%
    mutate(
      variable = recode(variable,
                        ezb_leitzins = "EZB-Leitzins",
                        hicp_inflation = "HICP-Inflation",
                        inflationsabweichung = "Inflationsabweichung (π-π*)",
                        output_gap = "Output-Gap")
    )
}

gt_stats <- descr_tbl(df, num_vars) %>%
  gt(rowname_col = "variable") %>%
  fmt_number(columns = everything(),
             decimals = 2) %>%
  tab_header(
    title = md("**Deskriptive Statistik**"),
    subtitle = paste0(format(min(df$datum), "%b %Y"), " – ",
                      format(max(df$datum), "%b %Y"))
  ) %>%
  cols_label(
    n = "n", mean = "Mittelwert", sd = "Std.-Abw.",
    min = "Min", q25 = "Q1", median = "Median", q75 = "Q3", max = "Max"
  )

gt_stats

# ---- 4) Korrelationen (gt + Heatmap) ----
cor_mat <- df %>%
  select(all_of(num_vars)) %>%
  cor(use = "pairwise.complete.obs")

# gt-Tabelle der Korrelationen
gt_cor <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  gt(rowname_col = "Variable") %>%
  fmt_number(everything(), decimals = 2) %>%
  tab_header(title = md("**Korrelationsmatrix**"))

gt_cor

# Heatmap (visuell)
cor_long <- as_tibble(cor_mat, rownames = "x") %>%
  pivot_longer(-x, names_to = "y", values_to = "corr")

p_cor <- ggplot(cor_long, aes(x, y, fill = corr)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-1, 1), name = "ρ") +
  labs(x = NULL, y = NULL, title = "Korrelations-Heatmap") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_cor

# ---- 5) Zeitreihenplots ----
# 5a) Leitzins & Inflation (mit Zielinflation 2%)
p_ts <- ggplot(df, aes(datum)) +
  geom_line(aes(y = ezb_leitzins, color = "EZB-Leitzins"), linewidth = 0.9) +
  geom_line(aes(y = hicp_inflation, color = "HICP-Inflation"), linewidth = 0.9, alpha = 0.9) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  scale_color_manual(values = c("EZB-Leitzins" = "#2C7FB8", "HICP-Inflation" = "#D95F0E")) +
  labs(title = "EZB-Leitzins und HICP-Inflation",
       x = NULL, y = "Prozent", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p_ts

# 5b) Inflationsabweichung & Output-Gap
p_ts_gap <- df %>%
  select(datum, inflationsabweichung, output_gap) %>%
  pivot_longer(-datum, names_to = "variable", values_to = "wert") %>%
  mutate(variable = recode(variable,
                           inflationsabweichung = "Inflationsabweichung (π-π*)",
                           output_gap = "Output-Gap")) %>%
  ggplot(aes(datum, wert, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("#7570B3", "#1B9E77")) +
  labs(title = "Abweichungen im Zeitverlauf", x = NULL, y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p_ts_gap
# inflation gap und ouput gap verhalten sich synchron - das spricht für den 
# gleichen faktor von 0.5 den taylor vorschlägt

# 5c) 12-Monats-Gleitmittel für Inflation (Glättung)
df <- df %>%
  arrange(datum) %>%
  mutate(
    hicp_ma12 = zoo::rollmean(hicp_inflation, k = 12, fill = NA, align = "right")
  )

p_infl_ma <- ggplot(df, aes(datum)) +
  geom_line(aes(y = hicp_inflation), alpha = 0.3) +
  geom_line(aes(y = hicp_ma12), linewidth = 1) +
  labs(title = "HICP-Inflation: Monatswerte & 12M-Gleitmittel",
       x = NULL, y = "Prozent") +
  theme_minimal(base_size = 12)

p_infl_ma

# ---- 6) Verteilungen (Histogramme + Dichte) ----
plot_hist <- function(v, v_lab, fill_col = "#3182BD") {
  ggplot(df, aes(x = .data[[v]])) +
    geom_histogram(bins = 30, fill = fill_col, alpha = 0.8) +
    geom_density(linewidth = 1) +
    labs(x = v_lab, y = "Häufigkeit") +
    theme_minimal(base_size = 12)
}

p_h1 <- plot_hist("ezb_leitzins", "EZB-Leitzins", "#2C7FB8")
p_h2 <- plot_hist("hicp_inflation", "HICP-Inflation", "#D95F0E")
p_h3 <- plot_hist("inflationsabweichung", "Inflationsabweichung (π-π*)", "#FDB863")
p_h4 <- plot_hist("output_gap", "Output-Gap", "#1B9E77")

p_hists <- (p_h1 | p_h2) / (p_h3 | p_h4) +
  plot_annotation(title = "Verteilungen der Kernvariablen")
p_hists
# ---- 7) Scatter mit Trend (Zins vs. Gaps) ----
p_sc1 <- ggplot(df, aes(inflationsabweichung, ezb_leitzins)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Inflationsabweichung (π-π*)", y = "EZB-Leitzins",
       title = "Zins vs. Inflationsabweichung") +
  theme_minimal(base_size = 12)

p_sc2 <- ggplot(df, aes(output_gap, ezb_leitzins)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Output-Gap", y = "EZB-Leitzins",
       title = "Zins vs. Output-Gap") +
  theme_minimal(base_size = 12)

p_scatter <- p_sc1 | p_sc2

p_scatter

# ---- 8) Saisonale Muster (Inflation nach Monat) ----
p_season <- ggplot(df, aes(monat, hicp_inflation)) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(title = "Saisonale Verteilung: HICP-Inflation nach Monat",
       x = NULL, y = "Prozent") +
  theme_minimal(base_size = 12)

p_season

# ---- 9) Plots anzeigen & optional speichern ----
# Anzeigen
print(p_ts)
print(p_ts_gap)
print(p_infl_ma)
print(p_hists)
print(p_scatter)
print(p_cor)
print(p_season)

# Optional: Speichern (Ordner "figures" wird erstellt, falls nicht vorhanden)
dir.create("figures", showWarnings = FALSE)
ggsave("figures/zeitreihe_leitzins_inflation.png", p_ts, width = 10, height = 5, dpi = 300)
ggsave("figures/zeitreihe_gaps.png", p_ts_gap, width = 10, height = 5, dpi = 300)
ggsave("figures/hicp_ma12.png", p_infl_ma, width = 10, height = 5, dpi = 300)
ggsave("figures/verteilungen.png", p_hists, width = 12, height = 8, dpi = 300)
ggsave("figures/scatter_zins_gaps.png", p_scatter, width = 10, height = 5, dpi = 300)
ggsave("figures/cor_heatmap.png", p_cor, width = 6, height = 5, dpi = 300)
ggsave("figures/saisonal_inflation.png", p_season, width = 8, height = 5, dpi = 300)

# ---- 10) Bonus: Kennzahlen nach Subperioden (gt) ----
# (z.B. Pre-Krise, Niedrigzinsphase, Post-2020 – Grenzen kannst du anpassen)
cuts <- as.Date(c("2008-01-01","2014-01-01","2020-01-01"))
df_periods <- df %>%
  mutate(
    periode = case_when(
      datum <  cuts[1] ~ "Bis 2007",
      datum >= cuts[1] & datum < cuts[2] ~ "2008–2013",
      datum >= cuts[2] & datum < cuts[3] ~ "2014–2019",
      datum >= cuts[3] ~ "Seit 2020",
      TRUE ~ NA_character_
    )
  )

gt_stats_periods <- df_periods %>%
  group_by(periode) %>%
  summarise(
    across(all_of(num_vars),
           list(mean = ~mean(., na.rm = TRUE),
                sd   = ~sd(., na.rm = TRUE)),
           .names = "{.col}__{.fn}")
  ) %>%
  pivot_longer(-periode, names_to = c("variable",".value"), names_sep = "__") %>%
  mutate(variable = recode(variable,
                           ezb_leitzins = "EZB-Leitzins",
                           hicp_inflation = "HICP-Inflation",
                           inflationsabweichung = "Inflationsabw.",
                           output_gap = "Output-Gap")) %>%
  arrange(match(periode, c("Bis 2007","2008–2013","2014–2019","Seit 2020"))) %>%
  gt(groupname_col = "periode") %>%
  fmt_number(columns = c(mean, sd), decimals = 2) %>%
  cols_label(mean = "Mittelwert", sd = "Std.-Abw.",
             variable = "Variable") %>%
  tab_header(title = md("**Kennzahlen nach Subperioden**"))

gt_stats_periods
