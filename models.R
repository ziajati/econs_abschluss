# ===== Setup =====
library(zoo)
library(lmtest)     # bptest, coeftest
library(sandwich)   # robuste SE
library(car)        # VIF
library(tseries)    # Jarque-Bera
library(gt)
library(broom)
library(stats)    # acf, pacf, Box.test
library(tseries)   # für den ADF-Test
library(dynlm)
library(dplyr)
library(ggplot2)
library(tidyverse)


# 1) Daten laden & Index sauber setzen
df <- read.csv("taylor_rule_interpolated.csv", stringsAsFactors = FALSE)
df$datum <- as.yearmon(df$datum, format = "%Y-%m")

# 2) OLS: Basismodell (statische Taylor-Regel)
model_ols <- lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = df)
sm = summary(model_ols)
print(sm)

# Koefﬁzienten-Tabelle vorbereiten
coef_df <- as.data.frame(sm$coefficients)
coef_df$Variable <- rownames(coef_df)
coef_df$Interpretation <- ifelse(coef_df$`Pr(>|t|)` <= 0.05, "signifikant", "nicht signifikant")

# Spaltenreihenfolge anpassen
coef_df <- coef_df[, c("Variable", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "Interpretation")]

# gt-Tabelle erzeugen
gt(coef_df) %>%
  tab_header(title = "Regressionskoeffizienten Modell 1") %>%
  cols_label(
    Variable = "Variable",
    Estimate = "Schätzwert",
    `Std. Error` = "Standardfehler",
    `t value` = "t-Wert",
    `Pr(>|t|)` = "p-Wert",
    Interpretation = "Interpretation (α = 0,05)"
  ) %>%
  fmt_number(
    columns = c(Estimate, `Std. Error`, `t value`, `Pr(>|t|)`),
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )


# Einzelwerte berechnen
r2 <- sm$r.squared
adj_r2 <- sm$adj.r.squared
sigma <- sm$sigma
df_model <- sm$df[1]
df_resid <- sm$df[2]
fstat <- sm$fstatistic
f_p <- 1 - pf(fstat["value"], df1 = fstat["numdf"], df2 = fstat["dendf"])

# DataFrame für Zusammenfassung
info_df <- data.frame(
  Kennzahl = c("R²", "Adjustiertes R²", "F-Statistik", "p-Wert (F)", "df Modell", "df Residuen", "sigma"),
  Wert = c(round(r2, 4), round(adj_r2, 4), round(fstat["value"], 4), round(f_p, 4), df_model, df_resid, round(sigma, 4))
)

# gt-Tabelle
gt(info_df) %>%
  tab_header(title = "Modellkennzahlen Modell 1") %>%
  fmt_number(columns = "Wert", decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )


# 2) OLS: Basismodell Spline (statische Taylor-Regel)
model_ols_spline <- lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = df_spline)
sm = summary(model_ols_spline)
print(sm)

# Koefﬁzienten-Tabelle vorbereiten
coef_df <- as.data.frame(sm$coefficients)
coef_df$Variable <- rownames(coef_df)
coef_df$Interpretation <- ifelse(coef_df$`Pr(>|t|)` <= 0.05, "signifikant", "nicht signifikant")

# Spaltenreihenfolge anpassen
coef_df <- coef_df[, c("Variable", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "Interpretation")]

# gt-Tabelle erzeugen
gt(coef_df) %>%
  tab_header(title = "Regressionskoeffizienten Modell 2") %>%
  cols_label(
    Variable = "Variable",
    Estimate = "Schätzwert",
    `Std. Error` = "Standardfehler",
    `t value` = "t-Wert",
    `Pr(>|t|)` = "p-Wert",
    Interpretation = "Interpretation (α = 0,05)"
  ) %>%
  fmt_number(
    columns = c(Estimate, `Std. Error`, `t value`, `Pr(>|t|)`),
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

# Einzelwerte berechnen
r2 <- sm$r.squared
adj_r2 <- sm$adj.r.squared
sigma <- sm$sigma
df_model <- sm$df[1]
df_resid <- sm$df[2]
fstat <- sm$fstatistic
f_p <- 1 - pf(fstat["value"], df1 = fstat["numdf"], df2 = fstat["dendf"])

# DataFrame für Zusammenfassung
info_df <- data.frame(
  Kennzahl = c("R²", "Adjustiertes R²", "F-Statistik", "p-Wert (F)", "df Modell", "df Residuen", "sigma"),
  Wert = c(round(r2, 4), round(adj_r2, 4), round(fstat["value"], 4), round(f_p, 4), df_model, df_resid, round(sigma, 4))
)

# gt-Tabelle
gt(info_df) %>%
  tab_header(title = "Modellkennzahlen Modell 2") %>%
  fmt_number(columns = "Wert", decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

# ===== Standard-Checks =====

## A) Heteroskedastizität
bp <- bptest(model_ols)  # Breusch-Pagan
bp

# Goldfeld-Quandt Test
gq_test <- gqtest(model_ols, order.by = df$inflationsabweichung)
print(gq_test)

# White Test (robuste Version des BP-Tests)
white_test <- bptest(model_ols, ~ inflationsabweichung * output_gap + I(inflationsabweichung^2) + I(output_gap^2), data = df)
print(white_test)

# Robuste (HC1) SE anzeigen – nur zur Einordnung, falls BP/White signifikant
coeftest(model_ols, vcov = vcovHC(model_ols, type = "HC1"))

## B) Normalität der Residuen
res <- resid(model_ols)
jarque.bera.test(res)       # JB-Test
# X-squared = 0.44472, df = 2, p-value = 0.8006
# p=0.80 → H0 (Normalität) nicht verworfen → Residuen sehen nach JB normal aus.

# Q-Q-Plot (schnell visuell)
qqnorm(res); qqline(res)
# Deutliche Abweichungen in den Rändern (dicke Tails, leichte S-Kurve), im Zentrum recht gut.

## C) Multikollinearität
vif(model_ols)              # > 10 = stark, ~5 = auffällig
# inflationsabweichung           output_gap 
# 1.1736               1.1736 
# keine multikolinearität

vif_data <- data.frame(Variable = names(vif_values), VIF = as.numeric(vif_values))
ggplot(vif_data, aes(x = Variable, y = VIF)) +
  geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
  geom_hline(yintercept = 5,  linetype = "dashed", color = "orange", alpha = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 1, y = 4.7, label = "unauffällig", hjust = 0, vjust = 1, size = 3) +
  annotate("text", x = 1, y = 9.7, label = "kritisch",    hjust = 0, vjust = 1, size = 3) +
  labs(title = "Variance Inflation Factors (Multikollinearität)", y = "VIF", x = "Regressor") +
  theme_minimal()

# 1) Durbin–Watson (Lag-1-Autokorrelation)
dwtest(model_ols)
# DW = 0.02 (!!), p < 0.001 → extrem starke positive Autokorrelation.

# 2) Breusch–Godfrey (Autokorrelation)
bgtest(model_ols, order = 3)
# LM test = 272.94, df = 4, p-value < 2.2e-16 -> Autokorrelation
bgtest(model_ols, order = 6)
# LM test = 273.01, df = 8, p-value < 2.2e-16
bgtest(model_ols, order = 12)
# LM test = 273.07, df = 12, p-value < 2.2e-16

# 3) ACF/PACF der Residuen (Grafiken)
res <- resid(model_ols)
par(mfrow = c(1,2))
acf(res, main = "ACF Residuen (OLS)")
pacf(res, main = "PACF Residuen (OLS)")
par(mfrow = c(1,1))
# sehr langsamer abfall ACF Residiuen
# PACF hoher Ausschlag bei 1
# Muster passt zu AR(1) Modell

# 4) Ljung–Box auf mehrere Lags (Globaltest auf „keine Autokorrelation“)
Box.test(res, lag = 12, type = "Ljung-Box")


# ===== Daten laden und formatieren für dynamische Analyse =====
df <- read.csv("taylor_rule_interpolated.csv", stringsAsFactors = FALSE)
df_ts <- zoo(df[, c("ezb_leitzins","inflationsabweichung","output_gap")], order.by = df$datum)

# ===== Stationarität testen (ADF-Test) =====
adf.test(na.omit(df_ts$ezb_leitzins))          # Leitzins
adf.test(na.omit(df_ts$inflationsabweichung))  # Inflationsabweichung
adf.test(na.omit(df_ts$output_gap))            # Output-Gap


# Taylor in Levels + AR(1) auf der abhängigen
model_dyn1 <- dynlm(ezb_leitzins ~ inflationsabweichung + output_gap +
                         L(ezb_leitzins, 1),
                       data = df_ts)
# Summary-Objekt
sm <- summary(model_dyn1)

# --- Koefﬁzienten-Tabelle vorbereiten (dein Format) ---
coef_df <- as.data.frame(sm$coefficients)
coef_df$Variable <- rownames(coef_df)
coef_df$Interpretation <- ifelse(coef_df$`Pr(>|t|)` <= 0.05, "signifikant", "nicht signifikant")

# Optional: schönere Variablennamen
pretty_names <- c(
  "(Intercept)" = "Konstante",
  "inflationsabweichung" = "Inflationslücke",
  "output_gap" = "Output-Gap",
  "L(ezb_leitzins, 1)" = "Lag(1) Leitzins"
)
coef_df$Variable <- ifelse(coef_df$Variable %in% names(pretty_names),
                           pretty_names[coef_df$Variable],
                           coef_df$Variable)

# Spaltenreihenfolge anpassen
coef_df <- coef_df[, c("Variable", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "Interpretation")]

# --- gt-Tabelle erzeugen ---
gt_tbl <- gt(coef_df) %>%
  tab_header(title = "Regressionskoeffizienten Modell 1") %>%
  cols_label(
    Variable = "Variable",
    Estimate = "Schätzwert",
    `Std. Error` = "Standardfehler",
    `t value` = "t-Wert",
    `Pr(>|t|)` = "p-Wert",
    Interpretation = "Interpretation (α = 0,05)"
  ) %>%
  fmt_number(
    columns = c(Estimate, `Std. Error`, `t value`),
    decimals = 2
  ) %>%
  # p-Werte: sinnvoll formatiert (sehr kleine p als wissenschaftliche Notation)
  fmt(
    columns = `Pr(>|t|)`,
    fns = function(x) ifelse(x < 0.0001, formatC(x, format = "e", digits = 2),
                             formatC(x, format = "f", digits = 4))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

# Tabelle anzeigen
gt_tbl


model_dyn2 <- dynlm(ezb_leitzins ~ inflationsabweichung + output_gap + L(ezb_leitzins, 2), data = df_ts)
model_dyn3 <- dynlm(ezb_leitzins ~ inflationsabweichung + output_gap + L(ezb_leitzins, 3), data = df_ts)
model_dyn4 <- dynlm(ezb_leitzins ~ inflationsabweichung + output_gap + L(ezb_leitzins, 4), data = df_ts)

# Modelle vergleichen
AIC(model_dyn1, model_dyn2, model_dyn3, model_dyn4)
BIC(model_dyn1, model_dyn2, model_dyn3, model_dyn4)

crit <- tribble(
  ~Model,       ~df, ~AIC,        ~BIC,
  "model_dyn1",  5,  658.9748,    677.0768,
  "model_dyn2",  5,  812.0603,    830.1441,
  "model_dyn3",  5,  836.3615,    854.4271,
  "model_dyn4",  5,  844.9518,    862.9991
) %>%
  mutate(
    dAIC = AIC - min(AIC),
    dBIC = BIC - min(BIC)
  )

# ---- 1) Kompakte gt-Tabelle (mit Hervorhebung des besten) ----
gt_tbl <- crit %>%
  select(Model, df, AIC, BIC, dAIC, dBIC) %>%
  gt() %>%
  tab_header(title = "Modellvergleich (AIC/BIC)") %>%
  cols_label(Model = "Modell", df = "df", AIC = "AIC", BIC = "BIC",
             dAIC = "ΔAIC", dBIC = "ΔBIC") %>%
  fmt_number(columns = c(AIC, BIC, dAIC, dBIC), decimals = 2) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_body(rows = AIC == min(AIC), columns = AIC),
      cells_body(rows = BIC == min(BIC), columns = BIC),
      cells_column_labels(everything())
    )
  )
gt_tbl


# 2. Residuenanalyse
bgtest(model_dyn1, order = 6)           # Autokorrelation
Box.test(resid(model_dyn1), lag = 12, type = "Ljung-Box")

jarque.bera.test(resid(model_dyn1))     # Normalität
# Q-Q-Plot der Residuen
res <- resid(model_dyn1)

qqnorm(res, main = "Q-Q-Plot der Residuen")
qqline(res, col = "red", lwd = 2)

# Goldfeld-Quandt Test (Heteroskedastizität)
mf <- model.frame(model_dyn1)  # enthält genau die im Modell verwendeten Zeilen
gq_test <- gqtest(model_dyn1, order.by = mf$inflationsabweichung)
gq_test
# GQ = 0.99634, df1 = 134, df2 = 134, p-value = 0.5084 -> keine Heterokedastizität 

# Breusch-Pagan
bp <- bptest(model_dyn1)
bp
# P = 34.582, df = 3, p-value = 1.493e-07 -> Heterokedastizität 
# 3. Multikollinearität
v <- vif(model_dyn1)

df_vif <- data.frame(
  variable = names(v),
  VIF = as.numeric(v),
  row.names = NULL
) %>%
  arrange(VIF) %>%
  mutate(variable = factor(variable, levels = variable))

ggplot(df_vif, aes(x = variable, y = VIF)) +
  geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange", alpha = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 1, y = 4.7, label = "good", hjust = 0, vjust = 1, size = 3) +
  annotate("text", x = 1, y = 9.7, label = "acceptable", hjust = 0, vjust = 1, size = 3) +
  labs(title = "Variance Inflation Factors (Multicollinearity)",
       y = "VIF", x = "Regressor") +
  theme_minimal()

# 5. Robuste Standardfehler
coeftest(model_dyn1, vcov = NeweyWest(model_dyn1, lag = 6, prewhite = FALSE))

fit <- model_dyn1   # dein bereits geschätztes dynamisches Modell

# Robuste (HAC) Kovarianz
V  <- vcovHAC(fit)
cf <- coef(fit)
se <- sqrt(diag(V))
df <- df.residual(fit)

# Spaltennamen im Modell
name_intercept <- "(Intercept)"
name_pi        <- "inflationsabweichung"
name_gap       <- "output_gap"

# Benchmarks (Ur-Taylor 1993)
b0_benchmark   <- 2.0   # i* = r* + π*; hier 2.0 (an deine Arbeit angepasst)
bpi_benchmark  <- 0.5
bgap_benchmark <- 0.5

# t- und p-Werte (zweiseitig; Taylor-Prinzip optional einseitig)
t_i   <- (cf[name_intercept] - b0_benchmark) / se[name_intercept]
p_i   <- 2 * pt(abs(t_i), df = df, lower.tail = FALSE)

t_pi  <- (cf[name_pi] - bpi_benchmark) / se[name_pi]
p_pi  <- 2 * pt(abs(t_pi), df = df, lower.tail = FALSE)

t_gap <- (cf[name_gap] - bgap_benchmark) / se[name_gap]
p_gap <- 2 * pt(abs(t_gap), df = df, lower.tail = FALSE)

# (Optional) Taylor-Prinzip: H0: beta_pi <= 1 vs H1: > 1
t_princ <- (cf[name_pi] - 1) / se[name_pi]
p_princ <- pt(t_princ, df = df, lower.tail = FALSE)

# Differenz-Tabelle (saubere Spaltennamen)
diff_tbl <- data.frame(
  Test      = c(
    sprintf("Intercept i* = %.1f", b0_benchmark),
    sprintf("β_π = %.1f", bpi_benchmark),
    sprintf("β_y = %.1f", bgap_benchmark),
    "Taylor-Prinzip: β_π > 1"
  ),
  Schaetzer = c(cf[name_intercept], cf[name_pi], cf[name_gap], cf[name_pi]),
  Benchmark = c(b0_benchmark, bpi_benchmark, bgap_benchmark, 1),
  Differenz = c(cf[name_intercept]-b0_benchmark,
                cf[name_pi]-bpi_benchmark,
                cf[name_gap]-bgap_benchmark,
                cf[name_pi]-1),
  SE_HAC    = c(se[name_intercept], se[name_pi], se[name_gap], se[name_pi]),
  t_Wert    = c(t_i, t_pi, t_gap, t_princ),
  p_Wert    = c(p_i, p_pi, p_gap, p_princ),
  Testart   = c("zweiseitig", "zweiseitig", "zweiseitig", "einseitig")
)

gt(diff_tbl) |>
  tab_header(title = "Abweichungen zu Taylor (1993)") |>
  cols_label(
    Schaetzer = "Schätzer",
    SE_HAC    = "Std.-fehler (HAC)",
    t_Wert    = "t-Wert",
    p_Wert    = "p-Wert"
  ) |>
  fmt_number(
    columns = c("Schaetzer","Benchmark","Differenz","SE_HAC","t_Wert","p_Wert"),
    decimals = 3
  )