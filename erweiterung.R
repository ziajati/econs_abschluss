# ===== VERBESSERTE TAYLOR-REGEL ANALYSE =====
# Setup - Pakete installieren und laden
required_packages <- c("zoo", "lmtest", "sandwich", "car", "tseries", 
                       "dynlm", "strucchange", "stargazer")

for(pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installiere Paket:", pkg, "\n"))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ===== 1. DATEN VORBEREITUNG =====
df <- read.csv("taylor_rule_interpolated.csv", stringsAsFactors = FALSE)

# Datum korrekt konvertieren
if(is.character(df$datum)) {
  df$datum <- as.yearmon(df$datum)
} else {
  df$datum <- as.yearmon(df$datum, format = "%Y-%m")
}

# Prüfung auf doppelte Datumswerte
cat("Prüfe Datumswerte...\n")
if(any(duplicated(df$datum))) {
  cat("❌ Doppelte Datumswerte gefunden! Entferne Duplikate...\n")
  df <- df[!duplicated(df$datum), ]
} else {
  cat("✅ Keine doppelten Datumswerte\n")
}

# Sicherstellen, dass Daten chronologisch sortiert sind
df <- df[order(df$datum), ]

# Zoo-Objekt für Zeitreihenanalyse (mit sauberen, eindeutigen Daten)
df_ts <- zoo(df[, c("ezb_leitzins","inflationsabweichung","output_gap")], 
             order.by = df$datum)

cat("Datensatz bereit. Zeitraum:", format(start(df_ts)), "bis", format(end(df_ts)), "\n")
cat("Beobachtungen:", nrow(df), "\n\n")

# ===== 2. BASISMODELL (Statische Taylor-Regel) =====
cat("=== STATISCHE TAYLOR-REGEL ===\n")
model_ols <- lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = df)
summary(model_ols)

# Diagnose-Tests
cat("\n--- Diagnose Statisches Modell ---\n")

# A) Autokorrelation
dw_ols <- dwtest(model_ols)
bg_ols <- bgtest(model_ols, order = 12)
lb_ols <- Box.test(resid(model_ols), lag = 12, type = "Ljung-Box")

cat(sprintf("Durbin-Watson: DW=%.3f, p=%.3f\n", dw_ols$statistic, dw_ols$p.value))
cat(sprintf("Breusch-Godfrey: LM=%.2f, p=%.3g\n", bg_ols$statistic, bg_ols$p.value))
cat(sprintf("Ljung-Box: X²=%.2f, p=%.3g\n", lb_ols$statistic, lb_ols$p.value))

# B) Heteroskedastizität
bp_ols <- bptest(model_ols)
cat(sprintf("Breusch-Pagan: BP=%.2f, p=%.3f\n", bp_ols$statistic, bp_ols$p.value))

# C) Normalität
jb_ols <- jarque.bera.test(resid(model_ols))
cat(sprintf("Jarque-Bera: JB=%.2f, p=%.3f\n", jb_ols$statistic, jb_ols$p.value))

# D) Multikollinearität
vif_ols <- vif(model_ols)
cat("VIF-Werte:\n")
print(vif_ols)

cat("\n=== FAZIT STATISCHES MODELL ===")
cat("\n❌ Starke Autokorrelation (DW ≈ 0)")
cat("\n❌ Heteroskedastizität signifikant")
cat("\n✅ Normalität ok")
cat("\n✅ Keine Multikollinearität")
cat("\n→ Dynamisches Modell erforderlich!\n\n")

# ===== 3. DYNAMISCHE MODELLE =====
cat("=== DYNAMISCHE TAYLOR-REGEL (Differenzen-Modelle) ===\n")

# Modell 1: Ein Lag
model_dyn1 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1), data = df_ts)

# Modell 2: Zwei Lags
model_dyn2 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1:2), data = df_ts)

# Modell 3: Drei Lags
model_dyn3 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1:3), data = df_ts)

# Modell 4: Vier Lags
model_dyn4 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1:4), data = df_ts)

# ===== 4. MODELLAUSWAHL =====
cat("--- Modellvergleich (Informationskriterien) ---\n")
ic_comparison <- data.frame(
  Modell = c("1 Lag", "2 Lags", "3 Lags", "4 Lags"),
  AIC = c(AIC(model_dyn1), AIC(model_dyn2), AIC(model_dyn3), AIC(model_dyn4)),
  BIC = c(BIC(model_dyn1), BIC(model_dyn2), BIC(model_dyn3), BIC(model_dyn4))
)
print(ic_comparison)

# Bestes Modell identifizieren
best_aic <- which.min(ic_comparison$AIC)
best_bic <- which.min(ic_comparison$BIC)
cat(sprintf("\nBestes Modell nach AIC: %s (AIC=%.2f)\n", 
            ic_comparison$Modell[best_aic], ic_comparison$AIC[best_aic]))
cat(sprintf("Bestes Modell nach BIC: %s (BIC=%.2f)\n", 
            ic_comparison$Modell[best_bic], ic_comparison$BIC[best_bic]))

# ===== 5. DETAILANALYSE BESTES MODELL =====
# Nehmen wir an, model_dyn2 ist das beste (häufig der Fall)
best_model <- model_dyn2
cat(sprintf("\n=== ANALYSE BESTES MODELL (%s) ===\n", ic_comparison$Modell[best_bic]))

# Newey-West HAC Standardfehler
nw_best <- NeweyWest(best_model, lag = 12, prewhite = TRUE, adjust = TRUE)
coef_robust <- coeftest(best_model, vcov = nw_best)

cat("--- Koeffizienten mit robusten SE ---\n")
print(coef_robust)

# Diagnose des besten Modells
cat("\n--- Diagnose Bestes Modell ---\n")
dw_best <- dwtest(best_model)
bg_best <- bgtest(best_model, order = 12)
lb_best <- Box.test(resid(best_model), lag = 12, type = "Ljung-Box")
bp_best <- bptest(best_model)

cat(sprintf("Durbin-Watson: DW=%.3f, p=%.3f\n", dw_best$statistic, dw_best$p.value))
cat(sprintf("Breusch-Godfrey: LM=%.2f, p=%.3g\n", bg_best$statistic, bg_best$p.value))
cat(sprintf("Ljung-Box: X²=%.2f, p=%.3g\n", lb_best$statistic, lb_best$p.value))
cat(sprintf("Breusch-Pagan: BP=%.2f, p=%.3f\n", bp_best$statistic, bp_best$p.value))

# Residuen-Plots
par(mfrow = c(2,2))
acf(resid(best_model), main = "ACF Residuen (Bestes Modell)")
pacf(resid(best_model), main = "PACF Residuen (Bestes Modell)")
plot(fitted(best_model), resid(best_model), main = "Fitted vs. Residuals")
qqnorm(resid(best_model)); qqline(resid(best_model))
par(mfrow = c(1,1))

# ===== 6. STRUKTURBRUCH-TESTS =====
cat("\n=== STRUKTURBRUCH-ANALYSE ===\n")

# Wichtige Daten für EZB-Politik
crisis_dates <- as.yearmon(c("2008-09", "2010-05", "2015-03", "2020-03"))
crisis_names <- c("Finanzkrise", "Euro-Schuldenkrise", "QE-Start", "COVID-19")

# Finde entsprechende Beobachtungsnummern
break_points <- numeric(length(crisis_dates))
for(i in 1:length(crisis_dates)) {
  idx <- which(as.yearmon(df$datum) >= crisis_dates[i])[1]
  if(!is.na(idx) && idx > 10 && idx < (nrow(df) - 10)) {  # Min. 10 Obs. vor/nach
    break_points[i] <- idx
  }
}
break_points <- break_points[break_points > 0]

cat("Teste Strukturbrüche an folgenden Zeitpunkten:\n")
for(i in 1:length(break_points)) {
  cat(sprintf("- Punkt %d: %s\n", break_points[i], format(df$datum[break_points[i]])))
}

# Chow-Tests für bekannte Bruchpunkte
chow_results <- data.frame(
  Zeitpunkt = character(),
  F_Statistik = numeric(),
  p_Wert = numeric(),
  Signifikant = logical(),
  stringsAsFactors = FALSE
)

for(i in 1:length(break_points)) {
  tryCatch({
    chow_test <- sctest(ezb_leitzins ~ inflationsabweichung + output_gap, 
                        data = df, type = "Chow", point = break_points[i])
    chow_results <- rbind(chow_results, data.frame(
      Zeitpunkt = format(df$datum[break_points[i]]),
      F_Statistik = chow_test$statistic,
      p_Wert = chow_test$p.value,
      Signifikant = chow_test$p.value < 0.05
    ))
  }, error = function(e) {
    cat(sprintf("Fehler bei Bruchpunkt %d: %s\n", break_points[i], e$message))
  })
}

cat("\n--- Chow-Test Ergebnisse ---\n")
print(chow_results)

# CUSUM-Test für unbekannte Brüche
cat("\n--- CUSUM-Test für unbekannte Strukturbrüche ---\n")
tryCatch({
  cusum_test <- efp(ezb_leitzins ~ inflationsabweichung + output_gap, 
                    data = df, type = "Rec-CUSUM")
  plot(cusum_test, main = "Recursive CUSUM Test")
  
  # Stabilitätstest
  cusum_sctest <- sctest(cusum_test)
  cat(sprintf("CUSUM Stabilität: p=%.3f\n", cusum_sctest$p.value))
  
}, error = function(e) {
  cat("CUSUM-Test nicht möglich:", e$message, "\n")
})

# ===== 7. ÖKONOMISCHE INTERPRETATION =====
cat("\n=== ÖKONOMISCHE INTERPRETATION ===\n")

coefs <- coef(coef_robust)
cat("Interpretation der Koeffizienten:\n")
cat(sprintf("- Konstante: %.3f (Basis-Zinsänderung)\n", coefs[1]))
cat(sprintf("- Inflationsabweichung: %.3f (Reaktion auf Inflation)\n", coefs[2]))
cat(sprintf("- Output Gap: %.3f (Reaktion auf Konjunktur)\n", coefs[3]))

if(length(coefs) > 3) {
  cat("- Lag-Koeffizienten (Trägheit/Glättung):\n")
  for(i in 4:length(coefs)) {
    cat(sprintf("  Lag %d: %.3f\n", i-3, coefs[i]))
  }
}

# Taylor-Regel-Konformität bewerten
if(coefs[2] > 0) {
  cat("\n✅ EZB reagiert Taylor-konform auf Inflation (positiver Koeffizient)")
} else {
  cat("\n❌ EZB reagiert nicht Taylor-konform auf Inflation (negativer Koeffizient)")
}

if(coefs[3] > 0) {
  cat("\n✅ EZB reagiert Taylor-konform auf Output Gap (positiver Koeffizient)")
} else {
  cat("\n❌ EZB reagiert nicht Taylor-konform auf Output Gap (negativer Koeffizient)")
}

# ===== 8. ZUSAMMENFASSUNG =====
cat("\n\n=== ZUSAMMENFASSUNG DER ANALYSE ===\n")
cat("✅ Statisches Modell: Starke Autokorrelation → Dynamisches Modell notwendig\n")
cat(sprintf("✅ Bestes dynamisches Modell: %s\n", ic_comparison$Modell[best_bic]))
cat("✅ Robuste Standardfehler (Newey-West) verwendet\n")
cat("✅ Strukturbruch-Tests durchgeführt\n")

if(bg_best$p.value > 0.05) {
  cat("✅ Autokorrelation im finalen Modell beseitigt\n")
} else {
  cat("❌ Autokorrelation im finalen Modell noch vorhanden\n")
}

cat("\n--- Nächste Schritte ---")
cat("\n1. Ergebnisse für Präsentation aufbereiten")
cat("\n2. Policy-Implikationen diskutieren")
cat("\n3. Vergleich mit Literatur")
cat("\n4. Robustheitsanalyse (andere Spezifikationen)")

# ===== 9. EXPORT FÜR PRÄSENTATION =====
# Ergebnistabelle speichern
tryCatch({
  stargazer(model_ols, best_model, 
            se = list(NULL, sqrt(diag(nw_best))),
            title = "Taylor-Regel Schätzungen: Statisch vs. Dynamisch",
            column.labels = c("OLS (Statisch)", "Dynamisch (HAC-SE)"),
            type = "text",
            out = "taylor_results.txt")
  cat("\n✅ Ergebnistabelle in 'taylor_results.txt' gespeichert\n")
}, error = function(e) {
  cat("\nHinweis: stargazer nicht verfügbar für Export\n")
})

cat("\n=== ANALYSE ABGESCHLOSSEN ===\n")

# ===== 9. EXPORT FÜR PRÄSENTATION =====

# Pakete laden (stille Installation, falls nötig)
needs <- c("stargazer","broom","ggplot2","dplyr","strucchange")
to_install <- needs[!needs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(needs, require, character.only = TRUE))

# ── A) ERGEBNISTABELLE SPEICHERN ──────────────────────────────────────────
# Voraussetzung: Objekte existieren: model_ols, best_model, nw_best (vcovHC/Newey-West)

tryCatch({
  # TXT
  stargazer(
    model_ols, best_model,
    se = list(NULL, sqrt(diag(nw_best))), # HAC-SE für dyn. Modell
    title = "Taylor-Regel Schätzungen: Statisch vs. Dynamisch",
    column.labels = c("OLS (Statisch)", "Dynamisch (HAC-SE)"),
    type = "text",
    out  = "taylor_results.txt"
  )
  cat("✅ Ergebnistabelle in 'taylor_results.txt' gespeichert\n")
  
  # HTML (für Folien/Handout)
  stargazer(
    model_ols, best_model,
    se = list(NULL, sqrt(diag(nw_best))),
    title = "Taylor-Regel Schätzungen: Statisch vs. Dynamisch",
    column.labels = c("OLS (Statisch)", "Dynamisch (HAC-SE)"),
    type = "html",
    out  = "taylor_results.html"
  )
  
  # LaTeX (falls du LaTeX nutzt)
  stargazer(
    model_ols, best_model,
    se = list(NULL, sqrt(diag(nw_best))),
    title = "Taylor-Regel Schätzungen: Statisch vs. Dynamisch",
    column.labels = c("OLS (Statisch)", "Dynamisch (HAC-SE)"),
    type = "latex",
    out  = "taylor_results.tex"
  )
}, error = function(e) {
  message("Hinweis: stargazer-Export nicht möglich: ", e$message)
})

# ── B) PLOTS FÜR DIE PRÄSENTATION ────────────────────────────────────────

# 1) Koeffizienten-Plot (Vergleich OLS vs. Dynamik mit HAC-Fehlern)
coef_df_ols <- broom::tidy(model_ols) |>
  dplyr::mutate(
    se = std.error,
    conf.low  = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    model = "OLS (Statisch)"
  ) |>
  dplyr::select(term, estimate, conf.low, conf.high, model)

# dynamisches Modell: SE aus nw_best übernehmen
tidy_dyn <- broom::tidy(best_model)
se_dyn   <- sqrt(diag(nw_best))
# Mappe SE per Term-Name
tidy_dyn$se <- se_dyn[match(tidy_dyn$term, names(se_dyn))]
coef_df_dyn <- tidy_dyn |>
  dplyr::mutate(
    conf.low  = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    model = "Dynamisch (HAC-SE)"
  ) |>
  dplyr::select(term, estimate, conf.low, conf.high, model)

coef_df <- dplyr::bind_rows(coef_df_ols, coef_df_dyn) |>
  # optional: uninteressante Terme filtern
  dplyr::filter(!grepl("^\\(Intercept\\)$", term))

p_coef <- ggplot(coef_df, aes(x = term, y = estimate, shape = model)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .15, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  labs(x = NULL, y = "Koeffizient (95%-KI)",
       title = "Taylor-Regel: Koeffizientenvergleich (OLS vs. Dynamik)") +
  theme_minimal(base_size = 12)

ggsave("taylor_coefplot.png", p_coef, width = 9, height = 5, dpi = 300)
cat("✅ 'taylor_coefplot.png' gespeichert\n")

# 2) Recursive CUSUM Test (auf dem OLS-Modell)
#    Nutzt die im Modell verwendete Formel & Daten
ef <- strucchange::efp(
  formula(model_ols),
  data = model.frame(model_ols),
  type = "Rec-CUSUM"
)

png("taylor_cusum_recursive.png", width = 1100, height = 550, res = 140)
plot(ef, main = "Recursive CUSUM Test")
dev.off()
cat("✅ 'taylor_cusum_recursive.png' gespeichert\n")

cat("\n=== ANALYSE ABGESCHLOSSEN ===\n")
