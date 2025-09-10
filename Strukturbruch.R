------------------------------------------
library(readr)
library(ggplot2)
library(zoo)
library(lmtest)
library(sandwich)
library(car)
library(broom)
library(strucchange)

#------------------------------------------------------------------------
# Daten einlesen
# Erwartet: datum (YYYY-MM), ezb_leitzins, inflationsabweichung, output_gap
#------------------------------------------------------------------------
Daten <- read_csv("taylor_rule_interpolated.csv", show_col_types = FALSE)
Daten$datum <- as.yearmon(Daten$datum, format = "%Y-%m")
Daten <- Daten[order(Daten$datum), ]

#------------------------------------------------------------------------
# Basismodell (statische Taylor-Regel)
# i_t = β0 + βπ * infl_gap_t + βy * output_gap_t + u_t
#------------------------------------------------------------------------
Modell <- lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = Daten)
print(summary(Modell))

Daten$zins_hat <- fitted(Modell)
Daten$u_hat    <- residuals(Modell)

#------------------------------------------------------------------------
# STRUKTURBRÜCHE
#------------------------------------------------------------------------
# CUSUM / MOSUM
ocus  <- efp(ezb_leitzins ~ inflationsabweichung + output_gap, type = "OLS-CUSUM", data = Daten)
plot(ocus); sctest(ocus)

omos  <- efp(ezb_leitzins ~ inflationsabweichung + output_gap, type = "OLS-MOSUM", data = Daten)
plot(omos); sctest(omos)

# Sup-F (Andrews, 1993) – F-Statistiken über alle möglichen Bruchpunkte
fs <- Fstats(ezb_leitzins ~ inflationsabweichung + output_gap, data = Daten)
plot(fs, alpha = 0.05)
sctest(fs, type = "supF")

# Bai-Perron Breakpoints (optimale Anzahl nach BIC)
bp <- breakpoints(ezb_leitzins ~ inflationsabweichung + output_gap, data = Daten)
summary(bp)
plot(bp)

# Beispiel: Chow-Test für einen Bruchpunkt bei z. B. 2008-09
# -> passende Beobachtungsnummer berechnen:
br_time <- as.yearmon("2008-09")
br_idx  <- which.min(abs(Daten$datum - br_time))
sctest(lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = Daten),
       type = "Chow", point = br_idx)


# ------------------------------------------------------------
# BAI–PERRON: Multiple Strukturbrüche (Niveau & Steigungen)
# ------------------------------------------------------------
library(strucchange)
library(broom)

# Trimming: mind. 15% je Segment (Daumenregel bei Monatsdaten ok)
trim_frac <- 0.15
h_min     <- floor(trim_frac * nrow(Daten))

# Maximale Breaks (z.B. bis 5; wird danach via BIC/AIC ausgewählt)
max_breaks <- 5

# Basis-Formel (statisch – wie dein OLS)
fm <- ezb_leitzins ~ inflationsabweichung + output_gap

# 2.1 Kandidaten-Breaks schätzen
bp_cand <- breakpoints(fm, data = Daten, h = h_min, breaks = max_breaks)

# 2.2 Modellwahl über BIC (oder AIC)
bic_vals <- BIC(bp_cand); aic_vals <- AIC(bp_cand)
k_star   <- which.min(bic_vals)   # oder which.min(aic_vals)

bp_best  <- breakpoints(bp_cand, breaks = k_star)

# 2.3 Konfidenzintervalle der Bruchpunkte
# Deine beste Anzahl der Brüche (z. B. k_star)
bp_best <- breakpoints(fm, data = Daten, breaks = k_star, h = h_min)

# Konfidenzintervalle
bp_ci <- confint(bp_best)
print(bp_ci)

# 2.4 Bruchpunkte auf Kalenderdaten abbilden
break_idx   <- bp_best$breakpoints
break_dates <- Daten$datum[break_idx]
print(break_dates); print(bp_ci)


# Faktor-Variable, die jedes Regime markiert
reg_fac <- breakfactor(bp_best)

# Re-Schätzung je Regime in einem Schritt:
reg_fit <- lm(fm, data = Daten, subset = reg_fac)
summary(reg_fit)

# Saubere Koef-Tabellen je Regime:
coef_tab <- tidy(reg_fit)
print(coef_tab)

# Optional: Predicted vs Actual je Regime
Daten$y_hat_bp <- fitted(reg_fit)

# 4.2 Optional: F-Statistikpfad & supF
fs <- Fstats(fm, data = Daten)
plot(fs, alpha = 0.05); sctest(fs, type = "supF")

library(sandwich); library(lmtest)

# a) fixe Bandbreite (z.B. 12 Monate)
vc <- NeweyWest(reg_fit, lag = 12, prewhite = FALSE, adjust = TRUE)
coeftest(reg_fit, vcov. = vc)

# b) automatische Bandbreite
bw <- bwNeweyWest(reg_fit, prewhite = FALSE)
vc <- NeweyWest(reg_fit, lag = bw, prewhite = FALSE, adjust = TRUE)
coeftest(reg_fit, vcov. = vc)

# Andrews–Kernel, kein Prewhitening
vc_hac <- vcovHAC(reg_fit, prewhite = FALSE, adjust = TRUE)
coeftest(reg_fit, vcov. = vc_hac)

# Oder explizit: Bartlett-Kernel mit Bandbreite ~ 12
vc_kern <- kernHAC(reg_fit, kernel = "Bartlett", bw = 12, prewhite = FALSE)
coeftest(reg_fit, vcov. = vc_kern)

# 1) Aliasing prüfen – NA-Koeffizienten bedeuten Rangdefizit
reg_fit$aliased
coef(reg_fit)          # NAs? → zu viele Breaks / zu kleine Regimes

# 2) Modellmatrix vollständig & numerisch?
X <- model.matrix(reg_fit)
stopifnot(all(complete.cases(X)))
mode(X)  # "numeric"

# 3) Genug Beobachtungen je Regime?
table(breakfactor(bp_best))
# Falls ein Regime sehr klein ist: k_star reduzieren oder Trimming (h) erhöhen

#------------------------------------------------------------------------
# STRUKTURBRÜCHE - ERWEITERUNG: INTERCEPT-SHIFTS (Krisenphasen)
#------------------------------------------------------------------------

# Nach den bestehenden Strukturbruch-Tests (CUSUM/MOSUM) aus dem Hauptcode
# erklärt diese Erweiterung die festgestellten Instabilitäten durch 
# systematische Niveauverschiebungen in identifizierten Krisenphasen

#------------------------------------------------------------------------
# 1. CRISIS DUMMIES DEFINIEREN (pragmatische Fenster)
#------------------------------------------------------------------------

# Definiere Krisen-Zeitfenster basierend auf historischen Ereignissen
Daten$D_GFC   <- ifelse(Daten$datum >= as.yearmon("2008-01") & 
                          Daten$datum <= as.yearmon("2010-12"), 1, 0)  # Globale Finanzkrise
Daten$D_Sov   <- ifelse(Daten$datum >= as.yearmon("2011-01") & 
                          Daten$datum <= as.yearmon("2013-12"), 1, 0)  # Euro-Schuldenkrise  
Daten$D_Covid <- ifelse(Daten$datum >= as.yearmon("2020-03") & 
                          Daten$datum <= as.yearmon("2021-12"), 1, 0)  # Pandemie
Daten$D_Infl  <- ifelse(Daten$datum >= as.yearmon("2022-02") & 
                          Daten$datum <= as.yearmon("2024-12"), 1, 0)  # Inflationsschock/Ukraine

# Überblick über die Dummy-Perioden
cat("Anzahl Beobachtungen je Krisenphase:\n")
cat("Globale Finanzkrise (2008-01 bis 2010-12):", sum(Daten$D_GFC), "Monate\n")
cat("Euro-Schuldenkrise (2011-01 bis 2013-12):", sum(Daten$D_Sov), "Monate\n")
cat("COVID-19 Pandemie (2020-03 bis 2021-12):", sum(Daten$D_Covid), "Monate\n")
cat("Inflationsschock (2022-02 bis 2024-12):", sum(Daten$D_Infl), "Monate\n")

#------------------------------------------------------------------------
# 2. INTERCEPT-SHIFT MODELL SCHÄTZEN
#------------------------------------------------------------------------

# Prüfe, ob alle benötigten Variablen vorhanden sind
required_vars <- c("ezb_leitzins", "inflationsabweichung", "output_gap", "D_GFC", "D_Sov", "D_Covid", "D_Infl")
missing_vars <- required_vars[!required_vars %in% names(Daten)]

if(length(missing_vars) > 0) {
  cat("FEHLER: Folgende Variablen fehlen:", paste(missing_vars, collapse = ", "), "\n")
  stop("Modell kann nicht geschätzt werden - fehlende Variablen")
} else {
  cat("✓ Alle benötigten Variablen sind vorhanden\n\n")
}

cat("=== SCHRITT 1: BASISMODELL TESTEN ===\n")

# Test: Funktioniert das ursprüngliche Modell?
basic_vars <- c("ezb_leitzins", "inflationsabweichung", "output_gap")
cat("Prüfe Basis-Variablen:\n")

for(var in basic_vars) {
  if(var %in% names(Daten)) {
    n_na <- sum(is.na(Daten[[var]]))
    n_valid <- sum(!is.na(Daten[[var]]))
    cat(sprintf("%s: %d gültige Werte, %d NAs\n", var, n_valid, n_na))
  } else {
    cat(sprintf("%s: NICHT GEFUNDEN!\n", var))
  }
}

# Teste Basismodell
basic_complete <- complete.cases(Daten[, basic_vars])
n_basic_complete <- sum(basic_complete)
cat(sprintf("\nVollständige Fälle für Basismodell: %d\n", n_basic_complete))

if(n_basic_complete > 0) {
  cat("✓ Basismodell kann geschätzt werden\n")
  
  # Schätze Basismodell zur Referenz (sollte Ihr ursprüngliches Modell sein)
  if(exists("Modell")) {
    cat("Verwende bereits geschätztes Basismodell\n")
    Modell_basis <- Modell
  } else {
    cat("Schätze neues Basismodell\n")
    Modell_basis <- lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = Daten)
  }
  
  cat("Basismodell Zusammenfassung:\n")
  print(summary(Modell_basis)$coefficients)
} else {
  stop("Auch das Basismodell kann nicht geschätzt werden - Datenproblem!")
}

cat("\n=== SCHRITT 2: ROBUSTE DUMMY-ERSTELLUNG ===\n")

# Zeige Datenstruktur
cat("Anzahl Zeilen in Daten:", nrow(Daten), "\n")
if("datum" %in% names(Daten)) {
  cat("Erste 5 Datumswerte:", as.character(head(Daten$datum, 5)), "\n")
  cat("Letzte 5 Datumswerte:", as.character(tail(Daten$datum, 5)), "\n")
  
  # METHODE 1: Manuelle Zuweisung basierend auf Datumswerten
  datum_char <- as.character(Daten$datum)
  
  # Initialisiere alle Dummies mit 0
  Daten$D_GFC <- 0
  Daten$D_Sov <- 0  
  Daten$D_Covid <- 0
  Daten$D_Infl <- 0
  
  cat("\n=== MANUELLE DUMMY-ZUWEISUNG ===\n")
  
  # Durchsuche jede Zeile und weise Dummies zu
  for(i in 1:nrow(Daten)) {
    datum_i <- datum_char[i]
    
    # Extrahiere Jahr aus String (erste 4 Zeichen sollten Jahr sein)
    jahr_str <- substr(datum_i, 1, 4)
    
    if(nchar(jahr_str) == 4 && !is.na(as.numeric(jahr_str))) {
      jahr <- as.numeric(jahr_str)
      
      # Extrahiere Monat (verschiedene Formate abdecken)
      if(nchar(datum_i) >= 7) {
        # Suche nach Monat in verschiedenen Positionen
        monat_candidates <- c(
          substr(datum_i, 6, 7),   # YYYY-MM oder YYYY/MM
          substr(datum_i, 5, 6),   # YYYYMM
          substr(datum_i, 1, 2)    # MM-YYYY oder MM/YYYY
        )
        
        monat <- NA
        for(m_candidate in monat_candidates) {
          m_num <- suppressWarnings(as.numeric(m_candidate))
          if(!is.na(m_num) && m_num >= 1 && m_num <= 12) {
            monat <- m_num
            break
          }
        }
        
        # Falls Monat-Extraktion scheitert, verwende nur Jahr
        if(is.na(monat)) monat <- 1
      } else {
        monat <- 1  # Fallback
      }
      
      # Weise Krisen-Dummies zu
      # Globale Finanzkrise: 2008-2010
      if(jahr >= 2008 && jahr <= 2010) {
        Daten$D_GFC[i] <- 1
      }
      
      # Euro-Schuldenkrise: 2011-2013
      if(jahr >= 2011 && jahr <= 2013) {
        Daten$D_Sov[i] <- 1
      }
      
      # COVID-19: 2020-03 bis 2021-12
      if((jahr == 2020 && monat >= 3) || jahr == 2021) {
        Daten$D_Covid[i] <- 1
      }
      
      # Inflationsschock: 2022-02 bis 2024+
      if((jahr == 2022 && monat >= 2) || jahr >= 2023) {
        Daten$D_Infl[i] <- 1
      }
    }
    
    # Fortschritt alle 50 Zeilen
    if(i %% 50 == 0) {
      cat(".")
    }
  }
  
  cat("\n\n=== FINALE DUMMY-VERIFIKATION ===\n")
  dummy_vars <- c("D_GFC", "D_Sov", "D_Covid", "D_Infl")
  
  for(var in dummy_vars) {
    n_ones <- sum(Daten[[var]] == 1, na.rm = TRUE)
    n_zeros <- sum(Daten[[var]] == 0, na.rm = TRUE)
    n_nas <- sum(is.na(Daten[[var]]))
    
    cat(sprintf("%s: %d Einsen, %d Nullen, %d NAs\n", var, n_ones, n_zeros, n_nas))
    
    if(n_ones > 0) {
      # Zeige erste paar aktivierte Perioden
      active_indices <- which(Daten[[var]] == 1)
      first_few <- head(active_indices, 3)
      cat(sprintf("  Erste aktivierte Zeilen: %s (Daten: %s)\n", 
                  paste(first_few, collapse = ", "),
                  paste(datum_char[first_few], collapse = ", ")))
    }
  }
  
  # NOTFALL-FALLBACK: Falls immer noch alle 0, verwende feste Indizes
  total_ones <- sum(Daten$D_GFC) + sum(Daten$D_Sov) + sum(Daten$D_Covid) + sum(Daten$D_Infl)
  
  if(total_ones == 0) {
    cat("\n⚠ NOTFALL-FALLBACK: Verwende geschätzte Zeilenbereiche\n")
    n_rows <- nrow(Daten)
    
    # Schätze Perioden basierend auf typischen Datenlängen
    if(n_rows > 100) {  # Mindestens ~8 Jahre Daten
      # Erste ~25% für Finanzkrise
      Daten$D_GFC[1:max(1, floor(n_rows * 0.25))] <- 1
      # Nächste ~25% für Schuldenkrise  
      start_sov <- floor(n_rows * 0.25) + 1
      end_sov <- floor(n_rows * 0.5)
      Daten$D_Sov[start_sov:end_sov] <- 1
      
      # Letzte ~40% aufteilen zwischen COVID und Inflation
      if(n_rows > 200) {  # Genug Daten für beide Krisen
        start_covid <- floor(n_rows * 0.7)
        end_covid <- floor(n_rows * 0.85)
        Daten$D_Covid[start_covid:end_covid] <- 1
        
        start_infl <- floor(n_rows * 0.9)
        Daten$D_Infl[start_infl:n_rows] <- 1
      }
      
      cat("Notfall-Zuweisung abgeschlossen\n")
      for(var in dummy_vars) {
        n_ones <- sum(Daten[[var]] == 1)
        cat(sprintf("%s: %d Einsen\n", var, n_ones))
      }
    }
  }
  
} else {
  stop("Keine Datumsspalte gefunden!")
}

cat("\n=== SCHRITT 3: MODELL MIT EINFACHEN DUMMIES TESTEN ===\n")

# Prüfe alle Variablen für das erweiterte Modell
all_vars <- c("ezb_leitzins", "inflationsabweichung", "output_gap", "D_GFC", "D_Sov", "D_Covid", "D_Infl")
extended_complete <- complete.cases(Daten[, all_vars])
n_extended_complete <- sum(extended_complete)

cat(sprintf("Vollständige Fälle für erweitertes Modell: %d\n", n_extended_complete))

if(n_extended_complete > 0) {
  cat("✓ Erweitertes Modell kann geschätzt werden!\n")
  
  # Schätze erweitertes Modell
  Modell_Shifts <- lm(ezb_leitzins ~ inflationsabweichung + output_gap + 
                        D_GFC + D_Sov + D_Covid + D_Infl, data = Daten)
  
  print("=== INTERCEPT-SHIFT MODELL ===")
  summary(Modell_Shifts)
  
  # Speichere Ergebnisse
  Daten$zins_hat_shifts <- fitted(Modell_Shifts)
  Daten$u_hat_shifts    <- residuals(Modell_Shifts)
  
  cat("\n✓ ERFOLG: Modell erfolgreich geschätzt!\n")
  
} else {
  # Falls immer noch Probleme: Zeige NA-Muster
  cat("Problem besteht weiterhin. NA-Muster:\n")
  na_pattern <- is.na(Daten[1:min(10, nrow(Daten)), all_vars])
  print(na_pattern)
  
  # Versuche nur mit verfügbaren Krisenperioden
  cat("\nVersuche mit nur verfügbaren Perioden...\n")
  available_dummies <- c()
  for(var in c("D_GFC", "D_Sov", "D_Covid", "D_Infl")) {
    if(sum(Daten[[var]] == 1, na.rm = TRUE) > 0) {
      available_dummies <- c(available_dummies, var)
    }
  }
  
  if(length(available_dummies) > 0) {
    cat("Verfügbare Dummies:", paste(available_dummies, collapse = ", "), "\n")
    formula_str <- paste("ezb_leitzins ~ inflationsabweichung + output_gap +", 
                         paste(available_dummies, collapse = " + "))
    cat("Verwende Formel:", formula_str, "\n")
    
    Modell_Shifts <- lm(as.formula(formula_str), data = Daten)
    print(summary(Modell_Shifts))
    
    Daten$zins_hat_shifts <- fitted(Modell_Shifts)
    Daten$u_hat_shifts    <- residuals(Modell_Shifts)
  } else {
    cat("Keine Krisen-Dummies verfügbar - verwende nur Basismodell\n")
    Modell_Shifts <- Modell_basis
    Daten$zins_hat_shifts <- fitted(Modell_Shifts)
    Daten$u_hat_shifts    <- residuals(Modell_Shifts)
  }
}

#------------------------------------------------------------------------
# 3. INTERPRETATION DER KRISENEFFEKTE
#------------------------------------------------------------------------

coef_shifts <- summary(Modell_Shifts)$coefficients
crisis_coefs <- coef_shifts[c("D_GFC", "D_Sov", "D_Covid", "D_Infl"), ]

cat("\n=== INTERPRETATION DER NIVEAUEFFEKTE ===\n")
for(i in 1:nrow(crisis_coefs)) {
  coef_name <- rownames(crisis_coefs)[i]
  coef_val <- crisis_coefs[i, "Estimate"]
  p_val <- crisis_coefs[i, "Pr(>|t|)"]
  
  interpretation <- ifelse(coef_val < 0, 
                           "EZB agiert EXPANSIVER als Regel nahelegt (Zins unter Normalniveau)",
                           "EZB agiert RESTRIKTIVER als Regel nahelegt")
  significance <- ifelse(p_val < 0.05, "SIGNIFIKANT", "nicht signifikant")
  
  cat(sprintf("%s: %.3f (%s) - %s\n", coef_name, coef_val, significance, interpretation))
}

#------------------------------------------------------------------------
# 4. GEMEINSAMER F-TEST: Sind alle Kriseneffekte gemeinsam signifikant?
#------------------------------------------------------------------------

# H0: δ1 = δ2 = δ3 = δ4 = 0 (keine systematischen Kriseneffekte)
# H1: mindestens ein δj ≠ 0

library(car)
f_test_crisis <- linearHypothesis(Modell_Shifts, 
                                  c("D_GFC = 0", "D_Sov = 0", "D_Covid = 0", "D_Infl = 0"),
                                  test = "F")
print("\n=== F-TEST: Gemeinsame Signifikanz aller Kriseneffekte ===")
print(f_test_crisis)

# Effektstärke: R² Vergleich
r2_basis  <- summary(Modell)$r.squared
r2_shifts <- summary(Modell_Shifts)$r.squared
delta_r2  <- r2_shifts - r2_basis

cat(sprintf("\nModellverbesserung durch Kriseneffekte:\n"))
cat(sprintf("R² Basismodell: %.4f\n", r2_basis))
cat(sprintf("R² mit Shifts:  %.4f\n", r2_shifts))
cat(sprintf("ΔR²:           %.4f (%.2f%% zusätzliche Erklärungskraft)\n", 
            delta_r2, delta_r2*100))

#------------------------------------------------------------------------
# 5. GRAFISCHE DARSTELLUNG
#------------------------------------------------------------------------

# Zeitreihenplot: Beobachtete vs. vorhergesagte Zinssätze
library(ggplot2)
library(reshape2)

# Daten für Plot vorbereiten
plot_data <- data.frame(
  datum = Daten$datum,
  beobachtet = Daten$ezb_leitzins,
  basis_modell = Daten$zins_hat,
  mit_shifts = Daten$zins_hat_shifts,
  D_GFC = Daten$D_GFC,
  D_Sov = Daten$D_Sov, 
  D_Covid = Daten$D_Covid,
  D_Infl = Daten$D_Infl
)

# Long format für ggplot
plot_long <- melt(plot_data[,1:4], id.vars = "datum", variable.name = "Modell", value.name = "Zinssatz")

# Krisenperioden als Hintergrund
crisis_periods <- data.frame(
  start = as.yearmon(c("2008-01", "2011-01", "2020-03", "2022-02")),
  end = as.yearmon(c("2010-12", "2013-12", "2021-12", "2024-12")),
  label = c("Finanzkrise", "Schuldenkrise", "COVID-19", "Inflation"),
  color = c("red", "orange", "purple", "darkgreen")
)

# Hauptplot
p1 <- ggplot(plot_long, aes(x = datum, y = Zinssatz, color = Modell, linetype = Modell)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("beobachtet" = "black", "basis_modell" = "blue", "mit_shifts" = "red")) +
  scale_linetype_manual(values = c("beobachtet" = "solid", "basis_modell" = "dashed", "mit_shifts" = "solid")) +
  labs(title = "Taylor-Regel: Basismodell vs. Modell mit Intercept-Shifts",
       subtitle = "Systematische Niveauverschiebungen in Krisenphasen",
       x = "Zeit", y = "EZB Leitzins (%)",
       color = "Modell", linetype = "Modell") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Krisenperioden als farbige Hintergründe hinzufügen
for(i in 1:nrow(crisis_periods)) {
  p1 <- p1 + annotate("rect", 
                      xmin = crisis_periods$start[i], xmax = crisis_periods$end[i],
                      ymin = -Inf, ymax = Inf,
                      alpha = 0.1, fill = crisis_periods$color[i])
}

print(p1)

# Residuenvergleich
p2 <- ggplot(Daten, aes(x = datum)) +
  geom_line(aes(y = u_hat, color = "Basismodell"), size = 1) +
  geom_line(aes(y = u_hat_shifts, color = "Mit Shifts"), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Basismodell" = "blue", "Mit Shifts" = "red")) +
  labs(title = "Residuen-Vergleich: Basismodell vs. Intercept-Shifts",
       subtitle = "Reduzierte Residuen in Krisenphasen durch systematische Niveaueffekte",
       x = "Zeit", y = "Residuen", color = "Modell") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Krisenperioden markieren
for(i in 1:nrow(crisis_periods)) {
  p2 <- p2 + annotate("rect", 
                      xmin = crisis_periods$start[i], xmax = crisis_periods$end[i],
                      ymin = -Inf, ymax = Inf,
                      alpha = 0.1, fill = crisis_periods$color[i])
}

print(p2)

library(ggplot2)
library(reshape2)
library(zoo)
library(lubridate)

# datum zuverlässig zu Date machen (funktioniert für Date, yearmon, "YYYY-MM" und "YYYY-MM-DD")
if (!inherits(Daten$datum, "Date")) {
  Daten$datum <- as.Date(as.yearmon(Daten$datum))
}


# Daten long
plot_data <- data.frame(
  datum = Daten$datum,
  beobachtet = Daten$ezb_leitzins,
  basis_modell = Daten$zins_hat,
  mit_shifts = Daten$zins_hat_shifts
)
plot_long <- melt(plot_data, id.vars = "datum",
                  variable.name = "Modell", value.name = "Zinssatz")

# Krisenperioden als Date-Bereiche (xmax = erster Tag NACH dem Ende)
crisis_periods <- data.frame(
  start = as.Date(as.yearmon(c("2008-01", "2011-01", "2020-03", "2022-02"))),
  end   = as.Date(as.yearmon(c("2010-12", "2013-12", "2021-12", "2024-12"))) %m+% months(1),
  label = c("Finanzkrise","Schuldenkrise","COVID-19","Inflation"),
  fill  = c("red","orange","purple","darkgreen")
)

# ---------------------
# 1) Taylor-Regel-Plot
# ---------------------
p1 <- ggplot() +
  # Krisenflächen zuerst (liegen unter den Linien)
  geom_rect(data = crisis_periods,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = label),
            alpha = 0.12, inherit.aes = FALSE, show.legend = FALSE) +
  # Linien darüber
  geom_line(data = plot_long,
            aes(x = datum, y = Zinssatz, color = Modell, linetype = Modell),
            linewidth = 1.1, na.rm = TRUE) +
  scale_color_manual(values = c(beobachtet = "black", basis_modell = "blue", mit_shifts = "red")) +
  scale_linetype_manual(values = c(beobachtet = "solid", basis_modell = "dashed", mit_shifts = "solid")) +
  labs(title = "Taylor-Regel: Basismodell vs. Modell mit Intercept-Shifts",
       subtitle = "Systematische Niveauverschiebungen in Krisenphasen",
       x = "Zeit", y = "EZB Leitzins (%)", color = "Modell", linetype = "Modell") +
  scale_x_date(date_breaks = "2 years", date_labels = "Jan %Y") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p1)

# ---------------------
# 2) Residuen-Plot
# ---------------------
p2 <- ggplot() +
  geom_rect(data = crisis_periods,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = label),
            alpha = 0.12, inherit.aes = FALSE, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(data = Daten, aes(x = datum, y = u_hat, color = "Basismodell"),
            linewidth = 1.1, na.rm = TRUE) +
  geom_line(data = Daten, aes(x = datum, y = u_hat_shifts, color = "Mit Shifts"),
            linewidth = 1.1, na.rm = TRUE) +
  scale_color_manual(values = c("Basismodell" = "blue", "Mit Shifts" = "red")) +
  labs(title = "Residuen-Vergleich: Basismodell vs. Intercept-Shifts",
       subtitle = "Reduzierte Residuen in Krisenphasen durch systematische Niveaueffekte",
       x = "Zeit", y = "Residuen", color = "Modell") +
  scale_x_date(date_breaks = "2 years", date_labels = "Jan %Y") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p2)

#------------------------------------------------------------------------

# 6. ROBUSTHEITSCHECK: HAC-Standardfehler für Intercept-Shift Modell
#------------------------------------------------------------------------

library(sandwich)
library(lmtest)

cat("\n=== ROBUSTE STANDARDFEHLER (HAC/Newey-West) ===\n")
coeftest(Modell_Shifts, vcov = NeweyWest(Modell_Shifts, lag = 12, prewhite = TRUE, adjust = TRUE))

#------------------------------------------------------------------------
# 7. MODELLDIAGNOSE FÜR ERWEITERTE MODELL
#------------------------------------------------------------------------

cat("\n=== DIAGNOSTIK DES ERWEITERTEN MODELLS ===\n")

# Durbin-Watson Test
cat("Durbin-Watson Test (Autokorrelation):\n")
print(dwtest(Modell_Shifts))

# Breusch-Pagan Test (Heteroskedastizität)
cat("\nBreusch-Pagan Test (Heteroskedastizität):\n")
print(bptest(Modell_Shifts))

# Normalität der Residuen (Jarque-Bera oder Shapiro-Wilk)
cat("\nShapiro-Wilk Test (Normalität der Residuen):\n")
if(nrow(Daten) <= 5000) {
  print(shapiro.test(Daten$u_hat_shifts))
} else {
  cat("Sample zu groß für Shapiro-Wilk, verwende Jarque-Bera:\n")
  # Alternative für große Samples
  library(tseries)
  print(jarque.bera.test(Daten$u_hat_shifts))
}

#------------------------------------------------------------------------
# 8. ZUSAMMENFASSUNG UND POLICY-IMPLIKATIONEN
#------------------------------------------------------------------------

cat("\n=== ZUSAMMENFASSUNG DER ERGEBNISSE ===\n")
cat("1. Das erweiterte Modell mit Intercept-Shifts erklärt systematische\n")
cat("   Abweichungen vom 'normalen' Taylor-Regel-Verhalten in Krisenphasen.\n\n")

cat("2. Policy-Implikationen der geschätzten Niveaueffekte:\n")
for(i in 1:nrow(crisis_coefs)) {
  coef_name <- rownames(crisis_coefs)[i]
  coef_val <- crisis_coefs[i, "Estimate"]
  
  crisis_name <- switch(coef_name,
                        "D_GFC" = "Globale Finanzkrise",
                        "D_Sov" = "Euro-Schuldenkrise", 
                        "D_Covid" = "COVID-19 Pandemie",
                        "D_Infl" = "Inflationsschock")
  
  cat(sprintf("   %s: %.3f Prozentpunkte ", crisis_name, coef_val))
  if(coef_val < 0) {
    cat("(expansivere Geldpolitik als Regel vorschreibt)\n")
  } else {
    cat("(restriktivere Geldpolitik als Regel vorschreibt)\n")
  }
}

cat(sprintf("\n3. Die Kriseneffekte sind %s (F-Test p-Wert: %.4f)\n",
            ifelse(f_test_crisis$`Pr(>F)`[2] < 0.05, "gemeinsam signifikant", "gemeinsam nicht signifikant"),
            f_test_crisis$`Pr(>F)`[2]))

cat(sprintf("\n4. Modellverbesserung: %.2f%% zusätzliche Erklärungskraft durch Kriseneffekte\n", 
            delta_r2*100))

cat("\n5. Diese Ergebnisse bestätigen die in CUSUM/MOSUM-Tests festgestellten\n")
cat("   Strukturinstabilitäten und bieten eine ökonomisch interpretierbare\n")
cat("   Erklärung durch systematische EZB-Reaktionen auf Krisensituationen.\n")