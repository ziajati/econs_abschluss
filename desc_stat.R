# ==============================================================================
# EZB Geldpolitik und Taylor-Regel: Interpolation und Deskriptive Analyse
# ==============================================================================

# Pakete laden und installieren falls nötig
required_packages <- c("readr", "dplyr", "tidyr", "lubridate", "ggplot2", 
                       "gridExtra", "corrplot", "forecast", "zoo", 
                       "moments", "tseries", "psych")

# Pakete installieren und laden
for(pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installiere Paket:", pkg, "\n"))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# 1. DATEN EINLESEN UND VORBEREITEN
# ==============================================================================

# CSV-Datei einlesen (Semikolon als Trennzeichen)
data <- read_delim("data/taylor_rule.csv", delim = ";", locale = locale(decimal_mark = ","))

# Datenstruktur anzeigen
str(data)
head(data)

# Datum in Date-Format konvertieren
data$date <- as.Date(paste0(data$date, "-01"))

# Variablen umbenennen für bessere Lesbarkeit und Datentypen korrigieren
colnames(data)[colnames(data) == "ezb_leitzins"] <- "interest_rate"
colnames(data)[colnames(data) == "HICP"] <- "inflation"

# Numerische Konvertierung sicherstellen
data$interest_rate <- as.numeric(data$interest_rate)
data$inflation <- as.numeric(data$inflation)
data$output_gap <- as.numeric(data$output_gap)

# TAYLOR-REGEL VARIABLEN: Inflationsabweichung berechnen
data$inflation_target <- 2.0  # EZB Zielinflation
data$inflation_gap <- data$inflation - data$inflation_target

# Leere Spalten entfernen
data <- data[, c("date", "interest_rate", "inflation", "inflation_target", "inflation_gap", "output_gap")]

cat("Datenstruktur nach Bereinigung (mit Taylor-Regel Variablen):\n")
str(data)
cat(sprintf("\nZielinflation gesetzt auf: %.1f%%\n", data$inflation_target[1]))

# ==============================================================================
# 2. OUTPUT GAP INTERPOLATION
# ==============================================================================

cat("=== OUTPUT GAP INTERPOLATION ===\n")

# Verfügbare Output Gap Werte anzeigen
output_gap_available <- data[!is.na(data$output_gap), c("date", "output_gap")]
cat("Verfügbare Output Gap Werte:\n")
print(output_gap_available)

# Verschiedene Interpolationsmethoden

# Methode 1: Lineare Interpolation
data$output_gap_linear <- na.approx(data$output_gap, na.rm = FALSE)

# Methode 2: Spline-Interpolation
data$output_gap_spline <- na.spline(data$output_gap, na.rm = FALSE)

# Methode 3: Kubische Interpolation mit zoo
data$output_gap_cubic <- na.fill(data$output_gap, "extend")
data$output_gap_cubic <- na.approx(data$output_gap, method = "linear", na.rm = FALSE)

# Für die Hauptanalyse verwenden wir die lineare Interpolation
data$output_gap_interpolated <- data$output_gap_linear

# Interpolationsergebnisse visualisieren
p_interpolation <- ggplot(data, aes(x = date)) +
  geom_point(aes(y = output_gap), color = "red", size = 3, alpha = 0.7, na.rm = TRUE) +
  geom_line(aes(y = output_gap_linear), color = "blue", linetype = "solid", size = 0.8) +
  geom_line(aes(y = output_gap_spline), color = "green", linetype = "dashed", size = 0.8) +
  labs(title = "Output Gap Interpolation: Vergleich der Methoden",
       subtitle = "Rote Punkte: Original-Daten, Blaue Linie: Linear, Grüne Linie: Spline",
       x = "Datum", y = "Output Gap (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

print(p_interpolation)

# ==============================================================================
# 3. DESKRIPTIVE STATISTIKEN
# ==============================================================================

cat("\n=== DESKRIPTIVE STATISTIKEN ===\n")

# Hauptvariablen für Taylor-Regel Analyse
main_vars <- c("interest_rate", "inflation_gap", "output_gap_interpolated")
analysis_data <- data[, c("date", "interest_rate", "inflation_gap", "output_gap_interpolated")]

cat("Taylor-Regel Variablen für Analyse:\n")
cat("- interest_rate: EZB Leitzins\n")
cat("- inflation_gap: Inflationsabweichung vom Ziel (Inflation - 2%)\n")
cat("- output_gap_interpolated: Output Gap (interpoliert)\n")

# Deskriptive Statistiken berechnen
desc_stats_list <- list()
for(var in main_vars) {
  desc_stats_list[[var]] <- describe(analysis_data[[var]])
}

# Zu einem DataFrame kombinieren
desc_stats <- do.call(rbind, desc_stats_list)
desc_stats <- round(desc_stats, 4)
rownames(desc_stats) <- c("EZB Leitzins", "Inflationsabweichung", "Output Gap")

print("Deskriptive Statistiken (Taylor-Regel Variablen):")
print(desc_stats)

# Zusammenfassende Statistiken berechnen
cat("\nBerechne zusammenfassende Statistiken für Taylor-Regel Variablen...\n")

# Prüfen der Datentypen
cat("Datentypen:\n")
for(var in main_vars) {
  cat(sprintf("%s: %s (Klasse: %s)\n", var, typeof(analysis_data[[var]]), 
              paste(class(analysis_data[[var]]), collapse = ", ")))
}

# Sicherstellen, dass alle Variablen numerisch sind
for(var in main_vars) {
  analysis_data[[var]] <- as.numeric(analysis_data[[var]])
}

# Zusammenfassende Statistiken einzeln berechnen
summary_stats <- data.frame(
  Variable = c("EZB Leitzins", "Inflationsabweichung", "Output Gap"),
  stringsAsFactors = FALSE
)

# Statistiken für jede Variable berechnen
for(i in 1:length(main_vars)) {
  var_data <- analysis_data[[main_vars[i]]]
  var_data <- var_data[!is.na(var_data)]  # NA-Werte entfernen
  
  summary_stats$Min[i] <- min(var_data)
  summary_stats$Q1[i] <- quantile(var_data, 0.25)
  summary_stats$Median[i] <- median(var_data)
  summary_stats$Mean[i] <- mean(var_data)
  summary_stats$Q3[i] <- quantile(var_data, 0.75)
  summary_stats$Max[i] <- max(var_data)
  summary_stats$SD[i] <- sd(var_data)
}

# Auf 4 Dezimalstellen runden
summary_stats[,2:8] <- round(summary_stats[,2:8], 4)

print("Zusammenfassende Statistiken (Taylor-Regel Variablen):")
print(summary_stats)

# ==============================================================================
# 4. ZEITREIHENVISUALISIERUNG
# ==============================================================================

# Einzelne Zeitreihen plotten (Taylor-Regel spezifisch)
p1 <- ggplot(data, aes(x = date, y = interest_rate)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 0.5) +
  labs(title = "EZB Leitzins", x = "Datum", y = "Zinssatz (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(data, aes(x = date, y = inflation_gap)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Inflationsabweichung (π - π*)", 
       subtitle = "Zielinflation: 2.0%",
       x = "Datum", y = "Abweichung (Prozentpunkte)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p3 <- ggplot(data, aes(x = date, y = output_gap_interpolated)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "darkgreen", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(title = "Output Gap (interpoliert)", x = "Datum", y = "Output Gap (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Zusätzlich: Original-Inflation vs. Ziel
p4 <- ggplot(data, aes(x = date)) +
  geom_line(aes(y = inflation), color = "red", size = 1) +
  geom_line(aes(y = inflation_target), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Inflation vs. Zielinflation", 
       subtitle = "Durchgezogene Linie: Tatsächliche Inflation, Gestrichelt: Ziel (2%)",
       x = "Datum", y = "Inflation (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Alle vier Grafiken kombinieren
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)

# Normalisierte Darstellung aller Taylor-Regel Variablen in einem Plot
data_normalized <- data %>%
  select(date, interest_rate, inflation_gap, output_gap_interpolated) %>%
  mutate(
    interest_rate_norm = scale(interest_rate)[,1],
    inflation_gap_norm = scale(inflation_gap)[,1],
    output_gap_norm = scale(output_gap_interpolated)[,1]
  ) %>%
  pivot_longer(cols = ends_with("_norm"), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = case_when(
    Variable == "interest_rate_norm" ~ "EZB Leitzins",
    Variable == "inflation_gap_norm" ~ "Inflationsabweichung",
    Variable == "output_gap_norm" ~ "Output Gap"
  ))

p_combined <- ggplot(data_normalized, aes(x = date, y = Value, color = Variable)) +
  geom_line(size = 1) +
  labs(title = "Taylor-Regel Variablen: Normalisierte Zeitreihen",
       subtitle = "Alle Variablen standardisiert (Mittelwert=0, SD=1)",
       x = "Datum", y = "Standardisierte Werte",
       color = "Variable") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom")

print(p_combined)

# ==============================================================================
# 5. KORRELATIONSANALYSE
# ==============================================================================

cat("\n=== KORRELATIONSANALYSE ===\n")

# Korrelationsmatrix berechnen
cor_matrix <- cor(analysis_data[, main_vars], use = "complete.obs")
print("Korrelationsmatrix:")
print(round(cor_matrix, 4))

# Korrelationsmatrix visualisieren
corrplot(cor_matrix, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.8, tl.col = "black",
         addCoef.col = "black", number.cex = 0.8,
         title = "Korrelationsmatrix der Hauptvariablen",
         mar = c(0,0,2,0))

# ==============================================================================
# 6. STATIONARITÄTSTEST VORBEREITUNG
# ==============================================================================

cat("\n=== STATIONARITÄTSEIGENSCHAFTEN (Überblick) ===\n")

# ADF-Test für alle Variablen (Augmented Dickey-Fuller Test)
library(tseries)

variables <- main_vars
adf_results <- data.frame(
  Variable = character(),
  ADF_Statistic = numeric(),
  P_Value = numeric(),
  Stationary = character(),
  stringsAsFactors = FALSE
)

for(var in variables) {
  adf_test <- adf.test(analysis_data[[var]], alternative = "stationary")
  adf_results <- rbind(adf_results, data.frame(
    Variable = var,
    ADF_Statistic = adf_test$statistic,
    P_Value = adf_test$p.value,
    Stationary = ifelse(adf_test$p.value < 0.05, "Ja", "Nein")
  ))
}

print("ADF-Test Ergebnisse (Stationarität):")
print(adf_results)

# ==============================================================================
# 7. STRUKTURBRUCH-ANALYSE VORBEREITUNG
# ==============================================================================

cat("\n=== MÖGLICHE STRUKTURBRÜCHE (visuelle Inspektion) ===\n")

# Wichtige historische Ereignisse markieren
crisis_dates <- data.frame(
  date = as.Date(c("2008-09-01", "2010-05-01", "2015-03-01", "2020-03-01", "2022-02-01")),
  event = c("Finanzkrise", "Euro-Schuldenkrise", "QE-Start", "COVID-19", "Ukraine-Krieg"),
  color = c("red", "orange", "blue", "purple", "darkred")
)

# Plot mit markierten Ereignissen
p_breaks <- ggplot(data, aes(x = date)) +
  geom_line(aes(y = interest_rate), color = "blue", size = 1) +
  geom_vline(data = crisis_dates, aes(xintercept = date, color = event), 
             linetype = "dashed", size = 0.8, alpha = 0.7) +
  labs(title = "EZB Leitzins mit möglichen Strukturbrüchen",
       x = "Datum", y = "Zinssatz (%)",
       color = "Ereignis") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom")

print(p_breaks)

# ==============================================================================
# 8. DATEN EXPORTIEREN
# ==============================================================================

# Bereinigte Daten mit Taylor-Regel Variablen speichern
final_data <- data %>%
  mutate(
    # date-Spalte in yearmon umwandeln
    datum = as.yearmon(date, format = "%Y-%m")
  ) %>%
  select(datum, interest_rate, inflation, inflation_target, inflation_gap, output_gap_interpolated) %>%
  rename(
    ezb_leitzins = interest_rate,
    hicp_inflation = inflation,
    ziel_inflation = inflation_target,
    inflationsabweichung = inflation_gap,
    output_gap = output_gap_interpolated
  )

write.csv(final_data, "taylor_rule_interpolated.csv", row.names = FALSE)

cat("\n=== ZUSAMMENFASSUNG ===\n")
cat("✓ Output Gap erfolgreich interpoliert (lineare Methode)\n")
cat("✓ Inflationsabweichung berechnet (Inflation - 2% Ziel)\n")
cat("✓ Deskriptive Statistiken für Taylor-Regel Variablen berechnet\n")
cat("✓ Zeitreihenvisualisierungen erstellt\n")
cat("✓ Korrelationsanalyse durchgeführt\n")
cat("✓ Stationaritätstests vorbereitet\n")
cat("✓ Strukturbruch-Analyse vorbereitet\n")
cat("✓ Bereinigte Daten in 'taylor_rule_interpolated.csv' gespeichert\n")
cat("\nTaylor-Regel bereit für ökonometrische Modellierung!\n")
cat("i_t = r* + π_t + α(π_t - π*) + β(y_t) + ε_t\n")

# Datensatz-Informationen ausgeben
cat(sprintf("\nDatensatz: %d Beobachtungen von %s bis %s\n", 
            nrow(final_data), 
            format(min(final_data$datum), "%m/%Y"),
            format(max(final_data$datum), "%m/%Y")))
cat(sprintf("Zielinflation: %.1f%%\n", final_data$ziel_inflation[1]))

# Missing Values Check
missing_check <- final_data %>%
  summarise_all(~sum(is.na(.)))
cat("\nFehlende Werte nach Interpolation:\n")
print(missing_check)