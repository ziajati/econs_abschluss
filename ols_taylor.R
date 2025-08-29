# ===== Setup =====
library(zoo)
library(lmtest)     # bptest, coeftest
library(sandwich)   # robuste SE
library(car)        # VIF
library(tseries)    # Jarque-Bera

# 1) Daten laden & Index sauber setzen
df <- read.csv("taylor_rule_interpolated.csv", stringsAsFactors = FALSE)
df$datum <- as.yearmon(df$datum, format = "%Y-%m")

# 2) OLS: Basismodell (statische Taylor-Regel)
model_ols <- lm(ezb_leitzins ~ inflationsabweichung + output_gap, data = df)
summary(model_ols)

# ===== Standard-Checks =====

## A) Heteroskedastizität
bp <- bptest(model_ols)  # Breusch-Pagan
bp

# (optional White-Test-Variante mit Quadraten + Interaktion)
white <- bptest(model_ols, ~ inflationsabweichung + output_gap +
                  I(inflationsabweichung^2) + I(output_gap^2) +
                  inflationsabweichung:output_gap)
white

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

library(stats)    # acf, pacf, Box.test

# 1) Durbin–Watson (Lag-1-Autokorrelation)
dwtest(model_ols)
# DW = 0.02 (!!), p < 0.001 → extrem starke positive Autokorrelation.

# 2) Breusch–Godfrey (höhere Ordnungen; passe 'order' an, z.B. 4, 8, 12)
bgtest(model_ols, order = 4)
# LM test = 272.94, df = 4, p-value < 2.2e-16 -> Autokorrelation
bgtest(model_ols, order = 8)
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
# (z.B. 12 Lags für Monatsdaten ~ 1 Jahr)
Box.test(res, lag = 12, type = "Ljung-Box")

library(dynlm)
# ===== Daten laden & Zeitachse setzen =====
df <- read.csv("taylor_rule_interpolated.csv", stringsAsFactors = FALSE)

df_ts <- zoo(df[, c("ezb_leitzins","inflationsabweichung","output_gap")], order.by = df$datum)
# ===== Dynamisches Δ-Modell (Taylor-Regel in Änderungen) =====
# Δ i_t = a + b1*infl_gap_t + b2*output_gap_t + phi*Δ i_{t-1} + u_t
model_dyn <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                     L(diff(ezb_leitzins), 1),
                   data = df_ts)

# OLS-Output
summary(model_dyn)

# Robuste (HAC) Standardfehler für Zeitreihen – 12 Lags ~ 1 Jahr bei Monatsdaten
nw <- NeweyWest(model_dyn, lag = 12, prewhite = TRUE, adjust = TRUE)
coeftest(model_dyn, vcov = nw)

# ===== Kern-Diagnosen =====
# Autokorrelation
dwtest(model_dyn)
bgtest(model_dyn, order = 12)                      # bis 12 Lags
Box.test(resid(model_dyn), lag = 12, type = "Ljung-Box")

# Heteroskedastizität (Info, da wir ohnehin HAC-SE reporten)
bptest(model_dyn)

# Normalität (kurzer Check)
jarque.bera.test(na.omit(resid(model_dyn)))

# (Optional, hübsche Grafik)
par(mfrow = c(1,2))
acf(resid(model_dyn),  main = "ACF Residuen (Δ-Modell)")
pacf(resid(model_dyn), main = "PACF Residuen (Δ-Modell)")
par(mfrow = c(1,1))

# Mehrere Lags der Zinsänderung
model_dyn2 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1:2),
                    data = df_ts)

model_dyn3 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1:3),
                    data = df_ts)

model_dyn4 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                      L(diff(ezb_leitzins), 1:4),
                    data = df_ts)

# HAC-Standardfehler (Newey-West) für Zeitreihen
nw2 <- NeweyWest(model_dyn2, lag = 12, prewhite = TRUE, adjust = TRUE)
nw3 <- NeweyWest(model_dyn3, lag = 12, prewhite = TRUE, adjust = TRUE)
nw4 <- NeweyWest(model_dyn4, lag = 12, prewhite = TRUE, adjust = TRUE)

coeftest(model_dyn2, vcov = nw2)
coeftest(model_dyn3, vcov = nw3)
coeftest(model_dyn4, vcov = nw4)

# Diagnose: Autokorrelation
bgtest(model_dyn2, order = 12)
bgtest(model_dyn3, order = 12)
bgtest(model_dyn4, order = 12)

Box.test(resid(model_dyn2), lag = 12, type = "Ljung-Box")
Box.test(resid(model_dyn3), lag = 12, type = "Ljung-Box")
Box.test(resid(model_dyn4), lag = 12, type = "Ljung-Box")

# Modellvergleich: AIC/BIC
AIC(model_dyn, model_dyn2, model_dyn3, model_dyn4)
BIC(model_dyn, model_dyn2, model_dyn3, model_dyn4)

# Modell mit 12 Lags der Δi
model_dyn12 <- dynlm(diff(ezb_leitzins) ~ inflationsabweichung + output_gap +
                       L(diff(ezb_leitzins), 1:12),
                     data = df_ts)

# Robuste SE (Newey-West)
nw12 <- NeweyWest(model_dyn12, lag = 12, prewhite = TRUE, adjust = TRUE)
coeftest(model_dyn12, vcov = nw12)

# Autokorrelationstests
bgtest(model_dyn12, order = 12)
Box.test(resid(model_dyn12), lag = 12, type = "Ljung-Box")

# Modellvergleich
AIC(model_dyn4, model_dyn12)
BIC(model_dyn4, model_dyn12)

# 4) Residual-Diagnostik (ACF/PACF + Ljung–Box p-Wert) ---------------------
res_z <- zoo(residuals(model_dyn12), order.by = index(fit_d))   # gleiche Länge wie fitted
par(mfrow = c(1,2))
acf(res_z, main = "ACF Residuen (model_dyn12)")
pacf(res_z, main = "PACF Residuen (model_dyn12)")
par(mfrow = c(1,1))

lb <- Box.test(as.numeric(res_z), lag = 12, type = "Ljung-Box")
cat(sprintf("Ljung-Box(12): X^2=%.2f, p=%.3g\n", lb$statistic, lb$p.value))


