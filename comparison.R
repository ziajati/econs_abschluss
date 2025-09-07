# ===== Setup =====
library(zoo)
library(lmtest)     # bptest, coeftest
library(sandwich)   # robuste SE
library(car)        # VIF
library(tseries)    # Jarque-Bera
library(gt)
library(broom)

df <- read.csv("taylor_rule_interpolated.csv", stringsAsFactors = FALSE)
df$datum <- as.yearmon(df$datum)  

# Taylor-Zins nach Taylor (1993)
df$taylor_rate <- 2 + df$hicp_inflation + 0.5 * (df$inflationsabweichung) + 0.5 * df$output_gap

# Plot
library(ggplot2)
ggplot(df, aes(x = datum)) +
  geom_line(aes(y = ezb_leitzins, color = "EZB-Leitzins"), size = 1) +
  geom_line(aes(y = taylor_rate, color = "Taylor-Zins (geschätzt)"), size = 1, linetype = "dashed") +
  labs(title = "EZB-Leitzins vs. Taylor-Zins",
       x = "Zeit", y = "Zinssatz (%)") +
  scale_color_manual(values = c("EZB-Leitzins" = "blue", "Taylor-Zins (geschätzt)" = "red")) +
  theme_minimal()
