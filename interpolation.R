library(ggplot2)

# 1) CSV lesen
df <- read.csv("taylor_regel_daten.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, check.names = FALSE)

# 2) Erste 4 Spalten robust umbenennen (egal wie sie vorher hießen)
#    Erwartung: 1=JahrMonat, 2=EZB Leitzins, 3=HICP, 4=Output Gap
stopifnot(ncol(df) >= 4)
names(df)[1:4] <- c("DateYM", "EZB_Leitzins", "HICP", "OutputGap")

# 3) Datum bauen (aus "YYYY-MM" -> YYYY-MM-01) und numerische Spalten erzwingen
df$Datum <- as.Date(paste0(df$DateYM, "-01"))
df$EZB_Leitzins <- as.numeric(df$EZB_Leitzins)
df$HICP        <- as.numeric(df$HICP)
df$OutputGap   <- as.numeric(df$OutputGap)

# 4) Lineare Interpolation (inkl. Extrapolation an den Rändern)
known_x <- which(!is.na(df$OutputGap))
if (length(known_x) >= 2) {
  interp <- approx(
    x = known_x,
    y = df$OutputGap[known_x],
    xout = 1:nrow(df),
    method = "linear",
    rule = 2
  )
  df$OutputGap <- interp$y
}

# 5) Plots
ggplot(df, aes(x = Datum, y = EZB_Leitzins)) +
  geom_line() +
  labs(title = "EZB Leitzins", x = NULL, y = "%") +
  theme_minimal()

ggplot(df, aes(x = Datum, y = HICP)) +
  geom_line() +
  labs(title = "HICP", x = NULL, y = "% J/J") +
  theme_minimal()

ggplot(df, aes(x = Datum, y = OutputGap)) +
  geom_line() +
  labs(title = "Output Gap (linear interpoliert)", x = NULL, y = "Prozentpunkte") +
  theme_minimal()
