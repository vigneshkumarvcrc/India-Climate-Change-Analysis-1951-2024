library(dplyr)

out_dir <- "/Users/khkr/Documents/New project/CCPOP_final"
national <- read.csv(file.path(out_dir, "tables/national_annual_climate_1951_2024.csv"))
period_levels <- c(
  "1951-1959", "1960-1969", "1970-1979", "1980-1989",
  "1990-1999", "2000-2009", "2010-2019", "2020-2024"
)
national$Period <- factor(national$Period, levels = period_levels)

decadal <- national |>
  group_by(Period) |>
  summarise(
    Start_year = min(Year),
    End_year = max(Year),
    N_years = n(),
    Temp_SD_C = stats::sd(Tmean_C),
    Temp_CV_percent = 100 * stats::sd(Tmean_C) / (mean(Tmean_C) + 273.15),
    Rain_SD_mm = stats::sd(Rainfall_mm),
    Rain_CV_percent = 100 * stats::sd(Rainfall_mm) / mean(Rainfall_mm),
    Tmin_C = mean(Tmin_C),
    Tmean_C = mean(Tmean_C),
    Tmax_C = mean(Tmax_C),
    Rainfall_mm = mean(Rainfall_mm),
    Temp_Anomaly_C = mean(Temp_Anomaly_C),
    Rain_Anomaly_mm = mean(Rain_Anomaly_mm),
    .groups = "drop"
  ) |>
  mutate(
    Period = as.character(Period),
    Temp_Change_from_previous_C = c(NA, diff(Tmean_C)),
    Rain_Change_from_previous_mm = c(NA, diff(Rainfall_mm)),
    Temp_Change_from_1951_1959_C = Tmean_C - first(Tmean_C),
    Rain_Change_from_1951_1959_mm = Rainfall_mm - first(Rainfall_mm),
    Temp_Acceleration_C_per_decade2 = c(NA, NA, diff(diff(Tmean_C))),
    Rain_Acceleration_mm_per_decade2 = c(NA, NA, diff(diff(Rainfall_mm)))
  )

write.csv(decadal, file.path(out_dir, "tables/national_decadal_climate.csv"), row.names = FALSE)
cat("Decadal table repaired\n")
