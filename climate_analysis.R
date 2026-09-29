library(terra)
library(ggplot2)
library(dplyr)
library(scales)
library(viridisLite)

options(scipen = 999)
terraOptions(memfrac = 0.75, progress = 1)

base_dir <- "/Volumes/K1/climate_data"
out_dir <- "/Users/khkr/Documents/New project/CCPOP_final"
dirs <- c("tables", "charts", "maps", "rasters", "shapefiles", "manuscript", "qa")
invisible(lapply(file.path(out_dir, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

period_labels <- c(
  "1951-1959", "1960-1969", "1970-1979", "1980-1989",
  "1990-1999", "2000-2009", "2010-2019", "2020-2024"
)

read_sort_stack <- function(path) {
  r <- rast(path)
  yrs <- as.integer(sub("^X", "", names(r)))
  ord <- order(yrs)
  r <- r[[ord]]
  yrs <- yrs[ord]
  names(r) <- paste0("Y", yrs)
  if (!identical(yrs, 1951:2024)) {
    stop("Unexpected year coverage in ", path, ": ", paste(range(yrs), collapse = "-"))
  }
  list(r = r, years = yrs)
}

tmin_obj <- read_sort_stack(file.path(base_dir, "tmin/tmin_yearly_stack.tif"))
tmax_obj <- read_sort_stack(file.path(base_dir, "tmax/tmax_yearly_stack.tif"))
tmean_obj <- read_sort_stack(file.path(base_dir, "tmean/tmean_yearly_stack.tif"))
rain_obj <- read_sort_stack(file.path(base_dir, "precip/rain_yearly_stack.tif"))

years <- tmean_obj$years
tmin <- tmin_obj$r
tmax <- tmax_obj$r
tmean <- tmean_obj$r
rain <- rain_obj$r

india <- project(vect(file.path(base_dir, "india_boundary.shp")), crs(tmean))
states <- project(vect(file.path(base_dir, "Ind_st.shp")), crs(tmean))
districts <- project(vect(file.path(base_dir, "Ind_Dist.shp")), crs(tmean))
districts <- districts[toupper(trimws(districts$DIST)) != "ADJACENT COUNTRIES", ]

area <- cellSize(tmean[[1]], unit = "km")
area_mask <- ifel(!is.na(tmean[[1]]), area, NA)
denom <- global(area_mask, "sum", na.rm = TRUE)[1, 1]
area_weighted_series <- function(r) {
  as.numeric(global(r * area_mask, "sum", na.rm = TRUE)[, 1] / denom)
}

national <- data.frame(
  Year = years,
  Tmin_C = area_weighted_series(tmin),
  Tmean_C = area_weighted_series(tmean),
  Tmax_C = area_weighted_series(tmax),
  Rainfall_mm = area_weighted_series(rain)
)

baseline <- national$Year >= 1951 & national$Year <= 1980
national <- national |>
  mutate(
    Temp_Anomaly_C = Tmean_C - mean(Tmean_C[baseline]),
    Rain_Anomaly_mm = Rainfall_mm - mean(Rainfall_mm[baseline]),
    Temp_Annual_Change_C = c(NA, diff(Tmean_C)),
    Rain_Annual_Change_mm = c(NA, diff(Rainfall_mm))
  )

linear_summary <- function(y, unit) {
  fit <- lm(y ~ years)
  sm <- summary(fit)
  data.frame(
    Variable = unit,
    Start_year = min(years),
    End_year = max(years),
    Trend_per_year = unname(coef(fit)[2]),
    Trend_per_decade = unname(coef(fit)[2]) * 10,
    Total_fitted_change = unname(coef(fit)[2]) * (max(years) - min(years)),
    R_squared = sm$r.squared,
    P_value = sm$coefficients[2, 4]
  )
}

trend_summary <- bind_rows(
  linear_summary(national$Tmin_C, "Tmin (deg C)"),
  linear_summary(national$Tmean_C, "Tmean (deg C)"),
  linear_summary(national$Tmax_C, "Tmax (deg C)"),
  linear_summary(national$Rainfall_mm, "Rainfall (mm)")
)

decade_for_year <- function(y) {
  cut(
    y,
    breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2025),
    labels = period_labels,
    right = FALSE
  )
}
national$Period <- decade_for_year(national$Year)

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

rolling_slope <- function(y, window = 15) {
  n <- length(y)
  out <- rep(NA_real_, n)
  half <- floor(window / 2)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    if ((hi - lo + 1) >= 10) out[i] <- coef(lm(y[lo:hi] ~ years[lo:hi]))[2]
  }
  out
}

rolling <- data.frame(
  Year = years,
  Temperature_slope_C_per_decade = rolling_slope(national$Tmean_C) * 10,
  Rainfall_slope_mm_per_decade = rolling_slope(national$Rainfall_mm) * 10
)

temp_extreme_threshold <- quantile(national$Tmax_C[baseline], 0.90, na.rm = TRUE)
rain_extreme_threshold <- quantile(national$Rainfall_mm[baseline], 0.90, na.rm = TRUE)
national$Extreme_Warm_Year <- national$Tmax_C > temp_extreme_threshold
national$Extreme_Wet_Year <- national$Rainfall_mm > rain_extreme_threshold

temporal_risk <- national |>
  group_by(Period) |>
  summarise(
    Years = n(),
    Mean_temp_anomaly_C = mean(Temp_Anomaly_C),
    Mean_rain_anomaly_mm = mean(Rain_Anomaly_mm),
    Temp_variability_SD_C = stats::sd(Tmean_C),
    Rain_variability_CV_percent = 100 * stats::sd(Rainfall_mm) / mean(Rainfall_mm),
    Extreme_warm_years = sum(Extreme_Warm_Year),
    Extreme_wet_years = sum(Extreme_Wet_Year),
    .groups = "drop"
  ) |>
  mutate(
    Period = as.character(Period),
    Warm_frequency = Extreme_warm_years / Years,
    Wet_frequency = Extreme_wet_years / Years,
    Temporal_risk_index = rescale(pmax(Mean_temp_anomaly_C, 0)) +
      rescale(abs(Mean_rain_anomaly_mm)) +
      rescale(Temp_variability_SD_C) +
      rescale(Rain_variability_CV_percent) +
      rescale(Warm_frequency) +
      rescale(Wet_frequency),
    Temporal_risk_index = Temporal_risk_index / 6,
    Risk_class = cut(
      Temporal_risk_index,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c("Very low", "Low", "Moderate", "High", "Very high")
    )
  )

write.csv(national, file.path(out_dir, "tables/national_annual_climate_1951_2024.csv"), row.names = FALSE)
write.csv(decadal, file.path(out_dir, "tables/national_decadal_climate.csv"), row.names = FALSE)
write.csv(trend_summary, file.path(out_dir, "tables/national_trend_summary.csv"), row.names = FALSE)
write.csv(rolling, file.path(out_dir, "tables/national_rolling_trends.csv"), row.names = FALSE)
write.csv(temporal_risk, file.path(out_dir, "tables/temporal_climate_risk_by_period.csv"), row.names = FALSE)

slope_fun <- function(time) {
  force(time)
  function(x) {
    ok <- is.finite(x)
    if (sum(ok) < 30) return(NA_real_)
    coef(lm(x[ok] ~ time[ok]))[2]
  }
}

pvalue_fun <- function(time) {
  force(time)
  function(x) {
    ok <- is.finite(x)
    if (sum(ok) < 30) return(NA_real_)
    summary(lm(x[ok] ~ time[ok]))$coefficients[2, 4]
  }
}

spatial_metrics <- function(r, prefix, is_temperature = FALSE) {
  early_idx <- years <= 1989
  late_idx <- years >= 1990
  baseline_idx <- years <= 1980
  recent_idx <- years >= 1995
  first_period <- years <= 1959
  last_period <- years >= 2020

  trend <- app(r, slope_fun(years))
  trend_decade <- trend * 10
  pvalue <- app(r, pvalue_fun(years))
  early_slope <- app(r[[early_idx]], slope_fun(years[early_idx]))
  late_slope <- app(r[[late_idx]], slope_fun(years[late_idx]))
  acceleration <- (late_slope - early_slope) / (mean(years[late_idx]) - mean(years[early_idx]))
  variability_sd <- app(r, sd, na.rm = TRUE)
  mean_raster <- app(r, mean, na.rm = TRUE)
  variability_cv <- if (is_temperature) {
    variability_sd / (mean_raster + 273.15) * 100
  } else {
    variability_sd / mean_raster * 100
  }
  change <- app(r[[last_period]], mean, na.rm = TRUE) - app(r[[first_period]], mean, na.rm = TRUE)

  threshold <- app(r[[baseline_idx]], function(x) quantile(x, 0.90, na.rm = TRUE))
  baseline_freq <- app(r[[baseline_idx]], function(x) mean(x > quantile(x, 0.90, na.rm = TRUE), na.rm = TRUE))
  recent_exceed <- sum(r[[recent_idx]] > threshold, na.rm = TRUE) / sum(recent_idx)
  extreme_frequency_change_decade <- (recent_exceed - baseline_freq) * 10
  extreme_count <- sum(r > threshold, na.rm = TRUE)

  normalize_raster <- function(x, use_abs = FALSE, positive_only = FALSE) {
    if (use_abs) x <- abs(x)
    if (positive_only) x <- ifel(x < 0, 0, x)
    qs <- as.numeric(global(x, function(v) quantile(v, c(0.02, 0.98), na.rm = TRUE)))[1:2]
    y <- (x - qs[1]) / (qs[2] - qs[1])
    clamp(y, 0, 1)
  }

  if (is_temperature) {
    risk <- (
      normalize_raster(trend_decade, positive_only = TRUE) +
      normalize_raster(variability_sd) +
      normalize_raster(acceleration, positive_only = TRUE) +
      normalize_raster(extreme_frequency_change_decade, positive_only = TRUE)
    ) / 4
  } else {
    risk <- (
      normalize_raster(trend_decade, use_abs = TRUE) +
      normalize_raster(variability_cv) +
      normalize_raster(acceleration, use_abs = TRUE) +
      normalize_raster(extreme_frequency_change_decade, positive_only = TRUE)
    ) / 4
  }
  risk_class <- floor(clamp(risk, 0, 0.9999) * 5) + 1

  metrics <- c(
    trend, trend_decade, pvalue, acceleration, variability_sd, variability_cv,
    change, extreme_count, extreme_frequency_change_decade, risk, risk_class
  )
  names(metrics) <- c(
    paste0(prefix, "_trend_per_year"),
    paste0(prefix, "_trend_per_decade"),
    paste0(prefix, "_trend_pvalue"),
    paste0(prefix, "_acceleration_per_year2"),
    paste0(prefix, "_variability_sd"),
    paste0(prefix, "_variability_cv_percent"),
    paste0(prefix, "_change_1951_1959_to_2020_2024"),
    paste0(prefix, "_extreme_year_count"),
    paste0(prefix, "_extreme_frequency_change_per_decade"),
    paste0(prefix, "_risk_index"),
    paste0(prefix, "_risk_class")
  )
  metrics
}

temp_metrics <- spatial_metrics(tmean, "temp", TRUE)
rain_metrics <- spatial_metrics(rain, "rain", FALSE)

writeRaster(temp_metrics, file.path(out_dir, "rasters/temperature_spatial_metrics_1951_2024.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
writeRaster(rain_metrics, file.path(out_dir, "rasters/rainfall_spatial_metrics_1951_2024.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))

for (i in seq_len(nlyr(temp_metrics))) {
  writeRaster(temp_metrics[[i]], file.path(out_dir, "rasters", paste0(names(temp_metrics)[i], ".tif")),
              overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
}
for (i in seq_len(nlyr(rain_metrics))) {
  writeRaster(rain_metrics[[i]], file.path(out_dir, "rasters", paste0(names(rain_metrics)[i], ".tif")),
              overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
}

decadal_raster <- function(r) {
  layers <- lapply(period_labels, function(p) {
    idx <- national$Period == p
    app(r[[idx]], mean, na.rm = TRUE)
  })
  x <- rast(layers)
  names(x) <- gsub("-", "_", period_labels)
  x
}

temp_decadal_mean <- decadal_raster(tmean)
rain_decadal_mean <- decadal_raster(rain)
temp_decadal_change <- temp_decadal_mean - temp_decadal_mean[[1]]
rain_decadal_change <- rain_decadal_mean - rain_decadal_mean[[1]]
names(temp_decadal_change) <- paste0("temp_change_", gsub("-", "_", period_labels))
names(rain_decadal_change) <- paste0("rain_change_", gsub("-", "_", period_labels))

writeRaster(temp_decadal_mean, file.path(out_dir, "rasters/temperature_decadal_mean.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
writeRaster(rain_decadal_mean, file.path(out_dir, "rasters/rainfall_decadal_mean.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
writeRaster(temp_decadal_change, file.path(out_dir, "rasters/temperature_decadal_change_from_1951_1959.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
writeRaster(rain_decadal_change, file.path(out_dir, "rasters/rainfall_decadal_change_from_1951_1959.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))

extract_metrics <- function(v, id_fields) {
  temp_ex <- extract(temp_metrics, v, fun = mean, na.rm = TRUE)
  rain_ex <- extract(rain_metrics, v, fun = mean, na.rm = TRUE)
  attrs <- as.data.frame(v)[, id_fields, drop = FALSE]
  bind_cols(attrs, temp_ex[, -1, drop = FALSE], rain_ex[, -1, drop = FALSE])
}

state_table <- extract_metrics(states, c("State", "Region"))
district_table <- extract_metrics(districts, c("STATE", "DIST"))

state_table$combined_risk_index <- rowMeans(
  state_table[, c("temp_risk_index", "rain_risk_index")], na.rm = TRUE
)
district_table$combined_risk_index <- rowMeans(
  district_table[, c("temp_risk_index", "rain_risk_index")], na.rm = TRUE
)
classify_risk <- function(x) {
  cut(x, quantile(x, probs = seq(0, 1, 0.2), na.rm = TRUE),
      include.lowest = TRUE, labels = c("Very low", "Low", "Moderate", "High", "Very high"))
}
state_table$combined_risk_class <- classify_risk(state_table$combined_risk_index)
district_table$combined_risk_class <- classify_risk(district_table$combined_risk_index)

write.csv(state_table, file.path(out_dir, "tables/state_spatial_climate_metrics.csv"), row.names = FALSE)
write.csv(district_table, file.path(out_dir, "tables/district_spatial_climate_metrics.csv"), row.names = FALSE)

states_out <- states
districts_out <- districts
values(states_out) <- state_table
values(districts_out) <- district_table
writeVector(states_out, file.path(out_dir, "shapefiles/india_state_climate_metrics.shp"), overwrite = TRUE)
writeVector(districts_out, file.path(out_dir, "shapefiles/india_district_climate_metrics.shp"), overwrite = TRUE)

low_risk_districts <- district_table |>
  arrange(combined_risk_index) |>
  slice_head(n = 50)
high_risk_districts <- district_table |>
  arrange(desc(combined_risk_index)) |>
  slice_head(n = 50)
write.csv(low_risk_districts, file.path(out_dir, "tables/lowest_relative_risk_districts.csv"), row.names = FALSE)
write.csv(high_risk_districts, file.path(out_dir, "tables/highest_relative_risk_districts.csv"), row.names = FALSE)

theme_pub <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey35"),
    legend.position = "bottom"
  )
save_chart <- function(p, name, width = 11, height = 7) {
  ggsave(file.path(out_dir, "charts", paste0(name, ".png")), p,
         width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "charts", paste0(name, ".pdf")), p,
         width = width, height = height, device = "pdf", bg = "white")
}

p_temp_year <- ggplot(national, aes(Year)) +
  geom_line(aes(y = Tmin_C, color = "Minimum"), linewidth = 0.65) +
  geom_line(aes(y = Tmean_C, color = "Mean"), linewidth = 0.9) +
  geom_line(aes(y = Tmax_C, color = "Maximum"), linewidth = 0.65) +
  scale_color_manual(values = c(Minimum = "#2C7BB6", Mean = "#D7191C", Maximum = "#7F0000")) +
  labs(title = "India annual temperature series, 1951-2024",
       subtitle = "Area-weighted national means from chronologically sorted gridded layers",
       y = "Temperature (deg C)", color = NULL) + theme_pub
save_chart(p_temp_year, "01_temperature_yearly_min_mean_max")

p_temp_anom <- ggplot(national, aes(Year, Temp_Anomaly_C)) +
  geom_col(aes(fill = Temp_Anomaly_C >= 0), width = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
  labs(title = "Annual mean-temperature anomaly",
       subtitle = "Relative to the 1951-1980 baseline", y = "Anomaly (deg C)") + theme_pub
save_chart(p_temp_anom, "02_temperature_annual_anomaly")

p_temp_trend <- ggplot(national, aes(Year, Tmean_C)) +
  geom_line(color = "#3B3B3B", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#B2182B", fill = "#F4A6A6", linewidth = 1) +
  labs(title = "India mean-temperature trend, 1951-2024", y = "Mean temperature (deg C)") + theme_pub
save_chart(p_temp_trend, "03_temperature_trend_1951_2024")

p_temp_dec <- ggplot(decadal, aes(Period, Temp_Change_from_1951_1959_C)) +
  geom_col(aes(fill = Temp_Change_from_1951_1959_C), width = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, guide = "none") +
  labs(title = "Decadal warming relative to 1951-1959", x = NULL, y = "Temperature change (deg C)") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_chart(p_temp_dec, "04_decadal_warming_bar")

p_temp_acc <- ggplot(rolling, aes(Year, Temperature_slope_C_per_decade)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_line(color = "#B2182B", linewidth = 1) +
  labs(title = "Temperature trend acceleration and deceleration",
       subtitle = "Centered 15-year rolling linear slope", y = "Local trend (deg C/decade)") + theme_pub
save_chart(p_temp_acc, "05_temperature_acceleration_rolling_slope")

p_rain_trend <- ggplot(national, aes(Year, Rainfall_mm)) +
  geom_line(color = "#2166AC", linewidth = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "#08306B", fill = "#9ECAE1", linewidth = 1) +
  labs(title = "India annual rainfall trend, 1951-2024", y = "Annual rainfall (mm)") + theme_pub
save_chart(p_rain_trend, "06_rainfall_trend_1951_2024")

p_rain_anom <- ggplot(national, aes(Year, Rain_Anomaly_mm)) +
  geom_col(aes(fill = Rain_Anomaly_mm >= 0), width = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B35806"), guide = "none") +
  labs(title = "Annual rainfall anomaly",
       subtitle = "Relative to the 1951-1980 baseline", y = "Rainfall anomaly (mm)") + theme_pub
save_chart(p_rain_anom, "07_rainfall_annual_anomaly")

p_rain_dec <- ggplot(decadal, aes(Period, Rain_Change_from_1951_1959_mm)) +
  geom_col(aes(fill = Rain_Change_from_1951_1959_mm), width = 0.75) +
  scale_fill_gradient2(low = "#B35806", mid = "white", high = "#2166AC", midpoint = 0, guide = "none") +
  labs(title = "Decadal rainfall change relative to 1951-1959", x = NULL, y = "Rainfall change (mm)") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_chart(p_rain_dec, "08_decadal_rainfall_change_bar")

p_rain_acc <- ggplot(rolling, aes(Year, Rainfall_slope_mm_per_decade)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_line(color = "#2166AC", linewidth = 1) +
  labs(title = "Rainfall trend intensification and weakening",
       subtitle = "Centered 15-year rolling linear slope", y = "Local trend (mm/decade)") + theme_pub
save_chart(p_rain_acc, "09_rainfall_acceleration_rolling_slope")

p_risk <- temporal_risk |>
  ggplot(aes(Period, Temporal_risk_index, fill = Risk_class)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c(
    "Very low" = "#1A9850", "Low" = "#91CF60", "Moderate" = "#FEE08B",
    "High" = "#FC8D59", "Very high" = "#D73027"
  ), drop = FALSE) +
  labs(title = "Relative temporal climate-risk index by period",
       subtitle = "Composite of anomalies, variability, and extreme-year frequencies",
       x = NULL, y = "Relative risk index (0-1)", fill = "Class") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_chart(p_risk, "10_temporal_climate_risk")

map_png <- function(r, filename, title, unit, palette, zlim = NULL, categorical = FALSE) {
  png(file.path(out_dir, "maps", filename), width = 4200, height = 3300, res = 450, bg = "white")
  par(mar = c(2.5, 2.5, 4.2, 5.2), family = "sans")
  if (categorical) {
    plot(r, col = palette, axes = FALSE, main = title, plg = list(title = unit))
  } else {
    plot(r, col = palette, axes = FALSE, main = title, range = zlim, plg = list(title = unit))
  }
  plot(states, add = TRUE, border = adjustcolor("grey15", 0.55), lwd = 0.45)
  plot(india, add = TRUE, border = "black", lwd = 1)
  dev.off()
}

div_temp <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
div_rain <- colorRampPalette(c("#B35806", "#F7F7F7", "#2166AC"))(101)
seq_var <- viridis(101, option = "C")
risk_pal <- c("#1A9850", "#91CF60", "#FEE08B", "#FC8D59", "#D73027")

map_png(temp_metrics[["temp_trend_per_year"]], "temperature_trend_per_year.png",
        "Mean-temperature trend, 1951-2024", "deg C/year", div_temp)
map_png(temp_metrics[["temp_trend_per_decade"]], "temperature_trend_per_decade.png",
        "Mean-temperature trend, 1951-2024", "deg C/decade", div_temp)
map_png(temp_metrics[["temp_acceleration_per_year2"]], "temperature_acceleration.png",
        "Temperature trend acceleration", "deg C/year2", div_temp)
map_png(temp_metrics[["temp_variability_cv_percent"]], "temperature_variability_cv.png",
        "Interannual temperature coefficient of variation", "CV (%)", seq_var)
map_png(temp_metrics[["temp_change_1951_1959_to_2020_2024"]], "temperature_overall_change.png",
        "Temperature change: 2020-2024 minus 1951-1959", "deg C", div_temp)
map_png(temp_metrics[["temp_extreme_frequency_change_per_decade"]], "temperature_extreme_frequency_change.png",
        "Change in extreme warm-year frequency", "events/decade", div_temp)
map_png(temp_metrics[["temp_risk_class"]], "temperature_risk_class.png",
        "Relative temperature-risk class", "1 very low - 5 very high", risk_pal, categorical = TRUE)

map_png(rain_metrics[["rain_trend_per_year"]], "rainfall_trend_per_year.png",
        "Annual-rainfall trend, 1951-2024", "mm/year", div_rain)
map_png(rain_metrics[["rain_trend_per_decade"]], "rainfall_trend_per_decade.png",
        "Annual-rainfall trend, 1951-2024", "mm/decade", div_rain)
map_png(rain_metrics[["rain_acceleration_per_year2"]], "rainfall_acceleration.png",
        "Rainfall trend acceleration", "mm/year2", div_rain)
map_png(rain_metrics[["rain_variability_cv_percent"]], "rainfall_variability_cv.png",
        "Interannual rainfall coefficient of variation", "CV (%)", seq_var)
map_png(rain_metrics[["rain_change_1951_1959_to_2020_2024"]], "rainfall_overall_change.png",
        "Rainfall change: 2020-2024 minus 1951-1959", "mm", div_rain)
map_png(rain_metrics[["rain_extreme_frequency_change_per_decade"]], "rainfall_extreme_frequency_change.png",
        "Change in extreme wet-year frequency", "events/decade", div_rain)
map_png(rain_metrics[["rain_risk_class"]], "rainfall_risk_class.png",
        "Relative rainfall-risk class", "1 very low - 5 very high", risk_pal, categorical = TRUE)

for (i in seq_len(nlyr(temp_decadal_change))) {
  map_png(temp_decadal_change[[i]],
          paste0("temperature_change_", period_labels[i], ".png"),
          paste0("Temperature change relative to 1951-1959: ", period_labels[i]),
          "deg C", div_temp)
}
for (i in seq_len(nlyr(rain_decadal_change))) {
  map_png(rain_decadal_change[[i]],
          paste0("rainfall_change_", period_labels[i], ".png"),
          paste0("Rainfall change relative to 1951-1959: ", period_labels[i]),
          "mm", div_rain)
}

manifest <- data.frame(
  file = list.files(out_dir, recursive = TRUE, full.names = FALSE),
  stringsAsFactors = FALSE
)
manifest$bytes <- file.info(file.path(out_dir, manifest$file))$size
write.csv(manifest, file.path(out_dir, "output_manifest.csv"), row.names = FALSE)

cat("Climate publication analysis complete:", out_dir, "\n")
