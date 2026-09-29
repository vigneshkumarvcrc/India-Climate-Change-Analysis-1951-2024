library(terra)
library(dplyr)

base <- "/Volumes/K1/climate_data"
out <- "/Users/khkr/Documents/New project/CCPOP_final"

temp <- rast(file.path(out, "rasters/temperature_spatial_metrics_1951_2024.tif"))
rain <- rast(file.path(out, "rasters/rainfall_spatial_metrics_1951_2024.tif"))
area <- cellSize(temp[[1]], unit = "km")

area_share <- function(condition) {
  num <- global(ifel(condition, area, NA), "sum", na.rm = TRUE)[1, 1]
  den <- global(ifel(!is.na(condition), area, NA), "sum", na.rm = TRUE)[1, 1]
  100 * num / den
}

spatial_summary <- data.frame(
  Metric = c(
    "Area with warming trend",
    "Area with significant warming trend (p<0.05)",
    "Area with cooling trend",
    "Area with increasing rainfall trend",
    "Area with decreasing rainfall trend",
    "Area with significant rainfall trend (p<0.05)",
    "Area with high/very-high temperature risk",
    "Area with high/very-high rainfall risk"
  ),
  Percent_area = c(
    area_share(temp[["temp_trend_per_decade"]] > 0),
    area_share(temp[["temp_trend_per_decade"]] > 0 & temp[["temp_trend_pvalue"]] < 0.05),
    area_share(temp[["temp_trend_per_decade"]] < 0),
    area_share(rain[["rain_trend_per_decade"]] > 0),
    area_share(rain[["rain_trend_per_decade"]] < 0),
    area_share(rain[["rain_trend_pvalue"]] < 0.05),
    area_share(temp[["temp_risk_class"]] >= 4),
    area_share(rain[["rain_risk_class"]] >= 4)
  )
)
write.csv(spatial_summary, file.path(out, "tables/spatial_area_summary.csv"), row.names = FALSE)

states <- read.csv(file.path(out, "tables/state_spatial_climate_metrics.csv"))
districts <- read.csv(file.path(out, "tables/district_spatial_climate_metrics.csv"))

state_rankings <- bind_rows(
  states |> arrange(desc(temp_trend_per_decade)) |> slice_head(n = 10) |> mutate(Ranking = "Highest temperature trend"),
  states |> arrange(temp_trend_per_decade) |> slice_head(n = 10) |> mutate(Ranking = "Lowest temperature trend"),
  states |> arrange(desc(abs(rain_trend_per_decade))) |> slice_head(n = 10) |> mutate(Ranking = "Largest absolute rainfall trend"),
  states |> arrange(desc(combined_risk_index)) |> slice_head(n = 10) |> mutate(Ranking = "Highest combined risk"),
  states |> arrange(combined_risk_index) |> slice_head(n = 10) |> mutate(Ranking = "Lowest combined risk")
)
write.csv(state_rankings, file.path(out, "tables/state_rankings.csv"), row.names = FALSE)

district_rankings <- bind_rows(
  districts |> arrange(desc(temp_trend_per_decade)) |> slice_head(n = 20) |> mutate(Ranking = "Highest temperature trend"),
  districts |> arrange(temp_trend_per_decade) |> slice_head(n = 20) |> mutate(Ranking = "Lowest temperature trend"),
  districts |> arrange(desc(abs(rain_trend_per_decade))) |> slice_head(n = 20) |> mutate(Ranking = "Largest absolute rainfall trend"),
  districts |> arrange(desc(combined_risk_index)) |> slice_head(n = 20) |> mutate(Ranking = "Highest combined risk"),
  districts |> arrange(combined_risk_index) |> slice_head(n = 20) |> mutate(Ranking = "Lowest combined risk")
)
write.csv(district_rankings, file.path(out, "tables/district_rankings.csv"), row.names = FALSE)

print(spatial_summary)
print(states |> arrange(desc(combined_risk_index)) |> select(State, combined_risk_index, combined_risk_class) |> head(10))
print(states |> arrange(combined_risk_index) |> select(State, combined_risk_index, combined_risk_class) |> head(10))
