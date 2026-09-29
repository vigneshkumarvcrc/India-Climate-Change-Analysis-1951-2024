library(terra)

base <- "/Volumes/K1/climate_data"
out <- "/Users/khkr/Documents/New project/CCPOP_final"

temp <- rast(file.path(out, "rasters/temp_risk_index.tif"))
rain <- rast(file.path(out, "rasters/rain_risk_index.tif"))
combined <- (temp + rain) / 2
names(combined) <- "combined_climate_hazard_index"
combined_class <- floor(clamp(combined, 0, 0.9999) * 5) + 1
names(combined_class) <- "combined_climate_hazard_class"

writeRaster(combined, file.path(out, "rasters/combined_climate_hazard_index.tif"),
            overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
writeRaster(combined_class, file.path(out, "rasters/combined_climate_hazard_class.tif"),
            overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW"))

india <- project(vect(file.path(base, "india_boundary.shp")), crs(combined))
states <- project(vect(file.path(base, "Ind_st.shp")), crs(combined))
risk_pal <- c("#1A9850", "#91CF60", "#FEE08B", "#FC8D59", "#D73027")

png(file.path(out, "maps/combined_climate_hazard_class.png"),
    width = 4200, height = 3300, res = 450, bg = "white")
par(mar = c(2.5, 2.5, 4.2, 5.2), family = "sans")
plot(combined_class, col = risk_pal, axes = FALSE,
     main = "Relative combined climate-hazard class",
     plg = list(title = "1 very low - 5 very high"))
plot(states, add = TRUE, border = adjustcolor("grey15", 0.55), lwd = 0.45)
plot(india, add = TRUE, border = "black", lwd = 1)
dev.off()

cat("Combined climate-hazard map created\n")
