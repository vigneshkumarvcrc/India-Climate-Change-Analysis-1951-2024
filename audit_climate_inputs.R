library(terra)

base <- "/Volumes/K1/climate_data"
out <- "/Users/khkr/Documents/New project/CCPOP_work/audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rasters <- c(
  tmin = file.path(base, "tmin/tmin_yearly_stack.tif"),
  tmax = file.path(base, "tmax/tmax_yearly_stack.tif"),
  tmean = file.path(base, "tmean/tmean_yearly_stack.tif"),
  rain = file.path(base, "precip/rain_yearly_stack.tif")
)

raster_audit <- lapply(names(rasters), function(nm) {
  r <- rast(rasters[[nm]])
  g <- global(r, c("mean", "min", "max", "sd"), na.rm = TRUE)
  data.frame(
    variable = nm,
    layer = seq_len(nlyr(r)),
    layer_name = names(r),
    mean = g$mean,
    min = g$min,
    max = g$max,
    sd = g$sd,
    nrow = nrow(r),
    ncol = ncol(r),
    xmin = xmin(r),
    xmax = xmax(r),
    ymin = ymin(r),
    ymax = ymax(r),
    crs = crs(r, proj = TRUE)
  )
})
write.csv(do.call(rbind, raster_audit),
          file.path(out, "raster_layer_audit.csv"), row.names = FALSE)

vector_files <- c(
  india = file.path(base, "india_boundary.shp"),
  states = file.path(base, "Ind_st.shp"),
  districts = file.path(base, "Ind_Dist.shp")
)

vector_audit <- lapply(names(vector_files), function(nm) {
  v <- vect(vector_files[[nm]])
  data.frame(
    vector = nm,
    features = nrow(v),
    geomtype = geomtype(v),
    fields = paste(names(v), collapse = "|"),
    xmin = xmin(v),
    xmax = xmax(v),
    ymin = ymin(v),
    ymax = ymax(v),
    crs = crs(v, proj = TRUE)
  )
})
write.csv(do.call(rbind, vector_audit),
          file.path(out, "vector_audit.csv"), row.names = FALSE)

states <- vect(vector_files[["states"]])
districts <- vect(vector_files[["districts"]])
write.csv(as.data.frame(states), file.path(out, "state_attributes.csv"), row.names = FALSE)
write.csv(as.data.frame(districts), file.path(out, "district_attributes.csv"), row.names = FALSE)

cat("Audit complete\n")
