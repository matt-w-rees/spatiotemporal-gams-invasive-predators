
# load function to change clock time to radial time for overlap package
ClocktimeToTimeRad <- function(Clocktime,
                               timeformat = "%Y-%m-%d %H:%M:%S")
{
  DateTime2 <- strptime(as.character(Clocktime), format = timeformat, tz = "UTC")
  Time2     <- format(DateTime2, format = "%H:%M:%S", usetz = FALSE)
  Time.rad  <- (as.numeric(as.POSIXct(strptime(Time2, format = "%H:%M:%S", tz = "UTC"))) -
                  as.numeric(as.POSIXct(strptime("0", format = "%S", tz = "UTC")))) / 3600 * (pi/12)
  return(Time.rad)
}



# ADD SOLAR TIMES FOR DETECTIONS ------------------------------------------
# take a matrix of the coordinates
coords <- matrix(c(records$longitude, records$latitude), nrow = length(records$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 

# calculate the time of sunrise given the location and timing of each species detection 
sunriset_df <- as.data.frame(sunriset(coords, records$date_time, direction = "sunrise", POSIXct.out=TRUE))
# convert to radian time and add to the dataframe
sunrise_radian <- ClocktimeToTimeRad(sunriset_df$time) 
# calculate the time of sunset given the location and timing of each species detection 
sunsett_df <- as.data.frame(sunriset(coords, records$date_time, direction = "sunset", POSIXct.out=TRUE))
# convert to radian time and add to the dataframe
sunset_radian <- ClocktimeToTimeRad(sunsett_df$time) 
# save these as a separate matrix for the transtime function below 
anchors <- matrix(c(sunrise_radian, sunset_radian), nrow = length(sunset_radian))

# we also need to convert the detection time to radian time
det_radian <- ClocktimeToTimeRad(records$date_time)

# we can now calculate the average anchored time (time relative to both sunrise and sunset) for every detection
det_anchored <- transtime(det_radian, anchors)

# and we can convert radian time back to 24-hour time, where 06:00 represents sunrise and 18:00 sunset
records$hour_adj <- det_anchored * 3.819719 
head(records$det_anchored_24)

