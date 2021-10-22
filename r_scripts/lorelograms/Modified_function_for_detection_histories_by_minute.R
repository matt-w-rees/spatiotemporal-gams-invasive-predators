# Get occupancy detection histories by hours ------------------------------
library(camtrapR)

# modify function detectionHistory in package to get detection histories of 1-hour interval: modify lines 479, 483 (x2) and 549.

# Step 1: define helper functions
# for detectionHistory functions

checkCamOpColumnNames <- function(cameraOperationMatrix){
  camopTest <- try(as.Date(colnames(cameraOperationMatrix)), silent = TRUE)
  if(class(camopTest) == "try-error") stop(paste('could not interpret column names in camOp as Dates. Desired format is YYYY-MM-DD, e.g. "2016-12-31". First column name in your camera operation matrix is "', colnames(cameraOperationMatrix)[1], '"', sep = '' ), call. = FALSE)
}


createDateRangeTable <- function(cam.op,
                                 subset_species_tmp,
                                 buffer_tmp,
                                 stationCol_tmp,
                                 day1_tmp,
                                 occasionStartTime_tmp,
                                 maxNumberDays_tmp,
                                 timeZone_tmp)
{
  
  cam.tmp.min <- apply(cam.op, MARGIN = 1, function(X){min(which(!is.na(X)))})    # 1st day of each station
  cam.tmp.max <- apply(cam.op, MARGIN = 1, function(X){max(which(!is.na(X)))})    # last day of each station
  
  rec.tmp.min  <- aggregate(as.Date(subset_species_tmp$DateTime2, tz = timeZone_tmp),
                            list(subset_species_tmp[,stationCol_tmp]),
                            FUN = min)
  rec.tmp.max  <- aggregate(as.Date(subset_species_tmp$DateTime2, tz = timeZone_tmp),
                            list(subset_species_tmp[,stationCol_tmp]),
                            FUN = max)
  
  
  
  date_ranges <- data.frame(rec.min = rec.tmp.min[match(rownames(cam.op), rec.tmp.min[,1]), 2],       # first record
                            rec.max = rec.tmp.max[match(rownames(cam.op), rec.tmp.max[,1]), 2],       # last record
                            cam.min = as.POSIXct(colnames(cam.op)[cam.tmp.min], tz = timeZone_tmp),   # station setup date
                            cam.max = as.POSIXct(colnames(cam.op)[cam.tmp.max], tz = timeZone_tmp)    # station retrieval date
  )
  
  rownames(date_ranges) <- rownames(cam.op)
  
  # check if images were taken between setup and retrieval dates (Error if images outside station date range)
  if(any(date_ranges$rec.min < as.Date(date_ranges$cam.min, tz = timeZone_tmp), na.rm = TRUE)) stop(paste("record date outside camera operation date range: ",
                                                                                                          paste(rownames(date_ranges)[which(date_ranges$rec.min < as.Date(date_ranges$cam.min, tz = timeZone_tmp))], collapse = ", " )), call. = FALSE)
  if(any(date_ranges$rec.max > as.Date(date_ranges$cam.max, tz = timeZone_tmp), na.rm = TRUE)) stop(paste("record date outside camera operation date range: ",
                                                                                                          paste(rownames(date_ranges)[which(date_ranges$rec.max > as.Date(date_ranges$cam.max, tz = timeZone_tmp))], collapse = ", " )), call. = FALSE)
  
  # define when first occasion begins (to afterwards remove prior records in function    cleanSubsetSpecies)
  if(!hasArg(buffer_tmp)) buffer_tmp <- 0
  
  #if(day1_tmp == "station") {
  date_ranges$start_first_occasion <- date_ranges$cam.min + buffer_tmp * 3600 + occasionStartTime_tmp * 60     #each stations setup  + buffer + starttime MODIFIED FROM 86400 TO 3600 AND 3600 TO 60 
  #  } else {
  #    if(day1_tmp == "survey") {
  date_ranges$start_first_occasion_survey <- min(date_ranges$cam.min) + buffer_tmp * 3600 + occasionStartTime_tmp * 60    # first station's setup  + buffer + starttime
  #      } else {
  
  if(day1_tmp %in% c("survey", "station") == FALSE) {
    if(as.Date(day1_tmp, tz = timeZone_tmp) < min(as.Date(date_ranges$cam.min,  tz = timeZone_tmp))) stop(paste("day1 (", day1_tmp, ") is before the first station's setup date (",  min(as.Date(date_ranges$cam.min,  tz = timeZone_tmp)), ")", sep = ""))
    if(as.Date(day1_tmp, tz = timeZone_tmp) > max(as.Date(date_ranges$cam.max,  tz = timeZone_tmp))) stop(paste("day1 (", day1_tmp, ") is after the last station's retrieval date (",  max(as.Date(date_ranges$cam.max,  tz = timeZone_tmp)), ")", sep = ""))
    date_ranges$start_first_occasion <- as.POSIXlt(day1_tmp, tz = timeZone_tmp) + occasionStartTime_tmp * 60
  }
  
  
  
  # define when last occasion ends
  date_ranges$end_of_retrieval_day <- as.POSIXct(paste(date_ranges$cam.max, "23:59:59"), tz = timeZone_tmp, format = "%Y-%m-%d %H:%M:%S")    # end of retrieval day
  
  # if maxNumberDays is defined, find which is earlier: start + maxNumberDays or station retrieval?
  if(hasArg(maxNumberDays_tmp)) {
    
    if(day1_tmp %in% c("survey", "station") == FALSE){
      # count maximum number days from the beginning of each station's 1st occasion
      date_ranges$start_first_occasion_plus_maxNumberDays <- date_ranges$start_first_occasion_survey + (maxNumberDays_tmp * 3600) - 1   # -1 second ensures that last occasion does not spill into next day if occasionStartTime = 0
    } else {
      # count maximum number days from the beginning of survey's 1st occasion
      date_ranges$start_first_occasion_plus_maxNumberDays <- date_ranges$start_first_occasion + (maxNumberDays_tmp * 3600) - 1   # -1 second ensures that last occasion does not spill into next day if occasionStartTime = 0
    }
    
    for(xy in 1:nrow(date_ranges)){
      date_ranges$end_last_occasion[xy] <- min(date_ranges$end_of_retrieval_day[xy], date_ranges$start_first_occasion_plus_maxNumberDays[xy])   # use smaller value
    }
    attributes(date_ranges$end_last_occasion) <- attributes(date_ranges$start_first_occasion)   # assign the attributes: POSIX + time zone (to convert from numeric value back to date/time)
  } else {
    date_ranges$end_last_occasion <- date_ranges$end_of_retrieval_day
  }
  
  
  return(date_ranges)
  
}





adjustCameraOperationMatrix <- function(cam.op,
                                        date_ranges2,
                                        timeZone_tmp,
                                        day1_2
){
  
  
  
  if(any(date_ranges2$start_first_occasion > date_ranges2$end_last_occasion)){
    remove.these.stations <- which(date_ranges2$start_first_occasion > date_ranges2$end_last_occasion)
    if(length(remove.these.stations) == nrow(date_ranges2)) stop("In all stations, the occasions begin after retrieval. Choose a smaller buffer argument.")
    cam.op [remove.these.stations, ] <- NA
  }
  
  
  ##################################
  # set values before beginning of first occasion NA in camera operation matrix (so effort is 0 before): only relevant if buffer was used
  first_col_to_keep2  <- match(as.character(as.Date(date_ranges2$start_first_occasion, tz = timeZone_tmp)), colnames(cam.op))
  
  
  for(xxx in 1:length(first_col_to_keep2)){
    if(first_col_to_keep2[xxx] > 1){         # if it does not begin on 1st day of camera operation matrix
      cam.op[xxx, seq(1, first_col_to_keep2[xxx]-1)] <- NA
    }
  }
  
  # set values after end of last occasion NA in camera operation matrix (so effort is 0 afterwards)
  last_col_to_keep2  <- match(as.character(as.Date(date_ranges2$end_last_occasion, tz = timeZone_tmp)), colnames(cam.op))
  
  
  for(yyy in 1:length(last_col_to_keep2)){
    if(last_col_to_keep2[yyy]+1 < ncol(cam.op)){   # if it does not end on last day of camera operation matrix
      cam.op[yyy, seq(last_col_to_keep2[yyy]+1, ncol(cam.op))] <- NA
    }
  }
  
  
  ####################################
  # trim camera operation matrix (taking into account buffer, occasionStartTime, maxNumberDays)
  # based on data frame   date_ranges   computed by    createDateRangeTable
  
  if(day1_2 == "station") {    # 1st day of each station OR some specified date
    
    cam.tmp.min <- apply(cam.op, MARGIN = 1, function(X){min(which(!is.na(X)))})    # first occasion of each station
    cam.tmp.max <- apply(cam.op, MARGIN = 1, function(X){max(which(!is.na(X)))})    # last occasion of each station
    
    diff.days.tmp <- cam.tmp.max - cam.tmp.min
    
    cam.op2 <- matrix(NA,
                      nrow = nrow(cam.op),
                      ncol = max(diff.days.tmp)+1)
    
    # make all stations begin in 1st column
    for(l in 1:nrow(cam.op)){
      if(is.finite(diff.days.tmp[l])){
        cam.op2[l, 1:(diff.days.tmp[l]+1)] <- as.vector(cam.op[l,cam.tmp.min[l]:cam.tmp.max[l]])
      }
    }
    
    if(day1_2 == "station") {
      colnames(cam.op2) <- paste("day", 1:ncol(cam.op2), sep = "")
    }
    
    rownames(cam.op2) <- rownames(cam.op)
    
    cam.op <- cam.op2
    
    
  } else {
    
    
    # remove all columns of cam.op that were before beginning of 1st occasion
    first_col_to_keep <- match(as.character(min(as.Date(date_ranges2$start_first_occasion, tz = timeZone_tmp))), colnames(cam.op))
    if(!is.na(first_col_to_keep)){
      if(first_col_to_keep != 1){
        cam.op <- cam.op[,-seq(1, (first_col_to_keep-1))]
      }
    }
    
    # remove all columns of cam.op that were after end of last occasion / after retrieval of last camera
    last_col_to_keep  <- match(as.character(max(as.Date(date_ranges2$end_last_occasion, tz = timeZone_tmp))), colnames(cam.op))
    if(!is.na(last_col_to_keep)){
      if(last_col_to_keep != ncol(cam.op)){
        cam.op <- cam.op[,-seq((last_col_to_keep + 1), ncol(cam.op))]
      }
    }
    
  }
  return(cam.op)
}



cleanSubsetSpecies <- function(subset_species2 ,
                               stationCol2,
                               date_ranges2
){
  
  nrow_subset_species2 <- nrow(subset_species2)
  
  
  # remove records that were taken before beginning of first occasion (because of buffer, occasionStartTime, day1)
  corrected_start_time_by_record <- date_ranges2$start_first_occasion[match(subset_species2[,stationCol2], rownames(date_ranges2))]
  
  remove.these <- which(subset_species2$DateTime2 < corrected_start_time_by_record)
  if(length(remove.these) >= 1){
    subset_species2 <- subset_species2[-remove.these,]
    warning(paste(length(remove.these), "records out of", nrow_subset_species2, "were removed because they were taken within the buffer period, before day1 (if a date was specified), or before occasionStartTime on the 1st day"), call. = FALSE)
    if(nrow(subset_species2) == 0) stop("No more records left. The detection history would be empty.")
    rm(corrected_start_time_by_record, remove.these)
  }
  
  # remove records that were taken after end of last occasion (because of maxNumberDays)
  corrected_end_time_by_record <- date_ranges2$end_last_occasion[match(subset_species2[,stationCol2], rownames(date_ranges2))]
  
  remove.these2 <- which(subset_species2$DateTime2 > corrected_end_time_by_record)
  if(length(remove.these2) >= 1){
    subset_species2 <- subset_species2[-remove.these2,]
    warning(paste(length(remove.these2), "records out of", nrow_subset_species2, "were removed because they were taken after the end of the last occasion"), call. = FALSE)
    if(nrow(subset_species2) == 0) stop("No more records left. The detection history would be empty.")
    rm(corrected_end_time_by_record, remove.these2)
  }
  
  return(subset_species2)
}



calculateTrappingEffort <- function(cam.op,
                                    occasionLength2,
                                    scaleEffort2,
                                    includeEffort2,
                                    minActiveDaysPerOccasion2){
  
  ######################
  # calculate trapping effort by station and occasion
  
  if(occasionLength2 == 1){
    effort <- cam.op          # if occasionLength2 = 1 day, it is identical
  } else {
    effort <- matrix(NA, nrow = nrow(cam.op), ncol = ceiling(ncol(cam.op) / occasionLength2 ))
    
    index <- 1
    for(m in 1:ncol(effort)){    # for every occasion in the effort matrix
      # index for columns to aggregate
      if(index + occasionLength2 <= ncol(cam.op)){
        index.tmp <- index : (index + occasionLength2 - 1)
      } else {
        index.tmp <- index : ncol(cam.op)
      }
      
      # calculate effort as sum of active days per occasion
      effort[, m] <- apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = sum, na.rm = TRUE)
      # if full occasion NA in cam.op, make effort NA
      effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {sum(is.na(X))}) == length(index.tmp), NA, effort[,m])
      # if full occasion = 0 in cam.op, make effort NA
      effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {all(X == 0)}),   NA, effort[,m])
      # if full occasion is not 1 (i.e. all 0 or NA), set effort NA
      effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {all(X != 1)}),   NA, effort[,m])
      
      # set cells in effort matrix NA (according to input arguments)
      # this is later used to adjust the detection/non-detection matrix
      if(includeEffort2 == FALSE){
        if(hasArg(minActiveDaysPerOccasion2)){   # includeEffort = FALSE and minActiveDays is defined
          # if occasion has less active days than minActiveDays, set NA
          effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {sum(X, na.rm = TRUE) < minActiveDaysPerOccasion2}), NA, effort[,m])
        } else {                                 # includeEffort = FALSE and minActiveDays not defined
          # if any day NA in cam.op, make  occasion effort NA
          effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {any(is.na(X))}), NA, effort[,m])
          # if any day = 0 in cam.op, make occasion effort NA
          effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {any(X == 0)}),   NA, effort[,m])
          
          if(length(index.tmp) < occasionLength2){  # if occasion is shorter than occasionLength (i.e. last occasion), set NA.
            effort[, m] <- NA
          }
          
        }
      } else {
        if(hasArg(minActiveDaysPerOccasion2)){   # includeEffort = TRUE and minActiveDays is defined
          # if occasion has less actice days than minActiveDays, set NA
          effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {sum(X, na.rm = TRUE) < minActiveDaysPerOccasion2}), NA, effort[,m])
        } else {                                 # includeEffort = TRUE and minActiveDays not defined
          # if all days of occasion NA in cam.op, make  occasion effort NA
          #effort[, m] <- ifelse(apply(as.matrix(cam.op[,index.tmp]), MARGIN = 1, FUN = function(X) {all(is.na(X))}), NA, effort[,m])
          # if all days of occasion = 0 in cam.op, make occasion effort NA
          
        }
      }
      index <- index + occasionLength2
    }
    rm(index, index.tmp)
  }
  
  # scale effort, if required
  if(isTRUE(scaleEffort2)){
    if(occasionLength2 == 1) stop("cannot scale effort if occasionLength is 1")
    scale.eff.tmp <- scale(as.vector(effort))                       # scale effort (as a vector, not matrix)
    scale.eff.tmp.attr <- data.frame(effort.scaled.center = NA,     # prepare empty data frame
                                     effort.scaled.scale = NA)
    scale.eff.tmp.attr$effort.scaled.center[1] <- attr(scale.eff.tmp, which = "scaled:center")   # write scaling parameters to data frame
    scale.eff.tmp.attr$effort.scaled.scale [1] <- attr(scale.eff.tmp, which = "scaled:scale")
    effort <- matrix(scale.eff.tmp, nrow = nrow(effort), ncol = ncol(effort))                    # convert effort vector to matrix again
  }
  
  rownames(effort) <- rownames(cam.op)
  
  # return objects
  if(isTRUE(scaleEffort2)){
    return(list(effort, scale.eff.tmp.attr))
  } else {
    return(list(effort))
  }
}

# for all functions in which user specifies column names: error if spaces in column names
checkForSpacesInColumnNames <- function(...){
  
  z <- list(...)
  
  # if all arguments are of length 1, do
  if(all(sapply(z, FUN = length) == 1)){
    if(any(grepl(pattern = " ", x = unlist(z), fixed = TRUE))) stop("column names may not contain spaces: \n ",
                                                                    paste(names(z)[which(grepl(pattern = " ", x = unlist(z), fixed = TRUE))], "=",
                                                                          z[which(grepl(pattern = " ", x = unlist(z), fixed = TRUE))], collapse = "\n "),
                                                                    call. = FALSE)
  }
  # if the argument is of length >1, do
  if(any(sapply(z, FUN = length) > 1)){
    if(length(z) != 1) stop("this is a bug in 'checkForSpacesInColumnNames'. I'm sorry. Please report it.")
    if(any(grepl(pattern = " ", x = unlist(z[[1]]), fixed = TRUE))) stop("column names in '", names(z) ,"' may not contain spaces: \n ",
                                                                         paste(names(unlist(z))[which(grepl(pattern = " ", x = unlist(z), fixed = TRUE))], "=",
                                                                               z[[1]][which(grepl(pattern = " ", x = unlist(z), fixed = TRUE))], collapse = "\n "),
                                                                         call. = FALSE)
  }
}



# Function for detection histories by minute --------------------------------


# store modified function

myfunc_detectionHistory_minute <- function(recordTable,
                                         species,
                                         camOp,
                                         stationCol = "Station",
                                         speciesCol = "Species",
                                         recordDateTimeCol = "DateTimeOriginal",
                                         recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                         occasionLength,
                                         minActiveDaysPerOccasion,
                                         maxNumberDays,
                                         day1,
                                         buffer,
                                         includeEffort = TRUE,
                                         scaleEffort = FALSE,
                                         occasionStartTime = 0,
                                         datesAsOccasionNames = FALSE,
                                         timeZone,
                                         writecsv = FALSE,
                                         outDir)
{
  wd0 <- getwd()
  on.exit(setwd(wd0))
  #################
  # check input
  
  # check column names
  checkForSpacesInColumnNames(stationCol = stationCol, speciesCol = speciesCol, recordDateTimeCol = recordDateTimeCol)
  if(!is.data.frame(recordTable)) stop("recordTable must be a data frame", call. = FALSE)
  if(!stationCol %in% colnames(recordTable))  stop(paste('stationCol = "', stationCol, '" is not a column name in recordTable', sep = ''), call. = FALSE)
  if(!speciesCol %in% colnames(recordTable))  stop(paste('speciesCol = "', speciesCol, '" is not a column name in recordTable', sep = ''), call. = FALSE)
  if(!recordDateTimeCol %in% colnames(recordTable))  stop(paste('recordDateTimeCol = "', recordDateTimeCol,  '" is not a column name in recordTable', sep = ''), call. = FALSE)
  
  
  stopifnot(hasArg(species))
  stopifnot(is.character(species))
  stopifnot(length(species) == 1)
  
  stopifnot(hasArg(occasionLength))
  
  stopifnot(hasArg(recordTable))
  stopifnot(class(recordTable) == "data.frame")
  stopifnot(hasArg(camOp))
  
  stopifnot(hasArg(stationCol))
  stopifnot(length(stationCol) == 1)
  recordTable[,stationCol] <- as.character(recordTable[,stationCol])
  stopifnot(is.character(stationCol))
  
  stopifnot(hasArg(speciesCol))
  stopifnot(length(speciesCol) == 1)
  recordTable[,speciesCol] <- as.character(recordTable[,speciesCol])
  stopifnot(is.character(speciesCol))
  
  
  stopifnot(hasArg(recordDateTimeCol))
  stopifnot(length(recordDateTimeCol) == 1)
  recordTable[,recordDateTimeCol] <- as.character(recordTable[,recordDateTimeCol])   # make character to get rid of attributes. Will later assign time zone again
  stopifnot(is.character(recordDateTimeCol))
  
  
  if(hasArg(timeZone) == FALSE) {
    warning("timeZone is not specified. Assuming UTC", call. = FALSE)
    timeZone <- "UTC"
  }
  if(!is.element(timeZone , OlsonNames())){
    stop("timeZone must be an element of OlsonNames()")
  }
  stopifnot(is.logical(writecsv))
  
  occasionStartTime <- as.integer(round(occasionStartTime))
  if(occasionStartTime != 0 & !is.integer(occasionStartTime)) {stop ("occasionStartTime must be between 0 and 23")}
  if(occasionStartTime < 0 | occasionStartTime >= 24){stop ("occasionStartTime must be between 0 and 23")}
  
  occasionLength <- as.integer(round(occasionLength))
  if(occasionLength <= 0) stop("occasionLength must be a positive integer and not 0")
  if(occasionLength > ncol(camOp)/2) stop("occasionLength may not be greater than half the total number of days in camOp")
  stopifnot(is.numeric(occasionLength))
  
  if(hasArg(maxNumberDays)){
    maxNumberDays <- as.integer(maxNumberDays)
    if(maxNumberDays > ncol(camOp)) stop("maxNumberDays is larger than the number of columns of camOp")
    if(maxNumberDays < occasionLength) stop("maxNumberDays must be larger than or equal to occasionLength")
  }
  
  if(hasArg(buffer)) {
    stopifnot(is.numeric(buffer))
    buffer <- round(buffer)
    stopifnot(buffer >= 1)
  }
  
  stopifnot(c(speciesCol, recordDateTimeCol, stationCol) %in% colnames(recordTable))
  
  if(species %in% recordTable[,speciesCol] == FALSE) stop("species is not in speciesCol of recordTable")
  
  if(writecsv == TRUE){
    if(!file.exists(outDir)){stop("outDir does not exist")}
  }
  
  if(includeEffort == TRUE){
    if(!hasArg(scaleEffort)) stop("scaleEffort must be defined if includeEffort is TRUE")
    if(class(scaleEffort) != "logical") stop("scaleEffort must be logical (TRUE or FALSE)")
  } else {scaleEffort <- FALSE}
  
  if(hasArg(minActiveDaysPerOccasion)){
    stopifnot(is.numeric(minActiveDaysPerOccasion))
    stopifnot(minActiveDaysPerOccasion <= occasionLength)
  }
  
  #############
  # bring date, time, station ids into shape
  
  subset_species           <- subset(recordTable, recordTable[,speciesCol] == species)
  subset_species$DateTime2 <- as.POSIXlt(subset_species[,recordDateTimeCol], tz = timeZone, format = recordDateTimeFormat)
  
  # check consistency of argument day1
  stopifnot(class(day1) == "character")
  if(day1 == "survey") {day1switch <- 1} else {
    if(day1 == "station") {day1switch <- 2} else {
      try(date.test <- as.Date(day1), silent = TRUE)
      if(class(date.test) != "Date") stop('could not interpret argument day1: can only be "station", "survey" or a specific date (e.g. "2015-12-31")')
      if(hasArg(buffer)) stop("if buffer is defined, day1 can only be 'survey' or 'station'")
      suppressWarnings(rm(date.test))
      day1switch <- 3
    }
  }
  
  if("POSIXlt" %in% class(subset_species$DateTime2) == FALSE) stop("couldn't interpret recordDateTimeCol of recordTable using specified recordDateTimeFormat")
  if(any(is.na(subset_species$DateTime2))) stop("at least 1 entry in recordDateTimeCol of recordTable could not be interpreted using recordDateTimeFormat")
  
  checkCamOpColumnNames (cameraOperationMatrix = camOp)
  
  #MODIFIED FUNCTION: REPEATED EACH COLUMN OF camOp 24*60=1440 TIMES (1 PER MINUTE)
  camOp_minute <- matrix(, nrow=nrow(camOp), ncol=0)
  rownames(camOp_minute) <- rownames(camOp)
  names_col <- NULL
  for (i in 1:ncol(camOp)) {
    nm <- rep(colnames(camOp[i]), times=1440)
    names_col <- c(names_col, nm)
  }
  for(i in 1:ncol(camOp)) {
    Op <- matrix(rep(camOp[,i], time = 1440, byrow=FALSE), ncol=1440, nrow=nrow(camOp))
    camOp_minute <- cbind(camOp_minute, Op)
  }
  colnames(camOp_minute) <- names_col
  camOp <- camOp_minute
  #END MODIFIED CODE
  
  
  cam.op.worked0 <- as.matrix(camOp)
  
  if(all(as.character(unique(subset_species[,stationCol])) %in% rownames(cam.op.worked0)) == FALSE){
    (stop("Not all values of stationCol in recordTable are matched by rownames of camOp"))
  }
  
  ################################################
  # compute date range of stations and records
  arg.list0 <- list(cam.op = cam.op.worked0, subset_species_tmp = subset_species, stationCol_tmp = stationCol, day1_tmp = day1, occasionStartTime_tmp = occasionStartTime, timeZone_tmp = timeZone)
  
  if(hasArg(maxNumberDays))  arg.list0 <- c(arg.list0,   maxNumberDays_tmp = maxNumberDays)
  if(hasArg(buffer))   arg.list0 <- c(arg.list0, buffer_tmp =  buffer)
  
  date_ranges <- do.call(createDateRangeTable, arg.list0)
  
  rm(arg.list0)
  
  #######################
  # adjust camera operation matrix
  
  cam.op.worked <- adjustCameraOperationMatrix(cam.op = cam.op.worked0, date_ranges2 = date_ranges, timeZone_tmp = timeZone, day1_2 = day1)
  
  # append occasionStartTime (if != 0) to column names for output table
  if(occasionStartTime != 0){
    colnames(cam.op.worked) <- paste(colnames(cam.op.worked), "+", occasionStartTime, "h", sep = "")
  }
  
  ######################
  # calculate trapping effort by station and occasion
  arg.list0 <- list(cam.op = cam.op.worked, occasionLength2 = occasionLength, scaleEffort2 = scaleEffort, includeEffort2 = includeEffort)
  
  if(hasArg(minActiveDaysPerOccasion))  arg.list0 <- c(arg.list0, minActiveDaysPerOccasion2 = minActiveDaysPerOccasion)
  
  effort.tmp <- do.call(calculateTrappingEffort, arg.list0)
  
  rm(arg.list0)
  
  effort <- effort.tmp[[1]]
  if(isTRUE(scaleEffort))  scale.eff.tmp.attr <- effort.tmp[[2]]
  
  ###################
  # remove records that fall into buffer period or were taken after maxNumberDays
  
  subset_species <- cleanSubsetSpecies(subset_species2 = subset_species, stationCol2 = stationCol, date_ranges2 = date_ranges)
  
  ############
  #  define the 1st day of the effective survey period.
  if(day1 %in% c("survey")){
    time2 <- date_ranges$start_first_occasion_survey[match(subset_species[,stationCol], rownames(date_ranges))]
  } else {
    time2 <- date_ranges$start_first_occasion[match(subset_species[,stationCol], rownames(date_ranges))]
  }
  
  # calculate the occasion each record belongs into from the time difference between records and beginning of the first occasion
  subset_species$occasion <- as.numeric(ceiling((difftime(time1  = subset_species$DateTime2,
                                                          time2 =  time2,
                                                          units = "secs",
                                                          tz = timeZone)
                                                 / (occasionLength * 60)))) #CHANGED FROM 86400 TO 60: FROM DAY-INT TO MINUTE-INT
  
   if(max(subset_species$occasion) > ncol(effort)) {stop("encountered a bug. I'm Sorry. Please report it.")}
  
  ############
  # make detection history
  
  if(occasionLength == 1){                    # occasion length = 1
    record.hist <- cam.op.worked
    record.hist <- ifelse(record.hist == 0, NA, record.hist)    # if cameras not operational, set NA
    record.hist <- ifelse(record.hist >= 1, 0,  record.hist)    # if cameras operational, set to 0
    
    occasions.by.station <- tapply(X = subset_species$occasion, INDEX = subset_species[,stationCol], FUN = unique, simplify = FALSE)
    
    # fill detection matrix with 1 in appropriate cells
    for(xyz in which(sapply(occasions.by.station, FUN = function(x){!is.null(x)}))){
      if(any(occasions.by.station[[xyz]] < 0)) stop("this is a bug in the function (line 188). Please report it.", call. = FALSE)
      if(any(occasions.by.station[[xyz]] > ncol(record.hist))) stop("this is a bug in the function (line 189). Please report it.", call. = FALSE)
      record.hist[match(names(occasions.by.station)[xyz], rownames(record.hist)), occasions.by.station[[xyz]]] <- 1
    }
    record.hist[is.na(cam.op.worked)] <- NA   # remove the records that were taken when cams were NA (redundant with above:   # remove records taken after day1 + maxNumberDays)
    
    rm(occasions.by.station, xyz)
    
  } else {                                    # occasion length > 1
    record.hist <- effort                     # start with effort (including NAs)
    record.hist <- ifelse(!is.na(record.hist), 0, record.hist)     # set all cells 0 (unless they are NA)
    rownames(record.hist) <- rownames(cam.op.worked)
    
    occasions.by.station <- tapply(X = subset_species$occasion, INDEX = subset_species[,stationCol], FUN = unique, simplify = FALSE)
    
    # fill detection matrix with "1" in appropriate cells
    for(xyz in which(sapply(occasions.by.station, FUN = function(x){!is.null(x)}))){
      record.hist[match(names(occasions.by.station)[xyz], rownames(record.hist)), occasions.by.station[[xyz]]] <- 1
    }
    record.hist[is.na(effort)] <- NA     # just to make sure NA stays NA
  }
  
  # assign row names to output
  row.names(record.hist) <- row.names(effort) <- rownames(cam.op.worked)
  
  
  # assign column names
  if(isTRUE(datesAsOccasionNames)){
    
    seq.tmp <- seq(from = 1, by = occasionLength, length.out = ncol(record.hist))
    if(occasionStartTime != 0){
      colnames.tmp <- paste(colnames(cam.op.worked)[seq.tmp],
                            colnames(cam.op.worked)[seq.tmp + occasionLength], sep = "_")
    } else {
      colnames.tmp <- paste(colnames(cam.op.worked)[seq.tmp],
                            colnames(cam.op.worked)[seq.tmp + occasionLength - 1], sep = "_")
    }
    # make sure the last occasion ends on the last day
    
    if(day1 == "station"){
      colnames.tmp[length(colnames.tmp)] <- paste(colnames(cam.op.worked)[max(seq.tmp)],
                                                  colnames(cam.op.worked)[ncol(cam.op.worked)],
                                                  sep = "_")
    } else {   # if day1 = "survey" or date
      colnames.tmp[length(colnames.tmp)] <- paste(colnames(cam.op.worked)[max(seq.tmp)],
                                                  max(as.Date(date_ranges$end_last_occasion)),
                                                  sep = "_")
    }
    
    colnames(record.hist) <- colnames(effort) <- colnames.tmp
    
  }  else {
    colnames(record.hist) <- colnames(effort) <-  paste("o", seq(1,ncol(record.hist), by = 1), sep = "")
  }
  
  
  ################################################
  # save output as table
  if(day1switch == 1) day1string <- "_first_day_from_survey"
  if(day1switch == 2) day1string <- "_first_day_by_station"
  if(day1switch == 3) day1string <- paste("_first_day", day1, sep = "_")
  
  effortstring <- ifelse(isTRUE(includeEffort), "with_effort__", "no_effort__")
  maxNumberDaysstring <- ifelse(hasArg(maxNumberDays), paste("max",maxNumberDays,"days_", sep = ""), "")
  if(isTRUE(includeEffort)){
    scaleEffortstring <- ifelse(isTRUE(scaleEffort), "scaled_", "not_scaled_")
  } else {
    scaleEffortstring <- ""
  }
  
  # create names for the csv files
  outtable.name <- paste(species, "__record_history__", effortstring,
                         occasionLength, "_days_per_occasion_",
                         maxNumberDaysstring,
                         "_occasionStart", occasionStartTime,"h_",
                         day1string, "__",
                         Sys.Date(),
                         ".csv", sep = "")
  
  outtable.name.effort <- paste(species, "__effort__",
                                scaleEffortstring,
                                occasionLength, "_days_per_occasion_",
                                maxNumberDaysstring,
                                "_occasionStart", occasionStartTime,"h_",
                                day1string, "__",
                                Sys.Date(),
                                ".csv", sep = "")
  
  outtable.name.effort.scale <- paste(species, "__effort_scaling_parameters__",
                                      occasionLength, "_days_per_occasion_",
                                      maxNumberDaysstring,
                                      "_occasionStart", occasionStartTime,"h_",
                                      day1string, "__",
                                      Sys.Date(),
                                      ".csv", sep = "")
  
  if(isTRUE(writecsv)){
    setwd(outDir)
    write.csv(record.hist, file = outtable.name)
    if(isTRUE(includeEffort)){
      write.csv(effort, file = outtable.name.effort)
      if(hasArg(scaleEffort)){
        if(scaleEffort == TRUE)  write.csv(scale.eff.tmp.attr, file = outtable.name.effort.scale)
      }
    }
  }
  
  if(isTRUE(includeEffort)){
    if(scaleEffort == TRUE){
      return(list(detection_history = record.hist,
                  effort = effort,
                  effort_scaling_parameters = scale.eff.tmp.attr))
    } else {
      return(list(detection_history = record.hist,
                  effort = effort))
    }
  } else {
    return(list(detection_history = record.hist))
  }
}




# Apply to real data ------------------------------------------------------

# read data: 

# read recordTable file as created using function recordTable in the camtrapR package (see https://cran.r-project.org/web/packages/camtrapR/vignettes/DataExtraction.html#camera-operation)
recordData <- read.csv("") 

# read camera operation matrix as created using function cameraOperation in the camtrapR package (see https://cran.r-project.org/web/packages/camtrapR/vignettes/DataExtraction.html#camera-operation)
camop_problem <- read.csv("", header=TRUE, check.names = FALSE) 
row.names(camop_problem) <- camop_problem[,1]
camop_problem <- camop_problem[,-1] #to remove first column with station name 

# run modified function
DetHist_Minute <- myfunc_detectionHistory_Minute(recordTable = recordData,
                                             species = "GreyFox", #change species name
                                             camOp = camop_problem,
                                             stationCol =  "Station",
                                             speciesCol = "Species",
                                             recordDateTimeCol = "DateTimeOriginal", 
                                             recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                             occasionLength = 1, 
                                             includeEffort = FALSE,
                                             scaleEffort = FALSE, 
                                             day1 = "survey", 
                                             timeZone="US/Central") #adjust time zone
write.csv(x = DetHist_Minute$detection_history, file = "output/DetHist_Minute.csv" )
