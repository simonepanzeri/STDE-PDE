## HELPER FUNCTIONS FOR TIME CONVERSION ----------------------------------------
generateHour <- function() {
  
  h <- correctHour(as.character(sample(seq(from = 0, to = 23, by = 1), 1)))
  m <- correctHour(as.character(sample(seq(from = 0, to = 59, by = 1), 1)))
  
  return (paste0(h, ":", m))
}

correctHour <- function(x){
  
  if(nchar(x) == 1){
    x <- paste0("0", x)
  }
  
  return (x)
}

## DATASET ---------------------------------------------------------------------
readDataset <- function(district = "Rome", years = NULL){
  
  foldername <- paste0(getwd(), "/utils/data/", district, "/")

  if(!file.exists(paste0(foldername, district, "_data_pedestrians.RData")) ||
     !file.exists(paste0(foldername, district, "_data_vehicles.RData")) ||
     !file.exists(paste0(foldername, district, "_dataset.RData"))){
    
    variables <- c("ID", "Longitude", "Latitude", "AccidentType", "DateHour",
                   "Year", "Month", "Day", "DayType", "Hour")
    
    data_pedestrians <- data.frame(matrix(ncol = length(variables), nrow = 1))
    names(data_pedestrians) <- variables
    data_pedestrians$DateHour <- as.POSIXct(NA,"")
    
    data_vehicles <- data.frame(matrix(ncol = length(variables), nrow = 1))
    names(data_vehicles) <- variables
    data_vehicles$DateHour <- as.POSIXct(NA,"")
    
    dataset <- data.frame(matrix(ncol = length(variables), nrow = 1))
    names(dataset) <- variables
    dataset$DateHour <- as.POSIXct(NA,"")
    
    if(district == "Rome"){
      
      if(is.null(years)) years <- c(2006, 2022)
      if(length(years == 1)) years <- rep(years, 2)
      
      pb <- txtProgressBar(min = 1, max = (years[2] - years[1] + 1)*12, style = 3, char = "=")
        
      for(y in seq(from = years[1], to = years[2], by = 1)){
        
        for(m in 1:12){
          
          if(nchar(as.character(m)) == 1){
            m <- paste0("0", m)
          }
          
          filename <- paste0(foldername, y, "/road_accidents_", district, "_", y, "_", m, ".csv")
          
          if(file.exists(filename)){
            
            # reading the data
            data <- read.csv(filename, header = TRUE, sep = ";", quote = "")
            
            if(sum(names(data) == "Longitude")){
              
              names(data)[which(names(data) == "Longitude")] <- "Longitudine"
              
            }
            
            if(sum(names(data) == "Latitude")){
              
              names(data)[which(names(data) == "Latitude")] <- "Latitudine"
              
            }
              
            # selecting only the variables of interest
            data <- data %>%
              select("Protocollo", "DataOraIncidente", "Longitudine", "Latitudine", "NaturaIncidente")
            
            # reassigning variable names
            names(data) <- c("ID", "DateHour", "Longitude", "Latitude", "AccidentType")
            
            # keeping only one row per road accident
            data <- data %>% distinct(ID, .keep_all = TRUE)
            
            data$Longitude <- ifelse(data$Longitude == "", NA, as.numeric(sub(",", ".", data$Longitude)))
            data$Latitude <- ifelse(data$Latitude == "", NA, as.numeric(sub(",", ".", data$Latitude)))
            
            # discarding NAs (in longitudes and latitudes)
            data <- na.omit(data)
            
            if(nrow(data) > 0){
              
              # correcting coordinates
              for(i in 1:nrow(data)){
                
                n <- floor(log10(abs(data$Longitude[i]))) + 1
                data$Longitude[i] <- data$Longitude[i] / 10^(n-2)
                
                n <- floor(log10(abs(data$Latitude[i]))) + 1
                data$Latitude[i] <- data$Latitude[i] / 10^(n-2)
                
              }
              
              # converting time
              tmp <- data$DateHour
              data <- data[, -which(names(data) == "DateHour")]
              
              DateHour <- as.POSIXct(NA, "")
              for(i in 1:nrow(data)){
                
                d <- as.POSIXct(tmp[i], format = ifelse(nchar(as.character(tmp[i])) == 16,
                                                        "%d/%m/%Y %H:%M",
                                                        "%d/%m/%Y %H:%M:%S"))
                
                if(!is.na(d)){
                  
                  DateHour[i] <- d
                  
                } else {
                  
                  DateHour[i] <- as.POSIXct(paste0(tmp[i], " ", generateHour()),
                                            format = ifelse(nchar(as.character(paste0(tmp[i], " ", generateHour()))) == 16,
                                                            "%d/%m/%Y %H:%M",
                                                            "%d/%m/%Y %H:%M:%S"))
                  
                }
              }
              
              data$DateHour <- as.POSIXct(NA,"")
              data$DateHour <- DateHour
              
              rm(DateHour, tmp, d)
              
              # saving the current locale
              original_locale <- Sys.getlocale("LC_TIME")
              
              # setting the locale to English
              invisible(capture.output(Sys.setlocale("LC_TIME", "en_US.UTF-8")))
              
              # getting Day, DayType, Month, Year
              data$Year <- as.integer(format(data$DateHour, "%Y"))
              data$Month <- months(data$DateHour, abbreviate = TRUE)
              data$Day <- weekdays(data$DateHour, abbreviate = TRUE)
              data$DayType <- ifelse(data$Day != "Sat" & data$Day != "Sun", "weekday", "weekend")
              data$Hour <- format(data$DateHour, "%H:%M")
              
              # restoring the original locale
              invisible(capture.output(Sys.setlocale("LC_TIME", original_locale)))
              
              dataset <- rbind(dataset, data)
              data_pedestrians <- rbind(data_pedestrians, data[which(data$AccidentType == "Investimento di pedone"),])
              data_vehicles <- rbind(data_vehicles, data[which(data$AccidentType != "Investimento di pedone"),])
            }
          } else {
            
            m <- as.integer(m)
            setTxtProgressBar(pb, (y-years[1])*12+m)
            next
            
          }
          
          m <- as.integer(m)

          setTxtProgressBar(pb, (y-years[1])*12+m)

        }
      } 
    } else if(district == "Bergamo"){
  
      filename <- paste0(foldername, "road_accidents_", district, ".csv")
      
      if(file.exists(filename)){
        
        # reading the data
        data <- read.csv(filename, header = TRUE, sep = ",")
        
        # selecting only the variables of interest
        data <- data %>%
          select("Protocollo", "Data", "Ora", "NaturaIncidente", "Localizzazione")
        
        # reassigning variable names
        names(data) <- c("ID", "Date", "Hour", "AccidentType", "Coordinates")
        
        # keeping only one row per road accident
        data <- data %>% distinct(ID, .keep_all = TRUE)
        
        # discarding NAs (in longitudes and latitudes)
        data <- na.omit(data)
        
        data$Longitude <- unlist(lapply(strsplit(gsub("[(),]", "", data$Coordinates), "\\s* "), '[', 2))
        data$Latitude <- unlist(lapply(strsplit(gsub("[(),]", "", data$Coordinates), "\\s* "), '[', 1))
        data <- data[, -which(names(data) == "Coordinates")]
        
        # converting time
        tmp <- paste(data$Date, data$Hour)
        data$DateHour <- as.POSIXct(NA, "")
        data$DateHour <- as.POSIXct(tmp, format = "%d/%m/%Y %H:%M")
  
        rm(tmp)
  
        # saving the current locale
        original_locale <- Sys.getlocale("LC_TIME")
        
        # setting the locale to English
        invisible(capture.output(Sys.setlocale("LC_TIME", "en_US.UTF-8")))
        
        # getting Day, DayType, Month, Year
        data$Year <- as.integer(format(data$DateHour, "%Y"))
        data$Month <- months(data$DateHour, abbreviate = TRUE)
        data$Day <- weekdays(data$DateHour, abbreviate = TRUE)
        data$DayType <- ifelse(data$Day != "Sat" & data$Day != "Sun", "weekday", "weekend")
        data$Hour <- format(data$DateHour, "%H:%M")
        
        # restoring the original locale
        invisible(capture.output(Sys.setlocale("LC_TIME", original_locale)))
        
        data <- data[,variables]
        
        dataset <- rbind(dataset, data)
        data_pedestrians <- rbind(data_pedestrians, data[which(data$AccidentType == "Investimento di pedone"),])
        data_vehicles <- rbind(data_vehicles, data[which(data$AccidentType != "Investimento di pedone"),])
        
      }
    }
    
    dataset <- dataset[-1,]
    data_pedestrians <- data_pedestrians[-1,]
    data_vehicles <- data_vehicles[-1,]
    
    save(dataset, file = paste0(foldername, district, "_dataset.RData"))
    save(data_pedestrians, file = paste0(foldername, district, "_data_pedestrians.RData"))
    save(data_vehicles, file = paste0(foldername, district, "_data_vehicles.RData"))
    
  } else {
    
    load(paste0(foldername, district, "_dataset.RData"))
    load(paste0(foldername, district, "_data_pedestrians.RData"))
    load(paste0(foldername, district, "_data_vehicles.RData"))
    
  }
  
  return (list(dataset, data_pedestrians, data_vehicles))
  
}

loadDataset <- function(district = "Rome", years = NULL, mesh){
  
  foldername <- paste0(getwd(), "/utils/data/", district, "/")
  
  if(!file.exists(paste0(foldername, district, "_data.RData"))){
    
    dataset <- readDataset(district, years)
    data <- dataset[[1]]
    # data_pedestrians <- dataset[[2]]
    # data_vehicles <- dataset[[3]]
    
    idx_holidays <- as.numeric(filterHolidays(data = data, from = min(data$DateHour), to = max(data$DateHour)))
    data_holidays <- data[idx_holidays,]
    data_workdays <- data[-idx_holidays,]
    
    data_workdays <- data_workdays[which(data_workdays$Year >= years[1] & data_workdays$Year <= years[2]),]
    data_holidays <- data_holidays[which(data_holidays$Year >= years[1] & data_holidays$Year <= years[2]),]
    
    # Locations of Road Accidents
    locations_workdays <- cbind(data_workdays$Longitude, data_workdays$Latitude)
    locations_workdays <- fdaPDE::projection.points.1.5D(mesh, locations_workdays)
    locations_workdays <- data.frame(longitude = locations_workdays[,1], latitude = locations_workdays[,2])
    locations_holidays <- cbind(data_holidays$Longitude, data_holidays$Latitude)
    locations_holidays <- fdaPDE::projection.points.1.5D(mesh, locations_holidays)
    locations_holidays <- data.frame(longitude = locations_holidays[,1], latitude = locations_holidays[,2])
    
    # Times of Road Accidents
    hours <- as.numeric(format(strptime(data_workdays$Hour, format = "%H:%M"), "%H")) * 60
    minutes <- as.numeric(format(strptime(data_workdays$Hour, format = "%H:%M"), "%M"))
    times_workdays <- (hours + minutes) / (24 * 60)
    
    hours <- as.numeric(format(strptime(data_holidays$Hour, format = "%H:%M"), "%H")) * 60
    minutes <- as.numeric(format(strptime(data_holidays$Hour, format = "%H:%M"), "%M"))
    times_holidays <- (hours + minutes) / (24 * 60)
    
    data_workdays <- data.frame(lon = locations_workdays[,1],
                               lat = locations_workdays[,2],
                               times = times_workdays)
    data_holidays <- data.frame(lon = locations_holidays[,1],
                               lat = locations_holidays[,2],
                               times = times_holidays)
    
    save(data_workdays, data_holidays, file = paste0(foldername, district, "_data.RData"))
    
  } else {
    
    load(paste0(foldername, district, "_data.RData"))
      
  }
    
  return (list(data_workdays, data_holidays))
  
}

filterHolidays <- function(data, from, to){
  
  from_year <- as.integer(format(from, "%Y"))
  to_year <- as.integer(format(to, "%Y"))
  
  holidayIT <- holiday(year = from_year:to_year, Holiday = "EasterSunday")
  holidayIT <- c(holidayIT, holiday(year = from_year:to_year, Holiday = "EasterMonday"))  
  
  for(y in from_year:to_year){
    holidayIT <- c(holidayIT, as.timeDate(paste0(y, "-01-01")),
                              as.timeDate(paste0(y, "-01-06")),
                              as.timeDate(paste0(y, "-04-25")),
                              as.timeDate(paste0(y, "-05-01")),
                              as.timeDate(paste0(y, "-06-02")),
                              as.timeDate(paste0(y, "-08-15")),
                              as.timeDate(paste0(y, "-11-01")),
                              as.timeDate(paste0(y, "-12-08")),
                              as.timeDate(paste0(y, "-12-25")),
                              as.timeDate(paste0(y, "-12-26")))
  }
  
  return(which(isHoliday(as.timeDate(data$DateHour), holidays = holidayIT)))
    
}


# ---
# folder.img <- paste0(foldername,"imgs/")
# if(!dir.exists(folder.img))
#   dir.create(folder.img)
# 
# incidenti_bbox <- data.frame(x.min = (min(mesh$nodes[,1]) + 0.1 * sd(mesh$nodes[,1])), 
#                              x.max = (max(mesh$nodes[,1]) - 0.1 * sd(mesh$nodes[,1])),
#                              y.min = (min(mesh$nodes[,2]) + 0.1 * sd(mesh$nodes[,2])), 
#                              y.max = (max(mesh$nodes[,2]) - 0.1 * sd(mesh$nodes[,2])))
# 
# crashes_MI <- load_data()
# 
# mesi <- c("Jan", "Feb", "Mar", "Apr", "May", "June", "July",
#           "Aug", "Sep", "Oct", "Nov", "Dec")
# 
# num_data <- matrix(0, nrow=12,ncol=1)
# LOCATIONS <- list()
# IDXS <- list()
# for(i in 1:12){
#   idxs =  intersect(which(crashes_MI$MESE_INCIDENTE == i),
#                     intersect(
#                       intersect(
#                         which(crashes_MI$longitudine >= incidenti_bbox$x.min),
#                         which(crashes_MI$longitudine <= incidenti_bbox$x.max)), 
#                       intersect(
#                         which(crashes_MI$latitudine >= incidenti_bbox$y.min), 
#                         which(crashes_MI$latitudine <= incidenti_bbox$y.max))
#                     ))
#   IDXS[[i]] <- idxs
#   
#   locations <-  data.frame(lon = crashes_MI$longitudine[idxs], 
#                            lat = crashes_MI$latitudine[idxs])
#   
#   LOCATIONS[[i]] <- locations
#   
#   locations <- sf::st_as_sf(locations, 
#                             coords = c("lon", "lat"),
#                             crs = 4269) 
#   num_data[i] <- nrow(locations)
#   {  
#     milano_incidenti_plot <- ggplot() +
#       geom_sf(data = st_as_sf(filtered,"edges"), 
#               inherit.aes = FALSE, 
#               color = "black", 
#               size = .0005,
#               alpha = .6) +
#       geom_sf(data = locations, 
#               inherit.aes = FALSE, 
#               color = "red3", 
#               size = 0.75,
#               alpha = 1.) +
#       geom_rect( aes(xmin=incidenti_bbox$x.min, xmax=incidenti_bbox$x.max, 
#                      ymin=incidenti_bbox$y.min, ymax=incidenti_bbox$y.max), 
#                  color="red", fill="transparent")+
#       theme_bw() +
#       labs(title = paste("car crashes (", mesi[i], " 2018)",sep=""),
#            x = "Longitudine",
#            y = "Latitudine")+
#       theme(plot.title = element_text(hjust = 0.5))
#   }
#   
#   ggsave(filename = paste0(folder.img, district,"_crashes_", i,"_.pdf"), 
#          plot=milano_incidenti_plot,
#          width=6, height=6, units = "in")
#   
# }
# 
# locs <- matrix(nrow=0,ncol=2)
# idxs <- rep(0,length=0)
# crashes <- matrix(nrow=0, ncol=ncol(crashes_MI))
# for(i in 1:12){
#   idxs <- c(idxs, IDXS[[i]])
#   crashes <- rbind(crashes, crashes_MI[IDXS[[i]],])
#   locs <- rbind(locs, LOCATIONS[[i]])
# }
# 
# crashes$ID <- idxs
# save(LOCATIONS, IDXS, crashes, file=paste0(foldername,"data_raw.RData"))
# 
# # ------------------------------------------------------------------------------
# # projecting 
# locs <- fdaPDE::projection.points.1.5D(mesh, locs)
# crashes$longitudine <- locs[,1]
# crashes$latitudine <- locs[,2]
# 
# LOCS <- data.frame(lon = locs[,1], lat = locs[,2])
# LOCS <- sf::st_as_sf(LOCS, coords = c("lon", "lat"), crs = 4269) 
# {  
#   incidenti_plot <- ggplot() +
#     geom_sf(data = st_as_sf(filtered,"edges"), 
#             inherit.aes = FALSE, 
#             color = "black", 
#             size = .0005,
#             alpha = .6) +
#     geom_sf(data = LOCS, 
#             inherit.aes = FALSE, 
#             color = "red3", 
#             size = 0.75,
#             alpha = 1.) +
#     geom_rect( aes(xmin=incidenti_bbox$x.min, xmax=incidenti_bbox$x.max, 
#                    ymin=incidenti_bbox$y.min, ymax=incidenti_bbox$y.max), 
#                color="red", fill="transparent")+
#     theme_bw() +
#     labs(title = paste0(district," crashes (2018)",sep=""),
#          x = "Longitudine",
#          y = "Latitudine")+
#     theme(plot.title = element_text(hjust = 0.5))
# }
# 
# ggsave(filename = paste0(folder.img, district, "_crashes_2018.pdf"), 
#        plot=incidenti_plot,
#        width=6, height=6, units = "in")
