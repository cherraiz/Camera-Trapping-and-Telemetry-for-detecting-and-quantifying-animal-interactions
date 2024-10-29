telm.contacts <- function(data, ID = "ID",x = "x", y = "y", species = "sp",
                          date_time = "date_time", spatial_threshold, 
                          temporal_threshold) {
  
  stopifnot(lubridate::is.POSIXct(data$date_time[[1]]))
 
   data = data %>%
    rename(ID = ID, species = species, date_time = date_time, x = x, y = y) %>%
    data.frame()
  
  data$ID_loc<-c(1:dim(data)[1])
  dist = c()
  tdiff = c()
  datetime = c()
  ID_depositor = c()
  ID_acquirer = c()
  sp_depositor = c()
  sp_acquirer = c()
  loc_x = c()
  loc_y = c()
  
  for (i in 1:length(unique(data$ID_loc))) {
    ID.i = unique(data$ID_loc)[i]
    data.i = data[data$ID_loc == ID.i,]
    x.i = data.i$x
    y.i = data.i$y
    fh_ref = data.i$date_time
    Indiv = data.i$ID
    ID_ref = data.i$ID_loc
    sp.i = data.i$species
    
    data.j = data[data$ID != Indiv,]
    data.j$hours = difftime(data.j$date_time, fh_ref, units = "hours")
    data.j = data.j[data.j$hours <= temporal_threshold & data.j$hours >= 0,]
    data.j$dist.i = sqrt((x.i - data.j$x)^2 + (y.i - data.j$y)^2)
    data.j = data.j[data.j$dist.i <= spatial_threshold, ]
    
    if (dim(data.j)[1] > 0) {
      for (j in 1:dim(data.j)[1]) {
        x.j = data.j$x[j]
        y.j = data.j$y[j]
        tdiff.ij = data.j$hours[j]
        dist.ij = data.j$dist.i[j]
        Indiv_par = data.j$ID[j]
        ID_par = data.j$ID_loc[j]
        sp.j = data.j$species[j]
        
        cat("\014") 
        cat(paste0(round(i / length(unique(data$ID_loc)), digits=3) * 100, 
                   '% completed...')) 
        Sys.sleep(.001)
        if (i == length(unique(data$ID_loc))) cat(': Done')
        
        xmean = ((abs(x.i - x.j)) / 2) + min(x.i, x.j)
        ymean = ((abs(y.i - y.j)) / 2) + min(y.i, y.j)
        
        loc_x = c(loc_x, xmean)
        loc_y = c(loc_y, ymean)
        dist = c(dist, dist.ij)
        tdiff = c(tdiff, tdiff.ij)
        datetime = c(datetime, fh_ref)
        ID_depositor = c(ID_depositor, Indiv)
        ID_acquirer = c(ID_acquirer, Indiv_par)
        sp_depositor = c(sp_depositor, sp.i)
        sp_acquirer = c(sp_acquirer, sp.j)
      }
    }
  }
  
  contacts <- data.frame(
    ID_depositor = ID_depositor,
    sp_depositor = sp_depositor,
    ID_acquirer = ID_acquirer,
    sp_acquirer = sp_acquirer,
    DateTime_deposition = as.POSIXct(datetime, origin = "1970-01-01",
                                     tz="Europe/Paris"),
    tdiff = as.difftime(tdiff, units = "hours"),
    Distance = dist,
    x = loc_x,
    y = loc_y)
  
  return(contacts)
}

ct.contacts <- function(data, CT = "id", species = "sp",
                        date_time = "date_time", temporal_threshold,
                        temp_gr_threshold) {
 
  # Assign group ID based on time difference and compute visit duration
  data. <- data %>% 
    arrange(!!as.name(CT), !!as.name(species), !!as.name(date_time)) %>%
    group_by(!!as.name(species), !!as.name(CT)) %>%
    mutate(group_id = 
             cumsum(abs(as.numeric(!!as.name(date_time) - 
                                     lag(!!as.name(date_time),
                                         default = min(!!as.name(date_time))), 
                                   units = "mins")) >= temp_gr_threshold) + 1) %>% 
    mutate(group_id = paste0(!!as.name(CT),!!as.name(species),"_",group_id)) %>% 
    group_by(group_id) %>%  mutate(
      duration = difftime(max(!!as.name(date_time)), min(!!as.name(date_time)),
                          units = "mins")) %>%
    rename(CT = CT, species = species, date_time = date_time) %>%
    data.frame()
  
  stopifnot(lubridate::is.POSIXct(data.$date_time[[1]]))
  
  ct_ID <- c()
  tdiff = c()
  datetime = c()
  sp_depositor = c()
  sp_acquirer = c()
  duration_dep = c()
  duration_acq = c()
  
  for (i in 1:length(unique(data.$group_id))) {
    ID.i = unique(data.$group_id)[i]
    data.i = data.[data.$group_id == ID.i,]
    ct.i = data.i$CT[1]
    fh_ref = data.i$date_time
    sp.i = data.i$species[1]
    dur.i = data.i$duration[1]
    
    data.j = data.[data.$group_id != ID.i,]
    data.j = data.j[data.j$CT == ct.i,]
    
    if (dim(data.j)[1] > 0) {
      data.j = data.j %>% group_by(group_id) %>% mutate(
        tdiff.min = difftime(min(date_time), max(fh_ref), units = "hours")) %>%
        mutate(
          tdiff.max = difftime(max(date_time), min(fh_ref), units = "hours")) %>%
        dplyr::select(-date_time)
      data.j$tdiff.min = ifelse(data.j$tdiff.min < 0 & 
                                  data.j$tdiff.max >= 0, 0, data.j$tdiff.min)
      data.j = data.j[data.j$tdiff.min <= temporal_threshold & 
                        data.j$tdiff.min >= 0,]
      data.j = data.j[!duplicated(data.j$group_id),] %>% data.frame()
    }
    
    if (dim(data.j)[1] > 0) {
      for (j in 1:length(unique(data.j$group_id))) {
        tdiff.ij = data.j$tdiff.min[j]
        sp.j = data.j$species[j]
        dur.j = data.j$duration[j]
        
        cat("\014") 
        cat(paste0(round(i / length(unique(data.$group_id)) * 100), 
                   '% completed...')) 
        Sys.sleep(0.001)
        if (i == length(unique(data.$group_id))) cat(': Done')
        
        ct_ID = c(ct_ID, ct.i)
        tdiff = c(tdiff, tdiff.ij)
        datetime = c(datetime, min(fh_ref))
        sp_depositor = c(sp_depositor, sp.i)
        sp_acquirer = c(sp_acquirer, sp.j)
        duration_dep = c(duration_dep, dur.i)
        duration_acq = c(duration_acq, dur.j)
      }
    }
  }
  
  contacts <- data.frame(
    ct_ID = ct_ID,
    sp_depositor = sp_depositor,
    duration_dep = as.difftime(duration_dep, units = "mins"),
    sp_acquirer = sp_acquirer,
    duration_acq = as.difftime(duration_acq, units = "mins"),
    DateTime_dep = as.POSIXct(datetime, origin = "1970-01-01",
                              tz="Europe/Paris"),
    tdiff = as.difftime(tdiff, units = "hours"))
  
  return(contacts)
}
