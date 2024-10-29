library(terra)
library(sf)
library(tidyverse)
library(ggplot2)
library(data.table)
library(ggpubr)
library(viridis)


# 1. Simulation ----
# Load traj simulatior function
source("~/traj_simulator_function.R")

# 1.1. Water points ----
# Any simulated points in any simulated frame area would work
res.lis <- list(no.agr = terra::vect(st_transform(st_read(
  'GIS/waterpoints_simulation_grid.shp'), crs = 32629)),
  med.agr = terra::vect(st_transform(st_read(
    'GIS/waterpoints_simulation_medium_agr.shp'), crs = 32629)),
  high.agr = terra::vect(st_transform(st_read(
    'GIS/waterpoints_simulation_agrreggated.shp'), crs = 32629)))

# 1.2. Defining resource needs
theta.A <- c(-0.87, -1.22, -1.56) # 2.5, 3.5, and 4.5 daily visits
theta.B <- c(-0.52, -0.87, -1.22) # 1.5, 2.5, and 3.5 daily visits
theta.C <- c(-0.35, -0.69, -1.04) # 1, 2, and 3 daily visits

# 1.3. Looping through resource needs
res.needs <- c("low", "medium", "high")

counter = 0
for (m in res.lis[3]) {
  counter = counter +1
  resource.dist = distance.rast(m, crs = "EPSG:32629", 
                                area = st_transform(st_read('GIS/frame.shp'), 
                                                    crs = 32629), resolution = 1)

  for (i in 3) {
    A = species(HRkm2 = 4, meanT = 70, sdT = 35, meanF = 8, sdF = 4,
                rhoT = 0.7, rhoF = 0.2, pTF = 0.07, pFT = 0.03, 
                theta = c(theta.A[i], 100), energy = 50)
    B = species(HRkm2 = 6, meanT = 40, sdT = 20, meanF = 4, sdF = 2,
                rhoT = 0.8, rhoF = 0.3, pTF = 0.07, pFT = 0.06, 
                theta = c(theta.B[i], 100), energy = 50)
    C = species(HRkm2 = 8, meanT = 50, sdT = 25, meanF = 6, sdF = 3,
                rhoT = 0.9, rhoF = 0.5, pTF = 0.1, pFT = 0.07, 
                theta = c(theta.B[i], 100), energy = 50)
    empty.traj = data.frame(matrix(nrow = 0, ncol = 8))
    colnames(empty.traj) = c("ID","x", "y", "behavior", "angle", "distance", "energy", 
                             "step")
    Atraj = empty.traj
    Btraj = empty.traj
    Ctraj = empty.traj
    
    for (j in 1:50) {
      value = NA
      while(is.na(value) | value > 500){  # Less than 500 m from resource
        init.coords = spatSample(resource.dist, size=1, xy=TRUE)
        value = init.coords[3]}
      simu.traj =  sim.traj.rech(sp = A, env = m, n = 8640, 
                                 init.coords = init.coords[1:2], 
                                 resolution = 1, 
                                 res.dist = 5)
      simu.traj$ID = paste0("A", j)
      Atraj = rbind(Atraj, simu.traj)
    }
    
    for (j in 1:50) {
      value = NA
      while(is.na(value) | value > 500){
        init.coords = spatSample(resource.dist, size=1, xy=TRUE)
        value = init.coords[3]}
      
      simu.traj =  sim.traj.rech(sp = B, env = m, n = 8640, 
                                 init.coords = init.coords[1:2], 
                                 resolution = 1, 
                                 res.dist = 5)
      simu.traj$ID = paste0("B", j)
      Btraj = rbind(Btraj,simu.traj)
    }
    
    for (j in 1:50) {
      value = NA
      while(is.na(value) | value > 500){
        init.coords = spatSample(resource.dist, size=1, xy=TRUE)
        value = init.coords[3]}
      simu.traj =  sim.traj.rech(sp = C, env = m, n = 8640, 
                                 init.coords = init.coords[1:2], 
                                 resolution = 1, 
                                 res.dist = 5)
      simu.traj$ID = paste0("C", j)
      Ctraj = rbind(Ctraj, simu.traj)
    }
    
    Atraj$sp = "A"
    Btraj$sp = "B"
    Ctraj$sp = "C"
    
    traj = rbind(Atraj, rbind(Btraj, Ctraj))
    agr = ifelse(counter == 1, "no.agr.", ifelse(counter == 2, "medium.agr.", 
                                                 "high.agr."))
    assign(paste0(agr,res.needs[i],".need"), value = traj)
    rm(list = c("value", "Atraj", "Btraj", "Ctraj", "traj", "A",
                "B", "C", "empty.traj", "init.coords"))
  }
}

sce.1 <- no.agr.low.need
sce.2 <- no.agr.medium.need
sce.3 <- no.agr.high.need
sce.4 <- medium.agr.low.need
sce.5 <- medium.agr.medium.need
sce.6 <- medium.agr.high.need
sce.7 <- high.agr.low.need
sce.8 <- high.agr.medium.need
sce.9 <- high.agr.high.need


# 2. Estimating contacts ----
# 2.1. Preparing data ----
# Adding time to the simulated data and deleting first 3 steps
for (i in paste0("sce.", 1:9)){
  scetime <- get(i) %>% 
    mutate(datetime = as.POSIXct(step*5*60, origin = torigin, 
                                 tz = Sys.timezone(location = TRUE))) %>% 
    group_by(ID) %>% slice(3:n())
  assign(i, scetime)
}

# 2.2. Locations to linestrings ----
for (k in paste0("sce.", 1:9)){
  copy <- get(k)
  line0 = sfheaders::sf_linestring(copy[1:2,], x = "x", y = "y", 
                                   keep = TRUE)[0,]
  for (j in unique(copy$ID)) {
    df = copy[copy$ID == j,]
    line = sfheaders::sf_linestring(df, x = "x", y = "y", keep = TRUE)[0,]
    for (i in 1:(nrow(df)-1)) {
      line <- bind_rows(line, sfheaders::sf_linestring(obj = df[i:(i+1),],
                                                       x = "x", y = "y",
                                                       keep = TRUE))
    }
    cat("\014") 
    cat(paste0(round(match(j, unique(copy$ID)) / 
                       length(unique(copy$ID)) * 100), 
               '% completed...')) 
    Sys.sleep(0.001)
    if (match(j, unique(copy$ID)) == 
        length(unique(copy$ID))) cat(': Done')
    line0 = bind_rows(line0, line)
  }
  assign(paste0("line.", k), value = line0)
}

# 2.3. Loading resource points as sf ----
study.area <- st_transform(st_read("starting_area.shp"), crs = 32629)
res.no.agr <- st_transform(st_read('waterpoints_simulation_grid.shp'), 
                           crs = 32629)
res.med.agr <- st_transform(st_read('waterpoints_simulation_medium_agr.shp'), 
                            crs = 32629)
res.high.agr <- st_transform(st_read('waterpoints_simulation_agrreggated.shp'), 
                             crs = 32629)

# 2.4. Obtaining telemetry data ----
# Randomly collaring 10 individuals of each species
for (m in 1:10){
  telm.effort <- c(48, 24, 12, 6, 4, 2)
  for (i in paste0("sce.", 1:9)){
    telm <- get(i) %>% select(sp, ID) %>% unique() %>% group_by(sp) %>% 
      slice_sample(n = 10)  # 10 individuals of each species
    telm <- get(i)[get(i)$ID %in% telm$ID,]
    
    # Resampling one every x steps (2 for 10 min, 6 for 30 min, etc.)
    for (j in 1:6){
      telm2 <- telm[seq(1, nrow(telm), telm.effort[j]),]
      assign(paste0("telm.", i, ".effort.", j, ".rep.", m), telm2)
    }
  }
}

# 2.5. Obtaining camera trap data ----
# Placing camera traps
for (m in 1:10){
  grid.effort <- c(7,11,13,17,22,33) # less than desired for adjusting the grid
  res.effort <- c(5,10,15,20,25,30)
  for (i in 1:6){
    # Grid
    ct.grid <- spatSample(vect(study.area), method = "regular", size = grid.effort[i]) %>%
      st_as_sf() %>% mutate(id = paste0("G", row_number())) %>% 
      mutate(type = "Grid") %>% cbind(data.frame(st_coordinates(.)))
    
    # Adding random monitoring for the resources
    ct.no.agr <- res.no.agr %>% 
      select(geometry) %>% slice_sample(n = res.effort[i]) %>% 
      mutate(id = paste0("A", row_number())) %>% mutate(type = "Agreggation") %>% 
      cbind(data.frame(st_coordinates(.))) %>% rbind(ct.grid)
    ct.med.agr <- res.med.agr %>% 
      select(geometry)  %>% slice_sample(n = res.effort[i]) %>% 
      mutate(id = paste0("A", row_number())) %>% mutate(type = "Agreggation") %>% 
      cbind(data.frame(st_coordinates(.))) %>% rbind(ct.grid)
    ct.high.agr <- res.high.agr %>% 
      select(geometry)  %>% slice_sample(n = res.effort[i]) %>% 
      mutate(id = paste0("A", row_number())) %>% mutate(type = "Agreggation") %>% 
      cbind(data.frame(st_coordinates(.))) %>% rbind(ct.grid)
    
    # Assign names
    assign(paste0("ct.no.agr.effort.", i, ".rep.", m), ct.no.agr)
    assign(paste0("ct.med.agr.effort.", i, ".rep.", m), ct.med.agr)
    assign(paste0("ct.high.agr.effort.", i, ".rep.", m), ct.high.agr)
  }
}

# Capturing animals' trajectories with the camera traps
for (i in ls(pattern = "^line.sce")) {
  traj <- sf::st_set_crs(get(i), 32629)
  sc <- as.numeric(gsub(".*?([0-9]+).*", "\\1", i))
  cts <- ifelse(sc <= 3, "ct.no.agr", 
                ifelse(sc >= 7, "ct.high.agr", "ct.med.agr"))
  cts <- ls(pattern = cts)
  for (j in cts){
    photos <- traj %>% select(sp, datetime, geometry) %>%
      st_intersection(st_buffer(get(j), dist = 5, endCapStyle = "SQUARE")) %>%
      data.frame() %>% select(-geometry) # 100m2 of detection
    
    assign(paste0("photo.", gsub("line.", "", i), ".",
                  sub("^([^.]+.[^.]+.[^.]+)\\.*", "\\2", j)), photos)
    
    cat("\014")
    cat(paste0(round(match(j, cts)/length(cts)*100), " % \n scenario ",
               match(i, ls(pattern = "^line.sce")), " of ",
               length(ls(pattern = "^line.sce"))))
  }
}

# 3. Contact estimation ----
# Load contact functions
source("Rscripts/contact_functions.R")

# 3.1. Telemetry contacts estimation ----
for (i in ls(pattern = "^telm.sce")[7:54]){
  contacts.telm <- telm.contacts(data = get(i), spatial_threshold = 25, ID = "ID",
                                 x = "x", y = "y", species = "sp",
                                 temporal_threshold = 48, date_time = "datetime")
  assign(paste0("contacts.", i), contacts.telm)
}

# 3.2. Camera trapping contacts estimation ----
for (i in ls(pattern = "^photo.sce")){
  contacts.ct <- ct.contacts(data = get(i), temporal_threshold = 48, CT = "id",
                             species = "sp", date_time = "datetime", 
                             temp_gr_threshold = 10)
  assign(paste0("contacts.", i), contacts.ct)
}


# 4. Approaches' comparison
# 4.1. Summarizing interactions by species pair ----
res <- data.frame(Interact = character(),
                  Scenario = integer(), 
                  Method = character(), 
                  Effort = integer(),
                  rep = integer(),
                  direct = numeric(),
                  indirect = numeric(),
                  total = numeric())
res.real <- data.frame(Scenario = integer(), 
                       direct.real = numeric(),
                       indirect.real = numeric(),
                       total.real = numeric())

int <- data.frame(Interact = c("AA", "AB", "AC", "BA", "BB", "BC", "CA", "CB", "CC"))

for(i in 1:9) {
  for (j in 1:6) {
    for (n in 1:10) {
      telm <- get(paste0("contacts.telm.sce.", i, ".effort.",j, 
                         ".rep.", n)) %>%
        mutate(Interact = paste0(sp_depositor, sp_acquirer)) %>%
        group_by(Interact) %>% 
        mutate(direct = sum(tdiff == 0)) %>% 
        mutate(indirect = sum(tdiff > 0)) %>%
        mutate(total = direct + indirect) %>%
        dplyr::select(Interact, direct, indirect, total) %>% 
        summarise(across(.cols = 1:3,.fns = max)) %>%
        merge(int, all.y = T) %>% 
        replace(is.na(.), 0)
      
      telm$Scenario = i
      telm$Method = "telm"
      telm$Effort = j
      telm$rep = n
      
      res <- rbind(res, telm)
      
      for (m in c("G", "A", "all")) {
        ct <- get(paste0("contacts.photo.sce.", i, ".effort.",j,
                         ".rep.", n)) %>%
          filter(if (m != "all") {gsub('[[:digit:]]+', '', ct_ID) == m} else TRUE) %>%
          mutate(Interact = paste0(sp_depositor, sp_acquirer)) %>%
          group_by(Interact) %>% 
          mutate(direct = sum(tdiff == 0)) %>% 
          mutate(indirect = sum(tdiff > 0)) %>%
          mutate(total = direct + indirect) %>%
          dplyr::select(Interact, direct, indirect, total) %>% 
          summarise(across(.cols = 1:3,.fns = max)) %>%
          merge(int, all.y = T) %>% 
          replace(is.na(.), 0)
        
        ct$Scenario = i
        ct$Method = paste0("ct.", m)
        ct$Effort = j
        ct$rep = n
        
        res <- rbind(res, ct)
      }
    }
    real <- get(ls(pattern = paste0("^contacts.real.sce.", i))) %>%
      mutate(Interact = paste0(sp_depositor, sp_acquirer)) %>%
      group_by(Interact) %>% 
      mutate(direct.real = sum(tdiff == 0)) %>% 
      mutate(indirect.real = sum(tdiff > 0)) %>%
      mutate(total.real = direct.real + indirect.real) %>%
      dplyr::select(Interact, direct.real, indirect.real, total.real) %>% 
      summarise(across(.cols = 1:3,.fns = max))
    
    real$Scenario = i
    
    res.real <- rbind(res.real, real)
  }
}

res <- merge(unique(res), unique(res.real))
res <- res[substr(res$Interact, start = 1, stop = 1) != 
             substr(res$Interact, start = 2, stop = 2),]

# 4.2. Approaches' goodness of fit ----
method <- data.frame(Method = character(), 
                     Effort = integer(),
                     directR2 = numeric(),
                     indirectR2 = numeric())

counter = 0
for(i in c("telm", "ct.G", "ct.A", "ct.all")){
  for (j in 1:6){
    counter = counter +1
    m.dir <- glm(formula =  direct.real ~ direct, 
                 data = res[res$Method == i & res$Effort == j,], 
                 family = poisson())
    m.indir <- glm(formula =  indirect.real ~ indirect, 
                   data = res[res$Method == i & res$Effort == j,],
                   family = poisson())
    
    method[counter, 1] <- i
    method[counter, 2] <- j
    method[counter, 3] <- with(summary(m.dir), 1 - deviance/null.deviance)
    method[counter, 4] <- with(summary(m.indir), 1 - deviance/null.deviance)
  }
}

method$Method <- factor(method$Method, levels = c("telm", "ct.G","ct.A", "ct.all"))
method <- method %>% 
  mutate(Method2 = recode(Method,
                          ct.G = "CT: Regular Grid",
                          telm = "GPS",
                          ct.A = "CT: Resource Monitoring",
                          ct.all = "CT: all"
  ))

ggplot(data = method, mapping = aes(x = Effort, y = directR2, 
                                    col = Method2, shape = Method2)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.3)  +
  scale_color_viridis(discrete = TRUE) + 
  scale_shape_manual(values = c(17, 16, 15, 18))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  scale_y_continuous(name = expression("Direct contacts estimation R"^2),  
                     limits = c(0,0.8),
                     breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)) +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title= element_blank())

ggplot(data = method, mapping = aes(x = Effort, y = indirectR2, 
                                    col = Method2, shape = Method2)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.3)  +
  scale_color_viridis(discrete = TRUE) + 
  scale_shape_manual(values = c(17, 16, 15, 18))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  scale_y_continuous(name = expression("Indirect contacts estimation R"^2), 
                     limits = c(0,0.8),
                     breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)) +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title= element_blank())


# 4.3. Approaches' missing pairs ----
sample <- res %>% 
  count(Scenario, Method, Effort, rep) %>% 
  group_by(Method, Effort) %>% 
  summarize()

missdirect <- res[res$direct == 0 & res$direct.real != 0,]
missindirect <- res[res$indirect ==0 & res$indirect.real != 0,]

missdirect <- missdirect %>% count(Scenario, Method, Effort, rep) %>% 
  merge(count(res, Scenario, Method, Effort, rep)[,1:4], all.y = TRUE) %>% 
  replace(is.na(.), 0) %>%
  group_by(Method, Effort) %>%
  summarize(median = median(n), max = max(n), min = min(n)) %>%
  merge(sample, all.y =  TRUE) %>%
  replace(is.na(.), 0)
missdirect$Method <- factor(missdirect$Method, levels = c("telm", "ct.G","ct.A", "ct.all"))
missdirect <- missdirect %>% 
  mutate(Method2 = recode(Method,
                          ct.G = "CT: Regular Grid",
                          telm = "GPS",
                          ct.A = "CT: Resource Monitoring",
                          ct.all = "CT: all"
  ))

missindirect <- missindirect %>% count(Scenario, Method, Effort, rep) %>% 
  merge(count(res, Scenario, Method, Effort, rep)[,1:4], all.y = TRUE) %>% 
  replace(is.na(.), 0) %>%
  group_by(Method, Effort) %>%
  summarize(median = median(n), max = max(n), min = min(n)) %>%
  merge(sample, all.y =  TRUE) %>%
  replace(is.na(.), 0)
missindirect$Method <- factor(missindirect$Method, levels = c("telm", "ct.G","ct.A", "ct.all"))
missindirect <- missindirect %>%  mutate(Method2 = recode(Method,
                                                          ct.G = "CT: Regular Grid",
                                                          telm = "GPS",
                                                          ct.A = "CT: Resource Monitoring",
                                                          ct.all = "CT: all"
))
rm(sample)

ggplot(missdirect, aes(x=Effort, y=median, group = Method2)) +
  geom_errorbar(aes(ymin=min, ymax=max, color=Method2), width=0.7, 
                position = position_dodge(0.9),alpha=0.9, size=1.7) + 
  geom_point(aes(fill=Method2), size=5, position = position_dodge(0.9), shape = 21) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  scale_y_continuous(name = "Missing pairs", breaks = c(0,1,2,3,4,5,6)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title= element_blank())


ggplot(missindirect, aes(x=Effort, y=median, group = Method2)) +
  geom_errorbar(aes(ymin=min, ymax=max, color=Method2), width=0.7, 
                position = position_dodge(.9),alpha=0.9, size=1.7) + 
  geom_point(aes(fill=Method2), size=5, position = position_dodge(0.9), shape = 21) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  scale_y_continuous(name = "Missing pairs", breaks = c(0,1,2,3,4,5,6), limits = c(0,6)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.title= element_blank())


# 5. Factors influencing performance ----
res$Method <- factor(res$Method, levels = c("telm", "ct.G","ct.A", "ct.all"))
res <- res %>% mutate(Aggregation = case_when(
  Scenario <= 3 ~ 1,
  Scenario >= 7 ~ 3,
  TRUE ~ 2)) %>%
  mutate(Frequency = case_when(
    Scenario == 1 | Scenario == 4 | Scenario == 7 ~ 1,
    Scenario == 2 | Scenario == 5 | Scenario == 8 ~ 2,
    Scenario == 3 | Scenario == 6 | Scenario == 9 ~ 3
  )) %>%  mutate(Method2 = dplyr::recode(Method,
                                         ct.G = "CT-RG",
                                         telm = "GPS",
                                         ct.A = "CT-RM",
                                         ct.all = "CT-All"
  ))

res$Aggregation <- ordered(res$Aggregation, levels = c(1,2,3))
res$Frequency <- ordered(res$Frequency, levels = c(1,2,3))
res$Effort <- ordered(res$Effort, levels = c(1,2,3,4,5,6))

effects.dir <- glm.nb(formula = direct ~ Method2 + Aggregation + Frequency  + Effort +
                        Aggregation * Method2 + Frequency * Method2 + Effort * Method2,
                      data = res)

summary(effects.dir)

effects.indir <- glm.nb(formula = indirect ~ Method2 + Aggregation + Frequency  + Effort +
                          Aggregation * Method2 + Frequency * Method2 + Effort * Method2,
                        data = res)
summary(effects.indir)

res <- res %>% mutate(preddir = predict.glm(effects.dir, type = "response"),
                      predindir = predict.glm(effects.indir, type = "response"))

ggarrange(ggplot(res,  aes(x = Method2, y = preddir, color = Effort)) +
            geom_boxplot(size = 1) + theme_minimal() + 
            scale_y_log10(labels = c(0.001, parse(text="10^-2"), parse(text="10^-1"), 
                                     parse(text="10^0"), parse(text="10^1"), 
                                     parse(text="10^2"), 1000)) +
            scale_color_viridis_d() + ylab("Predicted\ndirect interactions") +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                  legend.title=element_blank(), legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5, size = 10),
                  text = element_text(size=10),
                  axis.title.y = element_text(margin = margin(r = 10))) + 
            guides(color = guide_legend(nrow = 1),) +
            ggtitle("Sampling\neffort")
          ,
          ggplot(res,  aes(x = Method2, y = preddir, color = Aggregation)) +
            geom_boxplot(size = 1) + theme_minimal() + scale_y_log10() +
            scale_color_viridis_d(end = 0.5) +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                  axis.title.y=element_blank(), axis.text.y=element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 10),
                  text = element_text(size=10)) +
            ggtitle("Resource\nspatial aggregation")
          ,
          ggplot(res,  aes(x = Method2, y = preddir, color = Frequency)) +
            geom_boxplot(size = 1) + theme_minimal() + scale_y_log10() +
            scale_color_viridis_d(end = 0.5) +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                  axis.title.y=element_blank(), axis.text.y=element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 10),
                  text = element_text(size=10)) +
            ggtitle("Frequency of visits\nto the resource")
          ,
          ggplot(res,  aes(x = Method2, y = predindir, color = Effort)) +
            geom_boxplot(size = 1) + theme_minimal() + 
            scale_y_log10(labels = c(1, parse(text="10^1"), parse(text="10^2"), 
                                     parse(text="10^3"), parse(text="10^4"), 
                                     parse(text="10^5"), 1000000)) +
            ylab("Predicted\nindirect interactions") +
            theme(axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.title.y = element_text(margin = margin(r = 10)),
                  text = element_text(size=10))
          ,
          ggplot(res,  aes(x = Method2, y = predindir, color = Aggregation)) +
            geom_boxplot(size = 1) + theme_minimal() + 
            scale_y_log10() +
            scale_color_viridis_d(end = 0.5) +
            theme(axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.title.y=element_blank(), axis.text.y=element_blank(),
                  text = element_text(size=10))
          ,
          ggplot(res,  aes(x = Method2, y = predindir, color = Frequency)) +
            geom_boxplot(size = 1) + theme_minimal() + scale_y_log10() +
            scale_color_viridis_d(end = 0.5) +
            theme(axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.title.y=element_blank(), axis.text.y=element_blank(),
                  text = element_text(size=10))
          , nrow = 2, ncol = 3, common.legend = T, legend = "bottom",
          heights = c(1.03, 1), widths = c(1.2, 1, 1))


# 6. Doñana National Park analysis ----
# 6.1. Doñana contacts estimation ----
ctdn <- read.csv("databases/cameratrap_data.csv")

ctdn$date_time <- as.POSIXct(ctdn$date_time, 
                             tz = "Europe/Paris", 
                             format = "%Y-%m-%d %H:%M:%OS")

ctdn <- ctdn %>% 
  arrange(ct_ID, sp, date_time) %>%
  group_by(ct_ID, sp) %>%
  mutate(tdiff = round((c(0, diff(date_time)))/60, digits = 1))

breaks <- getJenksBreaks(ctdn$tdiff[ctdn$tdiff < 30], k = 5)
breaks <- c(2.5, 10, 20, 30)

hist(ctdn$tdiff[ctdn$tdiff < 30])

for (i in breaks){
  contact.dn <- ct.contacts(data = ctdn, 
                            species = "sp", 
                            CT = "ct_ID", 
                            date_time = "date_time", 
                            temporal_threshold = 9, 
                            temp_gr_threshold = i)
  assign(paste0("contact.dn.thres.", i), contact.dn)
  
}

ctdn$day <- format(as.Date(ctdn$date_time,format="%Y-%m-%d %H:%M:%OS"), 
                   format = "%Y-%m-%d")

# Taking species with at least 15 days of appearance in regular grid
# and discarding intra-specific contacts
ctdn <- ctdn %>% filter(!str_detect(ct_ID, "A"))
sp <- ctdn %>% group_by(sp, day) %>% count() %>% mutate(n = 1) %>% ungroup() %>%
  group_by(sp) %>% summarise(n = sum(n)) %>% filter(n >=15) %>% pull(sp)

contact.dn <- contact.dn.thres.10[contact.dn.thres.10$sp_depositor %in% sp &
                                    contact.dn.thres.10$sp_acquirer %in% sp &
                                    contact.dn.thres.10$sp_depositor !=
                                    contact.dn.thres.10$sp_acquirer,]

# Camera trapping operability
operability <- data.frame(t(read.csv("databases/Operability_dn.csv", sep = ";", 
                                     header = FALSE,check.names = FALSE)))
names(operability) <- operability[1,]
operability <- operability[-1,]
operability <- operability[, grepl("A", names(operability))] %>%
  mutate_if(is.character, as.numeric) %>%
  mutate(operability.a = rowSums(.,na.rm = TRUE)) %>% dplyr::select(operability.a) %>%
  bind_cols(operability)
operability <- operability[, !grepl("A", names(operability))] %>% dplyr::select(-1,-2) %>%
  mutate_if(is.character, as.numeric) %>%
  mutate(operability.g = rowSums(.,na.rm = TRUE)) %>% dplyr::select(operability.g) %>%
  bind_cols(operability) %>% rename(Date = 3) %>% 
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>% dplyr::select(1:3)

# Contacts by species pairs
dn <- data.frame(Method = character(),
                 Interact = character(),
                 direct = numeric(),
                 indirect = numeric(),
                 total = numeric())

for (j in c("", "A", "all")) {
  ct <- contact.dn %>%
    filter(if (j != "all") {gsub('[[:digit:]]+', '', ct_ID) == j} else TRUE) %>%
    mutate(Interact = paste0(sp_depositor, " - ", sp_acquirer)) %>%
    group_by(Interact) %>% 
    mutate(direct = sum(tdiff == 0)) %>% 
    mutate(indirect = sum(tdiff > 0)) %>%
    mutate(total = direct + indirect) %>%
    dplyr::select(Interact, direct, indirect, total) %>% 
    summarise(across(.cols = 1:3,.fns = max)) %>%
    replace(is.na(.), 0)
  
  ct$Method = ifelse(j == "", "Regular grid", ifelse(
    j == "A", "Resource monitoring", "All"))
  
  dn <- rbind(dn, ct)
}

spcomb <- crossing(sp1 = sp, sp2 = sp, Method =unique(dn$Method))
spcomb <- spcomb[spcomb$sp1 != spcomb$sp2,] %>% 
  mutate(Interact = paste0(sp1, " - ", sp2)) %>% dplyr::select(Interact, Method)
dn <- merge(dn, spcomb, all.y = TRUE) %>% replace(is.na(.), 0)

# 6.2. Doñana missing pairs ----
npairs <- dn %>% filter(Method != "All") %>% 
  mutate(n = case_when(direct == 0 ~ 0, direct != 0 ~1)) %>%
  group_by(Method) %>% summarise(n = sum(n)) %>% mutate(Type = "Direct")
npairs <- dn %>% filter(Method != "All") %>% 
  mutate(n = case_when(indirect == 0 ~ 0, indirect != 0 ~1)) %>%
  group_by(Method) %>% summarise(n = sum(n)) %>% 
  mutate(Type = "Indirect") %>%  bind_rows(npairs)

# Number of interacting pairs detected
# scales::show_col(scales::viridis_pal()(4))
# scales::viridis_pal()(4)
ggplot(npairs, aes(fill=Method, y=n, x=Type)) + 
  geom_bar(position="dodge", stat="identity") + theme_minimal() +
  theme(legend.title= element_blank(),
        axis.title.y = element_text(margin = margin( r = 10,))) +
  scale_y_continuous(breaks = seq(0, 60, by = 10)) +
  xlab("Type of interaction") + 
  ylab("Interacting pairs detected") +
  scale_fill_manual(values = c("#31688EFF", "#35B779FF"))


# 6.3. Number of interactions by camera and date ----
dn.daily <- data.frame(Method = character(),
                       Date = character(),
                       direct = numeric(),
                       indirect = numeric(),
                       total = numeric())

for (j in c("", "A", "all")) {
  ct <- contact.dn %>%
    filter(if (j != "all") {gsub('[[:digit:]]+', '', ct_ID) == j} else TRUE) %>%
    mutate(Date = as.Date(DateTime_dep)) %>%
    group_by(Date) %>% 
    mutate(direct = sum(tdiff == 0)) %>% 
    mutate(indirect = sum(tdiff > 0)) %>%
    mutate(total = direct + indirect) %>%
    dplyr::select(Date, direct, indirect, total) %>% 
    summarise(across(.cols = 1:3,.fns = max)) %>%
    replace(is.na(.), 0)
  
  ct$Method = ifelse(j == "", "Regular grid", ifelse(
    j == "A", "Resource monitoring", "All"))
  
  dn.daily <- rbind(dn.daily, ct)
}

# Days with less than 5 operating cameras by each method (first days) are removed
sel.operab <- operability %>% filter(operability.g >= 5, operability.a >= 5)
dn.daily <- merge(dn.daily, sel.operab)

# Number of interactions by camera
dn.daily.dir <- dn.daily %>% filter(Method != "All") %>%
  mutate(contact.cam = 
           case_when(Method == "Regular grid" ~ direct/operability.g,
                     Method == "Resource monitoring" ~ direct/operability.a),
         intype = "Direct interactions")

dn.daily.indir <- dn.daily %>% filter(Method != "All") %>%
  mutate(contact.cam = 
           case_when(Method == "Regular grid" ~ indirect/operability.g,
                     Method == "Resource monitoring" ~ indirect/operability.a),
         intype = "Indirect interactions")
dn.daily <- bind_rows(dn.daily.dir, dn.daily.indir)

# Plotting
ggplot(dn.daily, aes(color=Method, linetype = intype,
                     y=contact.cam, x=Date)) + 
  geom_line(stat="identity", linewidth = 1.1) +
  theme_minimal() +  
  theme(legend.title= element_blank(),
        axis.title.y = element_text(margin = margin( r = 10,))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c("#31688EFF", "#35B779FF")) + 
  ylab("Number of interactions detected by camera") +
  scale_y_log10(breaks=c(0, 1, 3, 10, 30, 100, 250, 500, 1000),
                labels = c("0","1", "3","10", "30", "100", "250","500","1000")) +
  scale_linetype_manual(values = c(22, "solid")) 

# 6.4. Clustering analysis ----
data.rg <- intpairs %>% filter(Method == "Regular grid") %>% 
  column_to_rownames(var = "Interact") %>% dplyr::select(total)
data.rm <- intpairs %>% filter(Method == "Resource monitoring") %>% 
  column_to_rownames(var = "Interact") %>% dplyr::select(total)

# Perform hierarchical clustering using ward.D linkage
hierm.rg <- hclust(dist(data.rg), method = "ward.D")
hierm.rm <- hclust(dist(data.rm), method = "ward.D")

# Plot the dendrogram
par(mfrow=c(1,2))
plot(hierm.rg, main = "Regular grid dendrogram")
rect.hclust(hierm.rg, h = sd(data.rg)*0.2, border = 2:10)
plot(hierm.rm, main = "Resource monitoring  dendrogram")
rect.hclust(hierm.rm, h = sd(data.rm)*0.2, border = 3:12)

library(dendextend)

colramp <- colorRampPalette(c("deepskyblue4",
                              "darkslategray3",
                              "darkseagreen3",
                              "olivedrab3",
                              "yellow3",
                              "gold1",
                              "sienna1",
                              "indianred3",
                              "indianred4"))
dend <- dendlist(
  as.dendrogram(hierm.rm) %>%
    rotate(order = intpairs %>% filter(Method == "Resource monitoring") %>% 
             arrange(total) %>% pull(Interact)) %>%
    set("branches_k_color", value = colramp(7), h = sd(data.rm$total)*0.5) %>%
    #    set("labels_col", value = colramp(10), h = sd(data.rm$total)*0.2) %>%
    set("leaves_pch", 19)  %>% 
    set("branches_lwd", 3) %>%
    set("nodes_cex", 0.7) %>% 
    set("labels_cex", 0.8) %>%
    set("rank_branches"),
  as.dendrogram(hierm.rg) %>%
    rotate(order = intpairs %>% filter(Method == "Regular grid") %>% 
             arrange(total) %>% pull(Interact)) %>%
    set("branches_k_color", value = colramp(4), h = sd(data.rg$total)*0.5) %>%
    #    set("labels_col", value = colramp(6), h = sd(data.rg$total)*0.2) %>%
    set("leaves_pch", 19)  %>% 
    set("branches_lwd", 3) %>%
    set("nodes_cex", 0.7) %>% 
    set("labels_cex", 0.8) %>%
    set("rank_branches")
)

cor.dendlist(dend)
tanglegram(dend, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = FALSE, 
           highlight_branches_lwd=TRUE, 
           margin_inner=7.8,
           lwd=2, axes = F,
           rank_branches = TRUE,
           edge.lwd = 2,
           columns_width = c(2,0.8,2),
           color_lines = get_leaves_branches_col(dend[[1]]),
           type = "r"
)


# Cutting height sensitivity analysis
# Coefficient in number of clusters
coef <- c()
for (i in 0.05*1:20){
  coef <- c(coef, max(cutree(hierm.rg, h = sd(data.rg$total)*i))/
              max(cutree(hierm.rm, h = sd(data.rm$total)*i)))
}

# mean and confidence interval 95%
mean(coef)
sd(coef)
1.96*sd(coef)/sqrt(length(coef))

