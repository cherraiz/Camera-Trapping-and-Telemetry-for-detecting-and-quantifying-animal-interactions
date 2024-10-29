species <- setClass("species", slots=c(HRkm2 = "numeric", meanT = "numeric", 
                                       sdT = "numeric", meanF = "numeric", 
                                       sdF = "numeric", rhoT = "numeric", 
                                       rhoF = "numeric", pTF = "numeric", 
                                       pFT = "numeric", theta = "numeric",
                                       energy = "numeric"))  

sim.traj.rech <- function (sp, env, n, init.coords, resolution, res.dist) {  
  
  # Reduce the environment to a random Home Range with at least one resource spot 
  environ = env[0]
  while (length(environ) == 0) {
    try({
      homerange = rand.HR(HRkm2 = sp@HRkm2, sides = 5, 
                          init.coords = as.numeric(init.coords), 
                          crs = terra::crs(env))
      area = terra::vect(st_as_sfc(sf::st_bbox(homerange)))
      environ = terra::mask(env, mask = terra::vect(homerange))
    }, silent = TRUE)
  }
  env = distance.rast(environ, area = area, resolution = resolution)
  env[env<res.dist] <- 0
  env = terra::mask(env, mask = terra::vect(homerange))

  init.coords = terra::spatSample(env, size = 1, xy = TRUE, na.rm = TRUE)[1:2]

  # Extracting parameters
  meanT=sp@meanT; sdT=sp@sdT; meanF=sp@meanF; sdF=sp@sdF; rhoT=sp@rhoT; 
  rhoF=sp@rhoF; pTF=sp@pTF; pFT=sp@pFT
  
  # Database with locations
  track = data.frame()
  track[1,1] = init.coords[[1]]  # Origin x coordinate
  track[1,2] = init.coords[[2]]  # Origin y coordinate
  track[1,3] = "T" # Starts in transit behaviour
  track[1,4] = 0   # Turning angle of the first step (0 rad by default)
  track[1,5] = 0   # Step length of the first step (0 m. by default)
  track[1,6] = sp@energy
  track[1,7] = 1
  colnames(track) = c("x", "y", "behavior", "angle", "distance", "energy", 
                      "step")
  
  for (step in 2:n) {
    
    track[step, 7] = step
    
    # Progress
    cat("\014") 
    cat(paste0(round(step / n, digits=3) * 100, 
               '% completed...')) 
    Sys.sleep(.001)
    if (step == n) cat(': Done')  
    
    lon_candidate = -9999  
    lat_candidate = -9999 
    resource.dist = as.numeric(terra::extract(env, 
                                              matrix(c(track[step-1,1],
                                                       track[step-1,2]),1,2)))
    
    # Recharge dynamics:
    # Recharge if at a resource (max energy = 100)
    # Else, loose energy
    if (resource.dist == 0 && track[step-1,3] == "D") {
      track[step, 6] = track[step-1, 6] + sp@theta[2]
    } else {
      track[step, 6] = track[step-1, 6] + sp@theta[1]
    }
    if (track[step,6] > 100) {
      track[step,6] = 100
    }
    
    # Behavior switching process:
    # If the animal was discharged but it recharged, changes to transit
    # If energy is low, changes to discharged
    # Else, transition probabilities between foraging and transit
    if (track[step-1,3] == "D") {
      if (track[step,6] > 90) {
        track[step,3] = "T"
      } else {
        track[step,3] = "D"
      }
    } else {
      if (track[step, 6] < 10) {
        track[step,3] = "D"
      } else {
        # Behavior transition probability between foraging and transit
        if (track[step-1,3] == "T"){
          if (runif(1, min = 0, max = 1) > pTF) {
            track[step,3] = "T"
          } else{
            track[step,3] = "F"
          }
        }
        
        if (track[step-1,3] == "F"){
          if (runif(1, min = 0, max = 1) > pFT) {
            track[step,3] = "F"
          } else{
            track[step,3] = "T"
          }
        }
      }
    }
    
    counter = 0
    # For Transit behavior
    if (track[step,3] == "T") {
      while ( is.na(terra::extract(env,matrix(c(lon_candidate,
                                                lat_candidate),1,2)))) {
        if(counter < 2) {
          angle = CircStats::rwrpcauchy(1, track[step-1,4], rhoT)
          dist = rgamma(1, shape = (meanT/sdT)^2, scale=sdT^2/meanT)
          track[step,4] = angle
          track[step,5] = dist
          new_coords = dagR::anglePoint(c(track[step-1,1], track[step-1,2]),
                                        angle, dist) 
          lon_candidate = new_coords[1]; lat_candidate<-new_coords[2]
          counter = counter +1 
        } else {
          angle = CircStats::rwrpcauchy(1, track[step-1,4], 0)
          dist = rgamma(1, shape = (meanT/sdT)^2, scale=sdT^2/meanT)
          track[step,4] = angle
          track[step,5] = dist
          new_coords = dagR::anglePoint(c(track[step-1,1], track[step-1,2]),
                                        angle, dist) 
          lon_candidate = new_coords[1]; lat_candidate<-new_coords[2]
          counter = counter +1 
        }
      }
    }
    
    # For Foraging behavior
    if (track[step,3] == "F") {
      while ( is.na(terra::extract(env, matrix(c(lon_candidate,
                                                 lat_candidate),1,2)))) {
        if (counter < 2) {
          angle = CircStats::rwrpcauchy(1, track[step-1,4], rhoF)
          dist = rgamma(1, shape = (meanF/sdF)^2, scale=sdF^2/meanF)
          track[step,4] = angle
          track[step,5] = dist
          new_coords = dagR::anglePoint(c(track[step-1,1], track[step-1,2]), 
                                        angle, dist)
          lon_candidate  = new_coords[1]; lat_candidate<-new_coords[2]
          counter = counter +1
        } else {
          angle = CircStats::rwrpcauchy(1, track[step-1,4], 0)
          dist = rgamma(1, shape = (meanF/sdF)^2, scale=sdF^2/meanF)
          track[step,4] = angle
          track[step,5] = dist
          new_coords = dagR::anglePoint(c(track[step-1,1], track[step-1,2]), 
                                        angle, dist)
          lon_candidate  = new_coords[1]; lat_candidate<-new_coords[2]
          counter = counter +1
        }
      } 
    }
    
    # For discharged behavior
    if (track[step,3] == "D") {
      # For selecting the direction to move
      # IDs of the adjacent cells
      neig = terra::adjacent(env,   
                             terra::cellFromXY(env, matrix(c(track[step-1,1],  
                                                             track[step-1,2]), 1,2)),   
                             directions=8, pairs=FALSE )  
      # Get difference between value and optimum for each cell
      options = data.frame()  
      for (i in 1:length(neig)){  
        options[i,1] = neig[i]  
        options[i,2] = 0 - env[neig[i]]  
      } 
      
      # Select better value (minimum difference)
      option = na.omit(options)[abs(na.omit(options$lyr.1)) == min(abs(
        na.omit(options$lyr.1))),1]  
      opt_coords = terra::xyFromCell(env,option)  
      
      # Moving towards resource
      while (is.na(terra::extract(env, matrix(c(lon_candidate,
                                                lat_candidate),1,2)))) { 
        dist = rgamma(1, shape = (meanT/sdT)^2, scale=sdT^2/meanT)
        if (dist >= resource.dist) {
          dist = resource.dist  # If distance higher than distance to water then arrives
        }
        angle = atan2(opt_coords[2]-track[step-1,2], 
                      opt_coords[1]- track[step-1,1])
        track[step,4] = angle
        track[step,5] = dist
        new_coords = c(track[step-1,1] + dist*cos(angle), 
                       track[step-1,2] + dist*sin(angle))
        lon_candidate  = new_coords[1]; lat_candidate<-new_coords[2]
      }  
    }
    
    track[step,1] = new_coords[1]  
    track[step,2] = new_coords[2]
  }  
  track$step <- as.numeric(track$step)
  return(track)  
}  

rand.HR <- function(sides, HRkm2, init.coords, crs) {
  # Mean radius
  area = HRkm2*10^6
  radius = sqrt((2*area)/(sides*sin((2*pi)/sides)))
  
  # Randomize the radii and angles
  radii = rnorm(sides, radius, radius/5)
  angles = runif(sides, min = 0, max = 2*pi)
  angles = sort(angles)
  
  # X and Y displacement
  points = list(x=NULL, y=NULL)
  points$x = cos(angles) * radii
  points$y = sin(angles) * radii
  
  # Current area of the polygon
  m = matrix(unlist(points), ncol=2)
  m = rbind(m, m[1,])
  current.area = 0.5 * (sum(m[1:sides,1]*m[2:(sides+1),2]) - 
                          sum(m[1:sides,2]*m[2:(sides+1),1]))
  
  # Area correction of the displacement and new coordinates
  points$x = init.coords[1] + cos(angles) * radii * sqrt(area/current.area)
  points$y = init.coords[2] + sin(angles) * radii * sqrt(area/current.area)
  
  # Transform to sf polygon and smoothing
  poly = data.frame(x=points$x, y=points$y) %>% 
    sf::st_as_sf(coords = c("x", "y"), crs = crs) %>% 
    summarise((geometry = sf::st_combine(geometry))) %>% 
    sf::st_cast("POLYGON") %>% 
    smoothr::smooth(method = "chaikin") %>% 
    sf:: st_make_valid() %>% 
    sf:: st_union()
  
  return(poly)
}

distance.rast <- function(vect, crs = terra::crs(area), area, resolution, ...) {
  r = terra::rast(area, crs = crs,
                  extent = terra::ext(area), resolution = resolution)
  dist = terra::distance(r, y = vect)
  return(dist)
}
