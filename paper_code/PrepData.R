library(tidyverse)
library(sf)
library(tmap)
library(maptools)
library(igraph)

epsg_wgs84 <- 4326 # GPS CRS (WGS 84), coordinate system
epsg_ma <- 2249 # NAD83, MA mainland

# [ get data and filter based on district and street type ]

boston <- st_read("./data/Zoning_Subdistricts/Zoning_Subdistricts.shp", quiet = TRUE) %>%
  st_transform(epsg_ma) %>%
  # filter out disconnected regions:
  filter(DISTRICT != "Allston/Brighton Neighborhood" &
           DISTRICT != "East Boston Neighborhood" &
           DISTRICT != "Boston Harbor" )

streets <- st_read("./data/Boston_Street_Segments/Boston_Segments.shp") %>%
  st_transform(epsg_ma) %>% filter(lengths(st_intersects(., boston)) > 0) %>%
  # remove bridges, highways, islands, etc:
  filter(ST_TYPE != "BRG" & ST_TYPE != "CSWY" & ST_TYPE != "DM" &
           ST_TYPE != "HWY" & ST_TYPE != "IS")

# [ generate primal and dual graphs ]
dist_coords <- function (x, y) sqrt(sum((x - y) ^ 2))

shape <- streets %>% dplyr::select(OBJECTID, SHAPElen, NBHD_L) 
seg_id <- shape$OBJECTID # use this to attach information to the edges (vertex information in the line graph)
seg_totaldist <- shape$SHAPElen

# get endpoints
n <- nrow(shape)
seg_from <- matrix(nrow = n, ncol = 2)
seg_to <- matrix(nrow = n, ncol = 2)
for (i in 1:n) {
  cc <- st_coordinates(shape[i, ])[, 1:2]
  seg_from[i,] <- cc[1, ]
  seg_to[i,] <- cc[nrow(cc), ]
}

# find distance threshold
seg_n <- numeric(n)
for (i in 1:n) seg_n[i] <- nrow(st_coordinates(shape[i, ])) #gets the coordinates along each street
m <- sum(seg_n) - n #number of coordinates for each street
seg_dist <- numeric(m)
k <- 1
for (i in 1:n) {
  cc <- st_coordinates(shape[i, ])[, 1:2]
  for (j in 2:seg_n[i]) {
    seg_dist[k] <- dist_coords(cc[j, ], cc[j - 1, ])
    k <- k + 1
  }
}
dist_thold <- min(seg_dist[seg_dist > 0])

# NOTE: cluster tolerance to integrate points is 0.001m ~ 1.57e-10 lat-lon dist
# Source: http://webhelp.esri.com/arcgisserver/9.3/Java/index.htm#geodatabases/topology_in_arcgis.htm

# remove loops
seg_enddist <- numeric(n)
for (i in 1:n) seg_enddist[i] <- dist_coords(seg_from[i, ], seg_to[i, ]) #total length of street
valid_idx <- which(seg_enddist >= dist_thold) #one street with length zero
seg_id <- seg_id[valid_idx]
seg_from <- seg_from[valid_idx, ]
seg_to <- seg_to[valid_idx, ]
seg_totaldist <- seg_totaldist[valid_idx]
n <- length(valid_idx)

# [ unique points ]
point_id <- c(1)
point_coords <- matrix(seg_from[1, ], ncol = 2)
k <- 2
seg_fromid <- numeric(n)
seg_toid <- numeric(n)
for (i in 1:n) {
  # from:
  d <- apply(point_coords, 1, dist_coords, seg_from[i, ])
  if (min(d) < dist_thold) # already in point list?
    idx <- point_id[which.min(d)]
  else { # new point
    point_id <- c(point_id, k)
    point_coords <- rbind(point_coords, seg_from[i, ])
    idx <- k
    k <- k + 1
  }
  seg_fromid[i] <- idx
  # to:
  d <- apply(point_coords, 1, dist_coords, seg_to[i, ])
  if (min(d) < dist_thold) # already in point list?
    idx <- point_id[which.min(d)]
  else { # new point
    point_id <- c(point_id, k)
    point_coords <- rbind(point_coords, seg_to[i, ])
    idx <- k
    k <- k + 1
  }
  seg_toid[i] <- idx
}

# optional
shape <- shape[valid_idx, ]
shape$fromid <- seg_fromid
shape$toid <- seg_toid


# [ primal graph ]
n <- length(point_id)
m <- length(seg_id)
g <- make_empty_graph(n = n, directed = FALSE) %>%
  add_edges(c(rbind(seg_fromid, seg_toid)))
V(g)$id <- point_id
V(g)$lon <- point_coords[, 1]
V(g)$lat <- point_coords[, 2]
E(g)$id <- seg_id
E(g)$dist <- seg_totaldist

cg <- igraph::components(g)
valid_idx <- which(cg$membership == 1)
g <- induced_subgraph(g, valid_idx)

lo <- layout.norm(as.matrix(cbind(V(g)$lon,V(g)$lat)))
plot(g, vertex.size = 2, vertex.label = NA, layout = lo) 


# [ (dual) line graph ]
nl <- ecount(g)
ml <- sum(choose(igraph::degree(g), 2))
gl_fromid <- integer(ml)
gl_toid <- integer(ml)
gl_inters <- integer(ml)
gl_dist <- numeric(ml)
k <- 1
for (i in 1:vcount(g)) {
  e <- incident(g, i)
  ne <- length(e)
  if (ne > 1) {
    for (j in 1:(ne - 1)) {
      for (l in (j + 1):ne) {
        gl_fromid[k] <- e[j]
        gl_toid[k] <- e[l]
        gl_inters[k] <- i
        gl_dist[k] <- (e[j]$dist + e[l]$dist) / 2.
        k <- k + 1
      }
    }
  }
}

gl <- make_empty_graph(n = nl, directed = FALSE) %>%
  add_edges(c(rbind(gl_fromid, gl_toid)))
V(gl)$id <- E(g)$id
V(gl)$lon <- (V(g)[head_of(g, E(g))]$lon + V(g)[tail_of(g, E(g))]$lon) / 2.
V(gl)$lat <- (V(g)[head_of(g, E(g))]$lat + V(g)[tail_of(g, E(g))]$lat) / 2.
E(gl)$id <- gl_inters
E(gl)$dist <- gl_dist

lo_gl <- layout.norm(as.matrix(cbind(V(gl)$lon,V(gl)$lat))) #at midpoint of street
plot(gl, vertex.size = 2, vertex.label = NA, layout = lo_gl) 

# [ Vertex Attribute Information ]

# [ Crime: Residential Burglary ]
# Source: https://data.boston.gov/dataset/crime-incident-reports-august-2015-to-date-source-new-system

crime1 <- read_csv("./data/CrimeReport.csv") %>%
  filter(OFFENSE_CODE == 520 | OFFENSE_CODE == 521 |
           OFFENSE_CODE == 522) %>% # residential burglaries
  filter(!is.na(Lat)) %>%
  st_as_sf(coords = c("Long", "Lat")) %>%
  st_set_crs(epsg_wgs84) %>% st_transform(epsg_ma) %>%
  filter(lengths(st_intersects(., boston)) > 0) %>%
  mutate(date = mdy_hm(OCCURRED_ON_DATE)) %>%
  select(-c("OCCURRED_ON_DATE")) 

# match formating of earlier crime report
crime1$INCIDENT_NUMBER <- str_remove(crime1$INCIDENT_NUMBER, pattern = "\\-.*")
crime1$INCIDENT_NUMBER<- as.double(str_sub(crime1$INCIDENT_NUMBER, 2, -1L))

crime2 <- read_csv("./data/CrimeReport2.csv") %>%
  filter(OFFENSE_CODE == 520 | OFFENSE_CODE == 521 |
           OFFENSE_CODE == 522) %>% # residential burglaries
  filter(!is.na(Lat)) %>%
  st_as_sf(coords = c("Long", "Lat")) %>%
  st_set_crs(epsg_wgs84) %>% st_transform(epsg_ma) %>%
  filter(lengths(st_intersects(., boston)) > 0)%>%
  mutate(date = mdy_hm(OCCURRED_ON_DATE)) %>%
  select(-c("OCCURRED_ON_DATE"))

crime <- rbind(crime1,crime2) %>%
  distinct(INCIDENT_NUMBER, .keep_all = TRUE) %>%
  arrange(date)

remove(crime1, crime2)

# [ Crime to Street ]
crime_st_id <- st_nearest_feature(crime, shape)
t <- table(crime_st_id)
Y <- integer(nrow(shape))
Y[as.integer(names(t))] <- t

seg <- data.frame(seg_id, seg_totaldist, Y) #collect all information in seg dataframe 

shape <- shape %>%
  mutate(crime_count = Y)

# [ Wealth Tax Proxy ]
# based on residential parcel taxes and current set of intersections

# [ Parcels, set buffer of 10 ]
parcels <- st_read("./data/Parcels2016Data/output.shp") %>% st_transform(epsg_ma) %>%
  filter(lengths(st_intersects(., boston)) > 0) %>% st_buffer(10) 

res_parcels <- parcels %>% # residential parcels
  filter(PTYPE < 200 | PTYPE == 903 | PTYPE == 907 | PTYPE == 908) %>% # dorms
  st_buffer(10)

coll_parcels <- parcels %>% # college parcels
  filter(PTYPE == 123 | PTYPE == 122 | PTYPE == 977) %>%
  st_buffer(10)

STREET_BUFFER <- 30 

street_buffer <- shape %>% 
  st_buffer(STREET_BUFFER) #buffer around street

street_gross_tax <- st_intersection(res_parcels, street_buffer) %>%
  left_join(as_tibble(res_parcels %>% select(OBJECTID, geometry)), by = "OBJECTID") %>%
  mutate(gross_tax = as.numeric(st_area(geometry.x) / st_area(geometry.y) * AV_TOTAL)) %>%
  select(OBJECTID.1, gross_tax) %>%
  rename(OBJECTID = OBJECTID.1) %>%
  group_by(OBJECTID) %>% summarize(gross_tax = sum(gross_tax))

shape <- shape %>%
  left_join(as_tibble(street_gross_tax %>% select(OBJECTID, gross_tax))) %>%
  select(-geometry.x) %>% # remove polygons
  mutate(gross_tax = ifelse(is.na(gross_tax), 0, gross_tax))

seg$gross_tax <- shape$gross_tax

# [ Paper Figure 1]

#Visualizations from the LINE GRAPH perspective
vert <- data.frame(V(gl)$id)
names(vert) <- c("seg_id")
map <- left_join(vert, seg)

V(gl)$crime <- map$Y
V(gl)$gross_tax <- map$gross_tax

# Using map tools (export these to power point to edit)
#wealth
my_pal1 <- c('gray80',"#C7E9C0","#74C476" ,"#238B45","#00441B")
p_w <- tm_shape(shape) + tm_lines(col = "gross_tax", palette = my_pal1, size = .025,
                                  style = "fixed", border.alpha = .5, title.col = "Wealth Proxy",
                                  breaks = c(0, 100000, 600000, 1500000,  10000000, 500000000),
                                  labels = c("0-100", "100-600", "600-1500", "1500-10000", "10000+"), lwd= 3) +
  tm_legend( legend.text.size = 2, legend.title.size = 3, legend.width = 10) +
  tm_layout(frame = FALSE)


#crime
my_pal2 <- c('gray80','yellow','darkorange','red')
shape2 <- shape
names(shape2)[5] <- "Crime Count"

p_c <- tm_shape(shape2) + 
  tm_lines(col = "Crime Count", palette = my_pal2,  breaks = c(0, 1,2,5,30),labels = c("0", "1", "2-5", "5-30"), lwd = 3) +
  tm_legend(title = "Crime Count", legend.text.size = 2, legend.title.size = 3, legend.width = 10) +
  tm_layout(frame = FALSE)

# [ Paper Figure 2]

df <- data.frame(shape$crime_count, shape$gross_tax, shape$NBHD_L)
df$shape.NBHD_L[df$shape.NBHD_L == "na"] <- c("Boston")

df <- df %>%
  filter(shape.gross_tax < 6272010) 

g <- ggplot(df, aes(x = shape.gross_tax, y = shape.crime_count, color = shape.NBHD_L))+
  geom_jitter(size = 7, alpha = 0.7)+
  labs(y = "Residential Burglary", x = "Wealth Proxy", color = "Neighborhood") + 
  theme(
    axis.title = element_text(size = 20),  
    legend.text = element_text(size = 16)  
  )


# other preds to build dataframe for regression

# [ Police Distance ]
police <- st_read("./data/Boston_Police_Stations/Boston_Police_Stations.shp") %>%
  st_transform(epsg_ma) %>% filter(lengths(st_intersects(., boston)) > 0)

int_near_police <- st_nearest_feature(shape, police)

shape <- shape %>%
  mutate(pol_dist = st_distance(., police[int_near_police, ],
                                by_element = TRUE))

seg$pol_dist <- shape$pol_dist

# [ Zoning Subdistricts ]
sub_dist <- st_read("./data/Zoning_Subdistricts/Zoning_Subdistricts.shp") %>%
  st_transform(epsg_ma) %>% filter(lengths(st_intersects(., boston)) > 0)

shape2 <- shape %>%
  st_join(sub_dist %>% dplyr::select(SUBDISTRIC)) %>% 
  st_difference()

shape <- shape2 %>%
  distinct(OBJECTID, .keep_all = TRUE)

seg$district <- shape$SUBDISTRIC

# [ Population ]

# Sources:
# https://catalog.data.gov/dataset/tiger-line-shapefile-2016-state-massachusetts-current-census-tract-state-based
# https://docs.digital.mass.gov/dataset/massgis-data-2010-us-census-environmental-justice-populations
# https://docs.digital.mass.gov/dataset/massgis-data-datalayers-2010-us-census

# Note that the "POP100_RE" fields store the 100% population count as listed in
# the PL-94-171 demographics tables (the 'RE' in the field name stands for
# 'redistricting').
# The housing unit counts ("HU100_RE") were generated from the same PL-94-171 tables.
population <- st_read("./data/Population/CENSUS2010BLOCKGROUPS_POLY.shp") %>%
  st_transform(epsg_ma) %>% filter(lengths(st_intersects(., boston)) > 0)

street_pop <- population %>% st_intersection(street_buffer) %>%
  left_join(as_tibble(population %>% dplyr::select(GEOID10, geometry)), by = "GEOID10") %>%
  mutate(POP100_RE = as.numeric(st_area(geometry.x) / st_area(geometry.y) * POP100_RE)) %>%
  dplyr::select(OBJECTID, POP100_RE) %>%
  group_by(OBJECTID) %>% summarize(POP100_RE = sum(POP100_RE))

shape <- shape %>%
  left_join(as_tibble(street_pop) %>% dplyr::select(OBJECTID, POP100_RE))

seg$pop <- shape$POP100_RE

# [ College Distance ]
COLLEGE_BUFFER <- units::set_units(x = 300, value = m) # based on definition of school zone
street_near_coll <- st_nearest_feature(street_buffer, coll_parcels)

street_buffer <- street_buffer %>%
  mutate(col_dist = st_distance(., coll_parcels[street_near_coll, ],
                                by_element = TRUE),
         col_b = col_dist < COLLEGE_BUFFER)

shape$col_b <- street_buffer$col_b

seg$col_b <- street_buffer$col_b

# [ Residential Buffer: Which streets should we marginalize out of the network?] #
keep <- street_buffer %>%
  filter(lengths(st_intersects(., res_parcels)) > 0)

# [ Data is ready - lets marginalize out streets we don't want]

phi <- similarity_phi(dist = E(gl)$dist,
                      median_similarity = 0.8) 

sim <- exp(-E(gl)$dist /
             (max(E(gl)$dist) * phi))

#marginalize out nodes not in residential area:
rem <- symdiff(V(gl)$id, keep$OBJECTID) # These are the ids, not the location
rem2 <- which(V(gl)$id %in% rem == TRUE)

mL <- laplacian_marginal(gl, rem2, weights = sim)
e.mL <- eigen(mL, symmetric = TRUE)

m <- nrow(e.mL$vectors)
e.values <- e.mL$values[m:1] # reversing eigenvalues
e.vectors <- e.mL$vectors[, m:1] # reversing eigenvectors


# covariates
to_keep <- data.frame(intersect(V(gl)$id, keep$OBJECTID)) %>%
  rename(seg_id = "intersect.V.gl..id..keep.OBJECTID.")

covars <- seg %>% 
  right_join(to_keep) %>%
  rename(crime = Y)

y <- covars$crime

covars$pol_dist <- as.numeric(covars$pol_dist)
preds <- data.frame(intercept = 1,
                    log_tax = log(covars$gross_tax + 1),
                    log_police_dist = log(covars$pol_dist + 1),
                    log_population = log(covars$pop + 1),
                    college = covars$col_b)

# standardizing tax, police, population:
preds1 <- preds
for (i in 2:4) preds1[, i] <- scale(preds1[, i])
preds <- as.matrix(preds1)

# imputing one population value [population and street buffer did not intersect for one point]
preds$log_population[6468] <- 2.4

