# Packages ----
rm(list = ls())
library(openxlsx) # not used in arcticts.qmd
library(dplyr)
library(imputeTS)
library(zoo)
library(dtw)
library(dtwclust) # not used in arcticts.qmd
source("/Users/spencerblackwood/Documents/reu_project/code/fun_clucomp.R")
library(funtimes) # not used in arcticts.qmd
library(ggplot2)
library(dbscan)
library(leaflet)
library(mapview)
library(plot.matrix) # not used in arcticts.qmd



# Data ----


## Import ----
data <- read.xlsx('/Users/spencerblackwood/Documents/reu_project/dataraw/BenthicTimeSeries_UniqueID_0-5nm_DBO_063023_v3.xlsx', 3, cols = c(1:10, 12:13, 16, 24, 34)) %>%
  rename(Year = DataYear, Date = DataDate)
data$Year <- as.character(data$Year)
data$Date <- as.Date(data$Date, '%Y%m%d')
data$StationName_Standard <- gsub(' ', '', data$StationName_Standard)
data$DBOnme <- gsub(' ', '', data$DBOnme)
# check for 1-to-1 relationship between StationName_Standard and DBOnme
tmp1 <- data %>%
  group_by(DBOnme) %>%
  summarise(name_count = length(unique(StationName_Standard))) %>%
  filter(name_count > 1)
tmp2 <- data %>%
  group_by(StationName_Standard) %>%
  summarise(name_count = length(unique(DBOnme))) %>%
  filter(name_count > 1)
dim(tmp1)[1] == 0 & dim(tmp2)[1] == 0 # should be true!


## Narrow Pt. 1 ----
# filter out data outside of July 2010-2019 and data from DBO5
data_full <- data %>%
  mutate(Month = format(Date, '%m')) %>%
  filter(Month == '07' & Year %in% paste0("20", as.character(10:19))) %>%
  filter(!grepl('DBO5', DBOnme))


## Aggregate ----
data_avg <- data_full %>%
  group_by(DBOnme, Year) %>%
  summarise(Temp_Avg = mean(Temp, na.rm = TRUE),
            Salinity_Avg = mean(Salinity, na.rm = TRUE),
            SedChla_Avg = mean(sedchla, na.rm = TRUE),
            BiomGC_Avg = mean(BiomGC, na.rm = TRUE),
            PhiGTE5_Avg = mean(phigte5, na.rm = TRUE))
Vars <- 3:ncol(data_avg)
for (varnum in Vars) {
  data_avg[[varnum]][is.nan(data_avg[[varnum]])] <- NA
}


## Narrow Pt. 2 ----
# remove stations with fewer than 'thr' observations
thr <- 5
stlst <- character()
Stations <- sort(unique(data_avg$DBOnme))
for (st in Stations) {
  tmp1 <- data_avg[data_avg$DBOnme == st, Vars]
  tmp2 <- apply(tmp1, 2, function(x) sum(!is.na(x)))
  if (all(tmp2 >= thr)) {
    stlst <- c(stlst, st)
  }
}
data_avg <- data_avg %>% filter(DBOnme %in% stlst)


## Widen ----
Times <- sort(unique(data_avg$Year))
for (varnum in Vars) {
  
  # create a table with only non-missing observations of the variable 'varnum'
  tmp1 <- data_avg[!is.na(data_avg[[varnum]]), c(1:2, varnum)]
  
  # create the wide matrix with values of the variable 'varnum'
  matx <- matrix(NA, nrow = length(Times), ncol = length(stlst), dimnames = list(Times, stlst))
  for (st in stlst) {
    tmp2 <- tmp1 %>% filter(DBOnme == st)
    matx[tmp2$Year, st] <- tmp2[[3]]
  }
  
  # rename
  assign(paste0('mat', varnum - 2), matx)
}


## Impute ----
for (i in 1:5) {
  matx <- get(paste0('mat', i))
  
  # pattern 1: identify & impute
  impnames1 <- colnames(matx)[apply(matx, 2, function(x) all(is.na(x[2])) & all(!is.na(x[c(1, 3:10)])))]
  for (st in impnames1) {
    matx[, st] <- na_kalman(matx[, st], model = 'StructTS')
  }
  
  # pattern 2: identify & impute
  impnames2 <- colnames(matx)[apply(matx, 2, function(x) all(is.na(x[c(1, 3:4)])) & all(!is.na(x[c(2, 5:10)])))]
  for (st in impnames2) {
    matx[2:10, st] <- na_kalman(matx[2:10, st], model = 'StructTS')
  }
  
  # replace
  assign(paste0('mat', i), matx)
}



# DTW ----


## Window Function ----
sakoeChibaWindowSB <- function (iw, jw, window.size, query.size, reference.size, ...) {
  if (query.size - reference.size == 1) {       # if query and reference have same length
    return(abs(jw + 1 - iw) <= window.size)
  }
  else if (query.size - reference.size == 0) {  # if query is 1 shorter than reference
    return(abs(jw - iw) <= window.size)
  }
  else if (query.size - reference.size == -3) { # if query is 4 shorter than reference
    return(abs(jw - 1 - iw) <= window.size)
  }
  else {
    return(FALSE)
  }
}


## Dissimilarity Matrices (Normalized) ----
for (varnum in 1:5) {
  
  # create an empty matrix
  dx <- matrix(NA, nrow = length(stlst), ncol = length(stlst), dimnames = list(stlst, stlst))
  matx <- get(paste0('mat', varnum))
  
  # build the upper triangle
  for (i in 1:length(stlst)) {
    for (j in i:length(stlst)) {
      
      # assign the shorter series to be the query
      tsq <- scale(na.omit(zoo(matx[, i], rownames(matx))), scale = FALSE) # subtract the mean
      tsr <- scale(na.omit(zoo(matx[, j], rownames(matx))), scale = FALSE) # ""
      if (length(tsq) > length(tsr)) {
        tmp <- tsr
        tsr <- tsq
        tsq <- tmp
      }
      
      # for two series of equal length, record average normalized dissimilarity
      if (length(tsq) == length(tsr)) {
        dx[i, j] <- (dtw::dtw(tsq, tsr, step.pattern = asymmetric,
                              window.type = sakoeChibaWindowSB, window.size = 2,
                              keep.internals = TRUE, open.end = TRUE,
                              open.begin = TRUE)$normalizedDistance +
                       dtw::dtw(tsr, tsq, step.pattern = asymmetric,
                                window.type = sakoeChibaWindowSB, window.size = 2,
                                keep.internals = TRUE, open.end = TRUE,
                                open.begin = TRUE)$normalizedDistance) / 2
      }
      
      # for two series of unequal length
      else {
        
        # case 1: shorter series has length 9, missing observation at time point 1
        if (length(tsr) - length(tsq) == 1) {
          dx[i, j] <- dtw::dtw(tsq, tsr, step.pattern = asymmetric,
                               window.type = sakoeChibaWindowSB, window.size = 2,
                               keep.internals = TRUE, open.end = TRUE,
                               open.begin = TRUE)$normalizedDistance
        }
        
        # case 2: shorter series has length 5, missing observations at time points 1-3, 9-10 
        else {
          tsr <- tsr[paste0("20", as.character(11:19))] # remove 2010 observation from reference
          dx[i, j] <- dtw::dtw(tsq, tsr, step.pattern = asymmetric,
                               window.type = sakoeChibaWindowSB, window.size = 2,
                               keep.internals = TRUE, open.end = TRUE,
                               open.begin = TRUE)$normalizedDistance
        }
      }
    }
  }
  
  # fill in the matrix and rename
  dx[lower.tri(dx)] <- t(dx)[lower.tri(dx)]
  assign(paste0('d', varnum), dx)
  assign(paste0('dist', varnum), as.dist(dx))
}



# Clustering ----

results <- list()
for (varnum in 1:5) {
  
  distx <- get(paste0('dist', varnum))
  
  
  ## Agglomerative and Centroid-Based
  
  # obtain optimal 'k' and method
  # ideally would use Davies-Bouldin Index instead of WB Ratio
  clucomp_results <- clucomp(distx, k = 2:6,
                             methods = c('complete', 'PAM', 'single', 'ward.D2'), 
                             criteria = c('avg.silwidth', 'ch', 'wb.ratio'))
  opt_k <- clucomp_results$opt_k
  opt_method <- clucomp_results$opt_method
  
  # store and print resulting clusters
  if (opt_method == 'PAM') {
    opt_cl <- cluster::pam(distx, k = opt_k, diss = TRUE)
    opt_cl_labels <- opt_cl$clustering
  } else {
    opt_cl <- hclust(distx, method = opt_method)
    opt_cl_labels <- cutree(opt_cl, k = opt_k)
  }
  clusters <- lapply(1:opt_k, function(x) names(opt_cl_labels[opt_cl_labels == x]))
  
  
  ## DBSCAN
  
  # identify optimal 'eps' given minPts = 2
  (knndist <- sort(kNNdist(distx, k = 1)))
  kNNdistplot(distx, minPts = 2)
  # results from visual inspection:
  epsvec <- c(0.3, 0.085, 2, 2.5, 1.8)
  
  # store and print resulting clusters
  dbs <- dbscan(distx, eps = epsvec[varnum], minPts = 2)$cluster
  clusters_dbs <- lapply(1:max(dbs), function(x) names(distx)[dbs == x])
  clusters_dbs_noise <- names(distx)[dbs == 0] # noise
  
  
  ## Final Results
  clresults <- list(opt_k = opt_k,
                    opt_method = opt_method,
                    clusters = clusters,
                    clusters_dbs = clusters_dbs,
                    clusters_dbs_noise = clusters_dbs_noise)
  results <- c(results, list(clresults))
}



# Appendix ----


## Plots ----

### Map of Stations
# coordinates of DBO regions 1-5
rect_coords <- data.frame(lat = c(63.769, 61.946, 61.847, 63.663, NA,
                                  65.111, 65.097, 64.482, 64.496, NA,
                                  66.786, 68.609, 68.572, 66.752, NA,
                                  71.867, 71.721, 70.682, 70.820),
                          lng = c(-172.312, -172.187, -175.791, -176.147, NA,
                                  -170.492, -167.860, -167.908, -170.481, NA,
                                  -171.330, -171.419, -166.481, -166.756, NA,
                                  -164.255, -160.507, -160.994, -164.553))
coords <- data_full[, c('DBOnme', 'Latitude', 'Longitude')] %>%
  filter(DBOnme %in% data_avg$DBOnme)
map <- leaflet(coords) %>%
  addTiles() %>%
  addCircleMarkers(radius = 3, stroke = TRUE, weight = 1, color = '#000', fillColor = '#FFF', fillOpacity = 0.5) %>%
  addPolygons(lng = rect_coords$lng, lat = rect_coords$lat, color = '#444444',
              weight = 2, opacity = 1, fillOpacity = 0)
mapshot(map, file = paste0('/Users/spencerblackwood/Documents/reu_project/output/map.png'))

### Maps of Clusters
varnames <- gsub('_Avg', '', colnames(data_avg)[3:ncol(data_avg)])
colors <- c('#FCE621', '#36BFDE', '#FF6154', '#473D9C', '#61C430', '#C252B8')
for (i in 1:length(results)) {
  clvar <- results[[i]]
  
  # globular clusters
  coordstmp <- coords %>% mutate(Color = NA)
  for (cl in 1:length(clvar['clusters'][[1]])) {
    clx <- clvar['clusters'][[1]][[cl]]
    coordstmp[coordstmp$DBOnme %in% clx, 'Color'] <- colors[cl]
  }
  coordstmp <- coordstmp %>% filter(!is.na(Color))
  mapx <- leaflet(coordstmp) %>%
    addTiles() %>%
    addCircleMarkers(radius = 4, stroke = TRUE, weight = 1, color = '#000', fillColor = coordstmp$Color, fillOpacity = 0.5)
  mapshot(mapx, file = paste0('/Users/spencerblackwood/Documents/reu_project/output/map', varnames[i], '.png'))
  
  # density-based clusters
  coordstmp <- coords %>% mutate(ColorD = NA)
  for (cl in 1:length(clvar['clusters_dbs'][[1]])) {
    clx <- clvar['clusters_dbs'][[1]][[cl]]
    coordstmp[coordstmp$DBOnme %in% clx, 'ColorD'] <- colors[cl]
  }
  coordstmp[coordstmp$DBOnme %in% clvar['clusters_dbs_noise'][[1]], 'ColorD'] <- '#000'
  coordstmp <- coordstmp %>% filter(!is.na(ColorD))
  mapxd <- leaflet(coordstmp) %>%
    addTiles() %>%
    addCircleMarkers(radius = 4, stroke = TRUE, weight = 1, color = '#000', fillColor = coordstmp$ColorD, fillOpacity = 0.5)
  mapshot(mapxd, file = paste0('/Users/spencerblackwood/Documents/reu_project/output/map', varnames[i], 'DBS.png'))
  
  assign(paste0('map', i), mapx)
  assign(paste0('map', i, 'd'), mapxd)
}

### Plots of Time Series by Cluster

axnames <- c('Temperature (ºC)', 'Salinity', 'Chlorophyll (mg/m^2)',
             'Biomass (gC/m^2)', '% Silt & Clay')
# globular cluster plots
for (varnum in 1:length(results)) {
  clvar <- results[[varnum]]
  coordstmp <- coords %>% mutate(Color = NA)
  for (cl in 1:length(clvar['clusters'][[1]])) {
    clx <- clvar['clusters'][[1]][[cl]]
    coordstmp[coordstmp$DBOnme %in% clx, 'Color'] <- colors[cl]
  }
  coordstmp <- coordstmp %>% filter(!is.na(Color))
  for (cl in 1:length(clvar['clusters'][[1]])) {
    clx <- clvar['clusters'][[1]][[cl]]
    matx <- get(paste0('mat', varnum))
    minval <- min(matx, na.rm = T)
    maxval <- max(matx, na.rm = T)
    png(file = paste0('/Users/spencerblackwood/Documents/reu_project/output/plot', varnum, '_', cl, '.png'),
        width = 960, height = 960, pointsize = 24)
    par(mar = c(5, 5.3, 2, 1.7) + 0.1)
    plot(2010:2019, matx[,clx[1]], type = 'b', col = colors[cl], lwd = 4,
         xlab = 'Year', ylab = '', ylim = c(minval, maxval), cex.axis = 1.5,
         cex.lab = 1.5)
    title(ylab = axnames[varnum], line = 3.7, cex.lab = 1.5)
    if (length(clx) > 1) {
      for (i in 2:length(clx)) {
        lines(2010:2019, matx[,clx[i]], type = 'b', col = colors[cl], lwd = 4)
      }
    }
    dev.off()
  }
}

# DBSCAN cluster plots
for (varnum in 1:length(results)) {
  clvar <- results[[varnum]]
  coordstmp <- coords %>% mutate(Color = NA)
  for (cl in 1:length(clvar['clusters_dbs'][[1]])) {
    clx <- clvar['clusters_dbs'][[1]][[cl]]
    coordstmp[coordstmp$DBOnme %in% clx, 'Color'] <- colors[cl]
  }
  coordstmp[coordstmp$DBOnme %in% clvar['clusters_dbs_noise'][[1]], 'Color'] <- '#000'
  coordstmp <- coordstmp %>% filter(!is.na(Color))
  for (cl in 1:length(clvar['clusters_dbs'][[1]])) {
    clx <- clvar['clusters_dbs'][[1]][[cl]]
    matx <- get(paste0('mat', varnum))
    minval <- min(matx, na.rm = T)
    maxval <- max(matx, na.rm = T)
    png(file = paste0('/Users/spencerblackwood/Documents/reu_project/output/plotD', varnum, '_', cl, '.png'),
        width = 960, height = 960, pointsize = 24)
    par(mar = c(5, 5.3, 2, 1.7) + 0.1)
    plot(2010:2019, matx[,clx[1]], type = 'b', col = colors[cl], lwd = 4,
         xlab = 'Year', ylab = '', ylim = c(minval, maxval), cex.axis = 1.5,
         cex.lab = 1.5)
    title(ylab = axnames[varnum], line = 3.7, cex.lab = 1.5)
    for (i in 2:length(clx)) {
      lines(2010:2019, matx[,clx[i]], type = 'b', col = colors[cl], lwd = 4)
    }
    dev.off()
  }
  clx <- clvar['clusters_dbs_noise'][[1]]
  png(file = paste0('/Users/spencerblackwood/Documents/reu_project/output/plotD', varnum, '_n.png'),
      width = 960, height = 960, pointsize = 24)
  par(mar = c(5, 5.3, 2, 1.7) + 0.1)
  plot(2010:2019, matx[,clx[1]], type = 'b', col = "#000000", lwd = 4,
       xlab = 'Year', ylab = '', ylim = c(minval, maxval), cex.axis = 1.5,
       cex.lab = 1.5)
  title(ylab = axnames[varnum], line = 3.7, cex.lab = 1.5)
  for (i in 2:length(clx)) {
    lines(2010:2019, matx[,clx[i]], type = 'b', col = "#000000", lwd = 4)
  }
  dev.off()
}


# Plots for OSM Poster
colors <- c('#ffff33', '#984ea3', '#ff7f00', '#4daf4a', '#377eb8', '#e41a1c')
## Map
clvar <- results[[1]]
coordstmp <- coords %>% mutate(Color = NA)
for (cl in 1:length(clvar['clusters'][[1]])) {
  clx <- clvar['clusters'][[1]][[cl]]
  coordstmp[coordstmp$DBOnme %in% clx, 'Color'] <- colors[cl]
}
coordstmp <- coordstmp %>% filter(!is.na(Color))
mapx <- leaflet(coordstmp) %>%
  addTiles() %>%
  addCircleMarkers(radius = 4, stroke = TRUE, weight = 1, color = '#000', fillColor = coordstmp$Color, fillOpacity = 0.5) %>%
  addPolygons(lng = rect_coords$lng, lat = rect_coords$lat, color = '#444444',
              weight = 2, opacity = 1, fillOpacity = 0)
mapshot(mapx, file = paste0('/Users/spencerblackwood/Documents/reu_project/output/poster/mapTemperature.png'))
## Time Series Plots
clvar <- results[[1]]
coordstmp <- coords %>% mutate(Color = NA)
for (cl in 1:length(clvar['clusters'][[1]])) {
  clx <- clvar['clusters'][[1]][[cl]]
  coordstmp[coordstmp$DBOnme %in% clx, 'Color'] <- colors[cl]
}
coordstmp <- coordstmp %>% filter(!is.na(Color))
for (cl in 1:length(clvar['clusters'][[1]])) {
  clx <- clvar['clusters'][[1]][[cl]]
  matx <- mat1
  minval <- min(matx, na.rm = T)
  maxval <- max(matx, na.rm = T)
  png(file = paste0('/Users/spencerblackwood/Documents/reu_project/output/poster/plot1_', cl, '.png'),
      width = 960, height = 960, pointsize = 24)
  par(mar = c(5, 5.3, 2, 1.7) + 0.1, bg = '#aad3df')
  plot(2010:2019, matx[,clx[1]], type = 'b', col = colors[cl], lwd = 4,
       xlab = 'Year', ylab = '', ylim = c(minval, maxval), cex.axis = 1.5,
       cex.lab = 1.5)
  title(ylab = 'Bottom Water Temperature (ºC)', line = 3.7, cex.lab = 1.5)
  if (length(clx) > 1) {
    for (i in 2:length(clx)) {
      lines(2010:2019, matx[,clx[i]], type = 'b', col = colors[cl], lwd = 4)
    }
  }
  dev.off()
}



### DTW plot example

st1 <- 'DBO3.2'
st2 <- 'DBO3.4'
tsq <- scale(na.omit(zoo(mat1[, st1], rownames(mat1))), scale = FALSE) # subtract the mean
tsr <- scale(na.omit(zoo(mat1[, st2], rownames(mat1))), scale = FALSE) # ""
tmp <- dtw::dtw(tsq, tsr, step.pattern = asymmetric,
  window.type = sakoeChibaWindowSB, window.size = 2,
                keep.internals = TRUE, open.end = TRUE,
                open.begin = TRUE)
png(file = '/Users/spencerblackwood/Documents/reu_project/output/plotdtwex.png')
par(mar = c(5, 5, 2, 2) + 0.1)
dtwPlotTwoWay(tmp, lwd = 1.5, match.col = '#444444', cex.lab = 1.2)
legend('topleft', legend=c(st1, st2), col = 1:2, lty = 1:2)
dev.off()




### Dissimilarity Plot
par(mar = c(5, 4, 4, 4) + 0.1)
plot(d5, las = 2, ylab = '', 
     main = 'DTW Distance Matrix (Grain Size)', xlab = '', border = NA)


## Imputation Test ----

if (FALSE) {
  ## FINDING THE BEST IMPUTATION METHOD
  ## RESULT: Imputation by Structural Model & Kalman Smoothing
  
  # Validation Function
  test_impute <- function(vec, ind) {
    vecna <- vec
    vecna[ind] <- NA
    tmp1 <- na_interpolation(vecna, option = 'linear')
    e1 <- vec[ind] - tmp1[ind] # real minus imputed
    tmp2 <- na_interpolation(vecna, option = 'spline')
    e2 <- vec[ind] - tmp2[ind]
    tmp3 <- na_kalman(vecna, model = 'StructTS')
    e3 <- vec[ind] - tmp3[ind]
    tmp4 <- na_kalman(vecna, model = 'auto.arima')
    e4 <- vec[ind] - tmp4[ind]
    return (c(e1, e2, e3, e4))
  }
  apply(tmat, 2, function(x) test_impute(x, ind = 2))
  
  # Test 1: filling in the second value in series of 10
  # prepare the matrix
  tmat1 <- mat1[,c(1:9, 11, 13:14)]
  tmat2 <- mat2[,c(1:9, 11, 13:14)]
  tmat3 <- mat3[,c(1:9,13)]
  tmat4 <- mat4[,1:9]
  tmat5 <- mat5[,c(1:9,13)]
  tmat <- cbind(tmat1, tmat2, tmat3, tmat4, tmat5)
  # results
  results1 <- apply(tmat, 2, function(x) test_impute(x, ind = 2))
  results1.1 <- apply(results1[, 1:12], 1, function(x) mean(x^2))
  results1.2 <- apply(results1[, 13:24], 1, function(x) mean(x^2))
  results1.3 <- apply(results1[, 25:34], 1, function(x) mean(x^2))
  results1.4 <- apply(results1[, 35:43], 1, function(x) mean(x^2))
  results1.5 <- apply(results1[, 44:53], 1, function(x) mean(x^2))
  test1 <- rbind(results1.1, # 3
                 results1.2, # 3
                 results1.3, # 3
                 results1.4, # 3
                 results1.5) # 3
  apply(test1, 1, which.min) # results as a comment directly above
  
  # Test 2: filling in the second and third values in series of 9
  # prepare the matrix
  tmat1 <- mat1[2:10,1:16]
  tmat2 <- mat2[2:10,1:16]
  tmat3 <- mat3[2:10,1:16]
  tmat4 <- mat4[2:10,c(1:9,14:16)]
  tmat5 <- mat5[2:10,1:16]
  tmat <- cbind(tmat1, tmat2, tmat3, tmat4, tmat5)
  # results
  results2 <- apply(tmat, 2, function(x) test_impute(x, ind = 2:3))
  results2.1.3 <- apply(results2[c(1, 3, 5, 7), 1:16], 1, function(x) mean(x^2))
  results2.1.4 <- apply(results2[c(2, 4, 6, 8), 1:16], 1, function(x) mean(x^2))
  results2.2.3 <- apply(results2[c(1, 3, 5, 7), 17:32], 1, function(x) mean(x^2))
  results2.2.4 <- apply(results2[c(2, 4, 6, 8), 17:32], 1, function(x) mean(x^2))
  results2.3.3 <- apply(results2[c(1, 3, 5, 7), 33:48], 1, function(x) mean(x^2))
  results2.3.4 <- apply(results2[c(2, 4, 6, 8), 33:48], 1, function(x) mean(x^2))
  results2.4.3 <- apply(results2[c(1, 3, 5, 7), 49:60], 1, function(x) mean(x^2))
  results2.4.4 <- apply(results2[c(2, 4, 6, 8), 49:60], 1, function(x) mean(x^2))
  results2.5.3 <- apply(results2[c(1, 3, 5, 7), 61:76], 1, function(x) mean(x^2))
  results2.5.4 <- apply(results2[c(2, 4, 6, 8), 61:76], 1, function(x) mean(x^2))
  test2 <- rbind(results2.1.3, # 4
                 results2.1.4, # 1
                 results2.2.3, # 3
                 results2.2.4, # 3
                 results2.3.3, # 3
                 results2.3.4, # 3
                 results2.4.3, # 1
                 results2.4.4, # 3
                 results2.5.3, # 3
                 results2.5.4) # 3
  apply(test2, 1, which.min) # results as a comment directly above
}



