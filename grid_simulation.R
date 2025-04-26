library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)
library(readxl)
library(ggpubr)
library(grid)
library(gridExtra)
library(patchwork)
library(cowplot)
library(terra)
library(tidyterra)

### Treeline
trel <- st_read("Data/treeline_la/treeline_la.shp")

### Tundra
stars_cavm <- read_stars("Data/RasterCAVM/raster_cavm_v1.tif")

bfr  <- st_point(c(0, 0)) %>% st_sfc(crs = st_crs(trel)) %>%
  st_buffer(as.numeric(max(st_point(c(0, 0)) %>% st_sfc(crs = st_crs(trel)) %>% 
                             st_distance(trel %>% st_cast("POINT")) %>% suppressWarnings())) + 250000) %>%
  st_transform(st_crs(stars_cavm))

cavmRast <- stars_cavm %>% st_as_stars(downsample = 1) %>% st_crop(bfr) %>%
  st_apply(., 1:2, function(x) ifelse(x[1]<92, 1, ifelse(x[1]==92, NA, ifelse(x[1]==93, -1, 0))), FUTURE = T) %>% 
  setNames("cavm")

### Parameters
{
  kruseData <- read.table("Data/MigrationRates.csv", skip = 1, sep = ";")[,-6] %>%
    as_tibble() %>% setNames(c("Region", "Scenario", "Year", "Rate", "Range")) %>%
    mutate(Rate = as.numeric(Rate), Range = as.numeric(Range), Scenario = substring(Scenario, 1, 7))
  
  
  MigRates <- kruseData %>% group_split(Scenario) %>% lapply(function(x) {
    cbind(median(x$Rate, na.rm = T), mean(x$Range, na.rm = T))
  }) %>% Reduce("rbind", .) %>% as_tibble() %>% 
    mutate(Scenario = sapply(kruseData %>% group_split(Scenario), function(x) x$Scenario[1]), .before = "V1") %>%
    setNames(c("Scenario", "Median", "Range")) %>%
    mutate(pexp = 
        sapply(1:3, function(x) {
          kmSeq <- seq(0, 50, by = 0.05)
          pseq  <- seq(0.001, 0.8, by = 0.001)
          out <- sapply(pseq, function(y) {
            pexp(kmSeq, y, lower.tail = FALSE)[which.min(abs(kmSeq - Median[x]))]
          })
          pseq[which.min(abs(out-0.5))]
          }))
  
  # plot(tibble(kmSeq = seq(0, 50, by = 0.05), pext = NA), ylim = c(0,1), type = "n")
  # lapply(1:3, function(x) lines(seq(0, 50, by = 0.05), pexp(seq(0, 50, by = 0.05), MigRates$pexp[x], lower.tail = FALSE), 
                               # col = c("goldenrod1", "orange", "darkred")[x]))
  
  plot(tibble(kmSeq = seq(0, 50, by = 0.05), pext = NA), ylim = c(0,1), type = "n")
  lapply(1:3, function(x) lines(seq(0, 50, by = 0.05), 1-pnorm(seq(0, 50, by = 0.05), MigRates$Median[x], 2), 
                                col = c("goldenrod1", "orange", "darkred")[x]))
  }

path <- "~/Output/"
rcps <- c('2.6', '4.5', '8.5')
  
for(rcp in 1:3) {
  
  dir.create(glue::glue("{path}/{rcps[rcp]}"))
  
  treelDist <- c(cavmRast %>% rast(), (cavmRast %>% rast() %>% mutate(distance = ifelse(is.na(cavm) | cavm!=0, NA, 1)))[[2]] %>% distance()/1000)
  
  writeRaster(treelDist[[1]], glue::glue("{path}/{rcps[rcp]}/sim_2010.tif"))
  
  for(y in seq(2020,2300, by = 10)) {
    out <- app(treelDist, function(x, p1) ifelse(x[1]==1 & runif(1) > pnorm(x[2], p1, 2.5), 0, x[1]), p1 = MigRates$Median[rcp], cores = 5) %>%
      setNames("cavm")
    
    writeRaster(out, glue::glue("{path}/{rcps[rcp]}/sim_{y}.tif"))
    
    treelDist <- c(out, (out %>% mutate(distance = ifelse(is.na(cavm) | cavm!=0, NA, 1)))[[2]] %>% distance()/1000)
  }
}


### Images
map <-  rnaturalearth::ne_countries(scale = 50, returnclass = 'sf') %>% st_geometry() %>%
  st_transform(st_crs(bfr)) %>% st_intersection(bfr) %>% st_buffer(25) %>% st_union() %>% st_simplify(dTolerance = 30)


path <- '/Users/slisovsk/Documents/treelineMigration_rasts'
fls  <- tibble(path = path, fls = list.files(path, pattern = "*.tif", recursive = T)) %>%
  mutate(scenario = sapply(strsplit(fls, "/"), function(x) x[[1]]),
         year = as.numeric(sapply(strsplit(fls, "_"), function(x) gsub(".tif", "", x[[2]]))))

for(y in unique(fls$year)) {
  
  rasts <- lapply(fls %>% filter(year == y) %>% dplyr::select(path, fls) %>% apply(1, function(x) paste(x[1], x[2], sep = '/')),
                   function(x) read_stars(x) %>% setNames("treeline"))
             
  plot_list <- lapply(rasts, function(x) {
    pl <- ggplot() +
      geom_stars(data = x, mapping = aes(fill = as.factor(treeline)), downsample = 2, show.legend = FALSE) +
      scale_fill_manual(values = c("tan", "olivedrab", "white"), limits = as.factor(c(1, 0, -1)), na.value = "transparent") +
      geom_sf(data = map, fill = NA, linewidth = 0.25) +
      geom_sf(data = bfr, fill = NA, linewidth = 0.4) +
      coord_sf() +
      theme_void()
  
    return(pl)
  })
  
  out <- ggarrange(plotlist = plot_list, nrow = 1, ncol = 3, labels = c('Scenario: 2.6', '4.5', '8.5'),
                   font.label = list(size = 50))
  
  ggsave(glue::glue("~/Documents/treelineMigration_rasts/tmp/tmp_{y}.png"), out, 
         width = 1593/30, height = 451/30, bg = "white", limitsize = FALSE)        
 }


### Animation
library(magick)
library(magrittr)

list.files(path='~/Documents/treelineMigration_rasts/tmp', pattern = '*.png', full.names = TRUE) %>%
  image_read() %>%
  image_join() %>%
  image_animate(fps = 2) %>%
  image_write("~/Documents/treelineMigration_rasts/tmp/Treeline_25022025.gif")



