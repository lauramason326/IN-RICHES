library(ggmap)
library(ggplot2)
library(tidyverse)


Sites<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")

colnames(Sites)

map <- Sites %>%
  select(1,6:8,12,13)%>%
  na.omit(.)


glr<-get_googlemap(center=c(lon=-105.0844,lat=40.5853),zoom=5,maptype="terrain")
#fort collins: 40.5853° N, 105.0844° W

map$Long_fixed <- iconv(map$Long, from = "", to = "UTF-8", sub = "")
map$Long <- as.numeric(as.character(map$Long_fixed))
map$Lat <- as.numeric(as.character(map$Lat))

ggmap(glr) + 
  geom_point(
    data = map,
    aes(x = Long, y = Lat, fill = type, shape = dataset),
    size = 3, stroke = 0.5
  ) + 
  scale_shape_manual(values = c(21, 22, 23, 24)) + 
  scale_fill_manual(values = c(
    "#CFFC00", "#008FB3","#E56A54", "#D9782D", "#7E5475", 
    "#105456", "#FFC038", "#C8C372", "#006144", "#82C503"
  )) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),  # force filled shape in fill legend
    shape = guide_legend(override.aes = list(fill = "grey"))  # neutral fill in shape legend
  )
