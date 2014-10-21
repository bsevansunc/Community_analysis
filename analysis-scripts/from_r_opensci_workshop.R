# Code located at : bit.ly/ro-unc

# Reinhart-Rogoff scandal The Excel Depression ... read this!
# What is Markdown?
# Use Sweave!

# rmaps : for making interactive maps

library(rgbif)
library(dismo)
library(maptools)
library(plyr)



files <- list.files(paste(system.file(package="dismo"), '/ex', sep=''), 'grd', full.names=TRUE)
predictors <- stack(files)



data(wrld_simpl)

# The code below looks up the name in the gbif data base

namesdf <- name_lookup(query='bradypus*', rank="SPECIES", limit=60)$data
head(namesdf)



# Key values for each species (gbif naming code, gets rid of duplicates):

keys <- na.omit(unique(namesdf$nubKey))

# Search for the data (comes back as a list):

df <- occ_search(keys, georeferenced=TRUE, limit=100, return='data')

# Convert each individual list to a data frame:

df <- df[sapply(df, class) == "data.frame"]

# One data frame:

df2 <- ldply(df)

# Just lat and lon
df2 <- data.frame(lon=df2$longitude, lat=df2$latitude)

df2 = df2[df2$lon!=0,]

plot(predictors, 1)
plot(wrld_simpl, add=TRUE)
points(df2, col='blue', pch = 19, cex = .4)

#####

library(rgbif)
library(plyr)
library(doMC)



spplist <- c('Geothlypis trichas','Tiaris olivacea','Pterodroma axillaris','Calidris ferruginea','Pterodroma macroptera','Gallirallus australis','Falco cenchroides','Telespiza cantans','Oreomystis bairdi','Cistothorus palustris')



keys <- sapply(spplist, function(x) name_backbone(x, rank="species")$usageKey)
# remove NULLs
keys <- compact(keys)



countrynames <- as.character(isocodes$gbif_name)[1:25]



occ_by_countries <- function(spkey){
  occ_count_safe <- plyr::failwith(NULL, occ_count)
  tmp <- lapply(countrynames, function(x) occ_count_safe(spkey, country=x))
  names(tmp) <- countrynames
  tmp[grep("No enum", tmp)] <- NA
  tmp
}



registerDoMC(cores=4) 
out <- llply(keys, occ_by_countries)

# The above is equivalent to the below

## out <- lapply(keys, occ_by_countries)



names(out) <- spplist
df <- ldply(out, function(x){
  tmp <- ldply(x)
  names(tmp)[1] <- "country"
  tmp
})
df <- na.omit(df) # remove NAs (which were caused by errors in country names)



# Get only countries found in
df_foundin <- df[!df$V1==0,]



head(df_foundin); tail(df_foundin)

####

library(spocc)


out <- occ(query='Accipiter striatus', from='gbif')
out # prints summary of output data
out$gbif # GBIF data w/ metadata
out$ebird$data # empty
out$gbif$meta #  metadata, your query parameters, time the call executed, etc. 
out$gbif$data # just data


# inat is iNaturalist data (citizen scientist ... program?)
out <- occ(query='Accipiter striatus', from=c('gbif','inat'))
df <- occ2df(out)
head( df ); tail( df )

#### 

library(scales)
library(ggplot2)
library(doMC)
library(plyr)
library(rnoaa)



urls <- seaiceeurls(mo='Apr', pole='N')[1:12]



# registerDoMC(cores = 2)
# out <- llply(urls, noaa_seaice, storepath="~/", .parallel = TRUE)
out <- lapply(urls, noaa_seaice, storepath="~/")
names(out) <- seq(1979, 1990, 1)
df <- ldply(out)
head(df)



ggplot(df, aes(long, lat, group=group)) + 
  geom_polygon(fill = "steelblue") +
  theme_bw() +
  facet_wrap(~ .id)


##### Climate
library(rWBclimate)



gbr.dat.t <- get_ensemble_temp("GBR", "annualavg", 1900,2100)
gbr.dat.t <- subset(gbr.dat.t,gbr.dat.t$percentile == 50)

ggplot(gbr.dat.t,aes(x=fromYear,y=data,group=scenario, colour=scenario)) + 
  theme_bw(base_size=20) + 
  geom_point() + 
  geom_path() + 
  labs(x="Year", y="Annual Average Temperature in 20 year increments")



gbr.dat.p <- get_ensemble_precip("GBR", "annualavg", 1900,2100)
gbr.dat.p <- subset(gbr.dat.p,gbr.dat.p$percentile == 50)
ggplot(gbr.dat.p, aes(x=fromYear,y=data,group=scenario, colour=scenario)) + 
  theme_bw(base_size=20) + 
  geom_point() + 
  geom_path() + 
  labs(x="Year", y="Annual Average precipitation in mm")



country.list <- c("ISL", "FIN", "NOR", "SWE")
country.dat <- get_ensemble_stats(country.list, "mavg", "tmin_means")
####### Subset data Exclude A2 scenario
country.dat.b1 <- subset(country.dat, country.dat$scenario == "b1")
# choose just one percentile
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$percentile == 50)
# get just one year period
country.dat.b1 <- subset(country.dat.b1, country.dat.b1$fromYear == 2081)


ggplot(country.dat.b1, aes(x = month, y = data, group = locator, colour = locator)) + 
  geom_point() + geom_path() + ylab("Average monthly minimum temperature") + 
  theme_bw() + xlab("Month")

###

library(spocc)

red_tailed_hawk <- occ(query = "Buteo jamaicensis", from = c("gbif","ebird"),limit=35, ebirdopts = list(region='US'))

rt_hawk <-  occ2df(red_tailed_hawk)




#### Map with leaflet.js
mapleaflet(rt_hawk)

## put the map up as a gist

# Push map as a java script object to github
# geoJsn
mapgist(data = rt_hawk)

options(github.username = 'ropensciWks', github.password = 'ropensci1')


spp <- c('Accipiter gentilis','Accipiter striatus','Accipiter cooperii')
dat <- occ(query = spp, from = "gbif", gbifopts = list(georeferenced = TRUE))
dat <- fixnames(dat)
dat <- occ2df(dat)
mapgist(data = dat, color = c("#976AAE", "#6B944D", "#BD5945"))

###
library(spocc)
library(rgbif)
## Get some of our hawk data
### Prevent strings from becoming factors
options(stringsAsFactors = FALSE) # To ensure that character data remains as character data instead of factors
spp <- c('Accipiter gentilis','Accipiter striatus','Accipiter cooperii')
dat <- occ(query = spp, from = "gbif", gbifopts = list(georeferenced = TRUE))
dat <- fixnames(dat)
dat <- occ2df(dat)

map_dat <- stylegeojson(input = dat, var = "name", color = c("#8BA8D9", "#8BD99D", "#FFEF0D"), symbol="zoo")
mapgist(map_dat)

# See Mapbox for map symbol list

### Add common names
dat$common_name <- dat$name
dat$common_name[dat$name %in% 'Accipiter gentilis'] <- "Northern goshawk"
dat$common_name[dat$name %in% 'Accipiter striatus'] <- "Sharp-shinned hawk"
dat$common_name[dat$name %in% 'Accipiter cooperii'] <- "Coopers hawk"


map_dat <- stylegeojson(input = dat, var = "name", color = c("#8BA8D9", "#8BD99D", "#FFEF0D"), size="small")
mapgist(map_dat)




spp <- c('Accipiter gentilis','Accipiter striatus','Accipiter cooperii')
dat <- occ(query = spp, from = c("gbif","inat"), gbifopts = list(georeferenced = TRUE))
dat <- fixnames(dat)
dat <- occ2df(dat)

map_dat <- stylegeojson(input = dat, var_col = "name", color = c("#8BA8D9", "#8BD99D", "#FFEF0D"), var_sym = "prov", symbol=c("g","i"))
mapgist(map_dat)

###

### This is all the code to clean the data from source, but we'll use a cleaned version in class

## Grab our beer data  
options(stringsAsFactors = FALSE)

### Source page: https://opendata.socrata.com/Government/UtahBeerTaxmap/6i4w-nzeq
utahBT <- read.csv("http://opendata.socrata.com/api/views/6i4w-nzeq/rows.csv")
## Clean beer data

# delete header
utahBT <- utahBT[-1,]

### parse lat lon data
lat <- rep(NA,dim(utahBT)[1])
lon <- rep(NA,dim(utahBT)[1])
loc_text <- strsplit(utahBT$Location,",")
for(i in 1:length(loc_text)){
  if(length(loc_text[[i]]) > 2){
    lat[i] <- as.numeric(strsplit(loc_text[[i]][2],"(",fixed=T)[[1]][2])
    lon[i] <- as.numeric(strsplit(loc_text[[i]][3],")"))
  }
}


utahBT$lat <- lat
utahBT$lon <- lon
utahBT <- utahBT[,-3]
utahBT <- utahBT[!is.na(utahBT$lat),]
colnames(utahBT) <- c("City","Distribution","latitude","longitude")
utahBT$Distribution <- as.numeric(gsub("$","",utahBT$Distribution,fixed=T))
write.csv(utahBT,"data/utahbt.csv")


## To download data use this URL: https://raw.github.com/ropensci/workshops-unc-2014-02/master/data/utahbt.csv
###  Now we can add some categories, let's start with binning

# Lump into categories (using log due to large spread):

utahBT$taxLev <- cut(log(utahBT$Distribution),breaks=10,labels=F)

# Add color ramp palette:

utahBT$"marker-color" <- colorRampPalette(c("blue", "red"))(10)[utahBT$taxLev]
utahBT$"marker-size" <- rep("small",dim(utahBT)[1])
mapgist(utahBT)




# The below does not work on a local machine
togeojson(input = "data/abiemagn/abiemagn.shp", method = "local", destpath = paste(getwd(),"/",sep=""), outfilename = "abiesmagmap")
gist("abiesmagmap.geojson", description = "Abies magnifica polygons")


# For rMaps, I have to install directly from github

library(rMaps)
map <- Leaflet$new()
map$setView(c(44.4758, -73.2119), zoom = 15)
map$tileLayer(provider = 'OpenStreetMap.Mapnik')
map$marker(
  c(44.4758, -73.2119),
  bindPopup = 'This is where I met Brian'
)
map

# SPATIAL STATS IN R

library(rWBclimate)
### Create path to store kml's
dir.create("~/kmltmp")
options(kmlpath="~/kmltmp")

usmex <- c(273:284,328:365)
usmex.basin <- create_map_df(usmex)


## Download temperature data
temp.dat <- get_historical_temp(usmex, "decade" )
temp.dat <- subset(temp.dat,temp.dat$year == 2000 )


#create my climate map
usmex.map.df <- climate_map(usmex.basin,temp.dat,return_map=F)


## Grab some species occurrence data for the 8 tree species.

splist <- c("Acer saccharum",
            "Abies balsamea",
            "Arbutus texana",
            "Betula alleghaniensis",
            "Chilopsis linearis",
            "Conocarpus erectus",
            "Populus tremuloides",
            "Larix laricina")
out <- occ(query = splist, from="gbif",limit = 100)

## Now just create the base temperature map
usmex.map <- ggplot()+geom_polygon(data=usmex.map.df,aes(x=long,y=lat,group=group,fill=data,alpha=.8))+scale_fill_continuous("Average annual \n temp: 1990-2000",low="yellow",high="red")+ guides(alpha=F)+theme_bw()

out <- fixnames(out)
out_df <- occ2df(out)
## And overlay of gbif data
usmex.map + geom_point(data=out_df,aes(y=latitude,x=longitude,group=name,colour= name)) + xlim(-125,-59)+ylim(5,55)



library(devtools)
library(rWBclimate)
library(rgbif)
library(plyr)
library(ggplot2)
library(sp)


##### Converts decimal latitude to a spatial object with WGS84 coords
to_pts <- function(x,y){
  ### check valid coords
  xy <- data.frame(cbind(x,y))
  xy <- xy[(xy$x < 180) & (xy$x > -180),]
  xy <- xy[(xy$y < 90) & (xy$y > -90),]
  return(SpatialPoints(xy,proj4string=CRS("+proj=longlat +datum=WGS84")))
}


### Download map polygons
usmex <- c(273:284,328:365)
usmex.basin <- create_map_df(usmex)

## Download temperature data
temp.dat <- get_historical_temp(usmex, "decade" )
temp.dat <- subset(temp.dat,temp.dat$year == 2000 )

## Create a spatial polygon dataframe binding kml polygons to temperature data
temp_sdf <- kml_to_sp(usmex.basin,df=temp.dat)

### Let's plot it just to compare to our KML plots
spplot(temp_sdf,"data",col.regions = rainbow(100, start = 4/6, end = 1))

### Species name list
splist <- c("Acer saccharum",
            "Abies balsamea",
            "Arbutus texana",
            "Betula alleghaniensis",
            "Chilopsis linearis",
            "Conocarpus erectus",
            "Populus tremuloides",
            "Larix laricina")
## Common name list
cname <- c("Sugar Maple",
           "Balsam Fir",
           "Texas Madrone",
           "Yellow Birch",
           "Desert Willow",
           "Mangrove shrub",
           "Quaking Aspen",
           "American Larch"
)

## Set up the keys to download from GBIF
keys <- sapply(splist, function(x) name_backbone(name=x, kingdom='plants')$speciesKey, USE.NAMES=FALSE)

### This may take awhile, it's lots of data
mdat <- occ_search(taxonKey=keys, limit=200, return='data', georeferenced=TRUE)

### Convert gbif floating points to spatial polygons
newlist <- list()
for(i in 1:length(cname)){
  newlist[[cname[i]]] <- to_pts(mdat[[i]]$longitude,mdat[[i]]$latitude)
}


N <- vector()
spp <- vector()
dat <- vector()
coords <- matrix(NA,ncol=2,nrow=0)
### Get averages
for(i in 1:length(newlist)){
  tmp_t <- over(newlist[[i]],temp_sdf)$data  
  coords <- rbind(coords,cbind(newlist[[i]]$x[!is.na(tmp_t)],newlist[[i]]$y[!is.na(tmp_t)]))
  N[i] <- sum(!is.na(tmp_t))
  tmp_t <- tmp_t[!is.na(tmp_t)]
  
  dat <- c(dat,tmp_t)
  spp <- c(spp,rep(cname[i],length(tmp_t)))
  
}


coords <- data.frame(coords)
coords$spp <- as.factor(spp)
coords$dat <- dat
colnames(coords) <- c("Longitude","Latitude","Species","Temp")


usmex.map.df <- climate_map(usmex.basin,temp.dat,return_map=F)

usmex.map <- ggplot()+geom_polygon(data=usmex.map.df,aes(x=long,y=lat,group=group,fill=data,alpha=.8))+scale_fill_continuous("Average annual \n temp: 1990-2000",low="yellow",high="red")+ guides(alpha=F)+theme_bw()

usmex.map + geom_point(data=coords,aes(y=Latitude,x=Longitude,group=Species,colour= Species)) + xlim(-125,-59)+ylim(5,55)



ggplot(coords,aes(x=Species,y=Temp))+geom_boxplot()+coord_flip()

summary_data <- ddply(coords,.(Species), summarise,mlat = mean(Latitude),mtemp = mean(Temp),sdlat = sd(Latitude),sdtemp = sd(Temp))


ggplot(summary_data,aes(x=mlat,y=mtemp))+geom_point()+ geom_errorbar(aes(ymin = mtemp - sdtemp, ymax=mtemp+sdtemp,width=.5)) + geom_errorbarh(aes(xmin = mlat - sdlat, xmax=mlat+sdlat,width=.5)) + stat_smooth(method="lm",se=F)+xlab("Mean Latitude") + ylab("Mean Temperature (C)") + theme_bw()

ggplot(summary_data,aes(x=mlat,y=mtemp,label=Species))+geom_text()+xlab("Mean Latitude") + ylab("Mean Temperature (C)") + theme_bw()+ xlim(10,50)

