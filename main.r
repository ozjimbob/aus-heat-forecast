stime=Sys.time()
library(sp)
library(maptools)
library(raster)
library(mgcv)
library(gstat)

# Set because my desktop RStudio installation is set to use Internet 2 setting, which don't work.
setInternet2(use=F)

cat("Loading 95th percentile temperature raster...")
# Load raster of 95th percentile average daily temperature
temp95=raster("avt95.tif")
cat("Done.\n")

# Load australia shape
cat("Loading shapefile of Australia...")
aus = readShapePoly("aus.shp")
cat("Done.\n")

# Load in mean daily temperature for previous 30 days.  We want both maximum and minimum temperature and will average them.
cat("Determining current date and daily temperature raster base URLs...")
dnow=as.Date(Sys.time())
g=2:30
dnow=dnow-g
dnow=format(dnow,"%Y%m%d")
f_max_list=paste("http://www.bom.gov.au/web03/ncc/www/awap/temperature/maxave/daily/grid/0.05/history/nat/",dnow,dnow,".grid.Z",sep="")
f_min_list=paste("http://www.bom.gov.au/web03/ncc/www/awap/temperature/minave/daily/grid/0.05/history/nat/",dnow,dnow,".grid.Z",sep="")
cat("Done.\n")

# We will download the daily maximum temperature rasters for the previous month, and put them in a list.
cat("Downloading daily maximum temperature rasters...")
rstore=list()
for(idx in seq_along(f_max_list)){
  tempfile=paste("tempm",idx,".Z",sep="")
  download.file(f_max_list[idx],tempfile,mode="wb")
  cmd=paste("gunzip ",tempfile,sep="")
  system(cmd)
  rasterfile=paste("tempm",idx,sep="")
  thisrast=raster(rasterfile)
  readAll(thisrast)
  rstore[[idx]]=thisrast
}
# Turn the list into a brick
cat("Bricking...")
rstore=brick(rstore)
mean_max = mean(rstore)
rstore <- NULL
gc()
cat("Done.\n")

cat("Deleting rasters...")
# Delete the rasters and their zip files
unlink("temp*")
cat("Done.\n")

cat("Downloading daily minimum temperature rasters...")
# Do the same for the minimum temperature
rstore=list()
for(idx in seq_along(f_min_list)){
  tempfile=paste("tempx",idx,".Z",sep="")
  download.file(f_min_list[idx],tempfile,mode="wb")
  cmd=paste("gunzip ",tempfile,sep="")
  system(cmd)
  rasterfile=paste("tempx",idx,sep="")
  thisrast=raster(rasterfile)
  readAll(thisrast)
  rstore[[idx]]=thisrast
}
cat("Bricking...")
# Make them into a brick
rstore=brick(rstore)
mean_min = mean(rstore)
rstore <- NULL
gc()
cat("Done.\n")

cat("Deleting rasters...")
# Delete the temporary files
unlink("temp*")
cat("Done.\n")

cat("Calculating mean temperature raster...")
# Mean temperature is the average of maximum and minimum
mean_temp = (mean_max + mean_min) / 2
cat("Done.\n")

# Load station locations
cat("Loading station locations and converting to SpatialPoints...")
st_locs=read.csv("ftp://ftp2.bom.gov.au/anon/gen/fwo/IDY02126.dat",skip=1,col.names=c("bom_st","bom_name","met_code","lat","lon","ele","wmo","state"))
coordinates(st_locs)=c("lon","lat")
cat("Done.\n")

# Load numerical forecast data
#Run at 8am, will do for the following 3 days and produce the result prior to 9am (AEST Time)
cat("Downloading numerical forcast data...")
forcast=read.csv("ftp://ftp2.bom.gov.au/anon/gen/fwo/IDY02123.dat",skip=3,na.strings="-9999.0")
cat("Done.\n")
cat("Removing missing values and calculating mean temp...")
# Remove missing forcast values
forcast=subset(forcast,!is.na(forcast$per))
# Calculate mean temperature from maximum and minimum
forcast$temp = (forcast$amax + forcast$amin)/2
cat("Done.\n")

cat("Kriging forecast points to rasters...")
# We will store the rasters in a list, again
nowtemp=list()
for(dx in 1:3){
  # Extract the forecast for the relevent day in the future
  forcast_d = subset(forcast,per==dx)
  forcast_d$stn.7.=as.numeric(as.character(forcast_d$stn.7.))
  # For each station we have in the locational table, extract the temperature
  for(dx2 in seq_along(st_locs)){
    this.temp=subset(forcast_d,forcast_d$stn.7.==st_locs$bom_st[dx2])$temp[1]
    st_locs$temp[dx2]=this.temp
  }
  # Clip to Australia, we don't need outlying territories
  in_vec=coordinates(st_locs)[,1] > 112.5 & coordinates(st_locs)[,1] < 155.8 & coordinates(st_locs)[,2] > -44.5 & coordinates(st_locs)[,2]< 9.6 
  st_locs2=st_locs[in_vec,]
  st_locs2=subset(st_locs2,!is.na(st_locs2$temp))
  
  # Spatial krieg of temperature data
  xy <- data.frame(xyFromCell(mean_temp, 1:ncell(mean_temp)))
  v<-variogram(temp~1,st_locs2)
  m<-fit.variogram(v,vgm(1,"Sph",200,1))
  kr=gstat(NULL,"temp",temp~1,data=st_locs2,model=m)

  temp.op = mean_temp
  coordinates(xy)=c("x","y")
  kr.p=predict(kr,newdata=xy)
  values(temp.op) = as.vector(kr.p@data$temp.pred)
  nowtemp[[dx]]=temp.op
}
cat("Done.\n")

cat("Bricking forcast raters...")
# Brick the raster list and take the mean for the next 3 days
nowtemp = brick(nowtemp)
nowtemp = mean(nowtemp)
cat("Done.\n")

# EHF Calculation based on the long term, short term and forecast rasters.
cat("Calculating EHF raster...")
diff_ac=nowtemp-mean_temp
one_ac=nowtemp
values(one_ac)=1
diff_ac=max(one_ac,diff_ac)
diff_lt = nowtemp - temp95
EHF = diff_lt * diff_ac
cat("Done.\n")

# Mask the EHF
cat("Masking EHF...")
EHF= mask(EHF,aus)
cat("Done.\n")

# Produce output maps
# First, the raw EHF data
cat("Output EHF map to disk...")
png("output//EHF_map.png",1200,1200)
par(mai=c(0.6,0.6,1,0.6))
plot(EHF,main="Heat Index Forecast",col=colorRampPalette(c("blue","green","yellow","orange","red"))(100),horizontal=T)
plot(aus, add=T)
dev.off()
cat("Done.\n")

cat("Calculate positive EHF cells...")
# Now we identify cells with a positive EHF, ie. Heat wave
pos=EHF
values(pos)[values(pos<=0)]=0
values(pos)[values(pos>0)]=1
cat("Done.\n")

cat("Loading EHF Threshold raster...")
# Load in the heat wave threshold raster
thresh=raster("EHF_threshold.tif")
cat("Done.\n")

cat("Calcualte Threshold multipliers... ")
values(thresh)[values(thresh<1)]=1


# Produce a raster where cells not in heat wave are shown as zero
# Cells in heat wave are positive
# And the value of cells indicate exceedence of local EHF 85% threshold
zones=EHF/thresh
values(zones)[values(zones<=1)]=0
values(zones)[values(zones>1) & values(zones<=2)]=1
values(zones)[values(zones>2) & values(zones<=3)]=2
values(zones)[values(zones>3)]=3
values(zones)[values(pos)==0]=NA
cat("Done.\n")

# Output zones map
cat("Output threshold zone map to disk...")
png("output//EHF_zones.png",1200,1200)
par(mai=c(0.6,0.6,1,0.6))
plot(zones,main="Heat Wave Threshold",col=colorRampPalette(c("yellow","orange","red","purple"))(4),breaks=c(0,1,2,3,4),lab.breaks=c("None","Heatwave","Severe","Extreme","Extreme"),horizontal=T)
plot(aus, add=T)
dev.off()
cat("Done.\n")

# For our location table, join it to the raster to find which stations exceed
cat("Extract zone values for met station table...")
hw_stats_v=extract(zones,st_locs2)
st_locs2$hw=hw_stats_v
st_data=st_locs2@data
hw_stations=subset(st_data,!is.na(hw))
cat("Done.\n")


# Produce HTML output. First make an empty vector to hold the HTML
output=c()
# Now add text.
cat("Produce heatwave table...")
output=c(output,"<h1>Locations with heat wave forecast</h1>")
### Cycle through heatwave stations
statelist = unique(hw_stations$state)
for(thestate in statelist){
  output=c(output,paste("<h2>",thestate,"</h2>"),"<br />\n")
  state_sub=subset(hw_stations,state==thestate)
  loclist=paste(state_sub$bom_name,"<br />\n")
  output=c(output,loclist)
}
cat("Done.\n")

# Now for severe/extreme heat waves
cat("Produce extreme/severe heatwave table...")
output=c(output,"<h2>Locations with extreme/severe heatwaves forecast</h1>")
hw_stations=subset(hw_stations,hw>0)

### Cycle through heatwave stations
statelist = unique(hw_stations$state)
for(thestate in statelist){
  output=c(output,paste("<h2>",thestate,"</h2>"),"<br />\n")
  state_sub=subset(hw_stations,state==thestate)
  loclist=paste(state_sub$bom_name,"<br />\n")
  output=c(output,loclist)
}
cat("Done.\n")

# Produce html page
cat("Format HTML...")
output=c("<h1>EHF Forecast Values for Next 3 days.</h1><img src=\"EHF_map.png\" /><br />\n",output)
output=c("<h1>Heat wave category.</h1><img src=\"EHF_zones.png\" /><br />\n",output)
output=c("<html><head><title>Australian Heat Wave Forecast</title></head>",output)
output=c(output,"</body></html>")
write(output,"output\\index.html")
cat("Done.\n")
etime=Sys.time()
tdiff=etime-stime
cat(paste("Execution time: ",tdiff," ",attributes(h)$units," \n",sep=""))