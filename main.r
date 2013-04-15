library(sp)
library(maptools)
library(raster)
library(mgcv)
library(gstat)

# Load raster of 95th percentile average daily temperature
temp95=raster("avt95.tif")

# Load australia shape
aus = readShapePoly("aus.shp")
# Load in mean daily temperature for previous 30 days
dnow=as.Date(Sys.time())
g=2:30
dnow=dnow-g
dnow=format(dnow,"%Y%m%d")
f_max_list=paste("http://www.bom.gov.au/web03/ncc/www/awap/temperature/maxave/daily/grid/0.05/history/nat/",dnow,dnow,".grid.Z",sep="")
f_min_list=paste("http://www.bom.gov.au/web03/ncc/www/awap/temperature/minave/daily/grid/0.05/history/nat/",dnow,dnow,".grid.Z",sep="")

rstore=list()
for(idx in seq_along(f_max_list)){
  print(idx)
  tempfile=paste("tempm",idx,".Z",sep="")
  download.file(f_max_list[idx],tempfile)
  cmd=paste("gunzip ",tempfile,sep="")
  system(cmd)
  rasterfile=paste("tempm",idx,sep="")
  thisrast=raster(rasterfile)
  readAll(thisrast)
  rstore[[idx]]=thisrast
}

rstore=brick(rstore)
mean_max = mean(rstore)
rstore <- NULL
gc()

system("rm temp*")

rstore=list()
for(idx in seq_along(f_min_list)){
  print(idx)
  tempfile=paste("tempm",idx,".Z",sep="")
  download.file(f_min_list[idx],tempfile)
  cmd=paste("gunzip ",tempfile,sep="")
  system(cmd)
  rasterfile=paste("tempm",idx,sep="")
  thisrast=raster(rasterfile)
  readAll(thisrast)
  rstore[[idx]]=thisrast
}

rstore=brick(rstore)
mean_min = mean(rstore)
rstore <- NULL
gc()

system("rm temp*")

mean_temp = (mean_max + mean_min) / 2

# Load station location
st_locs=read.csv("ftp://ftp2.bom.gov.au/anon/gen/fwo/IDY02126.dat",skip=1,col.names=c("bom_st","bom_name","met_code","lat","lon","ele","wmo","state"))
coordinates(st_locs)=c("lon","lat")

# Load numerical forecast data
forcast=read.csv("ftp://ftp2.bom.gov.au/anon/gen/fwo/IDY02122.dat",skip=3,na.strings="-9999.0")
forcast=subset(forcast,!is.na(forcast$per))
forcast$temp = (forcast$amax + forcast$amin)/2
nowtemp=list()
for(dx in 1:3){
  forcast_d = subset(forcast,per==dx)
  forcast_d$stn.7.=as.numeric(as.character(forcast_d$stn.7.))
  for(dx2 in seq_along(st_locs)){
    this.temp=subset(forcast_d,forcast_d$stn.7.==st_locs$bom_st[dx2])$temp[1]
    print(st_locs$bom_st[dx2])
    print(this.temp)
    st_locs$temp[dx2]=this.temp
  }
  #st_locs$temp = forcast_d$temp
  in_vec=coordinates(st_locs)[,1] > 112.5 & coordinates(st_locs)[,1] < 155.8 & coordinates(st_locs)[,2] > -44.5 & coordinates(st_locs)[,2]< 9.6 
  st_locs2=st_locs[in_vec,]
  st_locs2=subset(st_locs2,!is.na(st_locs2$temp))
  
  xy <- data.frame(xyFromCell(mean_temp, 1:ncell(mean_temp)))
  v<-variogram(temp~1,st_locs2)
  m<-fit.variogram(v,vgm(1,"Sph",200,1))
  kr=gstat(NULL,"temp",temp~1,data=st_locs2,model=m)

  temp.op = mean_temp
  coordinates(xy)=c("x","y")
  #proj4string(xy)="+proj=longlat +datum=WGS84"
  kr.p=predict(kr,newdata=xy)
  #proj4string(temp.op)="+proj=longlat +datum=WGS84"
  values(temp.op) = as.vector(kr.p@data$temp.pred)
  nowtemp[[dx]]=temp.op
}

nowtemp = brick(nowtemp)
nowtemp = mean(nowtemp)

diff_ac=nowtemp-mean_temp
one_ac=nowtemp
values(one_ac)=1
diff_ac=max(one_ac,diff_ac)

diff_lt = nowtemp - temp95
EHF = diff_lt * diff_ac
EHF= mask(EHF,aus)
# Produce output maps
plot(EHF,main="Heat Index Forecast",col=colorRampPalette(c("blue","green","yellow","orange","red"))(100))
plot(aus, add=T)