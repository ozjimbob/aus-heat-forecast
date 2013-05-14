# Process for individual station
library(sp)
library(maptools)
library(raster)
library(mgcv)
library(gstat)

st_root="G:\\SILO_Data\\complete_stations\\"
o_path="G:\\heat\\"
files=list.files("G:\\SILO_Data\\complete_stations\\", pattern=".csv")
out_frame=data.frame(station=0,EHF_t=0,min_95=0,min_99=0,max_95=0,max_99=0,mean_95=0,mean_99=0)

for(the_file in seq_along(files)){
  print(the_file/length(files))
  station=substr(files[the_file],1,nchar(files[the_file])-4)
  st_path=paste(st_root,station,".csv",sep="")
  data=read.csv(st_path)
  
  data$t_mean = (data$t_max + data$t_min) / 2
  data$c_date=as.Date(as.character(data$c_date))
  st_95=as.Date("1960-01-01")
  en_95=as.Date("1989-12-31")
  sub_data=subset(data,data$c_date > st_95 & data$c_date < en_95)
  
  temp_95 = quantile(sub_data$t_mean,0.95)
  
  v_EHF=rep(0,length(data$t_mean))
  
  for(idx in 365:length(data$t_mean)){
    now_temp=mean(data$t_mean[idx:(idx+2)],na.rm=T)
    mean_temp = mean(data$t_mean[(idx-30):(idx-1)],na.rm=T)
    diff_ac=now_temp-mean_temp
    diff_ac=max(1,diff_ac)
    diff_lt = now_temp - temp_95
    EHF = diff_lt * diff_ac
    v_EHF[idx]=EHF
  }
  
  data$EHF=v_EHF
  
  pos=data$EHF[data$EHF>0]
  thr=quantile(pos,0.85)
  min_95=quantile(data$t_min,0.95)
  min_99=quantile(data$t_min,0.99)
  max_95=quantile(data$t_max,0.95)
  max_99=quantile(data$t_max,0.99)
  mean_95=quantile(data$t_mean,0.95)
  mean_99=quantile(data$t_mean,0.99)
  
  o_vec=c(station,thr,min_95,min_99,max_95,max_99,mean_95,mean_99)
  out_frame=rbind(out_frame,o_vec)
  this_data=data.frame(stnum=data$stnum,c_date=data$c_date,t_min=data$t_min,t_max=data$t_max,t_mean=data$t_mean,EHF=data$EHF)
  write.csv(this_data,paste(o_path,station,".csv",sep=""))
}

out_frame=subset(out_frame,station!="0")

out_frame$EHF_t=as.numeric(out_frame$EHF_t)
out_frame$min_95=as.numeric(out_frame$min_95)
out_frame$min_99=as.numeric(out_frame$min_99)
out_frame$max_95=as.numeric(out_frame$max_95)
out_frame$max_99=as.numeric(out_frame$max_99)
out_frame$mean_95=as.numeric(out_frame$mean_95)
out_frame$mean_99=as.numeric(out_frame$mean_99)
out_frame$station=as.numeric(out_frame$station)
write.csv(out_frame,"station_EHF_thresholds.csv")

st_locs=read.csv("ftp://ftp2.bom.gov.au/anon/gen/fwo/IDY02126.dat",skip=1,col.names=c("bom_st","bom_name","met_code","lat","lon","ele","wmo","state"))
coordinates(st_locs)=c("lon","lat")
for(dx2 in seq_along(st_locs)){
  this.temp=subset(out_frame,out_frame$station==st_locs$bom_st[dx2])$EHF_t[1]
  print(st_locs$bom_st[dx2])
  print(this.temp)
  st_locs$EHF[dx2]=this.temp
}

st_locs=subset(st_locs,!is.na(st_locs$EHF))
in_vec=coordinates(st_locs)[,1] > 112.5 & coordinates(st_locs)[,1] < 155.8 & coordinates(st_locs)[,2] > -44.5 & coordinates(st_locs)[,2]< 9.6 
st_locs2=st_locs[in_vec,]

temp95=raster("avt95.tif")

xy <- data.frame(xyFromCell(temp95, 1:ncell(temp95)))
v<-variogram(EHF~1,st_locs2)
m<-fit.variogram(v,vgm(1,"Exp",5))
kr=gstat(NULL,"EHF",EHF~1,data=st_locs2,model=m)

temp.op = temp95
coordinates(xy)=c("x","y")
kr.p=predict(kr,newdata=xy)
values(temp.op) = as.vector(kr.p@data$EHF.pred)

writeRaster(temp.op,"EHF_threshold.tif",format="GTiff")