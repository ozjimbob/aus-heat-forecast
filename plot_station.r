## Plot graph for one of a number of selected stations

stations=list("Hobart"=94029,
              "New Norfolk"=95052,
              "Launceston"=91237,
              "Burnie"=91009,
              "Swansea"=92038,
              "Queenstown"=97091,
              "Alice"=15540)

o_path="G:\\heat\\"
threshold_frame=read.csv("station_EHF_thresholds.csv")

get_station_data=function(select_station){
  station_number=as.numeric(stations[select_station])
  station_data=read.csv(paste(o_path,station_number,".csv",sep=""),stringsAsFactors=F)
  station_data$c_date=as.Date(station_data$c_date)
  station_data
}

get_threshold=function(select_station){
  station_number=as.numeric(stations[select_station])
  this_threshold=subset(threshold_frame,threshold_frame$station==station_number)
  this_threshold$EHF_t[1]
}

plot_station=function(select_station){
  data=get_station_data(select_station)
  threshold=get_threshold(select_station)
  station_name=names(stations[select_station])
  plot(data$EHF ~ data$c_date,type="l",main=station_name,xlab="Date",ylab="EHF")
  abline(h=threshold,col="green")
  abline(h=threshold*2,col="orange")
  abline(h=threshold*3,col="red")
}

google_plot_station=function(select_station){
  require(googleVis)
  data=get_station_data(select_station)
  threshold=get_threshold(select_station)
  station_name=names(stations[select_station])
  data=data[(length(data$EHF)-1000):(length(data$EHF)),]
  out=gvisAnnotatedTimeLine(data,datevar="c_date",numvar="EHF",date.format="%Y-%m-%d")
  out
}