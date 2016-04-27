

library(raster)
library(rgdal)
library(gdalUtils)


infile = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132AAA04_B4.TIF"
outfile = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132AAA04_B4_temp1.tif"
cmd = paste("gdal_translate -of GTiff -co profile=baseline -co tfw=yes", infile, outfile)
shell(cmd)
gdal_translate(-of GTiff -co profile=baseline -co tfw=yes geotiff.tif baseline.tif)

tfw_file = sub("tif","tfw",outfile)
lines = as.numeric(readLines(tfw_file))
length(lines)

angle = 15
rad = angle*(pi/180)
lines[1] = lines[1]*cos(rad)
lines[2] = lines[1]*sin(rad)
lines[3] = -lines[4]*sin(rad)
lines[4] = lines[4]*cos(rad)
lines = as.character(lines)
writeLines(lines, tfw_file)


outfile1 = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132AAA04_B4_temp2.tif"
cmd = paste("gdal_translate -of GTiff",outfile, outfile1)
shell(cmd)

################################################################################################################
################################################################################################################
################################################################################################################


proj = projection(raster(infile))
proj = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84"

wktfile = sub("archv.tif","wkt.txt", fixfile)
projcmd = paste("gdalsrsinfo -o wkt", infile)
proj = system(projcmd, intern = TRUE)
write(proj, wktfile)


gcps = paste("-gcp", 2745, 983, 649170, 4828770, "-gcp", 790, 2439, 524870, 4741410, "-gcp", 2277, 2792, 618090, 4720230) #

outfile2 = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132AAA04_B4_temp3.tif"
gdaltrans_cmd = paste("gdal_translate -of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", proj, gcps, infile, outfile2)

gdaltrans_cmd = paste("gdal_translate", gcps, infile, outfile2)
system(gdaltrans_cmd)


outfile3 = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132AAA04_B4_temp4.tif"
gdalwarp_cmd = paste("gdalwarp -of Gtiff -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", 60, 60, outfile2, outfile3) #fixfile   "-tps"  "-order 2", "-order 3" 
system(gdalwarp_cmd)




#gdal translate command
tempname = sub("archv", "temp", fixfile) #"K:/scenes/034032/images/1976/LM10360321976248_archv.tif" 
gdaltrans_cmd = paste("gdal_translate -of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", wktfile, fixfile, tempname, gcpstr)
system(gdaltrans_cmd)



library(raster)
fixfiles = c("J:/l1g_warp/L1G_targz/LM40450301983103AAA03/LM40450301983103_archv.tif",
             "J:/l1g_warp/L1G_targz/LM40450301992176AAA03/LM40450301992176_archv.tif"
             )
reffile = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132_archv.tif"
outdir = "J:/l1g_warp/imgs"

make_tiepoint_img = function(imgfile,outdir){
  img = brick(imgfile)
  rows = nrow(img)
  cols = ncol(img)
  newrows = ceiling(0.4 * rows)
  newcols = ceiling(0.4 * cols)
  print(paste("rows:",newrows))
  print(paste("cols:",newcols))
  newbname = sub("archv.tif","archv.png",basename(imgfile))
  outpng = file.path(outdir,newbname)
  png(outpng, width = newcols, height=newrows)
  plotRGB(img,r=4, g=2, b=1, ext=NULL, maxpixels=10000000, stretch=NULL) #"lin" or "hist"
  dev.off()
  return(data.frame(origFile = imgfile, imgFile = newbname, origRows = rows, origCols = cols, imgRows = newrows, imgCols = newcols))
}

loop_tie_points = function(reffile, reffiles, outdir){
  refInfo = make_tiepoint_img(reffile,outdir)
  for(i in 1:length(fixfiles)){
    if(i == 1){
      fixInfo = make_tiepoint_img(fixfiles[i],outdir)
    } else{
      fixInfo = rbind(fixInfo, make_tiepoint_img(fixfiles[i],outdir))
    }
  }
  return(list(refInfo=refInfo,fixInfo=fixInfo))
}

info = loop_tie_points(reffile,fixfiles,outdir)

list(ref=list(file="test",point=c(list("orig"))))




writeplot = function(file,plot,gmBounds,date,tcb,tcg,tcw,tca, first=F,last=F, end=F, finalPlot=F){
  chipstrip = paste('"imgs/plot_',plot,'_chipstrip.png",',sep="")
  start=paste('{"plotID": ',plot,',','"chipStrip":',chipstrip,'"LatLon":[',gmBounds[1,2],',',gmBounds[1,1],'],','"bounds":[',gmBounds[2,2],',',gmBounds[3,2],',',gmBounds[3,1],',',gmBounds[2,1],'],','"Values": [')
  int=', '
  if(end == T & finalPlot == F){end=']},'} else if(end == T & finalPlot == T){end=']}'}
  imginfo = paste(
    paste('{"Year": ', date,',',sep=""),
    paste(' "TCB": ', tcb,',',sep=""),
    paste(' "TCG": ', tcg,',',sep=""),
    paste(' "TCW": ', tcw,',',sep=""),
    paste(' "TCA": ', tca,sep=""),
    '}',
    sep="")
  if(first == T){
    write(start, file, append=TRUE)
    write(paste(imginfo,int,sep=""), file, append=TRUE)
  }
  if(first != T & last != T){write(paste(imginfo,int,sep=""), file, append=TRUE)}
  if(last == T){
    write(imginfo, file, append=TRUE)
    write(end, file, append=TRUE)
  }
}




















