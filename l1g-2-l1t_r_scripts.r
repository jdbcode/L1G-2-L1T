

library(raster)
library(rjson)

fixfiles = c("J:/l1g_warp/L1G_targz/LM40450301983103AAA03/LM40450301983103_archv.tif",
             "J:/l1g_warp/L1G_targz/LM40450301992176AAA03/LM40450301992176_archv.tif"
)
reffile = "J:/l1g_warp/LM50450301985132AAA04/LM50450301985132_archv.tif"
outdir = "C:/git_repos/L1G-2-L1T"

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

prepare_tie_point_images = function(reffile, fixfiles, outdir){
  info = loop_tie_points(reffile,fixfiles,outdir)
  js = toJSON(info)
  js = paste("info =",js)
  outfile = file.path(outdir,"imgInfo.js")
  write(js, outfile)
}


get_gcp = function(d,p){
  r = raster(d$ref$file)
  refX = xFromCol(r, d$ref$point[[p]]$col)
  refY = yFromRow(r, d$ref$point[[p]]$row)
  fixCol = floor(d$fix$point[[p]]$col)
  fixRow = floor(d$fix$point[[p]]$row)
  gcp = paste("-gcp", fixCol, fixRow, refX, refY)
  return(gcp)    
}


json = "C:/Users/braatenj/Downloads/data.json"



data = fromJSON(file = json)
length(data)
for(i in 1:length(data)){
  i=1
  d= data[[i]]
  fixFile = d$fix$file
  for(p in 1:length(d$fix$point)){
    gcp = get_gcp(d,p)
    if(p == 1){gcps = gcp} else {gcps = paste(gcps, gcp)}
  }
  print(gcps)
  
  wktfile = sub("archv.tif","wkt.txt", fixFile)
  projcmd = paste("gdalsrsinfo -o wkt", fixFile)
  proj = system(projcmd, intern = TRUE)
  write(proj, wktfile)
  
  tempname = sub("archv", "temp", fixFile) #"K:/scenes/034032/images/1976/LM10360321976248_archv.tif" 
  gdaltrans_cmd = paste("gdal_translate -of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", wktfile, fixFile, tempname, gcps)
  system(gdaltrans_cmd)
  
  tempname1 = sub("temp", "temp1", tempname)
  gdalwarp_cmd = paste("gdalwarp -of Gtiff -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", 60, 60, tempname, tempname1) #fixfile   "-tps"  "-order 2", "-order 3" 
  system(gdalwarp_cmd)
  
  
}

