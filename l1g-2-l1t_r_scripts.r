

library(raster)
library(rjson)

make_tiepoint_img = function(imgfile,outdir){
  img = brick(imgfile)
  rows = nrow(img)
  cols = ncol(img)
  newrows = ceiling(1 * rows) #0.4
  newcols = ceiling(1 * cols) #0.4
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


initial_warp = function(json){
  data = fromJSON(file = json)
  len = length(data)
  for(i in 1:len){
    print(paste("working on file: ",i,"/",len))
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
    
    outname = sub("archv.tif", "archv_L1Gwarp.tif", fixFile)
    gdalwarp_cmd = paste("gdalwarp -of Gtiff -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tps -tr", 60, 60, tempname, outname) #fixfile   "-tps"  "-order 2", "-order 3" 
    system(gdalwarp_cmd)
    
    unlink(tempname)
  }
}


l1g2l1t_unpack = function(file){
  
  bname = basename(file)
  outdir = file.path(dirname(file),substr(bname,1,nchar(bname)-7))  
  untar(file, exdir=outdir, tar="internal")
  
  allfiles = list.files(outdir, full.names=T)
  tiffiles = allfiles[grep("TIF",allfiles)]
  otherfiles = allfiles[grep("TIF",allfiles, invert=T)]
  filebase = basename(tiffiles[1])
  
  imgid = substr(filebase, 1, 16)
  name = paste(imgid, "_archv.tif", sep = "")
  
  tempstack = file.path(outdir,paste(imgid, "_tempstack.tif", sep = ""))
  finalstack = file.path(outdir, name)
  
  junk = substr(filebase, 17,21)
  baseotherfiles = basename(otherfiles)
  for(h in 1:length(otherfiles)){baseotherfiles[h] = sub(junk, "", baseotherfiles[h])}
  newotherfiles =  file.path(outdir, baseotherfiles, fsep = .Platform$file.sep)
  
  ref = raster(tiffiles[1])
  s = stack(tiffiles[1],tiffiles[2],tiffiles[3],tiffiles[4])
  img = as.array(s)    
  b1bads = img[,,1]>1 #!=0
  b2bads = img[,,2]>1 #!=0
  b3bads = img[,,3]>1 #!=0
  b4bads = img[,,4]>1 #!=0
  bads = b1bads*b2bads*b3bads*b4bads
  
  img[,,1] = img[,,1]*bads
  img[,,2] = img[,,2]*bads
  img[,,3] = img[,,3]*bads
  img[,,4] = img[,,4]*bads
  
  cloudtest = img[,,1] > 130
  percent = round((length(which(cloudtest == T))/length(which((img[,,1]>0) == T)))*100, digits=4)
  outfile = sub("archv.tif", "cloud_cover.txt", finalstack)
  write(percent,file=outfile)
  
  s = setValues(s,img)    
  origproj = projection(s)
  s = as(s, "SpatialGridDataFrame")       
  writeGDAL(s, finalstack, drivername = "GTiff", options="INTERLEAVE=BAND", type = "Byte", mvFlag = 0) #, drivername = "GTiff"
  unlink(tiffiles)
}


