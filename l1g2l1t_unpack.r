

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


files = list.files("J:/l1g_warp/L1G_targz", "tar.gz$", full.names = T)
len = length(files)
for(i in 1:len){
  print(paste(i,"/",len,sep=""))
  l1g2l1t_unpack(files[i])
}


















