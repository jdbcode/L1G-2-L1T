

library(raster)
library(rgdal)
library(gdalUtils)
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
  plotRGB(img,r=4, g=2, b=1, ext=NULL, maxpixels=10000000, stretch="lin", colNA="black") #"lin" or "hist" or NULL
  dev.off()
  return(data.frame(origFile = imgfile, imgFile = newbname, origRows = rows, origCols = cols, imgRows = newrows, imgCols = newcols))
}


loop_tie_points = function(reffile, reffiles, outdir){
  print("preparing reference image...")
  refInfo = make_tiepoint_img(reffile,outdir)
  len = length(fixfiles)
  print("preparing fix image(s)...")
  for(i in 1:len){
    print(paste("  ",i,"/",len,sep=""))
    if(i == 1){
      fixInfo = make_tiepoint_img(fixfiles[i],outdir)
    } else{
      fixInfo = rbind(fixInfo, make_tiepoint_img(fixfiles[i],outdir))
    }
  }
  return(list(refInfo=refInfo,fixInfo=fixInfo))
}


prepare_tie_point_images = function(reffile, fixfiles, outdir){
  dir.create(outdir, showWarnings = F, recursive = T)
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
    d = data[[i]]
    fixFile = d$fix$file
    if(d$process == 0){
      dname = dirname(fixFile)
      dnamenew = paste(dname,"_NO_L1G2L1T", sep="")
      file.rename(dname,dnamenew)
      outfile = file.path(dnamenew,sub("archv.tif", "no_l1g2l1t.txt",basename(fixFile)))
      write(d$note, outfile)
      next()
    }
    
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
    
    outname = sub("archv.tif", "archv_l1g_warp.tif", fixFile)
    gdalwarp_cmd = paste("gdalwarp -of Gtiff -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tps -tr", 60, 60, tempname, outname) #fixfile   "-tps"  "-order 2", "-order 3" 
    system(gdalwarp_cmd)
    
    origsize = file.info(fixFile)$size
    warpsize = file.info(outname)$size
    ratio = warpsize/origsize
    if(ratio > 1.5 | ratio < 0.5){
      dname = dirname(fixFile)
      dnamenew = paste(dname,"_NO_L1G2L1T", sep="")
      file.rename(dname,dnamenew)
      outfile = file.path(dnamenew,sub("archv.tif", "no_l1g2l1t.txt",basename(fixFile)))
      write("weird warp - selected points may be incorrect", outfile)
      file.path(dnamenew,basename(tempname))
      unlink(c(file.path(dnamenew,basename(tempname)),file.path(dnamenew,basename(outname))))
      next()
    }
    
    unlink(tempname)
  }
}


trim_na_rowcol = function(imgfile, outimg, maskfile, outmask){
  #trim the stack
  x = raster(imgfile, band=1)
  cres = 0.5*res(x)
  #crs <- projection(x) 
  y = x
  x = matrix(as.array(x),nrow=nrow(x),ncol=ncol(x))
  r.na = c.na <- c()
  for(i in 1:nrow(x)) r.na <- c(r.na, all(is.na(x[i,])))
  for(i in 1:ncol(x)) c.na <- c(c.na, all(is.na(x[,i])))
  r1 = 1 + which(diff(which(r.na))>1)[1] 
  r2 = nrow(x) -  which(diff(which(rev(r.na)))>1)[1]
  c1 = 1 + which(diff(which(c.na))>1)[1] 
  c2 = ncol(x) - which(diff(which(rev(c.na)))>1)[1]
  #x = x[r1:r2,c1:c2]
  
  #if there are no NA rows and cols, then set to default row and col start and end
  if(is.na(r1)){
    if(r.na[1] == T){r1 = 1+length(which(r.na))} else {r1 = 1}
  }
  if(is.na(r2)){
    if(rev(r.na)[1] == T){r2 = nrow(x) - length(which(r.na))} else {r2 = nrow(x)}
  }
  if(is.na(c1)){
    if(c.na[1] == T){c1 = 1+length(which(c.na))} else {c1 = 1}
  }
  if(is.na(c2)){
    if(rev(c.na)[1] == T){c2 = ncol(x) - length(which(c.na))} else {c2 = ncol(x)}
  }
  
  xs = xFromCol(y,col=c(c1,c2)) + c(-1,1)*cres[1]
  ys = yFromRow(y,row=c(r2,r1)) + c(-1,1)*cres[2]
  
  #write out the trimmed file
  gdal_translate(src_dataset=imgfile, dst_dataset=outimg, of="GTiff", co="INTERLEAVE=BAND", projwin=c(xs[1],ys[2],xs[2],ys[1]))
  if(file.exists(maskfile) == T){gdal_translate(src_dataset=maskfile, dst_dataset=outmask, of="GTiff", co="INTERLEAVE=BAND", projwin=c(xs[1],ys[2],xs[2],ys[1]))}
}


l1g2l1t_unpack = function(file, proj){
  bname = basename(file)
  outdir = file.path(dirname(file),substr(bname,1,nchar(bname)-7))  
  untar(file, exdir=outdir, tar="internal")
  
  
  allfiles = list.files(outdir, full.names=T)
  tiffiles = allfiles[grep("TIF",allfiles)]
  otherfiles = allfiles[grep("TIF",allfiles, invert=T)]
  filebase = basename(tiffiles[1])
  
  imgid = substr(filebase, 1, 16)
  
  tempstack = file.path(outdir,paste(imgid, "_tempstack.tif", sep = ""))
  projstack = sub("tempstack", "projstack", tempstack)
  finalstack = file.path(outdir, paste(imgid, "_archv.tif", sep = ""))
  outprojfile = sub("archv.tif", "proj.txt", finalstack)
  
  junk = substr(filebase, 17,21)
  baseotherfiles = basename(otherfiles)
  for(h in 1:length(otherfiles)){baseotherfiles[h] = sub(junk, "", baseotherfiles[h])}
  newotherfiles =  file.path(outdir, baseotherfiles, fsep = .Platform$file.sep)
  
  file.rename(otherfiles, newotherfiles)
  file.rename(file, file.path(outdir,basename(file)))
  
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
  writeGDAL(s, tempstack, drivername = "GTiff", options="INTERLEAVE=BAND", type = "Byte", mvFlag = 0)
  
  
  write(proj, outprojfile)
  gdalwarp(srcfile=tempstack, dstfile=projstack, 
           s_srs=origproj, t_srs=proj, of="GTiff",
           r="near", srcnodata=0, dstnodata=0, multi=T,
           tr=c(60,60), co="INTERLEAVE=BAND")
  
  trim_na_rowcol(projstack, finalstack, "null", "null") 
  
  deleteThese = c(tempstack, projstack, tiffiles)
  unlink(deleteThese)
}


l1g2l1t_warp = function(reffile, fixfile, mode){
  
  #set default parameters
  search=35 #27 
  if(mode == "rmse"){search = 27} #it should be pretty close, no need to look over a larger region - keep at 27 since that is what the thresholds were selected at.
  
  #scale image values to center on mean and 1 unit variance (global mean and variance)
  scaleit = function(matrix){
    stnrd = (matrix - (mean(matrix, na.rm = TRUE)))/(sd(matrix, na.rm = TRUE))
    return(stnrd)
  }
  
  #make a kernal around a given point
  make_kernal = function(img, point1, windowsize){
    radius = floor(windowsize/2)
    ccol = colFromX(img, point1[1])
    crow = rowFromY(img, point1[2])
    mincol = ccol-radius
    maxcol = ccol+radius
    minrow = crow+radius
    maxrow = crow-radius
    return(extent(c(mincol,maxcol,minrow,maxrow)))
  }
  
  calc_rmse = function(info,reso){
    xresid = (info[,"refx"]-info[,"fixx"])^2 #get the residuals of each x
    yresid = (info[,"refy"]-info[,"fixy"])^2 #get the residuals of each y
    r = (sqrt(xresid+yresid))/reso #get the rmse of each xy point
    x_rmse = sqrt(mean(xresid))/reso
    y_rmse = sqrt(mean(yresid))/reso
    total_rmse = sqrt((x_rmse^2)+(y_rmse^2)) #total rmse including all points
    rmse_info = list(x_rmse=x_rmse, y_rmse=y_rmse, total_rmse=total_rmse, r=r)
    return(rmse_info)
  }
  
  make_gdaltrans_cmd = function(info, wktfile, fixfile, tempname){
    info[,"refx"] = info[,"refx"]+(reso/2)
    info[,"refy"] = info[,"refy"]-(reso/2)
    fixcol = paste(info[,"fixcol"]) #fix col for tie point
    fixrow = paste(info[,"fixrow"]) #fix row for tie point
    refx = paste(round(info[,"refx"]))  #fix x tie point coord 
    refy = paste(round(info[,"refy"]))  #fix y tie point coord
    gcpstr = paste(" -gcp", fixcol, fixrow, refx, refy, collapse="")
    gdaltrans_cmd = paste("gdal_translate -of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", wktfile, fixfile, tempname, gcpstr)
    return(gdaltrans_cmd)
  }
  
  #read in the fix image
  fiximg = raster(fixfile, band=3) #load the fix image
  origfiximg = fiximg #save a copy of the original fix image
  fiximgb1 = raster(fixfile, band=1)
  
  #load the ref image
  refimg = raster(reffile, 3) 
  
  #make sure that the ref and fix img are croppd to eachother
  refimg = intersect(refimg, fiximg)
  fiximg = intersect(fiximg, refimg)
  
  #calculate similarity index input values from the fix image subset
  values(fiximg) = scaleit(values(fiximg))
  values(fiximg)[is.na(values(fiximg))] = 0
  
  #calculate similarity index input values from the ref image subset
  values(refimg) = scaleit(values(refimg))
  values(refimg)[is.na(values(refimg))] = 0
  refimgsqr = refimg^2
  
  #get the resolution  
  reso = xres(refimg)
  
  #adjust the window and search size so that they are odd numbers
  if (search %% 2 == 0){search = search+1}
  #radius = floor(window/2) #radius of the window in pixels
  #nrc = search+(radius*2) #the reference extent length to slide over
  
  #sample the the reference image, laying down a regular grid of points to check
  #s = sampleRegular(refimg, sample, cells=T)
  #s = sampleRandom(refimg, 10000, cells=T)
  s = sampleRegular(refimg, 10000, cells=T)
  s = s[,1]
  xy = xyFromCell(refimg,s) #[,1] #get the xy coordinates for each good point
  
  #filter points in fiximg that fall on clouds
  theseones = cellFromXY(fiximgb1, xy)   #get fiximg cell index for sample 
  theseones = na.omit(theseones)
  a = fiximgb1[theseones] #extract values for fiximg cell sample
  b = which(fiximgb1[theseones] < 100) # fiximgb1[theseones] != NA) #exclude points that don't meet criteria
  
  #if the number of sample points is less than 10 delete the image return
  if(length(b) < 10){
    delete_files(fixfile, 2)
    return(0)
  }
  
  #subset the sample
  n_subsamp = 1000
  if(length(b) < n_subsamp){n_subsamp = length(b)}
  subsamp = sample(b, n_subsamp)
  xy = xy[subsamp,] #subset the original sample
  s = s[subsamp] #subset the original sample
  rowcol = rowColFromCell(refimg, s)
  
  #make an array to hold all information collected
  info = cbind(c(0),xy, rowcol[,2], rowcol[,1], array(0, c(length(xy[,1]), 8)))
  cnames = c("point","refx","refy","refcol","refrow","fixx","fixy","fixcol","fixrow","nmax", "max","edgedist","decision")
  colnames(info) = cnames
  
  #start a pdf file to hold image chips and similarity surfaces
  #pdf_file = sub("archv_l1g_warp.tif", "ccc_surface.pdf",fixfile)
  #unlink(pdf_file) #delete the pdf if it exists
  #pdf(file=pdf_file, width=10, heigh=7) #size of the pdf page
  #par(mfrow=c(2,3)) #number of trajectories to place on a page (columns, rows)
  
  #iterate process of creating a similarity surface for each check point in the sample
  window_size = c(101,201,275)
  for(size in 1:3){
    if(mode != "rmse"){
      if(size == 1){pdf_file = sub("archv_l1g_warp.tif", "ccc_surface_100w.pdf",fixfile)}
      if(size == 2){pdf_file = sub("archv_l1g_warp.tif", "ccc_surface_200w.pdf",fixfile)}
      if(size == 3){pdf_file = sub("archv_l1g_warp.tif", "ccc_surface_275w.pdf",fixfile)}
      unlink(pdf_file) #delete the pdf if it exists
      pdf(file=pdf_file, width=10, heigh=7) #size of the pdf page
      par(mfrow=c(2,3)) #number of trajectories to place on a page (columns, rows)
    }
    window = window_size[size]
    
    #adjust the window and search size so that they are odd numbers
    if (window %% 2 == 0){window = window+1}
    radius = floor(window/2) #radius of the window in pixels
    nrc = search+(radius*2) #the reference extent length to slide over
    
    print(paste("n goods =",length(which(info[,"decision"] == 1))))
    for(point in 1:length(info[,1])){ 
      #print(point) #print the point so we know where we're at
      if(info[point,"decision"] == 1){
        print("already good, skipping...")
        next
      }
      if(size == 1){info[point,"point"] = point} #info[point,1] = point #put the point number into the info table
      
      #make a subset of the reference image for the fiximg chip to slide over
      a = make_kernal(refimg, info[point,2:3], nrc)
      test = c(a@ymax,a@ymin,a@xmin,a@xmax)
      
      if(sum(is.na(test)) > 0){next}
      if(sum(test < 0) > 0){next}
      if(a@ymax > nrow(refimg) | a@ymin > nrow(refimg)){next}
      if(a@xmax > ncol(refimg) | a@xmin > ncol(refimg)){next}
      ext=extent(refimg,a@ymax,a@ymin,a@xmin,a@xmax)
      refsub = crop(refimg, ext)
      
      #make subset of fiximg (fiximg chip)
      a = make_kernal(fiximg, info[point,2:3], window)
      test = c(a@ymax,a@ymin,a@xmin,a@xmax)
      if(sum(is.na(test)) > 0){next}
      if(sum(test < 0) > 0){next}
      if(a@ymax > nrow(fiximg) | a@ymin > nrow(fiximg)){next}
      if(a@xmax > ncol(fiximg) | a@xmin > ncol(fiximg)){next}
      ext=extent(fiximg,a@ymax,a@ymin,a@xmin,a@xmax)
      fixsub = crop(fiximg, ext)
      
      #create numerator
      tofix = matrix(values(fixsub),ncol=window,byrow = T)
      
      if(length(tofix) %% 2 == 0) {
        print("skipping")
        next
      }
      
      num = focal(refsub, w=tofix ,fun=sum)
      
      #get refimg denom
      a = make_kernal(refimgsqr, info[point,2:3], nrc)
      ext=extent(refimgsqr,a@ymax,a@ymin,a@xmin,a@xmax)
      refsubsqr = crop(refimgsqr, ext)
      sumrefsubsqr = focal(refsubsqr, w=matrix(1,window, window)) #get the summed product of the refsubimg
      sumfixsubsqr = sum(values(fixsub)^2) #fiximg standard only gets calcuated once
      denom = sqrt(sumfixsubsqr*sumrefsubsqr)
      
      #badone=0
      if(cellStats(num, stat="sum") + cellStats(denom, stat="sum") == 0){next} 
      
      ncc = num/denom
      buf = (nrow(ncc)-search)/2
      off1 = buf+1
      off2 = buf+search
      ext = extent(ncc,off1,off2,off1,off2)
      ncc = crop(ncc,ext)
      nccv = values(ncc)
      nccm = matrix(nccv, ncol=sqrt(length(nccv)), byrow=T)
      
      x = y = seq(1,ncol(nccm),1)
      
      good = which(values(ncc) == maxValue(ncc))[1] 
      #good = good[1]
      
      ####
      xoffsetcoord = xFromCell(ncc, good)
      yoffsetcoord = yFromCell(ncc, good)
      xoffset = xoffsetcoord - info[point,"refx"]
      yoffset = yoffsetcoord - info[point,"refy"]
      info[point,"fixx"] = xoffsetcoord-(xoffset*2)
      info[point,"fixy"] = yoffsetcoord-(yoffset*2)
      ####
      
      #get the row and column numbers for the fix image
      origfiximg_x = info[point,"fixx"]
      origfiximg_y = info[point,"fixy"]
      #a = cellFromXY(origfiximg, c(info[point,"refx"],info[point,"refy"])) #fiximg
      a = cellFromXY(origfiximg, c(origfiximg_x,origfiximg_y))
      fiximgrc = rowColFromCell(origfiximg, a)
      info[point,"fixcol"] = fiximgrc[2] #info[point,8] = fiximgrc[2]
      info[point,"fixrow"] = fiximgrc[1] #info[point,9] = fiximgrc[1]
      
      
      #############################################################
      #screen by outlier
      ccc = scaleit(nccm)  
      maxmat = max(ccc, na.rm = T)
      rowcol = which(ccc == maxmat, arr.ind=T)
      r = ccc[rowcol[1],] 
      c = ccc[,rowcol[2]] 
      rmax = which(r == maxmat) #max(r)
      cmax = which(c == maxmat) #max(c)
      
      dist1 = abs(c((rmax - 1), (rmax - length(r)), (cmax - 1), (cmax - length(c))))/(floor(search/2)) 
      dist2 = min(dist1)
      
      #place filtering values in info table
      info[point,"nmax"] = length(which(ccc == maxmat)) #2 #info[point,10] = length(which(ccc == maxmat)) #2
      info[point,"max"] = round(maxmat, digits=1) #3 #info[point,11] = round(maxmat, digits=1) #3
      info[point,"edgedist"] = round(dist2, digits=2) #6 #info[point,14] = round(dist2, digits=2) #6
      
      #decide what surfaces are good\bad
      bad = array(0,3)
      bad[1] = info[point,"nmax"] > 1 #number of max peaks eq 1
      bad[2] = info[point,"max"] < 3 #peak ge to 3 standard devs from mean
      bad[3] = info[point,"edgedist"] < 0.12 #peak distance from edge >= 0.12
      info[point,"decision"] = sum(bad) == 0 #7
      #if(badone == 1){info[point,"decision"] = 0}
      
      #filter plots that will crash the plotting because of weird data points (na, NaN, (-)Inf)
      bad1 = is.na(ccc)
      bad2 = is.infinite(ccc) 
      bad3 = is.nan(ccc)
      badindex = which(bad1 == T | bad2 == T | bad3 == T)
      ccc[badindex] = 0
      if(length(which(ccc == 0)) == length(ccc)){next}
      
      #plot the cross correlation surface
      #title = paste(point,info$nmax[point],info$max[point],info$edgedist[point], info$decision[point], sep = ",")
      if(info[point,"decision"] == 1){status = "accept"} else {status = "reject"}
      title = paste("Point:", point, status)
      print(title)
      if(mode != "rmse"){
        ccc = ccc[nrow(ccc):1,]
        persp(x, y, ccc, theta = 30, phi = 30, expand = 0.5, col = 8, main=title)
      }
    }
    
    if(mode != "rmse"){
      dev.off() #turn off the plotting device
    }
  }
  print(paste("n goods =",length(which(info[,"decision"] == 1))))
  
  #write all the point info to a file
  if(mode != "rmse"){
    info_file = sub("archv_l1g_warp.tif", "info_full.csv",fixfile)
    write.csv(info, file=info_file, row.names = F) 
  }
  
  #subset the points that passed the surface tests
  these = which(info[,"decision"] == 1)
  
  
  rmse_outfile = file.path(dirname(fixfile),paste(substr(basename(fixfile),1,16),"_rmse.Rdata",sep=""))
  rmse_info = list(calc_rmse = F, x_rmse=NA, y_rmse=NA, total_rmse=NA, info=info)  
  
  
  
  #if the number of sample points is less than 10 delete the image and return
  if(length(these) < 10){
    if(mode != "rmse"){
      dname = dirname(fixfile)
      newdname = paste(dname,"_NO_L1G2L1T", sep="")
      file.rename(dname,newdname)
      outfile = sub("archv_l1g_warp.tif", "no_l1g2l1t.txt",fixfile)
      write("cant find enought points", outfile)
    } else{
      save(rmse_info, rmse_outfile)
    }
    return()
  }

  
  info = info[these,]
  
  #filter points based on rmse contribution
  if(mode != "rmse"){
    rmse = calc_rmse(info,reso)
    r = rmse$r
    sdr = sd(r)
    meanr = mean(r)
    limit = meanr+sdr*2
    goods = which(r <= limit)
    n_outliers = nrow(info)-length(goods)
    info = info[goods,]
    print(paste("getting rid of:",n_outliers,"outliers"))
    print(paste("there are still:",nrow(info),"points"))
    #maxr = endit = 10
    #while(maxr >2 & endit != 0){
    #  rmse = calc_rmse(info,reso)
    #  if (rmse$total_rmse != 0){contr = rmse$r/rmse$total_rmse} else contr = rmse$r #error contribution of each point
    #  maxr = max(contr) #while loop controler
    #  b = which(contr < 2) #subset finder - is point 2 times or greater in contribution
    #  info = info[b,] #subset the info based on good rsme
    #  endit = sum(contr[b]) #while loop controler
    #}
  }
    
  #if the number of sample points is less than 10 delete the image and return
  if(nrow(info) < 10){
    if(mode != "rmse"){
      dname = dirname(fixfile)
      newdname = paste(dname,"_NO_L1G2L1T", sep="")
      file.rename(dname,newdname)
      outfile = sub("archv_l1g_warp.tif", "no_l1g2l1t.txt",fixfile)
      write("cant find enought points", outfile)
    } else{
      rmse_info$info = info
      save(rmse_info, rmse_outfile)
    }
    return()
  }
  
  #if this is an rmse run, then save the info and get out
  if(mode == "rmse"){
    rmse = calc_rmse(info,reso)
    
    rmse_info$calc_rmse = T
    rmse_info$x_rmse=rmse$x_rmse
    rmse_info$y_rmse=rmse$y_rmse
    rmse_info$total_rmse=rmse$total_rmse
    rmse_info$info=info
    
    save(rmse_info, file=rmse_outfile)
    return()
  }
  
  
  #write out the filtered points that will be used in the transformation
  info_file = sub("archv_l1g_warp.tif", "info_sub.csv",fixfile)
  write.csv(info, file=info_file, row.names = F)
  
  
  
  
  
  
  #make some output file names
  tempname = sub("archv_l1g_warp", "temp", fixfile) #"K:/scenes/034032/images/1976/LM10360321976248_archv_l1g_warp.tif"
  outfile = sub("archv_l1g_warp.tif", "archv_l1g2l1t.tif", fixfile)
  wktfile = sub("archv_l1g_warp.tif","wkt.txt", fixfile)
  gcpfile = sub("archv_l1g_warp.tif", "gcp.txt", fixfile)
  
  #write out a project file for gdal translate to use
  projcmd = paste("gdalsrsinfo -o wkt", fixfile)
  proj = system(projcmd, intern = TRUE)
  write(proj, wktfile)
  
  #get the initial gcp string made
  gdaltrans_cmd = make_gdaltrans_cmd(info, wktfile, fixfile, tempname)
  
  #make sure that the gcpstr command is small enough to make it through the command prompt
  nchar_gcpstr = nchar(gdaltrans_cmd)
  counter = 1
  subsamp = F
  while(nchar_gcpstr > 8191){
    if(counter == 1){subsamp = T}
    subinfo = info[sample(1:nrow(info),(nrow(info)-counter)),]
    gdaltrans_cmd = make_gdaltrans_cmd(subinfo, wktfile, fixfile, tempname)
    nchar_gcpstr = nchar(gdaltrans_cmd)
    counter = counter + 1
  }
  
  #write out subset sample if the subset was sampled
  if(subsamp == T){
    info_sub_file_ = sub("archv_l1g_warp.tif", "info_sub_samp.csv",fixfile)
    write.csv(subinfo, file=info_sub_file_, row.names = F)
  }
  
  #write out the gcp metadata file
  write(gdaltrans_cmd, file=gcpfile) #, append=T)
  
  #gdal translate command
  #gdaltrans_cmd = paste("gdal_translate -of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", wktfile, fixfile, tempname, gcpstr)
  system(gdaltrans_cmd)
  
  #gdal warp command
  gdalwarp_cmd = paste("gdalwarp -of Gtiff -tps -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", reso, reso, tempname, outfile) #fixfile   "-tps"  "-order 2", "-order 3" 
  system(gdalwarp_cmd)
  
  #delete the temp file
  unlink(list.files(dirname(fixfile), pattern = "temp", full.names = T))
}


#for testing
#reffile = "J:/l1g_warp/test/refimg/LM50450301985132AAA04/LM50450301985132_archv.tif"
#fixfile = "J:/l1g_warp/test/LM50450301985132_archv_l1g_warp.tif"
#shift(brick(reffile), x=360, y=360, format = "GTiff", filename=fixfile, datatype = "INT1U", options=c("COMPRESS=NONE"))
#l1g2l1t_warp(reffile, fixfile)

