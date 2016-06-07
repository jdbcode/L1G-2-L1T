


####################################################
#MSS L1G-2-L1T warp steps
####################################################




####################################################
### MSS unpacker

dir = "J:/l1g_warp/test" #provide a directory of MSS L1G tar.gz files 
proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

###
files = list.files(dir, "tar.gz$", full.names = T)
len = length(files)
for(i in 1:len){
  print(paste(i,"/",len,sep=""))
  l1g2l1t_unpack(files[i], proj)
}
####################################################





####################################################
### Prepare images for the HTML L1G-2-L1T Tie Point Selector application

reffile = "J:/l1g_warp/test/LM50450301985132AAA04/LM50450301985132_archv.tif" #provide a well registered, cloud-free reference file (full file path to "*archv.tif" file)
fixfiles = c(
  "J:/l1g_warp/test/LM40450301982276AAA03/LM40450301982276_archv.tif",
  "J:/l1g_warp/test/LM40450301982292AAA01/LM40450301982292_archv.tif",
  "J:/l1g_warp/test/LM40450301982356AAA04/LM40450301982356_archv.tif",
  "J:/l1g_warp/test/LM40450301983039XXX03/LM40450301983039_archv.tif",
  "J:/l1g_warp/test/LM40450301983055AAA03/LM40450301983055_archv.tif"
  
) #provide a vector of L1G image outputs from the MSS unpacker function (full file paths to "*archv.tif" files) - must correspond to the wrs and path/row of the reference image
outdir = "J:/l1g_warp/test/tiepoints" #provide a full directory path to an existing folder where you want the application images and javascript file to go.

prepare_tie_point_images(reffile, fixfiles, outdir)

#after this finishes running you need to copy the "l1g-2-l1t_tie_point_selector.html" file into the output folder you defined above - open the file and find tie points
#when you have finished, you'll be prompted to download a json file whose filename is the input to the next step
####################################################





####################################################
### Perform an initial warp of the L1G images

json = "J:/l1g_warp/test/wrs2_045030_tie_points.json"  #full file path to the json file that was downloaded on completion of tie point finding

initial_warp(json)



#.........left to do is send the results into the MSSwarp function, the files that come out of here ("*archv_L1Gwarp.tif") are pretty good, if teh points were selected precisely
####################################################


fixfiles = list.files("J:/l1g_warp/test", "l1g_warp.tif", full.names=T, recursive=T)

len = length(fixfiles)
for(i in 1:len){
  print(paste(i,"/",len,sep=""))
  l1g2l1t_warp(reffile, fixfiles[i], window=275, search=27, sample=1000, refstart=c(0,0), fixstart=c(0,0))
}








