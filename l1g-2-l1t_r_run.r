


####################################################
#MSS L1G-2-L1T warp steps
####################################################




####################################################
### MSS unpacker

dir = "J:/l1g_warp/test/refimg" #provide a directory of MSS L1G tar.gz files 
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

refdir = "J:/l1g_warp/test/refimg" #provide a directory path to a reference image outputs from the MSS unpacker function
fixdir =  "J:/l1g_warp/test/fiximg" #provide a direcory path of L1G image outputs from the MSS unpacker function - must correspond to the wrs and path/row of the reference image
outdir = "J:/l1g_warp/test/tiepoints" #provide a full directory path to an existing folder where you want the application images and javascript file to go.



reffile = list.files(refdir, "archv.tif", full.names=T, recursive=T)
fixfiles = list.files(fixdir, "archv.tif", full.names=T, recursive=T)

prepare_tie_point_images(reffile, fixfiles, outdir)

#after this finishes running you need to copy the "l1g-2-l1t_tie_point_selector.html" file into the output folder you defined above - open the file and find tie points
#when you have finished, you'll be prompted to download a json file whose filename is the input to the next step
####################################################





####################################################
### Perform an initial warp of the L1G images

json = "J:/l1g_warp/test/tiepoints/wrs2_045030_tie_points.json"  #full file path to the json file that was downloaded on completion of tie point finding

initial_warp(json)

####################################################


dir = "J:/l1g_warp/test/fiximg"
reffile = "J:/l1g_warp/test/refimg/LM50450301985132AAA04/LM50450301985132_archv.tif"


fixfiles = list.files(dir, "l1g_warp.tif", full.names=T, recursive=T)
len = length(fixfiles)
for(i in 1:len){
  print(paste(i,"/",len,sep=""))
  l1g2l1t_warp(reffile, fixfiles[i])
}








