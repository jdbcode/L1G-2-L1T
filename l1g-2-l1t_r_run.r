


####################################################
#MSS L1G-2-L1T warp steps
####################################################




####################################################
### MSS unpacker

dir = "J:/l1g_warp/L1G_targz" #provide a directory of MSS L1G tar.gz files 

###
files = list.files(dir, "tar.gz$", full.names = T)
len = length(files)
for(i in 1:len){
  print(paste(i,"/",len,sep=""))
  l1g2l1t_unpack(files[i])
}
####################################################





####################################################
### Prepare images for the HTML L1G-2-L1T Tie Point Selector application

reffile = #provide a well registered, cloud-free reference file (full file path to "*archv.tif" file)
fixfiles = c() #provide a vector of L1G image outputs from the MSS unpacker function (full file paths to "*archv.tif" files) - must correspond to the wrs and path/row of the reference image
outdir = #provide a full directory path to an existing folder where you want the application images and javascript file to go.

prepare_tie_point_images(reffile, fixfiles, outdir)

#after this finishes running you need to copy the "l1g-2-l1t_tie_point_selector.html" file into the output folder you defined above - open the file and find tie points
#when you have finished, you'll be prompted to download a json file whose filename is the input to the next step
####################################################





####################################################
### Perform an initial warp of the L1G images

json = #full file path to the json file that was downloaded on completion of tie point finding

initial_warp(json)

#.........left to do is send the results into the MSSwarp function, the files that come out of here ("*archv_L1Gwarp.tif") are pretty good, if teh points were selected precisely
####################################################











