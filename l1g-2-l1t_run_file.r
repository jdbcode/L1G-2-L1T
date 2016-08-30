



# Unpack MSS L1G images ---------------------------------------------------------------

dir = "J:/l1g_warp/test/refimg" #full path to a directory of MSS L1G tar.gz files
proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #provide a proj.4 projection string

#######################################################################################
l1g2l1t_unpack_l1g(dir, proj)
#######################################################################################




# Unpack MSS L1T reference image ------------------------------------------

file = "J:/lcms_l1g2l1t/refimg/" #full path to an MSS L1T reference image
proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #provide a proj.4 projection string

#######################################################################################
l1g2l1t_unpack(file, proj)
#######################################################################################




# prepare L1G-2-L1T Initializer app images --------------------------------

reffile = #full path to an MSS L1T reference *archv.tif file
fixdir =  "J:/l1g_warp/test/fiximg" #full direcory path of L1G image outputs from the MSS unpacker function - must correspond to the wrs and path/row of the reference image
outdir = "J:/l1g_warp/test/tiepoints" #provide a full directory path to a folder where you want the application images and javascript file to go.
  
#######################################################################################
run_prepare_tie_point_images(reffile, fixdir, outdir)
#######################################################################################




# run the warping procedure -----------------------------------------------

tpfile = "J:/single_resampe_no_rotation/wrs1_012028_tie_points.json" #define the full path to the tie-point file that was downloaded from the MSS L1G-2-L1T Initializer app

#######################################################################################
run_l1g2l1t_warp(tpfile) #run the wapring script script - it will loop through all images included in the tie-point file
#######################################################################################





