



# Unpack MSS L1G images ---------------------------------------------------------------

dir = "J:/l1g2l1t_test/wrs1_012028/fiximg" #full path to a directory of MSS L1G tar.gz files
proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #provide a proj.4 projection string

#######################################################################################
l1g2l1t_unpack_l1g(dir, proj)
#######################################################################################




# Unpack MSS L1T reference image ------------------------------------------

file = "J:/l1g2l1t_test/wrs1_012028/refimg/LM20120281977191GMD04.tar.gz" #full path to an MSS L1T reference image
proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #provide a proj.4 projection string

#######################################################################################
l1g2l1t_unpack(file, proj)
#######################################################################################




# Prepare L1G-2-L1T initializer app images --------------------------------

reffile = "J:/l1g2l1t_test/wrs1_012028/refimg/LM20120281977191GMD04/LM20120281977191_archv.tif" #full path to an MSS L1T reference *archv.tif file
fixdir =  "J:/l1g2l1t_test/wrs1_012028/fiximg" #full direcory path of L1G image outputs from the MSS unpacker function - must correspond to the wrs and path/row of the reference image
outdir = "J:/l1g2l1t_test/wrs1_012028/itp" #provide a full directory path to a folder where you want the application images and javascript file to go.
  
#######################################################################################
run_prepare_tie_point_images(reffile, fixdir, outdir)
#######################################################################################




# Use the L1G-2-L1T initializer app to find starting tie points ------------

# copy the "l1g-2-l1t_initializer.html" file into the outdir folder defined in the above section. Open the HTML file and it will
# hitch up to a JavaScript file and image files in that folder. Follow the app instructions that are provided when the program opens.




# Run the warping procedure -----------------------------------------------

tpfile = "J:/l1g2l1t_test/wrs1_012028/wrs1_012028_tie_points.json" #define the full path to the tie-point file that was downloaded from the MSS L1G-2-L1T Initializer app
method = "tps" #select a method to warp the L1G images - the options are: "tps" (thin plate spline), "order 1" (first order polynomial), or "order 2" (second order polynomial)

#######################################################################################
run_l1g2l1t_warp(tpfile, method="tps") #run the wapring script script - it will loop through all images included in the tie-point file
#######################################################################################




# Check RMSE --------------------------------------------------------------

rmse_file = "J:/l1g2l1t_test/wrs1_012028/fiximg/LM10120281976224PAC08/LM10120281976224_rmse.Rdata" #for a single image define the path its "*_rmse.Rdata" file

#######################################################################################
print_rmse(rmse_file)
#######################################################################################


