# MSS L1G-2-L1T

## About
The L1G-2-L1T program uses image-to-image geo-registration to translate Landsat MSS L1G images to L1T images. L1G images are not precisily georeferenced and therefore are of little practical use, which is unfortunate since many MSS images with this processing status are otherwise decent images. L1G-2-L1T using an automated tie-point finding algorithm to find tie points between a given L1G image and an L1T image from the same WRS path-row footprint. It then uses GDAL to warp the L1G image based on the tie points that were found. Since L1G images are often many kilometers off, an interactive web-based initializer app is used to identify one coincident point between the L1T reference image and L1G image as a place for the automated process to start. The program is written in R and takes advantage of the R raster package and GDAL. It allows for bulk warping of many images from the same WRS path-row footprint.

Example of the JavaScript-based geo-registration initialization app [here](http://jdbcode.github.io/L1G-2-L1T/).

## Setup
##### Computer
The program was developed and tested on machines running Windows 7, 64-bit, with >= 8 GB RAM. 

##### Software
To run the program you need:
- [R](https://www.r-project.org/) (a free software environment for statistical computing and graphics)
- [RStudio](https://www.rstudio.com/) (a convenient frontend interface to the R environment)
- [GDAL](http://www.gdal.org/) (a program for reading, writing, and manipulating geospatial data)

Download and install R and RStudio using instructions in links above. To install GDAL use this example method:

1. Go to http://www.gisinternals.com/release.php
2. Click on the Downloads link for the version that best matches your system (we use MSVC 2010 - x64)
3. Download the Generic installer for the GDAL core components
4. Run the installer
5. Include GDAL in your system's environmental variable PATH
    1. Open Windows Control Panel and select System
    2. Click on Advanced system settings
    3. Click the Environmental Variables... button
    4. Under System variables, scroll down to the Path variable and click on it to highlight it
    5. Click the edit button
    6. Get your cursor to the end of the line, add a semi-colon (;) and add the path to the GDAL installation location. Example: C:\GDAL (this may not actually be the location on your system)

##### R libraries
You also need a few R libraries. Open RStudio and run the following line in the command prompt to install the dependencies.

<pre>install.packages(c("raster", "rgdal", "gdalUtils", "rjson"))</pre>

##### Images
MSS images should be downloaded from [Earth Explorer](http://earthexplorer.usgs.gov/). When you receive the images they will be compressed as .tar.gz files. You should retain a copy of these original compressed files, since that is what L1G-2-L1T will work with. For a given WRS path-row footprint that you are working on, you'll need a single, nearly cloud-free L1T MSS "reference" image with low positional RMSE to which MSS L1G images will be spatially matched to. Identify this image and download it. Then identify any L1G images that you want to include in your project. Remember that they must be from the same WRS path-row as the L1T reference image. The MSS L1G images can be cloudy, but maybe no more than 50% cloud cover. Download the L1G images.

##### Directory structure
You can organize your files however you like, but the following this is a recommended structure. Start by creating a folder somewhere for the WRS path-row you are working on: 

: somewhere\WRS_012028

within this path-row you are basically working with three kinds of data:

1. an L1T reference image (ref)
2. L1G images to be fixed (fix)
3. RGB images that will be created for use in initializer app (initializer)

Given this, you should create three sub-directories in your path-row folder

: somewhere\WRS_012028\
        : ref 
        : fix
        : initializer

In the *ref* folder place the compressed L1T reference image you downloaded. In the *fix* folder place the compressed L1G images you downloaded, and the files that go in the *initializer* folder will be populated automatically, execpt for the initializer HTML file (more on that later).

**With everything installed, downloaded, and files organized, you are ready to starting running the procedures***

## Running
You'll need three files that provided here on GitHub. To get them, click on the *Clone or Download* button on the upper right portion of this repository page and select *Download ZIP*. Decompress the folder and have a look inside. There is a script file (l1g-2-l1t_scripts.r) that will be sourced in RStudio, a "run file" (l1g-2-l1t_run_file.r) that can be used as a template for running all the procedures, and an HTML file (l1g-2-l1t_initializer.html) which launches the initializer app.

1. Open the two *.r* files in RStudio (double click them).
2. You need to source the script file - in RStudio, activate the *l1g-2-l1t_scripts.r* file by clicking its tab in the editor window and then click the *source* button on the upper-right side of the editor. This will load all the functions into the R environment.
3. Activate the *l1g-2-l1t_run_file.r* file by clicking on its tab in the editor window. Here is where you will define a couple variables for each procedure and then copy and paste each section into the RStudio command prompt to run it. Use the provided examples and variable definitions to fill them out correctly. Each section should be completed sequentially

#### Basic section description
- **Unpack MSS L1G images**: provide a full directory path to the L1G *fix* folder and an PROJ.4 projection string (see *Defining and PROJ.4 projection* section below) as inputs to the *run_prepare_tie_point_images* function. It will decompress and stack all L1G images in the provided directory and write out the new files in the same directory, each image in its own sub-directory.
- **Unpack MSS L1T reference image**: provide a full directory path to a WRS path-row corresponding to the L1G images just defined above and an PROJ.4 projection string (see more on this below) as inputs to the *l1g2l1t_unpack* function. It will decompress and stack the L1T image in the provided directory and write out the new files in the same directory. Note that the PROJ.4 projection needs to be the same as the one defined in the *Unpack MSS L1G images* section.
- **Prepare L1G-2-L1T initializer app images**: Once the L1T and L1G images have been prepared, run this section to create data used in the *L1G-2-L1T initializer app*. You'll need to define the reference file that was unpacked (...ref/*archv.tif), the directory where the unpacked fix files are (.../fix), and the output directory where you want the *L1G-2-L1T initializer app* files to be written to.
- **Use the L1G-2-L1T initializer app to find starting tie points**: Copy and paste the *l1g-2-l1t_initializer.html* file that was included in the ZIP file you downloaded from this repository into the output folder you defined in the previous *Prepare L1G-2-L1T initializer app images* section. Double click the *l1g-2-l1t_initializer.html* file and it will open in your default web browser (only tested in FireFox and Chrome). It will hitch up to the PNG images and the JavaScript file in the output directory to populate the application. Follow the direction in the modal that is displayed upon opening the app.
- **Run the warping procedure**: Once you've finished finding initial tie points and have downloaded the JSON file from the *L1G-2-L1T initializer app* you are ready to run the warping procedure. Cut and paste the downloaded JSON file to the main path-row directory you are working on to keep it with everything else, then define its directory location as a variable input in the *Run the warping procedure* section of the *l1g-2-l1t_run_file*. Also define what method to use in warping the image, and then copy and paste the section into the RStudio command prompt to run it.
- **Check the RMSE of images**: an ancillary function is provided to print the RMSE of image files. Provide the full file path to a *_rmse.Rdata* file and run the section in the command prompt to print the x, y, and total RMSE for the image (filtered by 2 standard deviations from the mean offset). If you want to inspect these *_rmse.Rdata* files yourself, use the R `load()` function with the full file path as an input and then use `str(rmse_info)`  to see the list components.

## Files produced

- _archv.tif - stacked 4-band MSS image from USGS
- _archv_l1g2l1t.tif - L1G image that has been warped to match reference L1T image
- _ccc_surface_100w.pdf - cross correlation surfaces for the 100 x 100 pixel window size
- _ccc_surface_200w.pdf - cross correlation surfaces for the 200 x 200 pixel window size
- _ccc_surface_275w.pdf - cross correlation surfaces for the 275 x 275 pixel window size
- _cloud_cover.txt - an estimate of cloud cover based on green band brightness
- _gdal_trans_opts.txt - GCP points and options sent to the GDAL translate command
- _info_full.csv - a list of all the GCP points found by the algorithm 
- _info_sub.csv - a list of filtered GCP points by two standard deviations from mean offset
- _MTL.txt - image metadata file provided with the image from USGS 
- _proj.txt - PROJ.4 string definition of the used image projection
- _rmse.Rdata - an R list object that contains x, y, and total RMSE (not filtered by standard deviation), and all points used to calcuate RMSE.
- _wkt.txt - well-known text definition of the used image projection

## Defining and PROJ.4 projection 
 
Projections transform the curved, 3-dimensional shape of the Earth to a flat 2-dimensional surface, which simplifies its representation for display and computation. Landsat USGS Earth Explorer images are delivered in the Universal Transverse Mercator (UTM) projection (Polar Stereographic projection for scenes with a center latitude greater than or equal to -63.0 degrees) using the World Geodetic System (WGS) 84datum. This projection is fairly specific to small geographic regions, with different parameters for each UTM zone. This is not ideal for LLR projects, since there is high potential for including multiple scenes that may have different projections, which makes spectral harmonization and mosaicking difficult. 

In LLR we have you define the projection that will be applied to all images in the project. Do some research to figure out what projection is best for your region - what projection are other data you might be using in? It's important that whatever projection you use preserves area. If you use the image data to classify land cover and then summarize the area of each class in the image, you want to ensure that area is equivalent over the entire image. We recommend the Albers Equal Area Conic projection. The parameters can be adjusted to suit local to continental regions. 

When you run the LLR unpacking steps you will be prompted to define a projection to use by declaring a 
PROJ.4 string. You can find an extensive list of projections at <a href="http://www.spatialreference.org/" target="_blank">www.spatialreference.org</a>. Use the search option to find named projections and regions, or scroll through the lists until you find what you are looking for. Click on the projection link and then click on the *Proj4* link. Copy and paste the projection string into a text file and save it for reference and easy access when you begin running LLR.

**Here are a few PROJ.4 examples:**

North American Albers Equal Area Conic:
<pre>+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs</pre>

Africa Albers Equal Area Conic:
<pre>+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs</pre>

Europe Albers Equal area Conic:
<pre>+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs</pre>

China Albers Equal Area Conic:
<pre>+proj=aea +lat_1=27 +lat_2=45 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs</pre> 

South America Albers Equal Area Conic:
<pre>+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs</pre>

Canada Albers Equal Area Conic:
<pre>+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs</pre> 

NAD83 UTM zone 10N:
<pre>+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs</pre> 