# open-clusters

I am using GitHub to organize my MSci Project work, which focuses on the photometry of open clusters for age and distance estimation. The images for analysis were collected by myself and my colleagues using the on-site telescope at Royal Holloway, University of London.

The primary aim of this project is to explore the capabilities of the DAOPHOT _starFinder_ algorithm in identifying stars from these images. My code utilises DAOPHOT and functions from the _photutils_ and _Astropy_ libraries to record a star's RA/DEC coordinates and perform aperture photometry to obtain the instrumental BVR magnitudes. These magnitudes will be calibrated against established values from the SIMBAD online database. The ultimate goal is to produce colour-magnitude diagrams for a series of target clusters, such as M39, NGC 744, and NGC 7686. These will be compared to sets of isochrones to make estimates of the targets' age and distance from Earth.

The DAOPHOT algorithm is implemented in **DAOPHOT-starfinder.py**, which outputs Excel files containing the coordinates and magnitude estimates. To run the program on different source images, you must use the `.fits` file's cluster name and BVR filter information to update the variables in the script, as these details influence the naming of the output Excel files. The program also generates:
- Starfield images
- Background distribution plots using the Photutils _Background2D_ function
- Starfield plots with points overlaid to indicate the 35 brightest stars identified by DAOPHOT  

Plots are displayed in interactive mode using Matplotlib's _ion_ function. If you are testing the program, ensure your IDE or console supports interactive plotting.

**pair_calibration.py** is a program that compares the coordinates of stars identified by the DAOPHOT algorithm to those listed in SIMBAD for the same cluster. The goal is to match stars between the datasets and between the BV images, allowing the instrumental and established BV magnitudes to be paired for each star in the cluster. A set of equations must be solved to calibrate the instrumental magnitudes, which involves matrix manipulation and estimations of 4 free-parameters. The program **iminuit-param-est.py** utilises the _iminuit_ Python interface to compute the MLEs for each parameter, allowing the calibrated BV magnitudes to be computed. Once calibration is complete, I will begin producing colour-magnitude diagrams for the clusters. Comparisons will be made to isochrones using the Harvard University MIST isochrone interpolation software.

The example `.fits` files for NGC 7686 can be used to test the code. An additional selection of my images, an output Excel file of paired star magnitudes, and a text file of the parameter estimates are provided for reference. This information is provided for academic purposes as part of my MSci Project and is not intended for redistribution or commercial use.
