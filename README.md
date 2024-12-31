# open-clusters

I am using GitHub to organize my MSci Project work, which focuses on the photometry of open clusters for age and distance estimation. The images for analysis were collected by myself and my colleagues using the on-site telescope at Royal Holloway, University of London.

The primary aim of this project is to explore the capabilities of the DAOPHOT and IRAF _starFinder_ algorithms in identifying stars from these images. These algorithms provide each star's RA/DEC coordinates and instrumental BVR magnitudes. The instrumental magnitudes will be calibrated against established values from the SIMBAD online database. The ultimate goal is to produce colour-magnitude diagrams for a series of target clusters, such as M39 and NGC 7686.

The DAOPHOT and IRAF algorithms are implemented in **dao-iraf-starfinder.py**, which outputs Excel files containing the coordinates and magnitude estimates. To run the program on different source images, you must use the `.fits` file's cluster name and BVR filter information to update the variables in the script, as these details influence the naming of the output Excel files. The program also generates:
- Starfield images
- Background distribution plots using the Photutils _Background2D_ function
- Starfield plots with points overlaid to indicate the 15 brightest stars identified by each algorithm  

Plots are displayed in interactive mode using Matplotlib's _ion_ function. If you are testing the program, ensure your IDE or console supports interactive plotting.

**calibration.py** is an initial program (currently under development) that compares the coordinates of stars identified by the DAOPHOT algorithm to those listed in SIMBAD for the same cluster. The goal is to match stars between the datasets, allowing the instrumental and established BVR magnitudes to be paired and calibrated. Once calibration is complete, I will begin producing colour-magnitude diagrams for the clusters.

The example `.fits` file for NGC 7686 can be used to test the code. Additionally, a selection of my images is provided for reference.
