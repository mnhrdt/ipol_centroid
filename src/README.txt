DESCRIPTION
-----------

This is an implementation of the centroid method for the
correction of turbulence.

This program is part of an IPOL article "Implementation of the
Centroid Method for the Correction of Turbulence", by Enric
Meinhardt-Llopis and Mario Micheli, 2013.

This is version 2.0 of the program, published on december 2013.
The source code and eventual future modifications are available from github:

	https://github.com/mnhrdt/ipol_centroid/tree/master/src



COMPILATION
-----------

Verify that you have the libraries libpng and libtiff installed, and a C
compiler.  Run "make".  This shall produce two executables: "centroid" and
"combine"



USAGE
-----

The "combine" executable combines several images using a simple pointwise
formula (max, min, average or median).  For example:

	./combine avg frame_*.png -o average_frame.png
	./combine min frame_*.png -o min_frame.png
	./combine weisz frame_*.png -o median_frame.png

The "centroid" executable combines several images using the centroid method
(based on optical flows computed by a multiscale version of Horn-Schunck).  For
example:

	./centroid frame_%02d.jpg 0 10 frame_00.jpg centroid_from_0.png
	./centroid frame_%02d.jpg 0 10 frame_01.jpg centroid_from_1.png
	./centroid frame_%02d.jpg 0 10 frame_02.jpg centroid_from_2.png

Notice that the programs can read images in many formats (PNG, JPEG, TIFF, PPM,
etc.) but can only write PNG or TIFF images.




FILES
-----

README.txt:    this file
Makefile:      rules for building the program
centroid.c:    implementation of the centroid method
combine.c:     implementation of simple combinations (average, min, max, median)
bicubic.c:     auxiliary functions for computing bicubic interpolation
iio.c:         auxiliary functions for reading and writing images (not reviewed)
iio.h:         header for iio.c
optical_flow/: directory containing an optical flow implementation
test/:         directory containing a test input sequene
LICENSE:       license and copyright details
agpl-3.0.txt:  AGPL license text




AUTHORS
-------

* Source code and article:
	- Enric Meinhardt-Lopis <enric.meinhardt@cmla.ens-cachan.fr>
	- Mario Micheli <micheli@uw.edu>

* Part of the source code:
	- Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>




COPYRIGNT AND LICENSE
---------------------

See file LICENSE for details.
