CC = c99

LFLAGS = -ljpeg -ltiff -lpng -lm

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = -g
CFLAGS = -O3 -fopenmp

SRC = centroid.c iio.c optical_flow/pyramidal_horn_schunck.c optical_flow/smooth_inpainting.c optical_flow/optical_flow.c

all: centroid imcombine

centroid: $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LFLAGS) -o centroid

imcombine: imcombine.c iio.c
	$(CC) $(CFLAGS) imcombine.c iio.c $(LFLAGS) -o imcombine

clean:
	rm -f centroid imcombine

zip:
	tar -C .. -acvf ../centroid-2.tar.gz src --transform s,src,centroid-2,
