/**
  * warp a huge image by bilinear interpolation
  * zoom-out the huge image and the warped huge image
  * 
  * assumption : the input is of size z*wOut * z*hOut for some integer z
  * 
  * compilation
  * 	gcc-5 -fopenmp -O3 warp_ground_truth.c -lfftw3 -ltiff -ljpeg -lpng -o warp_ground_truth
  * 	gcc -std=c99 -fopenmp -O3 warp_ground_truth.c -lfftw3 -ltiff -ljpeg -lpng -o warp_ground_truth
  * usage
  * 	./warp_ground_truth [image.png] x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1 wOut hOut
  * 	to send the points x,y,z,t (of the zoomed-out image) to the points a,b,c,d
  */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.c"
#include "zoom_out.h"
#include "parameters.h"
#include "aux_fun.h"



int main(int argc,char *argv[]){

	if (argc != 20) {
		printf("usage :\n\t[image.png] x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1 wOut hOut"); 
		return 1;
	}


	
//read the input
	char *filename_in = argv[1];

	double x[2] = {strtod(argv[2],NULL), strtod(argv[3],NULL)},
        y[2] = {strtod(argv[4],NULL), strtod(argv[5],NULL)},
        z[2] = {strtod(argv[6],NULL), strtod(argv[7],NULL)},
        t[2] = {strtod(argv[8],NULL), strtod(argv[9],NULL)},
        a[2] = {strtod(argv[10],NULL), strtod(argv[11],NULL)},
        b[2] = {strtod(argv[12],NULL), strtod(argv[13],NULL)},
        c[2] = {strtod(argv[14],NULL), strtod(argv[15],NULL)},
        d[2] = {strtod(argv[16],NULL), strtod(argv[17],NULL)};
	double H[3][3];
	homography_from_eight_points(H,a,b,c,d,x,y,z,t);

	int wOut = strtol(argv[18],NULL,10);
	int hOut = strtol(argv[19],NULL,10);

	float *img_big;
	int w,h,pd;

	img_big = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (h/hOut != w/wOut) {printf("image size is incompatible with wOut and hOut.\n");exit(1);}

//zoom-out
	float *img_small = malloc(wOut*hOut*pd*sizeof(float));
	if(img_small==NULL){
		printf("@warp_groundtruth.c @main : img_small is too large");
		exit(1);
	}
	zoom_out(img_big,img_small,w,h,pd,wOut,hOut);
	iio_save_image_float_vec("img_small.png",img_small,wOut,hOut,pd);
	free(img_small);

//ground_truth
	float *img_big_warped = malloc(pd*w*h*sizeof(float));
	if(img_big_warped==NULL){
		printf("@warp_groundtruth.c @main : img_big_warped is too large");
		exit(1);
	}

	//naive warping
	double fzoom=(double)h/(double)hOut;
	double HH[3][3];
	HH[0][0]=H[0][0];
	HH[0][1]=H[0][1];
	HH[0][2]=H[0][2]*fzoom;
	HH[1][0]=H[1][0];
	HH[1][1]=H[1][1];
	HH[1][2]=H[1][2]*fzoom;
	HH[2][0]=H[2][0]/fzoom;
	HH[2][1]=H[2][1]/fzoom;
	HH[2][2]=H[2][2];
	for(int l=0;l<pd;l++)
		for (int j = 0; j < h; j++)
			for (int i = 0; i < w; i++){
				double p[2] ={i,j};
				
				apply_homography_1pt(p, HH, p);
				int idx = pd*(w * j + i)+l;
				img_big_warped[idx] = bilinear_interpolation_at(img_big, w, h, pd, p[0], p[1], l);
	}

	//output
	float *img_grd = malloc(pd*wOut*hOut*sizeof(float));
	if(img_grd==NULL){
		printf("@warp_groundtruth.c @main : img_grd is too large");
		exit(1);
	}
	zoom_out(img_big_warped,img_grd,w,h,pd,wOut,hOut);
	free(img_big_warped);
	iio_save_image_float_vec("img_grd.png",img_grd,wOut,hOut,pd);
    free(img_grd);

	return 0;
}
