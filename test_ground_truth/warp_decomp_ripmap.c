/**
  * warp an image by an homography using decomposition method and ripmap
  * 
  * compilation
  * 	gcc-5 -fopenmp -O3 warp_decomp_ripmap.c -lfftw3 -ltiff -ljpeg -lpng -o warp_decomp_ripmap
  * 	gcc -std=c99 -fopenmp -O3 warp_decomp_ripmap.c -lfftw3 -ltiff -ljpeg -lpng -o warp_decomp_ripmap
  * usage
  * 	./warp_decomp_ripmap [image.png] x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1 wOut hOut
  *		to send the points x,y,z,t to a,b,c,d
  */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.c"
#include "decomp.h"
#include "zoom_out.h"
#include "ripmap.h"
#include "parameters.h"
#include "aux_fun.h"



int main(int argc,char *argv[]){

	if (argc != 20) {
		printf("usage :\n\t[image.png] x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1 wOut hOut"); 
		return 1;
	}


	
//read the input
	char *filename_in = argv[1];
	float *img;
	int w,h,pd;
	img = iio_read_image_float_vec(filename_in, &w, &h, &pd);

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



//Ripmap
	float *img_rip = malloc(3*wOut*hOut*sizeof(float));
	if(img_rip==NULL){
        printf("@compare_with_truth.c @main : img_rip is too large");
        exit(1);
    }
	
	if(pd==3){
        apply_homo_ripmap(img,img_rip,w,h,wOut,hOut,H);
	}else{//suppose pd=1
        float *img3 = malloc(3*w*h*sizeof(float));
        for(int i=0;i<w*h;i++){
            for(int l = 0;l<3;l++){
                img3[3*i+l]=img[i];
            }
        }
        apply_homo_ripmap(img3,img_rip,w,h,wOut,hOut,H);
	}
	iio_save_image_float_vec("img_rip.png",img_rip,wOut,hOut,3);
    free(img_rip);



//Decomposition
	float *img_dec = malloc(wOut*hOut*pd*sizeof(float));
	if(img_dec==NULL){
        printf("@compare_with_truth.c @main : img_dec is too large");
        exit(1);
    }

	apply_homography_decomp(img,img_dec,w,h,pd,wOut,hOut,H);
	iio_save_image_float_vec("img_dec.png",img_dec,wOut,hOut,pd);
    free(img_dec);
	


	free(img);
	return 0;
}
