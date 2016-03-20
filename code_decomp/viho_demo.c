//gcc-5 -fopenmp -O3 viho.c -I/usr/local/include/libiomp -I/usr/X11/include -I/Users/Hhhh/ENS/Stage_L3_math/homographies/code/jpeg-6b -L/usr/X11/lib -lfftw3 -lX11 -L/usr/local/Cellar/libtiff/4.0.3 -ltiff -ljpeg -lpng


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "iio.c"
#include "decomp.h"


/**
  * this program warps an image by an homography
  * the warping is such that
  *
  * - maps to -
  * x   -->   a
  * y   -->   b
  * z   -->   c
  * t   -->   d
  */



int main(int argc,char *argv[]){
	if (argc != 20) {
		printf("usage :\n\t[image.png] x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1 wOut hOut\n");
		return 1;
	}



	//image input
	char *filename_in = argv[1];
	float *img;
	int w,h,pd;
	img = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	//homography input
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

	//size of the output
	int wOut = strtol(argv[18],NULL,10);
	int hOut = strtol(argv[19],NULL,10);
	float *img_f = malloc(3*wOut*hOut*sizeof(float));



    //resample (and time)
	clock_t debutcpu,fincpu;
	double debutreal,finreal;
	debutcpu = clock();
	debutreal = omp_get_wtime();
	// if(pd==3){
 //        apply_homography_decomp(img,img_f,w,h,wOut,hOut,H);
	// }else{//suppose pd=1
 //        float *img3 = malloc(3*w*h*sizeof(float));
 //        for(int i=0;i<w*h;i++){
 //            for(int l = 0;l<3;l++){
 //                img3[3*i+l]=img[i];
 //            }
 //        }
 //        apply_homography_decomp(img3,img_f,w,h,wOut,hOut,H);
	// }
	apply_homography_decomp(img,img_f,w,h,pd,wOut,hOut,H);

	fincpu = clock();
	finreal = omp_get_wtime();
	printf("cputime :%fs\ntime : %fs\n",(double)(fincpu-debutcpu)/CLOCKS_PER_SEC,(double)(finreal-debutreal));



	//output
	iio_save_image_float_vec("img_f.png",img_f,wOut,hOut,pd);
	free(img);
	free(img_f);



	return 0;
}
