/**
  * compute the L1 and L2 norm of the difference between some sub-images of to images
  * 
  * compilation
  * 	gcc-5 -fopenmp -O3 compare_subimage.c -lfftw3 -ltiff -ljpeg -lpng -o compare_subimage
  * 	gcc -std=c99 -fopenmp -O3 compare_subimage.c -lfftw3 -ltiff -ljpeg -lpng -o compare_subimage
  * usage
  * 	./compare_subimage [image1.png] [image2.png]
  *		to compare all the image
  *		or
  * 	./compare_subimage [image1.png] [image2.png] i1 j1 i2 j2
  *		to compare the subimage between pixels (i1,j1) and (i2,j2)
  */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.c"



int main(int argc,char *argv[]){

	if (argc != 3 && argc != 7) {
		printf("usage :\n\t[image1.png] [image2.png]\nor\n\t[image1.png] [image2.png] i1 j1 i2 j2"); 
		return 1;
	}


	
//read the input
	char *filename1 = argv[1];
	char *filename2 = argv[2];
	float *img1, *img2;
	int w1,h1,pd1, w2,h2,pd2;

	img1 = iio_read_image_float_vec(filename1, &w1, &h1, &pd1);
	img2 = iio_read_image_float_vec(filename2, &w2, &h2, &pd2);

	int i1,j1,i2,j2; //coordinates of two opposite vertices of a sub-image of the output
	switch(argc){
		case 6:
		i1 = strtol(argv[2],NULL,10);
		j1 = strtol(argv[3],NULL,10);
		i2 = strtol(argv[4],NULL,10);
		j2 = strtol(argv[5],NULL,10);
		break;
		default:
		i1 = 0;
		j1 = 0;
		i2 = w1-1;
		j2 = h1-1;
		break;
	}


//error on a sub-image
	int i_tl, j_tl, //coordinates of the top left pixel of the sub-image
		i_br, j_br; //coordinates of the bottom right pixel of the sub-image
	if(i1<i2){
		i_tl = i1;
		i_br = i2;
	}else{
		i_tl = i2;
		i_br = i1;
	}
	if(j1<j2){
		j_tl = j1;
		j_br = j2;
	}else{
		j_tl = j2;
		j_br = j1;
	}

	//the sub-image must be a subset of the output image
	if(i_tl<0){i_tl=0;}
	if(j_tl<0){j_tl=0;}
	if(i_br>w1-1){i_br=w1-1;}
	if(j_br>h1-1){j_br=h1-1;}

	float l1 = 0., l2 = 0.;
	float diff;
	int idx;

	for(int l=0 ; l<pd1 ; l++)
		for(int i=i_tl ; i<=i_br ; i++)
			for(int j=j_tl ; j<=j_br ; j++){
				idx = pd1*(i+w1*j) + l;

				diff = fabs(img1[idx]-img2[idx]);
				l1 += diff;
				l2 += pow(diff,2);
	}

	//normalize
	int N = (i_br-i_tl+1)*(j_br-j_tl+1)*pd1; //size of the sub-image
    l1 = l1/(float) N;
    l2 = sqrt(l2/(float) N);

    //print
    printf("l1-error : %f\nl2-error : %f\n",l1,l2);
	return 0;
}
