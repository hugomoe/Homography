#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "parameters.h"
#include "aux_fun.h"

/**
  * the program homo_box use the fourth integral image to interpol an image by a piecewise affine function
  * and to convolve this function with the filter g3 (convolution of 3 boxes of width D)
  */

//build the integral images
void build_fourth_int(float *img,double *Img,int wh){
/**
  * @param
  *     img : row/column to integrate (1 level of color)
  *     wh : size of img
  *     Img : output row/column
  *         will contain 1-, 2-, 3- and 4-integral of img on points 0..wh
  *         suppose img is piecewise constant
  */
	//initialize values of integral images at 0 on 0
	Img[0] = 0;
	Img[wh+1] = 0;
	Img[2*wh+2] = 0;
	Img[3*wh+3] = 0;

	int u;
	double q,l;
	for(u=1;u<wh+1;u++){
		q = (double) img[u-1];
		//first integral (indexed from 0 to wh)
		Img[u] = Img[u-1] +  q;
		//second integral (indexed from wh+1 to 2*wh+1)
		Img[u+wh+1] = Img[u+wh] + Img[u-1] + q/2;
		//third integral (indexed from 2*wh+2 to 3*wh+2)
		Img[u+2*wh+2] = Img[u+2*wh+1] + Img[wh+u] + Img[u-1]/2 + q/6;
		//fourth integral (indexed from 3*wh+3 to 4*wh+3)
		Img[u+3*wh+3] = Img[u+3*wh+2] + Img[2*wh+1+u] + Img[wh+u]/2 + Img[u-1]/6 + q/24;
	}
}

//evaluate integral images
double eval_fourth_int(float *img,double *Img,double xy,int wh){
/**
  * @param
  *     img : row/column to integrate (1 level of color)
  *     wh : size of img
  *     Img : row/column containing 1-, 2-, 3- and 4-integral of img on points 0..wh
  *     xy : point where to evaluate Img
  */
    double r,Q;

    if(xy>=wh){ //compute the polynomial of degree 3 if xy exceed the size of the image
        r = (double) (xy-wh);
        return Img[4*wh+3] + r*(Img[3*wh+2] + r*(Img[2*wh+1]/2 + r * Img[wh]/6));
    }
    else if(xy<0){ //zero if xy is negative
        return 0;
    }
    else{ //compute F4(floor(xy)) + polynomial of interpolation between floor(xy) and floor(xy)+1
        int xys = floor(xy);
        r = (double) (xy-xys);
        Q = (double) img[xys];
        return Img[3*wh+3+xys] + r*(Img[2*wh+2+xys] + r*(Img[wh+1+xys]/2 + r*(Img[xys]/6 + r*Q/24)));
    }
}



//convolve an image with g3 (the image is considered as piecewise affine)
float convolve_img(float *img,double *Img,double xy,float d,int wh){
/**
  * @param
  *     img : row/column to integrate (1 level of color)
  *     wh : size of img
  *     Img : row/column containing 1-, 2-, 3- and 4-integral of img on points 0..wh
  *     xy : point where to evaluate Img
  *     d : zoom factor (determining standard deviation of g3)
  *     D : size of the box functions defining g3
  */
    
    double d_aux = pow(FINAL_BLUR,2)*pow(d,2)-pow(INITIAL_BLUR,2);
    if(d_aux<0.0){d_aux = 0.0;}
    double D = sqrt(2.0*d_aux+1)-1; //wells formula
    // double D = sqrt(2.0*d_aux); //iterated linear filters formula
    // double D = 2.0*sqrt(d_aux); //old formula
    //limit standard deviation D to avoid numerical problems when convolving with too sharp kernel
    if (D<0.05) {D = 0.05;}

    if(D>=wh/2 && xy>=0 && xy<=wh){
        return (float) Img[wh]/wh; //mean of img
    }else{
        //needed points (where to evaluate the fourth integral image to compute convolution)
        double xy1,xy2,xy3,xy4,xy5,xy6,xy7,xy8;

        xy1 = (double) xy + 3.*D/2. + 1.;
        xy2 = (double) xy + 3.*D/2.;
        xy3 = (double) xy + D/2. + 1.;
        xy4 = (double) xy + D/2.;
        xy5 = (double) xy - D/2. + 1.;
        xy6 = (double) xy - D/2.;
        xy7 = (double) xy - 3.*D/2. + 1.;
        xy8 = (double) xy - 3.*D/2.;

        // a. = valeur de l'image 4-int en .
        // b. = valeur de la derivé discrète de l'image 4-int
        // c. = valeur de la derivé discrète seconde de l'image 4-int
        // d. = valeur de la derivé discrète troisième de l'image 4-int
        // e. = valeur de la derivé discrète quatrième de l'image 4-int = valeur de la convolution

        /*
         * aX : value of 4-integral image in xyX
         * bX : value of discrete derivative of 4-integral image
         * cX : value of second discrete derivative of 4-integral image
         * dX : value of third discrete derivative of 4-integral image
         * eX : value of fourth discrete derivative of 4-integral image (returned value)
         */
        double a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,c1,c2,c3,d1,d2,e1;

        //evaluate fourth integral
        a1 = eval_fourth_int(img,Img,xy1,wh);
        a2 = eval_fourth_int(img,Img,xy2,wh);
        a3 = eval_fourth_int(img,Img,xy3,wh);
        a4 = eval_fourth_int(img,Img,xy4,wh);
        a5 = eval_fourth_int(img,Img,xy5,wh);
        a6 = eval_fourth_int(img,Img,xy6,wh);
        a7 = eval_fourth_int(img,Img,xy7,wh);
        a8 = eval_fourth_int(img,Img,xy8,wh);

        /*
         * discrete derivative of size 1 (to consider a piecewise constant function as a piecewise
         * affine function, there is an implicit convolution with box of size 1)
         */
        b1 = a1 - a2;
        b2 = a3 - a4;
        b3 = a5 - a6;
        b4 = a7 - a8;

        //discrete derivative
        c1 = (b1 - b2)/D;
        c2 = (b2 - b3)/D;
        c3 = (b3 - b4)/D;

        //discrete derivative
        d1 = (c1 - c2)/D;
        d2 = (c2 - c3)/D;

        //discrete derivative
        e1 = (d1 - d2)/D;

        return (float) e1;
    }
}


int apply_unidirectional_homography(float *img,float *img_f,int w,int h,int pd, int w_f,int h_f,int mu,int nu,int mu_f,int nu_f,double H[3][3]){
/**
  * @param
  *     img, img_f : initial and final images
  *     w,h, w_f,h_f : dimensions of images
  *     mu,nu, mu_f,nu_f : position of images (in the plan)
  *     H : homography matrix such that H[0][1]=H[0][3]=H[2][1]=0
  *
  * since a pixel has a size of 1, its inverse image has a size d, with d the absolute value of the derivative
  * at that point
  * in this function, x and y represent real coordinates (taking into account the position of the image)
  * and i and j represent indexes in the image (seen as an array)
  *
  * Img (with uppercase) represent an integral image, while img is an image
  */
    int l;

    //w_aux,h_aux, mu_aux,nu_aux pour l'image intermÈdiaire img_aux
    //mu_aux,nu_aux,w_aux,h_aux are the position and size of the auxiliary image
    int w_aux = w_f; //the second step does not change x
    int h_aux = h; //the first step does not change y
    int mu_aux = mu_f; //the second step does not change x
    int nu_aux = nu; //the first step does not change y

    float *column = malloc(w*sizeof(float));
    double *Img = malloc(4*(w+1)*sizeof(double)); //successive integrals of column
    float *img_aux = malloc(w_aux*h_aux*sizeof(float));
	
    float *row_aux = malloc(h_aux*sizeof(float));
    double *Img_aux = malloc(4*(h_aux+1)*sizeof(double)); //successive integrals of row_aux
    float *img_aux2 = malloc(w_f*h_f*sizeof(float));

	for(l=0;l<pd;l++){
		//first step
	    for(int j=0;j<h_aux;j++){
            for(int i=0;i<w;i++){column[i] = img[pd*(i+j*w)+l];} //extract the j-th column
            build_fourth_int(column,Img,w);

			#pragma parallel for
			for(int i=0;i<w_aux;i++){
				float x = (float) (i+mu_aux);
        float d = fabs((H[0][0]*H[2][2]-H[2][0]*H[0][2])/pow(H[2][0]*x+H[2][2],2)); //derivative with respect to x
				x = (H[0][0]*x+H[0][2])/(H[2][0]*x+H[2][2]) - (float) mu; //apply the homography
				img_aux[i+j*w_aux] = convolve_img(column,Img,x,d,w);
			}
		}
    if(l==0){
    	iio_save_image_float_vec("temp2_color0.png",img_aux,w_aux,h_aux,1);
    }



		//second step
		for(int i=0;i<w_f;i++){
			for(int j=0;j<h_aux;j++){row_aux[j] = img_aux[i+j*w_aux];} //extract the i-th row
			build_fourth_int(row_aux,Img_aux,h_aux);

			float x =(float) (i+mu_f);
			float d = fabs(H[1][1]/(H[2][0]*x+H[2][2])); //derivative with respect to y (does not depend on y)

			#pragma parallel for
			for(int j=0;j<h_f;j++){
				float y = (float) (j+nu_f);
				y = (H[1][1]*y+H[1][2])/(H[2][0]*x+H[2][2]) - (float) nu_aux; //apply the homography
				img_aux2[i+j*w_f] = convolve_img(row_aux,Img_aux,y,d,h_aux);
			}
		}



		//output
		for(int i=0;i<w_f*h_f;i++){img_f[pd*i+l]=img_aux2[i];}
	}
	free(column);
	free(Img);
	free(img_aux);
	free(row_aux);
	free(Img_aux);
	free(img_aux2);
	return 0;
}
