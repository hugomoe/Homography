#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "affine.h"
#include "homo_box.h"
#include "parameters.h"
#include "aux_fun.h"


/**
  * how are indexed coefficients of matrices :
  *
  * for the homographies
  * 0,0   0,1   0,2
  * 1,0   1,1   1,2
  * 2,0   2,1   2,2
  *
  * for the affine maps
  * 0,0   0,1   0,2
  * 1,0   1,1   1,2
  */





// compute the smallest rectangle containing the result of an affine mapping
void smallest_rectangle(double A[2][3], int wIn, int hIn, int *muOut, int *nuOut, int *wOut, int *hOut){
/**
  * @param
  *     A an affine map
  *     [0,wIn]*[0,hIn] a rectangle of the plan
  *     [muOut,muOut+wOut]*[nuOut,nuOut+hOut] an output rectangle
  *
  * @return
  *     [muOut,muOut+wOut]*[nuOut,nuOut+hOut] is the smallest rectangle (with integer-coordinated vertices) containing A([0,wIn]*[0,hIn])
  *     (rectangle starting at (muOut,nuOut) with dimensions wOut*hOut)
  */

	//vertices of the rectangle [0,wIn]*[0,hIn]
	double xIn[4] = {0,wIn,wIn,0};
	double yIn[4] = {0,0,hIn,hIn};
	//vertices of the parallelogram A([0,wIn]*[0,hIn])
	double xOut[4];
	double yOut[4];
	for(int k=0;k<4;k++){
		xOut[k]=xIn[k]*A[0][0]+yIn[k]*A[0][1]+A[0][2];
		yOut[k]=xIn[k]*A[1][0]+yIn[k]*A[1][1]+A[1][2];
	}

    //extrema of the coordinates of the output rectangle's vertices
    double xOutMin_double = xOut[0];
    double yOutMin_double = yOut[0];
    double xOutMax_double = xOut[0];
    double yOutMax_double = yOut[0];
    for(int k=1;k<4;k++){
        xOutMin_double = fmin(xOutMin_double,xOut[k]);
        yOutMin_double = fmin(yOutMin_double,yOut[k]);
        xOutMax_double = fmax(xOutMax_double,xOut[k]);
        yOutMax_double = fmax(yOutMax_double,yOut[k]);
    }
    int xOutMin = floor(xOutMin_double);
    int yOutMin = floor(yOutMin_double);
    int xOutMax = ceil(xOutMax_double);
    int yOutMax = ceil(yOutMax_double);

    //output
    *muOut = xOutMin;
    *nuOut = yOutMin;
    *wOut = xOutMax - xOutMin;
    *hOut = yOutMax - yOutMin;
}



// decompose H in H = A H0 B
void decomp(double H[3][3],double A[2][3],double H0[3][3],double B[2][3]){
/**
  * @param
  *     H an input homography
  *		H0 an output homography
  *		A,B output affinities
  */
    double a=H[0][0], b=H[0][1], p=H[0][2], c=H[1][0], d=H[1][1], q=H[1][2], r=H[2][0], s=H[2][1], t=H[2][2];
    double t0, t1;
	//assume r,s != 0,0

    //an infinity of (t0,t1) are possible. Here is a simple one
    if(fabs(a*s-b*r)<fabs(c*s-d*r)){
        t0 = 0.;
        t1 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(c*s-d*r);
    }else{
        t0 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(a*s-b*r);;
        t1 = 0.;
    }



    //translation of (t0,t1) (H becomes T_(t0,t1)*H)
    a += t0*r;
    b += t0*s;
    p += t0*t;
    c += t1*r;
    d += t1*s;
    q += t1*t;



	double N_phi = sqrt(pow(r,2.)+pow(s,2.));
	B[0][0] = r/N_phi;
	B[0][1] = s/N_phi;
	B[0][2] = 0.;
	B[1][0] = -s/N_phi;
	B[1][1] = r/N_phi;
	B[1][2] = 0.;

	double N_psi = sqrt(pow(a*s-b*r,2.)+pow(c*s-d*r,2.));
	A[0][0] = (c*s-d*r)/N_psi;
	A[0][1] = (a*s-b*r)/N_psi;
	A[0][2] = -t0;
	A[1][0] = -(a*s-b*r)/N_psi;
	A[1][1] = (c*s-d*r)/N_psi;
	A[1][2] = -t1;

	H0[0][0] = -(a*d-b*c)*N_phi/N_psi;
	H0[0][1] = 0.;
	H0[0][2] = (p*(c*s-d*r)-q*(a*s-b*r))/N_psi;
	H0[1][0] = 0.;
	H0[1][1] = -N_psi/N_phi;
	H0[1][2] = (p*(a*s-b*r)+q*(c*s-d*r))/N_psi;
	H0[2][0] = N_phi;
	H0[2][1] = 0.;
	H0[2][2] = t;
}



void apply_homography_decomp(float *img,float *img_f,int w,int h,int pd,int w_f,int h_f,double H[3][3]){

/**
  * @param
  *     img : input image
  *     img_f : output image
  *     w,h : width and heights of the input image
  *		w_f,h_f : width and heights of the output image
  *     H : inverse homography, such that img_f(x,y)=img(H(x,y))
  */

	if(H[2][2]!=0 && H[2][0]/H[2][2]==0 && H[2][1]/H[2][2]==0){    //H is an affine map
		double A[2][3] = {
			H[0][0]/H[2][2], H[0][1]/H[2][2], H[0][2]/H[2][2],
			H[1][0]/H[2][2], H[1][1]/H[2][2], H[1][2]/H[2][2]};
		apply_affine_map(img,img_f,w,h,pd,w_f,h_f,A);
	}else{  //H is an homography
		double A[2][3];
		double H0[3][3];
		double B[2][3];

		decomp(H,A,H0,B);   //compute the decomposition H = A H0 B



		//declare the size (w,h) and the position (mu,nu) of every intermediate image
		/*
		 * For an image at the position (mu,nu), the coordinates of the (i,j)-th pixel are (i+mu,j+nu)
		 * rectangleX is the rectangle of the plan [muX,muX+wX]*[nuX,nuX+hX] containing imgX
		 * Thus, rectangle1 is the input image and rectangle_f is the output image
		 *
		 *
		 * img1 is the input image (mu1=nu1=0, w1=w, h1=h)
		 * img2 and img3 are the intermediate images
		 * img4 is the output image (mu4=nu4=0, w4=w_f, h4=h_f)
		 */
		int mu2,nu2,w2,h2,
			mu3,nu3,w3,h3;

		//rectangle2 = invA(rectangle1) where invA = A^(-1)
		double detA = A[0][0]*A[1][1]-A[1][0]*A[0][1];
		double invA[2][3] = {
			A[1][1]/detA,
			-A[0][1]/detA,
			(A[0][1]*A[1][2]-A[1][1]*A[0][2])/detA,
			-A[1][0]/detA,
			A[0][0]/detA,
			(-A[0][0]*A[1][2]+A[1][0]*A[0][2])/detA};
		smallest_rectangle(invA,w,h,&mu2,&nu2,&w2,&h2);
		if((w2-w)%2!=0){w2++;}  //w2 must have the same parity than w
		if((h2-h)%2!=0){h2++;}  //h2 must have the same parity than h

		//rectangle3 = B(rectangle4)
		smallest_rectangle(B,w_f,h_f,&mu3,&nu3,&w3,&h3);
		if((w3-w_f)%2!=0){w3++;}    //w3 must have the same parity than w_f
		if((h3-h_f)%2!=0){h3++;}    //h3 must have the same parity than h_f



		//the affinities must be corrected so that they fit the positions of the rectangles
		/*
		 * the exact formulas are
		 * A[0][2] = mu2*A[0][0] + nu2*A[0][1] + A[0][2] - mu1;
		 * A[1][2] = mu2*A[1][0] + nu2*A[1][1] + A[1][2] - nu1;
		 * B[0][2] = mu4*B[0][0] + nu4*B[0][1] + B[0][2] - mu3;
		 * B[1][2] = mu4*B[1][0] + nu4*B[1][1] + B[1][2] - nu3;
		 * but a lot of term are zero
		 */
		A[0][2] = mu2*A[0][0] + nu2*A[0][1] + A[0][2];
		A[1][2] = mu2*A[1][0] + nu2*A[1][1] + A[1][2];
		B[0][2] = - mu3;
		B[1][2] = - nu3;



		///Application of the decomposition

		float *img2 = malloc(pd*w2*h2*sizeof(float));
		apply_affine_map(img,img2,w,h,pd,w2,h2,A);
    iio_save_image_float_vec("temp1.png",img2,w2,h2,pd);

		float *img3 = malloc(pd*w3*h3*sizeof(float));
		apply_unidirectional_homography(img2,img3,w2,h2,pd,w3,h3,mu2,nu2,mu3,nu3,H0);
    iio_save_image_float_vec("temp3.png",img3,w3,h3,pd);
		free(img2);

		apply_affine_map(img3,img_f,w3,h3,pd,w_f,h_f,B);
    iio_save_image_float_vec("temp4.png",img_f,w_f,h_f,pd);
		free(img3);
	}


	double p[2];

	//truncate the output image, because it has been symmetrized to prevent ringing
	for(int i=0;i<w_f;i++){
		for(int j=0;j<h_f;j++){
			p[0]=i; p[1]=j;
			apply_homography_1pt(p,H,p);
			if(p[0]<0 || p[0]>w || p[1]<0 || p[1]>h){
				for(int l=0;l<pd;l++){img_f[(j*w_f+i)*pd+l]=0;}
			}
		}
	}

}
