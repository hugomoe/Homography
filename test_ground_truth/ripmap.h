#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"
#include "aux_fun.h"



/**
  * construct the RipMap (a float-typed image 4 times larger than the input image, pd=3)
  */




//gives the coordinates in the RipMap depending on the size of the rectangles (given by i,j)
int coord(int i, int j,
	int u, int v,
	int w, int h, int l){
	/*
	 * i,j give the size of the rectangle
	 * the largest rectangle representing the original image corresponds to i=j=0
	 *
	 * u,v are the coordinates in the original image
	 * 
	 * w is the width of the original image, h its height
	 * 
	 * l is the current channel of color
	 */
	int x = good_modulus(u,w/pow(2,i)),
		y = good_modulus(v,h/pow(2,j));
	int ir = 2*(1-1/pow(2,i))*w + x, //abscissa in the RipMap
		jr = 2*(1-1/pow(2,j))*h + y; //ordinate in the RipMap
	return (ir + 2*w*jr)*3+l;
}



//convolve horizontally by a gaussian kernel
float gaussian_filter_h(float *r,int i,int j,int u,int v,int w,int h,int l,float *g){
	float a,b;
	a = b = 0;
	for(int f=0;f<TAPSR;f++){a += r[coord(i-1,j,2*u+f-(TAPSR-1)/2,v,w,h,l)]*g[f];}
	for(int f=0;f<TAPSR;f++){b += r[coord(i-1,j,2*u+1+f-(TAPSR-1)/2,v,w,h,l)]*g[f];}
	return (a+b)/2;
}
//convolve vertically by a gaussian kernel
float gaussian_filter_v(float *r,int i,int j,int u,int v,int w,int h,int l,float *g){
	float a,b;
	a = b = 0;
	for(int f=0;f<TAPSR;f++){a += r[coord(i,j-1,u,2*v+f-(TAPSR-1)/2,w,h,l)]*g[f];}
	for(int f=0;f<TAPSR;f++){b += r[coord(i,j-1,u,2*v+1+f-(TAPSR-1)/2,w,h,l)]*g[f];}
	return (a+b)/2;
}

//build the RipMap
void build_ripmap(float *img,float *r,int w,int h,int pd){
	int l,ll,i,j,u,v,w1,h2;
	int logw = (int) log2(w);
	
	//1D gaussian kernel
	float gauss1D[TAPSR];
	for(int f=0;f<TAPSR;f++){
		gauss1D[f]=exp(-pow((TAPSR-1)/2-f,2)/(2*pow(SIG,2)))/(sqrt(2*PI)*SIG);
	}
	float total = 0;
	for(int f=0;f<TAPSR;f++){total += gauss1D[f];}
	for(int f=0;f<TAPSR;f++){gauss1D[f] = gauss1D[f]/total;}
	
	//copy the original image in the rectangle i,j=0,0
	for(u=0;u<w;u++){
		for(v=0;v<h;v++){
			for(l=0;l<3;l++){
				if(l>pd-1){ll=pd-1;}else{ll=l;}
				r[coord(0,0,u,v,w,h,l)]=img[(u+v*w)*pd+ll];
			}
		}
	}

	//compute the rectangles obtained by successive horizontal zooms-out
	for(i=1;pow(2,i)<=w;i++){
		w1=w/pow(2,i);
		for(u=0;u<w1;u++){
			for(v=0;v<h;v++){
				for(l=0;l<3;l++){    
					r[coord(i,0,u,v,w,h,l)] = gaussian_filter_h(r,i,0,u,v,w,h,l,gauss1D);
				}
			}
		}
	}

	//for each previous rectangle, compute the rectangles obtained by successive vertical zooms-out
	for(i=0;pow(2,i)<=w;i++){
		w1=w/pow(2,i);
		for(j=1;pow(2,j)<=h;j++){
			h2=h/pow(2,j);
			for(u=0;u<w1;u++){
				for(v=0;v<h2;v++){
					for(l=0;l<3;l++){
						r[coord(i,j,u,v,w,h,l)] = gaussian_filter_v(r,i,j,u,v,w,h,l,gauss1D);
					}
				}
			}
		}
	}
}



//the distance function is the smallest rectangle containing the pre-image (by the locally affine-approximated warping)
//the rectangle is of size dw = |du/dx|+|du/dy| and dh = |dv/dx| + |dv/dy|


// notations :
// c = D[9]x+D[10]y+D[11]
// du/dx = D[1]y+D[2]/c
// du/dy = D[5]x+D[6]/c
// dv/dx = D[3]y+D[4]/c
// dv/dy = D[7]x+D[8]/c


void precal_D(double H[3][3],double *D){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1];
	D[2]=a[0]*a[8]-a[6]*a[2];
	D[3]=a[3]*a[7]-a[6]*a[4];
	D[4]=a[3]*a[8]-a[6]*a[5];
	D[5]=-D[1];
	D[6]=a[1]*a[8]-a[7]*a[2];
	D[7]=-D[3];
	D[8]=a[4]*a[8]-a[7]*a[5];
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}



//bilinear interpolation on the RipMap
static float bilinear_ripmap(float *x, int w, int h,
		float p, float q, int l, int d1, int d2){
	//a rectangle represent several pixels ; thus the coordinates p and q must be shifted
	//the center of the rectangle must be the coordinates on which to evaluate the image
	float pp = p - 0.5*(pow(2,d1)-1)/pow(2,d1);
	float qq = q - 0.5*(pow(2,d2)-1)/pow(2,d2);
	int ip = floor(pp);
	int iq = floor(qq);
	float a = x[coord(d1,d2,ip,iq,w,h,l)];
	float b = x[coord(d1,d2,ip+1,iq,w,h,l)];
	float c = x[coord(d1,d2,ip,iq+1,w,h,l)];
	float dd = x[coord(d1,d2,ip+1,iq+1,w,h,l)];
	return evaluate_bilinear_cell(a, b, c, dd, pp-ip, qq-iq);
}


//la fonction principale. On a distingué beaucoup de cas suivant dw<=1, dw>=logw, etc...

static float ripmap_interpolation_at(float *r, int w, int h,
		float x, float y, int l,double d[2]){
		int logw = (int) log2(w);
		int logh = (int) log2(h);
		int dw = floor(log2(d[0])) + 1;
		int dh = floor(log2(d[1])) + 1;
		float a,b,c,dd,u,v;
		if((pow(2,dw))>w && (pow(2,dh)>h)){return r[coord(logw,logh,0,0,w,h,l)];}
		if(dw<=1 && dh<=1){return bilinear_ripmap(r,w,h,x,y,l,0,0);}
		if(pow(2,dw)>w){
			if(dh<=1){return bilinear_ripmap(r,w,h,0,y,l,logw,0);}
			{
				a=bilinear_ripmap(r,w,h,0,y/pow(2,dh-1),l,logw,dh-1);
				b=bilinear_ripmap(r,w,h,0,y/pow(2,dh-2),l,logw,dh-2);
				return (dh-log2(d[1]))*b + (log2(d[1])-dh+1)*a;
			}
		}
		if(pow(2,dh)>h){
			if(dw<=1){return bilinear_ripmap(r,w,h,x,0,l,0,logh);}
			{
				a=bilinear_ripmap(r,w,h,x/pow(2,dw-1),0,l,dw-1,logh);
				b=bilinear_ripmap(r,w,h,x/pow(2,dw-2),0,l,dw-2,logh);
				return (dw-log2(d[0]))*b + (log2(d[0])-dw+1)*a;
			}
		}
		if(dw<=1){
			a=bilinear_ripmap(r,w,h,x,y/pow(2,dh-1),l,0,dh-1);
			b=bilinear_ripmap(r,w,h,x,y/pow(2,dh-2),l,0,dh-2);
			return (dh-log2(d[1]))*b + (log2(d[1])-dh+1)*a;}
		if(dh<=1){
			a=bilinear_ripmap(r,w,h,x/pow(2,dw-1),y,l,dw-1,0);
			b=bilinear_ripmap(r,w,h,x/pow(2,dw-2),y,l,dw-2,0);
			return (dw-log2(d[0]))*b + (log2(d[0])-dw+1)*a;
		
		}
		a=bilinear_ripmap(r,w,h,x/pow(2,dw-1),y/pow(2,dh-1),l,dw-1,dh-1);
		b=bilinear_ripmap(r,w,h,x/pow(2,dw-2),y/pow(2,dh-1),l,dw-2,dh-1);
		c=bilinear_ripmap(r,w,h,x/pow(2,dw-1),y/pow(2,dh-2),l,dw-1,dh-2);
		dd=bilinear_ripmap(r,w,h,x/pow(2,dw-2),y/pow(2,dh-2),l,dw-2,dh-2);
		u=(pow(2,dw)-d[0])/pow(2,dw);
		v=(pow(2,dh)-d[1])/pow(2,dh);
		return u*v*dd + (1-u)*v*b + (1-v)*u*c + (1-u)*(1-v)*a;
	}
	
	
	
	
	
int apply_homo_ripmap(float *img,float *img_f,int w,int h,int w_out,int h_out,double H[3][3]){
	
	double D[12];
	precal_D(H,D);

//on va construire le ripmap de img.
	float *ripmap=malloc(12*w*w*sizeof(float));
	build_ripmap(img,ripmap,w,h,3);


	for (int j = 0; j < h_out; j++)
	for (int i = 0; i < w_out; i++)
	{
		double p[2] = {i, j};
		apply_homography_1pt(p, H, p);
		
		double z[2]={i,j};
		double a;
		a = pow(D[9]*z[0]+D[10]*z[1]+D[11],2);	
		double dudx,dudy,dvdx,dvdy;
		dudx = (D[1]*z[1]+D[2])/a;
		dudy = (D[5]*z[0]+D[6])/a;
		dvdx = (D[3]*z[1]+D[4])/a;
		dvdy = (D[7]*z[0]+D[8])/a;
	
		double det = fabs( dudy*dvdx - dvdy*dudx );
		double d[2];
		d[0] = fabs(dudx) + fabs(dudy);
		d[1] = fabs(dvdx) + fabs(dvdy);

		//on doit faire une copie de p, sinon il ne veut pas décaler...
		double coo[2];
		coo[0]=p[0]; coo[1]=p[1];
	
		//On decale l'origine du rectangle pour que le parallélogramme soit à l'intérieur	
		if(dudx<0){coo[0] += dudx;}
		if(dudy<0){coo[0] += dudy;}
		if(dvdx<0){coo[1] += dvdx;}
		if(dvdy<0){coo[1] += dvdy;}
		
		for (int l = 0; l < 3; l++)
		{
			int idx = l + 3 * (w_out * j + i);
			float v = ripmap_interpolation_at(ripmap, w, h, coo[0], coo[1], l, d);
			img_f[idx] = v;
		}
	}
	
	
	//truncate in the end to have non-periodic result
	double p[2];
	for(int i=0;i<w_out;i++){
		for(int j=0;j<h_out;j++){
			p[0]=i; p[1]=j;
			apply_homography_1pt(p,H,p);
			if(p[0]<0 || p[0]>w || p[1]<0 || p[1]>h ){
				for(int l=0;l<3;l++){img_f[(j*w_out+i)*3+l]=0;}
			}
		}
	}
	
	return 0;
		
}
