// On fait un zoom in par zero padding.


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <fftw3.h>

#include "aux_fun.h"


int zoom_in (float *img, int w,int h,int pd,int kw,int kh,float *img_final)
{ void fourierForward(float* in,float* reOut,float* imOut,unsigned int w,unsigned int h);
  void fourierBackward(float* reIn,float* imIn,float* out,unsigned int w,unsigned int h);
/**
  * img the input image, of size w*h*pd
  * img_final the output image, of size kw*kh*pd
  */



	float *refftimg = malloc(w*h*sizeof(float));
  if(refftimg==NULL){
        printf("@fft_zoom @zoom : refftimg is too large");
        exit(1);
    }
	float *imfftimg = malloc(w*h*sizeof(float));
  if(imfftimg==NULL){
        printf("@fft_zoom @zoom : imfftimg is too large");
        exit(1);
    }
	float *img_aux = malloc(w*h*sizeof(float));
  if(img_aux==NULL){
        printf("@fft_zoom @zoom : img_aux is too large");
        exit(1);
    }
	float *refftimg_final = malloc(kw*kh*sizeof(float));
  if(refftimg_final==NULL){
        printf("@fft_zoom @zoom : refftimg_final is too large");
        exit(1);
    }
	float *imfftimg_final =malloc(kw*kh*sizeof(float));
  if(imfftimg_final==NULL){
        printf("@fft_zoom @zoom : imfftimg_final is too large");
        exit(1);
    }
	float *img_final_aux = malloc(kw*kh*sizeof(float));
  if(img_final_aux==NULL){
        printf("@fft_zoom @zoom : img_final_aux is too large");
        exit(1);
    }
	
  for(int l=0;l<pd;l++){ 

    //extract a channel of color
    for(int i=0;i<w*h;i++){
      img_aux[i]=img[pd*i+l];
    }

    //Fourier transform
    fourierForward(img_aux,refftimg,imfftimg,w,h);

    //zero-padding
    for (int i=0;i<kw;i++)
      for (int j=0;j<kh;j++){  
        refftimg_final[i+j*kw]=0; 
        imfftimg_final[i+j*kw]=0;
    }
    for (int i=0;i<w;i++)
      for (int j=0;j<h;j++){
        int idx = (kw-w)/2+i+kw*((kh-h)/2+j);
        refftimg_final[idx]=refftimg[w*j+i];     
        imfftimg_final[idx]=imfftimg[w*j+i];
    }

    //inverse Fourier transform
    fourierBackward(refftimg_final, imfftimg_final, img_final_aux, kw,kh);

    //output
    for(int i=0;i<kw;i++)
      for(int j=0;j<kh;j++){
        img_final[(i+j*kw)*pd+l]= ((float)kw/(float)w)*((float)kh/(float)h)*img_final_aux[i+j*kw];
    }

  }

  free(refftimg);
  free(imfftimg);
  free(img_aux);
	free(refftimg_final);
	free(imfftimg_final);
	free(img_final_aux);
	
  return 0 ;
}





//Fourier transform
void fourierForward(float* in, //input image
                    float* reOut, //real part of FFT of the output
                    float* imOut, //imaginary part of FFT of the output
                    unsigned int w, //width
                    unsigned int h) //height
{
  fftw_complex* spatial_repr;
  fftw_complex* frequency_repr;
  unsigned int i;
  unsigned int j;
  fftw_plan plan;
  int x,y;

  spatial_repr= fftw_malloc(sizeof(fftw_complex)*w*h);
  if(spatial_repr==NULL){
        printf("@fft_zoom @fourierForward : spatial_repr is too large");
        exit(1);
  }
  frequency_repr= fftw_malloc(sizeof(fftw_complex)*w*h);
  if(frequency_repr==NULL){
        printf("@fft_zoom @fourierForward : frequency_repr is too large");
        exit(1);
  }

  for(i=0;i<w*h;i++){
    spatial_repr[i][0] = in[i];
    spatial_repr[i][1] =  0.0f;
  }

  //compute execution plane
  plan=fftw_plan_dft_2d(h, w, spatial_repr, frequency_repr, FFTW_FORWARD, FFTW_ESTIMATE);
  if(plan==NULL){
        printf("@fft_zoom @fourierForward : plan has not been created");
        exit(1);
  }

  //compute FFT
  fftw_execute(plan);

  //recentered output
  for(j=0;j<h;j++)
    for(i=0;i<w;i++){
      x=good_modulus(i+w/2,w);
      y=good_modulus(j+h/2,h);
      reOut[y*w+x]=frequency_repr[j*w+i][0];
      imOut[y*w+x]=frequency_repr[j*w+i][1];
  }

  fftw_destroy_plan(plan);
  fftw_free(spatial_repr);
  fftw_free(frequency_repr);
 
}
		
 



//inverse Fourier transform
void fourierBackward(float* reIn, //real part of FFT of the input 
                     float* imIn, //imaginary part of FFT of the input
                     float* out, //output image
                     unsigned int w, //width
                     unsigned int h) //height
{
  fftw_complex* spatial_repr;
  fftw_complex* frequency_repr;
  unsigned int i;
  unsigned int j;
  int x,y;
  fftw_plan plan;

  spatial_repr= fftw_malloc(sizeof(fftw_complex)*w*h);
  if(spatial_repr==NULL){
        printf("@fft_zoom @fourierBackward : spatial_repr is too large");
        exit(1);
  }
  frequency_repr= fftw_malloc(sizeof(fftw_complex)*w*h);
  if(frequency_repr==NULL){
        printf("@fft_zoom @fourierBackward : frequency_repr is too large");
        exit(1);
  }

  //uncentered input
  for(j=0;j<h;j++)
    for(i=0;i<w;i++){
      x=good_modulus(i-w/2,w);
      y=good_modulus(j-h/2,h);
      frequency_repr[j*w+i][0]=reIn[y*w+x];
      frequency_repr[j*w+i][1]=imIn[y*w+x];
  }
  
  //compute execution plane
  plan=fftw_plan_dft_2d(h, w, frequency_repr, spatial_repr, FFTW_BACKWARD, FFTW_ESTIMATE);
  if(plan==NULL){
        printf("@fft_zoom @fourierBackward : plan has not been created");
        exit(1);
  }

  //compute inverse FFT
  fftw_execute(plan);

  //normalize and convert complex data to real data
  for(i=0;i<w*h;i++){
    out[i]=spatial_repr[i][0]/(float) (w*h);
  }

  fftw_destroy_plan(plan);
  fftw_free(spatial_repr);
  fftw_free(frequency_repr);
}
