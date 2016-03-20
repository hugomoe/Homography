#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "iio.h"
#include "iio.c"

/**
  * necessite iio.c et iio.h dans le même dossier
  *
  * mode d'emploi : mettre tous les fichiers .pgm dans le dossier
  * appeler ./pgmTo [img1.pgm] [img2.pgm] [img3.pgm] etc
  * les fichiers seront convertis en pgm (pour changer le format, changer directement dans le code)
  *
  * compilation :
  * c99 -O3 -DNDEBUG pgmTo.c -lX11 -ljpeg -ltiff -lpng -o pgmTo
  */


static int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}

int main(int argc, char *argv[]){
	if (argc < 2) {
		fprintf(stderr, "usage:\n\t%s [image.extension]\n", *argv);
		return 1;
	}

    int w,h,pd;
	for(int i=1;i<argc;i++){
        //printf("argc = %d\n",argc);
        char *filename_in = argv[i];
        float *img_in = iio_read_image_float_vec(filename_in,&w,&h,&pd);

        char filename_out[FILENAME_MAX];
        sprintf(filename_out,filename_in);
        char *extension = strchr(filename_out,'.');
        while(extension!=NULL){
            extension[0] = '_';
            extension = strchr(filename_out,'.');
        }
        extension = strrchr(filename_out,'_');
        sprintf(extension,".png");

        iio_save_image_float_vec(filename_out,img_in,w,h,pd);

        free(img_in);
        //free(extension);
	}
	return 0;
}
