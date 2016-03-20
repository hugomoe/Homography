#ifndef PARAMETERS
#define PARAMETERS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//parameters for raised cosine-weighted sinc
#define SUPP 4 //the support of the filter is included in [-SUPP,SUPP]
#define PERIOD 1. //period of the raised cosine-weighted sinc
#define BETA 0.36 //roll-off factor of the raised cosine-weighted sinc
#define PREC 10 //the precomputation of the values of the filter will be at precision 2^(-PREC)

//parameters for the gaussian kernel in ripmap.h
#define TAPSR 5 //number of non zero values of the kernel (must be odd)
#define SIG 0.6 //variance

//assumed initial and final blurs in homo_box.h
#define INITIAL_BLUR 0.7
#define FINAL_BLUR 0.7

//precision of the computation
#define PREC_EQ 20 //an equality will be at precision 2^(-PREC_EQ)
#define PI 3.14159265358979323

#endif //PARAMETERS
