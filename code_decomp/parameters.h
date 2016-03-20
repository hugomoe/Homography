#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//the support of the filter is included in [-SUPP,SUPP]
#define SUPP 4

//parameters of the filter
#define PERIOD 1. //period of the raised cosine-weighted sinc
#define BETA 0.36 //roll-off factor of the raised cosine-weighted sinc

//precision of the computation
#define PREC 10 //the precomputation of the values of the filter will be at precision 2^(-PREC)
#define PREC2 20 //an equality will be at precision 2^(-PREC2)
#define PI 3.14159265358979323

//a relaxed equality on doubles
bool eq(double a,double b){if(a<b+pow(2,-PREC2) && a>b-pow(2,-PREC2)){return true;}{return false;}}
