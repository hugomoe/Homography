#ifndef AUX_FUN
#define AUX_FUN

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"

//a relaxed equality on doubles
bool eq(double a,double b){if(a<b+pow(2,-PREC_EQ) && a>b-pow(2,-PREC_EQ)){return true;}{return false;}}

//a good definition of n%p, even for huge and negative numbers
int good_modulus(int nn, int p)
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





//MANIPULATING 3x3 MATRIX AND HOMOGRAPHIES
//transport a point by an homography
void apply_homography_1pt(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}
// compute the inverse homography (inverse of a 3x3 matrix)
double invert_homography(double invH[3][3], double H[3][3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = (a[4]*a[8]-a[5]*a[7])/det;
	r[1] = (a[2]*a[7]-a[1]*a[8])/det;
	r[2] = (a[1]*a[5]-a[2]*a[4])/det;
	r[3] = (a[5]*a[6]-a[3]*a[8])/det;
	r[4] = (a[0]*a[8]-a[2]*a[6])/det;
	r[5] = (a[2]*a[3]-a[0]*a[5])/det;
	r[6] = (a[3]*a[7]-a[4]*a[6])/det;
	r[7] = (a[1]*a[6]-a[0]*a[7])/det;
	r[8] = (a[0]*a[4]-a[1]*a[3])/det;
	return det;
}

// C = AoB, composition of two homographies (product of 3x3 matrices)
void compose_homographies(double C[3][3], double A[3][3], double B[3][3]){
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 3; k++)
			C[i][j] += A[i][k] * B[k][j];
	}
}

// Find the homography that changes the canonical projective basis into the given four points (x, y, z, w)
void homography_from_four_points(double H[3][3], double x[2], double y[2], double z[2], double w[2]){
	// fix the degree of freedom (assuming the four points are finite)
	double t = 1;

	// translation coefficients
	double p = x[0];
	double q = x[1];

	// "core" 2x2 system
	double A = w[0] - z[0];
	double B = y[0] - z[0];
	double C = w[1] - z[1];
	double D = y[1] - z[1];
	double P = z[0] - y[0] - w[0] + p;
	double Q = z[1] - y[1] - w[1] + q;
	double DET = A * D - B * C;
	double r = (D * P - B * Q) / DET;
	double s = (A * Q - C * P) / DET;
	if (!isnormal(DET))
		fprintf(stderr, "denormal! DET = %g\n", DET);

	// solve the rest of the diagonal system
	double a = w[0] * ( 1 + r ) - p;
	double b = y[0] * ( 1 + s ) - p;
	double c = w[1] * ( 1 + r ) - q;
	double d = y[1] * ( 1 + s ) - q;

	// fill-in the output
	H[0][0] = a; H[0][1] = b; H[0][2] = p;
	H[1][0] = c; H[1][1] = d; H[1][2] = q;
	H[2][0] = r; H[2][1] = s; H[2][2] = t;
}

// Find the homography that moves the four points (x,y,z,t) to (a,b,c,d)
void homography_from_eight_points(double H[3][3], double x[2], double y[2], double z[2], double w[2], double a[2], double b[2], double c[2], double d[2]){
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}





//BILINEAR INTERPOLATION
//extrapolate by a constant value
float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// bilinear interpolation
float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = getsample_cons(x, w, h, pd, ip  , iq  , l);
	float b = getsample_cons(x, w, h, pd, ip+1, iq  , l);
	float c = getsample_cons(x, w, h, pd, ip  , iq+1, l);
	float d = getsample_cons(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}




#endif //AUX_FUN
