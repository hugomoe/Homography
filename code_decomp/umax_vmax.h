#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"



//intersect a segment with vertical edges of the square [-1,1]*[-1,1]
void avi(double M[2][10], int *n, double u1, double v1, double u2, double v2){
/**
  * @name
  *     add_vertical_intersections
  * @param
  *     M : array containing coordinates of n points
  *         the points are stocked between 0 and n-1
  *     P1=(u1,v1), P2=(u2,v2) : [P1,P2] is the segment to intersect
  *         eq(u1,u2)=false
  */
    //a et b are the intersection points between the straight lines
    double ua = 1.;
    double ub = -1.;
    double va = v1+(ua-u1)*(v2-v1)/(u2-u1);
    double vb = v1+(ub-u1)*(v2-v1)/(u2-u1);

    //add a and b to M iff they belong to the intersected segments
    if(fabs(va)<=1. && ((u1<=ua && ua<=u2) || (u2<=ua && ua<=u1))){
        M[0][*n] = ua;
        M[1][*n] = va;
        *n = *n+1;
    }
    if(fabs(vb)<=1. && ((u1<=ub && ub<=u2) || (u2<=ub && ub<=u1))){
        M[0][*n] = ub;
        M[1][*n] = vb;
        *n = *n+1;
    }
}

//intersect a segment with horizontal edges of the square [-1,1]*[-1,1]
void ahi(double M[2][10], int *n, double u1, double v1, double u2, double v2){
/**
  * @name
  *     add_horizontal_intersections
  * @param
  *     M : array containing coordinates of n points
  *         the points are stocked between 0 and n-1
  *     P1=(u1,v1), P2=(u2,v2) : [P1,P2] is the segment to intersect
  *         eq(v1,v2)=false
  */
    //a et b are the intersection points between the straight lines
    double va = 1.;
    double vb = -1.;
    double ua = u1+(va-v1)*(u2-u1)/(v2-v1);
    double ub = u1+(vb-v1)*(u2-u1)/(v2-v1);

    //add a and b to M iff they belong to the intersected segments
    if(fabs(ua)<=1. && ((v1<=va && va<=v2) || (v2<=va && va<=v1))){
        M[0][*n] = ua;
        M[1][*n] = va;
        *n = *n+1;
    }
    if(fabs(ub)<=1. && ((v1<=vb && vb<=v2) || (v2<=vb && vb<=v1))){
        M[0][*n] = ub;
        M[1][*n] = vb;
        *n = *n+1;
    }
}

int umax_vmax(double *u, double *v, double A[2][2]){
/**
  * @param
  *     u, v : output for u_max and v_max
  *     A : 2x2 matrix such that img_f(x) = img(Ax)
  * @return
  *     1 iff error
  */
    double detA = A[0][0]*A[1][1]-A[0][1]*A[1][0];
    //a is the inverse of the transposed matrix of A
    double a[2][2] = {A[0][0]/detA, -A[1][0]/detA, -A[0][1]/detA, A[1][1]/detA};
    //(u1,v1) = a(1,1)
    double u1 = a[0][0]+a[0][1], v1 = a[1][0]+a[1][1];
    //(u2,v2) = a(1,-1)
    double u2 = a[0][0]-a[0][1], v2 = a[1][0]-a[1][1];

    int n = 0; //index where to put the next found point
    double M[2][10];

    //if a point is already inside the square, add them before intersecting
    if(fabs(u1)<=1 && fabs(v1)<=1){
        M[0][n]=u1;
        M[1][n]=v1;
        n++;
    }
    if(fabs(u2)<=1 && fabs(v2)<=1){
        M[0][n]=u2;
        M[1][n]=v2;
        n++;
    }

    if(eq(u1,u2)){ //the two points are vertically align ; do not call avi(M,&n,u1,v1,u2,v2);
        if(eq(v1,v2)){ //the two points are the same point
            printf("@umax_vmax : A is not invertible\n");
            return 1;
        }else if(eq(v1,-v2)){ //a([-1,1]*[-1,1]) is a rectangle
            *u = fmin(fabs(u1),1.);
            *v = fmin(fabs(v1),1.);
            return 0;
        }else{
            ahi(M,&n,u1,v1,u2,v2);
            ahi(M,&n,u1,v1,-u2,-v2);
            avi(M,&n,u1,v1,-u2,-u2);
        }
    }else if(eq(u1,-u2)){ //same cases than before but with (u1,v1) and (-u2,-v2)
        if(eq(v1,-v2)){
            printf("@umax_vmax : A is not invertible\n");
            return 1;
        }else if(eq(v1,v2)){
            *u = fmin(fabs(u1),1.);
            *v = fmin(fabs(v1),1.);
            return 0;
        }else{
            ahi(M,&n,u1,v1,-u2,-v2);
            ahi(M,&n,u1,v1,u2,v2);
            avi(M,&n,u1,v1,u2,u2);
        }
    }else{ //general case
        avi(M,&n,u1,v1,u2,v2);
        avi(M,&n,u1,v1,-u2,-v2);
        if(eq(v1,v2)){
            ahi(M,&n,u1,v1,-u2,-v2);
        }else if(eq(v1,-v2)){
            ahi(M,&n,u1,v1,u2,v2);
        }else{
            ahi(M,&n,u1,v1,u2,v1);
            ahi(M,&n,u1,v1,-u2,-v2);
        }
    }
    if(n==0){ //if no point has been found, then u_max=v_max=1
        *u = 1;
        *v = 1;
        return 0;
    }else{ //u_max, v_max are the maximal coordinates of all found points
        *u=fabs(M[0][0]);
        *v=fabs(M[1][0]);
        int i;
        for(i = 0;i<n;i++){
            *u=fmax(*u,fabs(M[0][i]));
            *v=fmax(*v,fabs(M[1][i]));
        }
        return 0;
    }
}
