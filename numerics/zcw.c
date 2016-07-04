/** 

   \file zcw.c

   ZCW angle sets

   (c) 2003 Thomas Neff

*/

#include <math.h>

#include "zcw.h"


#define NZCWSETS2 9

#define NZCWSETS3 15

int nzcw2[NZCWSETS2] = {21, 34, 55, 89, 144, 233, 377, 616, 987};

int nzcw3[NZCWSETS3] = {1, 10, 21, 35, 44, 98, 142, 200, 301, 504, 1086, 2022, 
			3993, 7991, 15990};


const int g2[NZCWSETS2][2] = {{1,8},
			      {1,13},
			      {1,21},
			      {1,34},
			      {1,55},
			      {1,89},
			      {1,144},
			      {1,233},
			      {1,377}};

const int g3[NZCWSETS3][3] = {{1,1,1},
                              {1,2,4},
			      {1,3,8},
			      {1,11,16},
                              {1,14,20},
			      {1,16,44},
			      {1,14,37},
			      {1,64,94},
			      {1,36,92},
			      {1,77,130},
			      {1,284,446},
			      {1,602,917},
			      {1,439,1867},
			      {1,807,2846},
			      {1,1619,6011}};


int nangles2(int idx)
{
  return nzcw2[idx];
}

int nangles3(int idx)
{
  return nzcw3[idx];
}

void getangles2(int idx, int i, double* alpha, double* beta)
{
  *alpha = 2*M_PI*((i*g2[idx][1]) % nzcw2[idx])/nzcw2[idx];
  *beta  = acos(-1.0 + 2.0*i*g2[idx][0]/nzcw2[idx]);
}

void getangles3(int idx, int i, double* alpha, double* beta, double* gamma)
{
  *alpha = 2*M_PI*((i*g3[idx][1]) % nzcw3[idx])/nzcw3[idx];
  *beta  = acos(1.0 - 2.0*i*g3[idx][0]/nzcw3[idx]);
  *gamma = 2*M_PI*((i*g3[idx][2]) % nzcw3[idx])/nzcw3[idx];
}
