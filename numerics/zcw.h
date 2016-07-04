/**  

   \file zcw.h

   

   (c) 2003 Thomas Neff

*/


#ifndef _ZCW_H
#define _ZCW_H

int nangles2(int i);
int nangles3(int i);

void getangles2(int idx, int i, double* alpha, double* beta);
void getangles3(int idx, int i, double* alpha, double* beta, double* gamma);

#endif
