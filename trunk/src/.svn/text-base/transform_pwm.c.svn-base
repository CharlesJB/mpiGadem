#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "gadem.h"

double **alloc_double_double(int ,int );

// independent motif model
void standardize_pwm(double **pwm,int pwmLen) {

   register int i,j;
   double sum;

   for (i=0; i<pwmLen; i++) {
      sum=0.0; for (j=0; j<4; j++)  sum +=pwm[i][j]; 
      if (sum>0.01) { 
         for (j=0; j<4; j++) pwm[i][j] /=sum;  
         for (j=0; j<4; j++) pwm[i][j]=(pwm[i][j]+PSEUDO_COUNT)/(1.0+PSEUDO_COUNT*4);
      }
      else { 
         for (j=0; j<4; j++) pwm[i][j]=0.25; 
      }
   }
}

// independent motif model
void log_pwm(double **pwm,double **logpwm,int pwmLen){

   register int i,j;

   for (i=0; i<pwmLen; i++) {
      for (j=0; j<4; j++) logpwm[i][j]=log((pwm[i][j]+PSEUDO_COUNT)/(1.0+PSEUDO_COUNT*4));
   }
}

// converting log likelihood ratio matrix to an integer matrix
void log_ratio_to_int(double **pwm,int **ipwm,int pwmLen,double *bfreq){

   register int i,j;
   double cellMin,cellMax;
   double **logodds;

   logodds=alloc_double_double(pwmLen,4);

   for (i=0; i<pwmLen; i++) {
      for (j=0; j<4; j++) logodds[i][j]=log((pwm[i][j]+PSEUDO_COUNT)/(1.0+PSEUDO_COUNT*4)/bfreq[j]);
   }

   cellMin=logodds[0][0]; cellMax=logodds[0][0];
   for (i=0; i<pwmLen; i++) {
      for (j=0; j<4; j++) {
        if (logodds[i][j]<cellMin) cellMin=logodds[i][j];
        if (logodds[i][j]>cellMax) cellMax=logodds[i][j];
      }
   }

   if (cellMax-cellMin<0.01) {
      for (i=0; i<pwmLen;  i++) {
         for (j=0; j<4; j++) ipwm[i][j]=DOUBLE_TO_INT_SCALE/4;
      }
   }
   else {
      for (i=0; i<pwmLen;  i++) {
         for (j=0; j<4; j++) {
            ipwm[i][j]=(int)(DOUBLE_TO_INT_SCALE*((logodds[i][j]-cellMin)/(cellMax-cellMin)));
         }
      }
   }
   if (logodds[0]) { free(logodds[0]); logodds[0]=NULL; }
   if (logodds)    { free(logodds);    logodds=NULL;    }
}

