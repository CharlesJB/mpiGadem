
#include <stdlib.h>
#include <string.h>
#include "defines.h"

double *alloc_double(int );

double *frequency(int numSeq,char **seq,char **rseq,int *seqLen) {

   register int i,j;
   int bcount[4];
   double sum;
   double *freq;

   freq=alloc_double(4);

   for (j=0; j<4; j++) bcount[j]=0;

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]; j++) {
         switch (seq[i][j]) {
            case 'a': (bcount[0])++; break;
            case 'c': (bcount[1])++; break;
            case 'g': (bcount[2])++; break;
            case 't': (bcount[3])++; break;
            default: break;
         }
      }
      for (j=0; j<seqLen[i]; j++) {
         switch (rseq[i][j]) {
            case 'a': (bcount[0])++; break;
            case 'c': (bcount[1])++; break;
            case 'g': (bcount[2])++; break;
            case 't': (bcount[3])++; break;
            default: break;
         }
      }
   }
   for (j=0; j<4; j++) {
      if (bcount[j]==0) freq[j]=PSEUDO_COUNT;
      else              freq[j]=(double)bcount[j]; 
   }
   sum=0; for (j=0; j<4; j++) sum +=freq[j];
   for (j=0; j<4; j++) freq[j] /=sum;

   return (freq);
}

