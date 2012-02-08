
#include <stdlib.h>
#include "gadem.h"

typedef struct mask_repeat MASK;
struct mask_repeat {
   int start,end,id;
};

//mask sites that have been found
void mask_sites(int predictedSiteCn,char **seq,char **rseq,int *seqLen,Sites *site,int pwmLen) {

   register int i,m;
   int pos,seqID;

   for (i=0; i<predictedSiteCn; i++) {
      pos=site[i].pos; seqID=site[i].seq;
      if (site[i].rev=='0') {
         for (m=0; m<pwmLen; m++) seq[seqID][pos+m]='n';
         for (m=0; m<pwmLen; m++) rseq[seqID][seqLen[seqID]-pwmLen+m-pos]='n';
      }
      else {
         for (m=0; m<pwmLen; m++) rseq[seqID][pos+m]='n';
         for (m=0; m<pwmLen; m++) seq[seqID][seqLen[seqID]-pwmLen+m-pos]='n';
      }
   }
}

void mask_repetitive(int *geneID,char **seq,int numSeq,int *seqLen,char *fileName) {

   register int i,j,k,l;
   char **kmer,*s1;
   char *maskedFileName;
   int maxNumKmer,maxKmerLen,klen,numKmer,pos,id,cn;
   // MASK *mask;
   //FILE *fp;

   maxKmerLen=20;
   maxNumKmer=30;
   kmer=alloc_char_char(maxNumKmer,maxKmerLen+1);
   s1=alloc_char(maxKmerLen+1);

   numKmer=5; klen=8;
   strcpy(kmer[0],"aaaaaaaa");
   strcpy(kmer[1],"tttttttt");
   strcpy(kmer[2],"cacacaca");
   strcpy(kmer[3],"tgtgtgtg");
   strcpy(kmer[4],"tatatata");

   // mask=(MASK *)calloc(numSeq*5,sizeof(MASK));
   // if (!mask) { printf("calloc for mask failed!\n"); exit(0); }

   maskedFileName=alloc_char(500);
   id=-1;
   for (i=0; i<strlen(fileName); i++) {
      if (fileName[i]=='/') id=i; 
   }
   if (id==-1) strcpy(maskedFileName,fileName);
   else {
      for (k=0,i=id+1; i<strlen(fileName); i++,k++) maskedFileName[k]=fileName[i]; maskedFileName[k]='\0'; 
   }
   strcat(maskedFileName,".mask");
   //fp=fopen(maskedFileName,"w");

   cn=0;
   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]-klen+1; j++) {
         for (k=0; k<klen; k++) s1[k]=seq[i][j+k]; s1[k]='\0';
         for (l=0; l<numKmer; l++) {
            pos=0;
            while (strncmp(s1,kmer[l],klen)==0) {
               switch (l) {
                  case 0: pos=pos+1; break;
                  case 1: pos=pos+1; break;
                  case 2: pos=pos+2; break;
                  case 3: pos=pos+2; break;
                  case 4: pos=pos+2; break;
                  default: break;
               }
               // pos++;
               for (k=0; k<klen; k++) s1[k]=seq[i][j+k+pos]; s1[k]='\0';
            };
            if (pos!=0) {
              // fprintf(fp,"%d:%d-%d\t",geneID[i],j+1,j+pos+klen-1);
              // for (k=0; k<pos+klen-1; k++) fprintf(fp,"%c",seq[i][j+k]); fprintf(fp,"\n");
               for (k=0; k<pos+klen-1; k++) seq[i][j+k]='n';
            }
         }
      }
   }

   numKmer=10; klen=12;
   strcpy(kmer[0],"ggaggaggagga");
   strcpy(kmer[1],"gaggaggaggag");
   strcpy(kmer[2],"agaagaagaaga");
   strcpy(kmer[3],"ctcctcctcctc");
   strcpy(kmer[4],"tcctcctcctcc");
   strcpy(kmer[5],"tcttcttcttct");
   strcpy(kmer[6],"tagtagtagtag");
   strcpy(kmer[7],"aataataataat");
   strcpy(kmer[8],"attattattatt");
   strcpy(kmer[9],"ataataataata");

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]-klen+1; j++) {
         for (k=0; k<klen; k++) s1[k]=seq[i][j+k]; s1[k]='\0';
         for (l=0; l<numKmer; l++) {
            pos=0;
            while (strncmp(s1,kmer[l],klen)==0) {
               switch (l) {
                  case 0: pos=pos+3; break;  // three letter repeats
                  case 1: pos=pos+3; break;  // three letter repeats
                  case 2: pos=pos+3; break;  // three letter repeats
                  case 3: pos=pos+3; break;  // three letter repeats
                  case 4: pos=pos+3; break;  // three letter repeats
                  case 5: pos=pos+3; break;  // three letter repeats
                  case 6: pos=pos+3; break;  // three letter repeats
                  case 7: pos=pos+3; break;  // three letter repeats
                  case 8: pos=pos+3; break;  // three letter repeats
                  case 9: pos=pos+3; break;  // three letter repeats
                  default: break;
               }
               //pos++;
               for (k=0; k<klen; k++) s1[k]=seq[i][j+k+pos]; s1[k]='\0';
            };
            if (pos!=0) {
               //fprintf(fp,"%d:%d-%d\t",geneID[i],j+1,j+pos+klen-1);
               //for (k=0; k<pos+klen-1; k++) fprintf(fp,"%c",seq[i][j+k]); fprintf(fp,"\n");
               for (k=0; k<pos+klen-1; k++) seq[i][j+k]='n';
            }
         }
      }
   }

   numKmer=1; klen=15;
   strcpy(kmer[0],"cagcagcagcagcag");

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]-klen+1; j++) {
         for (k=0; k<klen; k++) s1[k]=seq[i][j+k]; s1[k]='\0';
         for (l=0; l<numKmer; l++) {
            pos=0;
            while (strncmp(s1,kmer[l],klen)==0) {
               switch (l) {
                  case 0: pos=pos+3; break;  // three letter repeats
                  default: break;
               }
               // pos++;
               for (k=0; k<klen; k++) s1[k]=seq[i][j+k+pos]; s1[k]='\0';
            };
            if (pos!=0) {
               //fprintf(fp,"%d:%d-%d\t",geneID[i],j+1,j+pos+klen-1);
               //for (k=0; k<pos+klen-1; k++) fprintf(fp,"%c",seq[i][j+k]); fprintf(fp,"\n");
               for (k=0; k<pos+klen-1; k++) seq[i][j+k]='n';
            }
         }
      }
   }
   /*--------------------------------------
   numKmer=1; klen=16;
   strcpy(kmer[0],"catatatacatatata");

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]-klen+1; j++) {
         for (k=0; k<klen; k++) s1[k]=seq[i][j+k]; s1[k]='\0';
         for (l=0; l<numKmer; l++) {
            pos=0;
            while (strncmp(s1,kmer[l],klen)==0) {
               switch (l) {
                  case 0: pos=pos+2; break;  
                  default: break; 
               }
               pos++;
               for (k=0; k<klen; k++) s1[k]=seq[i][j+k+pos]; s1[k]='\0';
            };
            if (pos!=0) { 
               for (k=0; k<pos+klen-1; k++) seq[i][j+k]='n';
            }
         }
      }
   }
   --------------------------------------*/
 
   /*
   maskedFileName=alloc_char(500);
   id=-1;
   for (i=0; i<strlen(fileName); i++) {
      if (fileName[i]=='/') id=i; 
   }
   if (id==-1) strcpy(maskedFileName,fileName);
   else {
      for (k=0,i=id+1; i<strlen(fileName); i++,k++) maskedFileName[k]=fileName[i]; maskedFileName[k]='\0'; 
   }
   strcat(maskedFileName,".mask");
   fp=fopen(maskedFileName,"w");

   for (i=0; i<cn; i++) {
      fprintf(fp,"%s\n",geneID[mask[i].id]);
      for (j=mask[i].start; j<mask[i].end; j++) {
         fprintf(fp,"%c",seq[mask[i].id][j]); 
      } fprintf(fp,"\n");
   }
   */
   /*-----------------------------------------------
   for (i=0; i<numSeq; i++) {
      fprintf(fp,"%s\n",geneID[i]);
      for (j=0; j<seqLen[i]; j++) {
         fprintf(fp,"%c",seq[i][j]);
         if ((j+1)%50==0) fprintf(fp,"\n"); 
      }
   }
   ------------------------------------------------*/
   //fclose(fp);

   if (kmer[0]) { free(kmer[0]); kmer[0]=NULL; }
   if (kmer)    { free(kmer);    kmer=NULL;    }
   if (s1)      { free(s1);      s1=NULL;      }
   if (maskedFileName)  { free(maskedFileName); maskedFileName=NULL; }
}

