#include <stdio.h>
#include <stdlib.h>
#include "stdlib.h"
#include <string.h>
#include <math.h>
#include "defines.h"
#include "gadem.h"
#include <Rinternals.h>
#include <Rdefines.h>

void print_bed(Sites *site,int nsites,char **geneID,int *seqLen,int pwmLen,int id) {

   FILE *f1;
   char *fileName,*s1,*chr;
   int s,e,e2,len,start;
   register int i,j,k;

   s1=alloc_char(20);
   chr=alloc_char(20);
   fileName=alloc_char(500);
   sprintf(fileName,"%d.bed",id);
   f1=fopen(fileName,"w");
   for (i=0; i<nsites; i++) {
      len=strlen(geneID[site[i].seq]);

      s=-1; e=-1;
      for (j=0; j<len-3; j++) {
         if (geneID[site[i].seq][j]=='c' && geneID[site[i].seq][j+1]=='h' && geneID[site[i].seq][j+2]=='r') {
            s=j; break;
         }
      }
      for (j=s; j<len; j++) {
         if (geneID[site[i].seq][j]==':') {
            e=j; break;
         }
      }
      if (s!=-1 && e!=-1) {
         for (k=0,j=s; j<e; j++,k++)  chr[k]=geneID[site[i].seq][j];
         chr[k]='\0';
      }
      else {
         printf("%s chr not found! %d %d\n",geneID[site[i].seq],s,e); exit(0);
      }
      e2=-1;
      for (j=e+1; j<len; j++) {
         if (geneID[site[i].seq][j]=='-') {
            e2=j; break;
         }
      }
      if (e2!=-1) {
         for (k=0,j=e+1; j<e2; j++,k++)  s1[k]=geneID[site[i].seq][j];
         s1[k]='\0';
         start=atoi(s1);
      }
      else {
         printf("start not found!\n"); exit(0);
      }

      if (site[i].rev=='0') {
         if (site[i].pos>=0) fprintf(f1,"%s\t%d\t%d\n",chr,site[i].pos+start,site[i].pos+pwmLen+start-1);
      }
      else {
         if (site[i].pos>=0) fprintf(f1,"%s\t%d\t%d\n",chr,seqLen[site[i].seq]-site[i].pos-pwmLen+start,seqLen[site[i].seq]-site[i].pos+start-1);
      }
   }
   fclose(f1);
   if (fileName) { free(fileName); fileName=NULL; }
   if (s1)       { free(s1);       s1=NULL;       }
}

void print_motif(Sites *site,int nsites,char **seq,char **rseq,int *seqLen,int pwmLen,int id,double **opwm) {

   //FILE *f1;
   char *fileName;
   register int i,j;

   fileName=alloc_char(500);
   sprintf(fileName,"%d.seq",id);
  //f1=fopen(fileName,"w");
   for (i=0; i<nsites; i++) {
      if (site[i].rev=='0') {
         if (site[i].pos<0) {
            //for (j=site[i].pos; j<0; j++) fprintf(f1,"x"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(seq[site[i].seq][j]) {
                  //case 'a': fprintf(f1,"a"); break;
                  //case 'c': fprintf(f1,"c"); break;
                  //case 'g': fprintf(f1,"g"); break;
                  //case 't': fprintf(f1,"t"); break;
                  //case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         else {
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(seq[site[i].seq][j]) {
                  //case 'a': fprintf(f1,"a"); break;
                  //case 'c': fprintf(f1,"c"); break;
                  //case 'g': fprintf(f1,"g"); break;
                  //case 't': fprintf(f1,"t"); break;
                  //case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            //for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(f1,"x"); 
         }
         //fprintf(f1,"\n");
      }
      else {
         if (site[i].pos<0) {
            //for (j=site[i].pos; j<0; j++) fprintf(f1,"x"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(rseq[site[i].seq][j]) {
               //   case 'a': fprintf(f1,"a"); break;
                //  case 'c': fprintf(f1,"c"); break;
                //  case 'g': fprintf(f1,"g"); break;
                //  case 't': fprintf(f1,"t"); break;
                //  case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         else {
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(rseq[site[i].seq][j]) {
                 // case 'a': fprintf(f1,"a"); break;
                //  case 'c': fprintf(f1,"c"); break;
                //  case 'g': fprintf(f1,"g"); break;
                //  case 't': fprintf(f1,"t"); break;
                //  case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            //for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(f1,"x"); 
         }
         //fprintf(f1,"\n");
      }
   }
   //fclose(f1);
   if (fileName) { free(fileName); fileName=NULL; }

   // print out individual observed PWM in gadem format
   /*-----------------------------------------------------------------------
      fileName=alloc_char(500);
      sprintf(fileName,"%d.mx",id);
      f1=fopen(fileName,"w");

      fprintf(f1,"4\t%d\n",pwmLen);
      for (i=0; i<4; i++) {
         for (j=0; j<pwmLen; j++) {
            if (j<pwmLen-1) fprintf(f1,"%5.4f\t",opwm[j][i]);
            else fprintf(f1,"%5.4f\n",opwm[j][i]);
         }
      }
      fclose(f1);
      if (fileName) { free(fileName); fileName=NULL; }
   }
   -----------------------------------------------------------------------*/
}

SEXP print_result_R(Sites *site,int nsites,int numSeq,char **seq,char **rseq,int *seqLen,
   double logev,double **opwm,int pwmLen,int id,char *sdyad,char *pwmConsensus,int numCycle,
   double pvalueCutoff,double maxpFactor,int *geneID) {

   register int i,j;
int cn[4];//maxHeaderLen;
   int *seqCn;

	SEXP PWM;
	SEXP seqConsencus;
	SEXP motifname;
	SEXP motifname2;
	SEXP returnData;
	SEXP LengthSequence;
	
	SEXP SequencesIdent;
	SEXP StrandIdent;
	SEXP AccessionIdent;
	SEXP PositionIdent;
	SEXP SeqIden;
	SEXP PValue;
	SEXP GADEMList;

	PROTECT(returnData=NEW_LIST(5));
	PROTECT(GADEMList=NEW_LIST(6));


	PROTECT(PWM=allocMatrix(REALSXP,4,pwmLen));
	PROTECT(seqConsencus=NEW_CHARACTER(1));
	PROTECT(motifname=NEW_INTEGER(1));
	PROTECT(motifname2=NEW_CHARACTER(1));
	PROTECT(SequencesIdent=NEW_CHARACTER(nsites));
	PROTECT(PositionIdent=NEW_INTEGER(nsites));
	PROTECT(SeqIden=NEW_INTEGER(nsites));
	PROTECT(StrandIdent=NEW_CHARACTER(nsites));
	PROTECT(AccessionIdent=NEW_INTEGER(nsites));
	PROTECT(PValue=NEW_NUMERIC(nsites));
	PROTECT(LengthSequence=NEW_INTEGER(nsites));

	int increment_sequence=0; 	
	int compt=0;
   seqCn=alloc_int(numSeq);

   //maxHeaderLen=min(maxHeaderLen,MAX_SEQ_HEADER);

   for (i=0; i<numSeq; i++) seqCn[i]=0;
   for (i=0; i<nsites; i++) seqCn[site[i].seq]++; 
  
   for (i=0; i<4; i++) cn[i]=0; 
   for (i=0; i<numSeq; i++) {
      if (seqCn[i]==0) cn[0]++; 
      if (seqCn[i]==1) cn[1]++; 
      if (seqCn[i]==2) cn[2]++; 
      if (seqCn[i]>2)  cn[3]++; 
   }
   if (seqCn) { free(seqCn); seqCn=NULL; }

   for (i=0; i<nsites; i++) {
	//SET_STRING_ELT(AccessionIdent,increment_sequence,mkChar(geneID[site[i].seq]));
	INTEGER(AccessionIdent)[increment_sequence]=(geneID[site[i].seq]);
		
      if (site[i].rev=='0') {
         if (site[i].pos<0) {
				char sequence_conca[100]="";
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(seq[site[i].seq][j]) {
                 case 'a': strcat(sequence_conca,"A");break;
				 case 'c': strcat(sequence_conca,"C");break;
				 case 'g': strcat(sequence_conca,"G");break;
				 case 't': strcat(sequence_conca,"T");break;
				 case 'n': strcat(sequence_conca,"N");break;
                  default: break;
               }
            }
         }
         else {
		char sequence_conca[100]="";
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(seq[site[i].seq][j]) {
                 case 'a': strcat(sequence_conca,"A");break;
				 case 'c': strcat(sequence_conca,"C");break;
				 case 'g': strcat(sequence_conca,"G");break;
				 case 't': strcat(sequence_conca,"T");break;
				 case 'n': strcat(sequence_conca,"N");break;
                  default: break;
               }
            }
			SET_STRING_ELT(SequencesIdent,increment_sequence,mkChar(sequence_conca));
         }
    
         // print flanking region
         for (j=site[i].pos+pwmLen; j<min(site[i].pos+pwmLen+FLANKING_BASES,seqLen[site[i].seq]); j++) 
     
			SET_STRING_ELT(StrandIdent,increment_sequence,mkChar("+"));
			INTEGER(SeqIden)[increment_sequence]=site[i].seq+1;
			INTEGER(PositionIdent)[increment_sequence]=site[i].pos+1;
			DOUBLE_DATA(PValue)[increment_sequence]=site[i].pvalue;
			increment_sequence=increment_sequence+1;
      }
      else {
  
       if (site[i].pos<0) {
			char sequence_conca[50]="";
            //for (j=site[i].pos; j<0; j++) fprintf(fq,"X"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(rseq[site[i].seq][j]) {
                 case 'a': strcat(sequence_conca,"A");break;
				 case 'c': strcat(sequence_conca,"C");break;
				 case 'g': strcat(sequence_conca,"G");break;
				 case 't': strcat(sequence_conca,"T");break;
				 case 'n': strcat(sequence_conca,"N");break;
                  default: break;
               }
            }
         }
         else {
			char sequence_conca[50]="";
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(rseq[site[i].seq][j]) {
                  case 'a': strcat(sequence_conca,"A");break;
				  case 'c': strcat(sequence_conca,"C");break;
				  case 'g': strcat(sequence_conca,"G");break;
				  case 't': strcat(sequence_conca,"T");break;
				  case 'n': strcat(sequence_conca,"N");break;
                  default: break;
               }
            }
			SET_STRING_ELT(SequencesIdent,increment_sequence,mkChar(sequence_conca));
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            //for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(fq,"X"); 
         }
         // print flanking region
         for (j=site[i].pos+pwmLen; j<min(site[i].pos+pwmLen+FLANKING_BASES,seqLen[site[i].seq]); j++) 
      		SET_STRING_ELT(StrandIdent,increment_sequence,mkChar("-"));
			INTEGER(SeqIden)[increment_sequence]=site[i].seq+1;
			INTEGER(PositionIdent)[increment_sequence]=seqLen[site[i].seq]-site[i].pos;
			DOUBLE_DATA(PValue)[increment_sequence]=site[i].pvalue;
			increment_sequence=increment_sequence+1;

      }
   }


for (int aa=0;aa<pwmLen;aa++)
			{
				for(int bb=0;bb<4;bb++)
				{
					NUMERIC_DATA(PWM)[compt]=opwm[aa][bb];
					compt++;
				}
			}


		SET_STRING_ELT(seqConsencus,0,mkChar(pwmConsensus));
		INTEGER(LengthSequence)[0]=125;
		INTEGER(motifname)[0]=id;

		 const char base[] = "m";
		      char filename [ FILENAME_MAX ];
		      int number = id;
		      sprintf(filename, "%s%d", base, number);

			SET_STRING_ELT(motifname2,0,mkChar(filename));

	SET_VECTOR_ELT(returnData,0,seqConsencus);
	SET_VECTOR_ELT(returnData,2,LengthSequence);
	SET_VECTOR_ELT(returnData,4,motifname2);
	SET_VECTOR_ELT(returnData,1,PWM);
	SET_VECTOR_ELT(GADEMList,0,SequencesIdent);
	SET_VECTOR_ELT(GADEMList,1,StrandIdent);
	SET_VECTOR_ELT(GADEMList,2,PositionIdent);
	SET_VECTOR_ELT(GADEMList,3,PValue);
	SET_VECTOR_ELT(GADEMList,4,AccessionIdent);
	SET_VECTOR_ELT(GADEMList,5,SeqIden);

	SET_VECTOR_ELT(returnData,3,GADEMList);

	UNPROTECT(13);
	return (returnData);

}



