#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gadem.h"
#include "defines.h"

//int read_pwm0(char *fileName,double **pwm0) {
int read_pwm0(SEXP RSpwm,double **pwm0,int lengthMatrix) {
	
	int numCol; //numRow;
//FILE *fp;

	numCol=(lengthMatrix/4);
	int compt=0;

	for (int aa=0;aa<(lengthMatrix/4);aa++)
	{
		for(int bb=0;bb<4;bb++)
		{
			pwm0[aa][bb]=NUMERIC_DATA(RSpwm)[compt];
			compt++;
		}
	}

	return (numCol);
}

