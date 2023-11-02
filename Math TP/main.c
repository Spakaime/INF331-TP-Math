#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrice.h"


int main(int argc, char const *argv[])
{
	
	int n,chx;	
	Matrice M;


	printf("\t\t\t#########################################################\n");
	printf("\t\t\t#  DIAGONALISATION D'UNE Matrice CARREE D'ORDRE 2 OU 3  #\n");
	printf("\t\t\t#  ---------------------------------------------------  #\n");
	printf("\t\t\t#########################################################\n");

	printf("\n\n");
	printf("\t\t\t 1- Matrice carre d'ordre 2\n");
	printf("\t\t\t 2- Matrice carre d'ordre 3\n");

	printf("\n\t\t\tQuelle est la matrice que voulez diagonaliser ? :");
	scanf("%d",&chx);
	printf("--------------------------------------------------\n");
	switch(chx){
		case 1:
			n = chx +1;
			M = insertVal(M, n);
			M = valPropre(M, n);
			M = diagonaliser(M,n);
			afficheMatD(M, n);

		break;

		case 2 :
			n = chx +1;
			M = insertVal(M, n);
			M = valPropre(M, n);
			M = diagonaliser(M,n);
			afficheMatD(M, n);

		break;	
	}



scanf("%c",chx);
	return 0;
}
