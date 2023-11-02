#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrice.h"
/*------------------fonction pour optimiser la valeur absolue----------------*/
float absl(float a){
	if (a<0)
		return -1*a;
	else
		return a;
}

			/*#######################################################*/


/*---------fonction pour le calcul de la puissance d'un mombre--------------*/

float powfPl(float a,float n){
	float c;
	if (a!=0)
	{
		c = absl(a)/a;
		return c*pow(absl(a),n);		
	}
	else 
	return 0;
}

					/*#######################################################*/


/*-------------fonction permettant de resoudres les equations de degree 2 ------------*/

POl solutionDgr2(float a, float b, float c){
	float det;
	POl p;
	det = pow(b,2) - 4*a*c;
	if (det>=0)
	{
		p.x1 = (-b -sqrt(det)/(2*a));
		p.x2 = (-b +sqrt(det)/(2*a));
		p.nbrR = 2;
	}
	else
		p.nbrR = 0;
	return p;		
}
					/*#######################################################*/

/*----------------fonction permettant de resoudre un polynome de degree  3 ------------------------*/


POl solutionDgr3(float a, float b, float c, float d){

	float r0,r1;
	POl p;

	r0 = pow(b,2) - 3*a*c;
	r1 = (9*a*b*c - 2*pow(b,3) - 27*pow(a,2)*d)/(2*sqrt(pow(absl(r0),3)));

	if (r0>0)
	{	
	
		/*>>>>>>>>l'equation admet 3 racines reels<<<<<<<<*/
		if (absl(r1)<=1)
		{	
			p.x1 = (2*sqrt(r0)*cos((acos(r1)/3))-b)/(3*a);
			p.x2 = (2*sqrt(r0)*cos((acos(r1)/3 - 2*M_PI/3))-b)/(3*a);
			p.x3 = (2*sqrt(r0)*cos((acos(r1)/3 + 2*M_PI/3))-b)/(3*a);
			p.nbrR = 3;
		}
		else if(absl(r1)>1){
		/*>>>>>>>>l'equation admet 1 racines reels<<<<<<<<*/	
			float a0;
			a0 = powfPl(absl(r1) + sqrt(pow(r1,2)-1),1.0/3.0) + powfPl(absl(r1) - sqrt(pow(r1,2)-1),1.0/3.0);
			p.x1 = (sqrt(r0)*absl(r1)/(3*a*r1))*a0 - b/(3*a);
			p.nbrR = 1;

		}

	}
	else if (r0 == 0)
	{
		// l'equation admet une seul solution reel double
		p.x1 = (-b + powfPl((powfPl(b,3) - 27*pow(a,2)*p.d),1.0/3.0))/(3*a);
		p.x2 = p.x1;
		p.x3 = p.x1;
		p.nbrR = 3;

	}
	else if (r0 < 0)
	{	
		// l'equation admet une unique racine reel

		float a1;
		a1 = powfPl((r1 + sqrt(pow(r1,2)+1)),1.0/3.0) + powfPl((r1 - sqrt(pow(r1,2)+1)),1.0/3.0);
		p.x1 = (sqrt(absl(r0))/(3*a))*a1 - b/(3*a);
		p.nbrR = 1;

	}

return p;

}

						/*#######################################################*/

/*------------------------fonction permettant de remplir une matrice-------------------------------------*/

Matrice insertVal(Matrice M, int n){
	int i,j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("Entrer la valeur [%d][%d] :",i,j );
			scanf("%f",&M.mat[i][j]);
		}
	}

return M;
}	


						/*#######################################################*/

/*-------------------fonction permettant de calculer les valeurs propres de la matrice--------------------------*/

Matrice valPropre(Matrice M, int n){

	if (n == 2)
	{
			M.p.a = 1;
			M.p.b = -1*(M.mat[0][0]+M.mat[1][1]);
			M.p.c = M.mat[0][0]*M.mat[0][0] - M.mat[0][1]*M.mat[1][0];
			M.p = solutionDgr2(M.p.a,M.p.b,M.p.c);
		if (M.p.nbrR == 2)
			{
				printf("cette matrice peut etre diagonaliser dans le corp des reels\n");

			}
		else
		printf("Cette matrice ne peut pas etre diagonaliser dans le corp des reels\n");		

	}
	else{
		M.p.a = -1;
		M.p.b = M.mat[0][0] + M.mat[1][1] + M.mat[2][2];

		M.p.c = +M.mat[2][0]*M.mat[0][2] + M.mat[1][2]*M.mat[2][1]
				+ M.mat[1][0]*M.mat[0][1] - M.mat[0][0]*M.mat[1][1] 
				- M.mat[0][0]*M.mat[2][2] - M.mat[1][1]*M.mat[2][2];

		M.p.d = - M.mat[0][2]*M.mat[2][0]*M.mat[1][1] 
				+ M.mat[0][0]*M.mat[1][1]*M.mat[2][2]	 
				+ M.mat[0][2]*M.mat[1][0]*M.mat[2][1] 
				- M.mat[0][0]*M.mat[1][2]*M.mat[2][1] 
				- M.mat[2][2]*M.mat[0][1]*M.mat[1][0]  
				+ M.mat[0][1]*M.mat[2][0]*M.mat[1][2]; 

		M.p = solutionDgr3(M.p.a,M.p.b,M.p.c,M.p.d);
		
		if(M.p.nbrR == 3){
			printf("cette matrice peut etre diagonaliser dans le corp des reels\n");
		}
		else
		printf("Cette matrice ne peut pas etre diagonaliser dans le corp des reels\n");						
	}

return M;
}											

						/*#######################################################*/


/*--------------------fonction permettant d'afficher la  matrice diagonale -------------------------------------*/

void afficheMatD(Matrice M, int n){
	printf("-----------------------------------\n");
	int j,i;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%3.f    ",M.matD[i][j] );
		
		}
		printf("\n");
	}	

}			


						/*#######################################################*/

/*--------------------------fonction permettant de diagonaliser la matrice------------------------------*/

Matrice diagonaliser(Matrice M, int n){

	M.valP[0] = M.p.x1;
	M.valP[1] = M.p.x2;
	M.valP[2] = M.p.x3; 
	int i,j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
				{
					M.matD[i][j] = M.valP[i];
				}
			else
				M.matD[i][j] = 0;	
		
		}
	}	

return M;
}									
