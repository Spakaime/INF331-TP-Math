#ifndef __Matrice___
#define __Matrice___ 

typedef struct POl
{
	float a;
	float b;
	float c;
	float d;
	float x1;
	float x2;
	float x3;
	int nbrR;

}POl;

typedef struct Matrice
{
	struct POl p;
	float mat[3][3];
	float mat2[3][3];
	float matD[3][3];
	float valP[3];

}Matrice;

float absl(float a);
float powfPl(float a,float n);
POl solutionDgr3(float a, float b, float c, float d);
POl solutionDgr2(float a, float b, float c);
Matrice insertVal(Matrice M, int n);
Matrice valPropre(Matrice M, int n);
void afficheMatD(Matrice M, int n);
Matrice diagonaliser(Matrice M, int n);
#endif