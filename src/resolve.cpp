// cas {{{
#define _1D
//#define _2D
// }}}


#include <math.h>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "resolve.h"

#define _M (9)                                      //Number of conservative variable
#define _TMAX (1)
#define _GAP (0)                                    // taille du recouvrement entre les work groupe
#define _TRANSBLOCK 1                               // size of the cached transposed block
#define _NBWORKSX (1<<7)                            // number of work-items in a work-group
#define _NBLOCKSX (1)                               // number of work-groups
#define _NX (_NBLOCKSX*(_NBWORKSX-2*_GAP)+2*_GAP)   // number of volume finite
#define _NXTRANSBLOCK ( (_NX%_TRANSBLOCK==0)? _NX : _NX+_TRANSBLOCK-_NX%_TRANSBLOCK )

#define _NBWORKSY (1<<7)
#define _NBLOCKSY (1)
#define _NY (_NBLOCKSY*(_NBWORKSY-2*_GAP)+2*_GAP)
#define _NYTRANSBLOCK ( (_NY%_TRANSBLOCK==0)? _NY : _NY+_TRANSBLOCK-_NY%_TRANSBLOCK )
#define _CFL (0.5)
#define _SPLIT (1)                                  // affiche 1 maille sur _SPLIT dans le fichier .msh

#define Min(a,b) (((a) < (b)) ? (a) : (b))
#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Abs(a) ((a) > (0) ? (a) : (-a))

#ifdef _1D
    #define _LONGUEURX (10)                         //Longueur du domaine suivant x
    #define _LONGUEURY (10)                         //Longueur du domaine suivant y
    #define _XMIN (-5)
    #define _XMAX (5)
    #define _YMIN (-5)
    #define _YMAX (5)
#endif
#ifdef _2D
    // Orzag Tang
    #define _LONGUEURX (6.2831853)                  //Longueur du domaine suivant x
    #define _LONGUEURY (6.2831853)                  //Longueur du domaine suivant y
    #define _XMIN (0)
    #define _XMAX (6.2831853)
    #define _YMIN (0)
    #define _YMAX (6.2831853)
#endif



//#define _GAM (2)
#define _GAM (1.666666666666)
#define _PI (3.14159265359)
#define _CH (5)


// c'est les vraies fonctions pour les vecteurs je pense très fort
void copy(real* A, real* B){
	int i;

	for(i = 0; i < 9; i ++){
		A[i] = B[i];
	}
}


void conservatives(real* Y, real* W){

	real gam = _GAM;

	copy(W,Y);
	int i;
	for(i = 1; i <= 4; i ++){
		if(i == 2){
			W[i] = Y[2]/(gam-1) + Y[0]*(Y[1]*Y[1]+Y[3]*Y[3]+Y[4]*Y[4])/2 + (Y[7]*Y[7]+Y[5]*Y[5]+Y[6]*Y[6])/2;
		}
		else{
			W[i] = W[0] * W[i];
		}

	}
}

void primitives(real* Y, real* W){

	real gam = _GAM;

	copy(Y,W);
	int i;
	for(i = 1; i <= 4; i ++){
		if(i == 2){
			Y[i] = (gam-1)*(W[i] - W[0]*(W[1]*W[1]+W[3]*W[3]/(Y[0]*Y[0])+W[4]*W[4]/(Y[0]*Y[0]))/2 + (Y[7]*Y[7]+Y[5]*Y[5]+Y[6]*Y[6]));
		}
		else{
			Y[i] = Y[i] / Y[0];
		}
	}
}

void flux(real* W, real* vn, real* flux){

    real gam = _GAM;
    real un,bn,E;

    real Y[_M];

    primitives(W,Y);

    real b = Y[7];

    un = Y[1]*vn[0]+Y[3]*vn[1]+Y[4]*vn[2];
    bn = Y[7]*vn[0]+Y[5]*vn[1]+Y[6]*vn[2];

    E = Y[2]/(gam-1) + Y[0]*(Y[1]*Y[1]+Y[3]*Y[3]+Y[4]*Y[4])/2 + (Y[7]*Y[7]+Y[5]*Y[5]+Y[6]*Y[6])/2;

    flux[0] = Y[0]*un;
    flux[1] = Y[0]*un*Y[1] + (Y[2] + (Y[7]*Y[7] + Y[5]*Y[5] + Y[6]*Y[6])/2)*vn[0] - bn*Y[7];
    flux[2] = (E + Y[2] + (Y[7]*Y[7] + Y[5]*Y[5] + Y[6]*Y[6])/2)*un - (Y[7]*Y[1] + Y[5]*Y[3] + Y[6]*Y[4])*bn;
    flux[3] = Y[0]*un*Y[3] + (Y[2] + (Y[7]*Y[7] + Y[5]*Y[5] + Y[6]*Y[6])/2)*vn[1] - bn*Y[5];
    flux[4] = Y[0]*un*Y[4] + (Y[2] + (Y[7]*Y[7] + Y[5]*Y[5] + Y[6]*Y[6])/2)*vn[2] - bn*Y[6];
    flux[5] = -bn*Y[3] + un*Y[5] + Y[8]*vn[1];
    flux[6] = -bn*Y[4] + un*Y[6] + Y[8]*vn[2];
    flux[7] = -bn*Y[1] + un*Y[7] + Y[8]*vn[0];
    flux[8] = _CH*_CH*bn;

}

void Rusanov(real* Wl, real* Wr, real* n, real* rus){

  // calcul du flux de Rusanov en connaissant le flux numerique

  real flux1[_M];
  real flux2[_M];

  flux(Wl, n, flux1);
  flux(Wr, n, flux2);

  int i;
  //printf("rus1\n");
  for(i =0; i < _M; i ++){
    rus[i] = ((flux1[i] + flux2[i])/2 - 3*(Wr[i] - Wl[i]));
    //printf("irus : %d\n", i);
  }
  //printf("rus2\n");

}

void TimesStepCPU1D(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M], real* dtt){
  //la fonction prend un dtt en argument mais ej sais aps a quoi il correspond
  // c'est peut être le dt que je définis après mais c'est a vérifier
  // en plus c'est un pointer

  //Nx nombre de mailles du maillage dans la direction x
  // ça vaut peut etre _NXTRANSBLOCK : a verifier
  // je le met a 1 pour que ça compile
  real Nx = _NXTRANSBLOCK;
  real dx = (_XMAX - _XMIN)/Nx;

  int i, j;

  //tab suivant
  real Wns[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];

  //tab special
  real Wi[_M];
  real W1[_M];
  real W2[_M];
  real *zero = (real*)calloc(_M, sizeof(real));

  //norx est le vecteur normal nx : (1, 0, 0)
  real norx[3];
  norx[0] = 1;
  norx[1] = 0;
  norx[2] = 0;

  // utile pour le calcul
  real flux1[_M];
  real flux2[_M];
  real coucou;


  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i += _M){
    //printf(" début %d\n",i);
    //on recupere les composant W

    for(j = 0; j < _M; j ++){
      Wi[j] = Wn1[i+j];
    }


    //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
    if(i > 0 && i < _NXTRANSBLOCK*_NYTRANSBLOCK*(_M - 1)){
      for(j = 0; j < _M; j ++){
        W1[j] = Wn1[(i - 1) + j];
      }
      for(j = 0; j < _M; j ++){
        W2[j] = Wn1[(i + 1) + j];
      }
      Rusanov(Wi, W2, norx, flux1);
      Rusanov(W1, Wi, norx, flux2);
    }

    // est ce que ça vaut 0 en dehors de la grille ?
    // j ai mis oui pour l instant mais cest a verifier
    else if(i == 0){
      for(j = 0; j < _M; j ++){
        W2[j] = Wn1[(i + 1) + j];
      }
      Rusanov(Wi, W2, norx, flux1);
      Rusanov(zero, Wi, norx, flux2);
    }

    else{
      for(j = 0; j < _M; j ++){
        W1[j] = Wn1[(i - 1) + j];
      }
      Rusanov(Wi, zero, norx, flux1);
      Rusanov(W1, Wi, norx, flux2);
    }

    //printf(" milieu %d %f\n",i);

    for(j = 0; j < _M; j ++){
      //printf("Wns[%d+%d] : %f \n", i,j, Wns[i+j]); // Wns[i+j] ne provoque pas de segfault
      coucou = Wn1[i+j] - (*dtt/dx)*(flux1[j] - flux2[j]);
      //printf("coucou : %f\n", coucou); // coucou ne provoque pas de segfault
      Wns[i+j] = coucou; // provoque un segfault
      //printf("  j  : %d; i+j : %d, max : %d \n",j,i+j, _NXTRANSBLOCK*_NYTRANSBLOCK*_M);
    }



  }

  // est ce que la fonction se rappel elle meme et affiche la grille
  // ou on l'appelle dans une boucle ?

  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wn1[i] = Wns[i];
  }
  free(zero);

}


void TimesStepCPU2D(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M], real* dtt){
  //

  //Nx nombre de mailles du maillage dans la direction x
  // ça vaut peut etre _NXTRANSBLOCK : a verifier
  // de meme pour Ny
  real Nx = _NXTRANSBLOCK;
  real dx = (_XMAX - _XMIN)/Nx;

  real Ny = _NYTRANSBLOCK;
  real dy = (_YMAX - _YMIN)/Ny;

  int i, j;
  //tab suivant
  real Wns[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];

  //tab special
  real W[_M];
  real Wi1[_M];
  real Wi2[_M];
  real Wj1[_M];
  real Wj2[_M];
  real *zero = (real*)calloc(_M, sizeof(real));

  //norx est le vecteur normal nx : (1, 0, 0)
  real norx[3];
  norx[0] = 1;
  norx[1] = 0;
  norx[2] = 0;

  // utile pour le calcul
  real flux1[_M];
  real flux2[_M];
  real flux3[_M];
  real flux4[_M];

  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i += _M){

    //on recupere les composant W
    for(j = 0; j < _M; j ++){
      W[j] = Wn1[i + j];
    }

    //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
    if(i != 0 && i != _NXTRANSBLOCK*_NYTRANSBLOCK*(_M - 1)){
      //printf("%d\n",i);
      for(j = 0; j < _M; j ++){
        Wi1[j] = Wn1[(i - 1) + j];
      }
      for(j = 0; j < _M; j ++){
        Wi2[j] = Wn1[(i + 1) + j];
      }
      Rusanov(W, Wi2, norx, flux1);
      Rusanov(Wi1, W, norx, flux2);
    }

    // est ce que ça vaut 0 en dehors de la grille ?
    // j ai mis oui pour l instant mais cest a verifier
    else if(i == 0){
      for(j = 0; j < _M; j ++){
        Wi2[j] = Wn1[(i + 1) + j];
      }
      Rusanov(W, Wi2, norx, flux1);
      Rusanov(zero, W, norx, flux2);
    }

    else{
      for(j = 0; j < _M; j ++){
        Wi1[j] = Wn1[(i - 1) + j];
      }
      Rusanov(W, zero, norx, flux1);
      Rusanov(Wi1, W, norx, flux2);
    }

    if(i < _NXTRANSBLOCK){
      for(j = 0; j < _M; j ++){
        Wj2[j] = Wn1[(i + _NXTRANSBLOCK) + j];
      }
      Rusanov(W, Wj2, norx, flux3);
      Rusanov(zero, W, norx, flux4);
    }
    else if(i > _NXTRANSBLOCK * (_NYTRANSBLOCK - 1)){
      for(j = 0; j < _M; j ++){
        Wj1[j] = Wn1[(i - _NXTRANSBLOCK) + j];
      }
      Rusanov(W, zero, norx, flux3);
      Rusanov(Wj1, W, norx, flux4);
    }
    else{
      for(j = 0; j < _M; j ++){
        Wj1[j] = Wn1[(i - _NXTRANSBLOCK) + j];
      }
      for(j = 0; j < _M; j ++){
        Wj2[j] = Wn1[(i + _NXTRANSBLOCK) + j];
      }
      Rusanov(W, Wj2, norx, flux3);
      Rusanov(Wj1, W, norx, flux4);
    }

    for(j = 0; j < _M; j ++){
      Wns[i + j] = (real)Wn1[i + j] - (real)(*dtt/dx)*((real)flux1[j] - (real)flux2[j]) - (real)(*dtt/dy)*((real)flux3[j] - (real)flux4[j]);
    }

  }

  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wn1[i] = Wns[i];
  }
  free(zero);

}
