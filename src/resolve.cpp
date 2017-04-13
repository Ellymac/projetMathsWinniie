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


void Wexact(real* x, real* y, real* W){

#ifdef _1D
    real YL[_M], YR[_M], WR[_M], WL[_M];

//    // Test de Choc fort
//    YL[0] = 3;
//    YL[1] = 1.3;
//    YL[3] = 0;
//    YL[4] = 0;
//    YL[2] = 3;
//    YL[5] = 1;
//    YL[6] = 1;
//    YL[7] = 1.5;
//    YL[8] = 0;
//
//    YR[0] = 1;
//    YR[1] = 1.3;
//    YR[3] = 0;
//    YR[4] = 0;
//    YR[2] = 1;
//    YR[5] = 0.0707372016677029;
//    YR[6] = 0.9974949866040544;
//    YR[7] = 1.5;
//    YR[8] = 0;

//    // Test de Brio et Wu
//    YL[0] = 1;
//    YL[1] = 0;
//    YL[3] = 0;
//    YL[4] = 0;
//    YL[2] = 1;
//    YL[5] = 1;
//    YL[6] = 0;
//    YL[7] = 0.75;
//    YL[8] = 0;
//
//    YR[0] = 0.125;
//    YR[1] = 0;
//    YR[3] = 0;
//    YR[4] = 0;
//    YR[2] = 0.1;
//    YR[5] = -1;
//    YR[6] = 0;
//    YR[7] = 0.75;
//    YR[8] = 0;

    //Test de Dai et Woodward
    YL[0] = 1.08;
    YL[1] = 1.2;
    YL[3] = 0.01;
    YL[4] = 0.5;
    YL[2] = 0.95;
    YL[5] = 1.0155412503859613165;
    YL[6] = 0.56418958354775628695;
    YL[7] = 1.1283791670955125739;
    YL[8] = 0;

    YR[0] = 1;
    YR[1] = 0;
    YR[3] = 0;
    YR[4] = 0;
    YR[2] = 1;
    YR[5] = 1.1283791670955125739;
    YR[6] = 0.56418958354775628695;
    YR[7] = 1.1283791670955125739;
    YR[8] = 0;


    conservatives(YL, WL);
    conservatives(YR, WR);


    if(*x < 0)
        for(int i=0; i<_M; i++){
            W[i] = WL[i];
        }
    else
        for(int i=0; i<_M; i++){
            W[i] = WR[i];
        }
#endif
#ifdef _2D
// Orzag-Tang
    real Y[_M];

    real gam = _GAM;

    Y[0] = gam*gam;
    Y[1] = -sin(*y);
    Y[2] = gam;
    Y[3] = sin(*x);
    Y[4] = 0.0;
    Y[5] = sin(2*(*x));
    Y[6] = 0.0;
    Y[7] = -sin(*y);
    Y[8] = 0.0;

    conservatives(Y, W);
    //printf("%f ",W[0]);
#endif

}

// c'est les vraies fonctions pour les vecteurs je pense très fort
void copy(real* A, real* B){
  int i;

  for(i = 0; i < 9; i ++){
    A[i] = B[i];
  }
}


void conservatives(real* Y, real* W){

  int i;
  W[0] = Y[0];
  for(i = 1; i < 5; i++){
    if(i == 2){
      W[2] = (_GAM-1)*(Y[2] + Y[0]*(Y[1]*Y[1] + Y[4]*Y[4] + Y[3]*Y[3])/2 + (Y[5]*Y[5] + Y[6]*Y[6] + Y[7]*Y[7])/2);
    }
    else{
      W[i] = Y[i]*Y[0];
    }
  }
  for(i = 5; i < _M; i++){
    W[i] = Y[i];
  }
}

void primitives(real* Y, real* W){

  int i;
  Y[0] = W[0];
  for(i = 1; i < 5; i++){
    if(i == 2){
      Y[2] = W[2]/(_GAM-1) + Y[0]*(Y[1]*Y[1] + Y[3]*Y[3] + Y[4]*Y[4])/2 + (W[5]*W[5] + W[6]*W[6] + W[7]*W[7])/2;
    }
    else{
      Y[i] = W[i]/W[0];
    }
  }
  for(i = 5; i < _M; i++){
    Y[i] = W[i];
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
/*
  printf("\n\n\n\n");
  for(int i = 0; i < _M; i ++){
    printf("Wl[%d] = %f \n", i, Wl[i]);
    printf("Wr[%d] = %f \n", i, Wr[i]);
  }
*/
  flux(Wl, n, flux1);
  flux(Wr, n, flux2);

  int i;
  //printf("rus1\n");
  for(i =0; i < _M; i ++){
    rus[i] = (0.5*(flux1[i] + flux2[i]) - _CH*(Wr[i] - Wl[i]));
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


  //tab suivant
  real Wns[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];

  real Wcopy[_NXTRANSBLOCK*_NYTRANSBLOCK*_M];
  for(int  i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wcopy[i] = Wn1[i];
  }

  //tab special
  real Wi[_M];
  real W1[_M];
  real W2[_M];
  real zero[_M];
  for(int i = 0; i < _M; i++){
    zero[i] = (real)0;
  }

  real min = _XMIN;
  real max = _XMAX;

  //norx est le vecteur normal nx : (1, 0, 0)
  real norx[3];
  norx[0] = 1;
  norx[1] = 0;
  norx[2] = 0;

  // utile pour le calcul
  real flux1[_M];
  real flux2[_M];


  for(int i=0; i<_NXTRANSBLOCK; i++){
    for(int j=0; j < _NYTRANSBLOCK; j++){
      for(int k=0;k<_M;k++){
        Wi[k] = Wcopy[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
      }
      //int b;
      //scanf("%d",&b);
      //printf("rho : %f\n", Wn1[0]);
      //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
      if(i > 0 && i < _NXTRANSBLOCK ){
        for(int k = 0; k < _M; k ++){
          W1[k] = Wcopy[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i - 1];
          W2[k] = Wcopy[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i + 1];
        }
        Rusanov(Wi, W2, norx, flux1);
        Rusanov(W1, Wi, norx, flux2);
      }

      // est ce que ça vaut 0 en dehors de la grille ?
      // j ai mis oui pour l instant mais cest a verifier
      else if(i == 0){
        for(int k = 0; k < _M; k ++){
          W2[k] = Wcopy[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i + 1];
        }
        Wexact(&min, zero, W1);
        Rusanov(Wi, W2, norx, flux1);
        Rusanov(zero, Wi, norx, flux2);
      }

      else{
        for(int k = 0; k < _M; k ++){
          W1[k] = Wcopy[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i - 1];
        }
        Wexact(&max, zero, W2);
        Rusanov(Wi, zero, norx, flux1);
        Rusanov(W1, Wi, norx, flux2);
      }
      for(int k = 0; k < _M; k ++){
        Wns[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i ] = Wns[i + j] = Wcopy[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i ] - (*dtt/dx)*((real)flux1[k] - flux2[k]);
        //printf("i : %d, j : %d, k : %d \n", i, j, k);
       // printf("valeur : %f, flux1 : %f, flux2 : %f, Wns : %f \n", Wns[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i ], flux1[k], flux2[k], Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i ]);
      }
    }
    //printf("\n");
    //printf("\n");

    //printf("NX : %i; NY : %i\n", _NXTRANSBLOCK, _NYTRANSBLOCK);
  }
  for(int i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wn1[i] = Wns[i];
  }

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

  for(int i=0; i<_NXTRANSBLOCK; i++){
    for(int j=0; j < _NYTRANSBLOCK; j++){
      for(int k=0;k<_M;k++){
        W[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
      }
      //printf("rho : %f\n", Wn1[0]);
      //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
      if(i > 0 && i < _NXTRANSBLOCK ){
        for(int k = 0; k < _M; k ++){
          Wi1[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i - 1];
          Wi2[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i + 1];
        }
        Rusanov(W, Wi2, norx, flux1);
        Rusanov(Wi1, W, norx, flux2);
      }

      // est ce que ça vaut 0 en dehors de la grille ?
      // j ai mis oui pour l instant mais cest a verifier
      else if(i == 0){
        for(int k = 0; k < _M; k ++){
          Wi2[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i + 1];
        }
        Rusanov(W, Wi2, norx, flux1);
        Rusanov(zero, W, norx, flux2);
      }

      else{
        for(int k = 0; k < _M; k ++){
          Wi1[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i - 1];
        }
        Rusanov(W, zero, norx, flux1);
        Rusanov(Wi1, W, norx, flux2);
      }

      if(j > 0 && j < _NYTRANSBLOCK ){
        for(int k = 0; k < _M; k ++){
          Wj1[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ (j - 1)*_NXTRANSBLOCK + i];
          Wj2[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ (j + 1)*_NXTRANSBLOCK + i];
        }
        Rusanov(W, Wj2, norx, flux1);
        Rusanov(Wj1, W, norx, flux2);
      }

      // est ce que ça vaut 0 en dehors de la grille ?
      // j ai mis oui pour l instant mais cest a verifier
      else if(i == 0){
        for(int k = 0; k < _M; k ++){
          Wj2[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ (j + 1)*_NXTRANSBLOCK + i];
        }
        Rusanov(W, Wj2, norx, flux1);
        Rusanov(zero, W, norx, flux2);
      }

      else{
        for(int k = 0; k < _M; k ++){
          Wj1[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ (j - 1)*_NXTRANSBLOCK + i];
        }
        Rusanov(W, zero, norx, flux1);
        Rusanov(Wj1, W, norx, flux2);
      }

      for(int k = 0; k < _M; k ++){
        Wns[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i] - (*dtt/dx)*(flux1[k] - flux2[k]) - (*dtt/dy)*(flux3[k] - flux4[k]);
      }
    }
    //printf("\n");
    //printf("\n");

    //printf("NX : %i; NY : %i\n", _NXTRANSBLOCK, _NYTRANSBLOCK);
  }
  for(int i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wn1[i] = Wns[i];
  }



  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wn1[i] = Wns[i];
  }

}
