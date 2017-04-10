void Rusanov(real* Wl, real* Wr, real* n, real* rus){

  // calcul du flux de Rusanov en connaissant le flux numerique

  real flux1[9];
  real flux2[9];

  flux(Wl, n, flux1);
  flux(Wr, n, flux2);

  int i;

  for(i =0; i < 9; i ++){
    rus[i] = ((flux1[i] + flux2[i])/2 - 3*(Wr[i] - Wl[i]));
  }

}

void TimesStepCPU1D(real** Wn1){
  //

  //Nx nombre de mailles du maillage dans la direction x
  // ça vaut peut etre _NXTRANSBLOCK : a verifier
  // je le met a 1 pour que ça compile
  real Nx = 1;
  real dx = (_XMAX - _XMIN)/Nx;

  real alpha = dx/6;
  real dt = alpha - _CFL;

  int i, j;

  //tab suivant
  real **Wns = (real**)malloc(_NXTRANSBLOCK*_NYTRANSBLOCK*_M);
  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
    Wns[i] = (real*)malloc(9);
  }


  //norx est le vecteur normal nx : (1, 0, 0)
  real norx[3];
  norx[0] = 1;
  norx[1] = 0;
  norx[2] = 0;

  // utile pour le calcul
  real flux1[9];
  real flux2[9];


  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){

    //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
    if(i > 0 && i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1 ){
      Rusanov(Wn1[i], Wn1[i+1], norx, flux1);
      Rusanov(Wn1[i-1], Wn1[i], norx, flux2);
    }

    // est ce que ça vaut 0 en dehors de la grille ?
    // j ai mis oui pour l instant mais cest a verifier
    else if(i == 0){
      Rusanov(Wn1[i], Wn1[i+1], norx, flux1);
      Rusanov(0, Wn1[i], norx, flux2);
    }

    else{
      Rusanov(Wn1[i], 0, norx, flux1);
      Rusanov(Wn1[i-1], Wn1[i], norx, flux2);
    }


    for(j = 0; j < 9; j ++){
      Wns[i][j] = (real)Wn1[i][j] - (real)(dt/dx)*((real)flux1[i] - (real)flux2[i]);
    }

  }

  // est ce que la fonction se rappel elle meme et affiche la grille
  // ou on l'appelle dans une boucle ?

}

