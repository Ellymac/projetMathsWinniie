void Rusanov(real* Wl, real* Wr, real* n, real* rus){

  // calcul du flux de Rusanov en connaissant le flux numerique

  real flux1[_M];
  real flux2[_M];

  flux(Wl, n, flux1);
  flux(Wr, n, flux2);

  int i;

  for(i =0; i < _M; i ++){
    rus[i] = ((flux1[i] + flux2[i])/2 - 3*(Wr[i] - Wl[i]));
  }

}

void TimesStepCPU1D(real Wn1[_NXTRANSBLOCK*_NYTRANSBLOCK*_M], real* dtt){
  //la fonction prend un dtt en argument mais ej sais aps a quoi il correspond
  // c'est peut être le dt que je définis après mais c'est a vérifier
  // en plus c'est un pointer

  //Nx nombre de mailles du maillage dans la direction x
  // ça vaut peut etre _NXTRANSBLOCK : a verifier
  // je le met a 1 pour que ça compile
  real Nx = 1;
  real dx = (_XMAX - _XMIN)/Nx;

  int i, j;

  //tab suivant
  real *Wns = (real*)malloc(_NXTRANSBLOCK*_NYTRANSBLOCK*_M);

  //tab special
  real *Wi = (real*)malloc(_M);
  real *W1 = (real*)malloc(_M);
  real *W2 = (real*)malloc(_M);

  //norx est le vecteur normal nx : (1, 0, 0)
  real norx[3];
  norx[0] = 1;
  norx[1] = 0;
  norx[2] = 0;

  // utile pour le calcul
  real flux1[_M];
  real flux2[_M];


  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i + _M){

    //on recupere les composant W
    for(j = 0; j < _M; j ++){
      Wi[j] = Wn1[_M * i + j];
    }

    //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
    if(i > 0 && i < _NXTRANSBLOCK*_NYTRANSBLOCK*(_M - 1)){
      for(j = 0; j < _M; j ++){
        W1[j] = Wn1[_M * (i - 1) + j];
      }
      for(j = 0; j < _M; j ++){
        W2[j] = Wn1[_M * (i + 1) + j];
      }
      Rusanov(Wi, W2, norx, flux1);
      Rusanov(W1, Wi, norx, flux2);
    }

    // est ce que ça vaut 0 en dehors de la grille ?
    // j ai mis oui pour l instant mais cest a verifier
    else if(i == 0){
      for(j = 0; j < _M; j ++){
        W2[j] = Wn1[_M * (i + 1) + j];
      }
      Rusanov(Wi, W2, norx, flux1);
      Rusanov(0, Wi, norx, flux2);
    }

    else{
      for(j = 0; j < _M; j ++){
        W1[j] = Wn1[_M * (i - 1) + j];
      }
      Rusanov(Wi, 0, norx, flux1);
      Rusanov(W1, Wi, norx, flux2);
    }


    for(j = 0; j < _M; j ++){
      Wns[i * _M + j] = (real)Wn1[i * _M + j] - (real)(dt/dx)*((real)flux1[j] - (real)flux2[j]);
    }

  }

  // est ce que la fonction se rappel elle meme et affiche la grille
  // ou on l'appelle dans une boucle ?
  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; I ++){
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
  real dy = (_YMAX - _YMIN)/Ny

  int i, j;

  //tab suivant
  real *Wns = (real*)malloc(_NXTRANSBLOCK*_NYTRANSBLOCK*_M);

  //tab special
  real *W = (real*)malloc(_M);
  real *Wi1 = (real*)malloc(_M);
  real *Wi2 = (real*)malloc(_M);
  real *Wj1 = (real*)malloc(_M);
  real *Wj2 = (real*)malloc(_M);

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


  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i + _M){

    //on recupere les composant W
    for(j = 0; j < _M; j ++){
      W[j] = Wn1[_M * i + j];
    }

    //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
    if(i > 0 && i < _NXTRANSBLOCK*_NYTRANSBLOCK*(_M - 1)){
      for(j = 0; j < _M; j ++){
        W1[j] = Wn1[_M * (i - 1) + j];
      }
      for(j = 0; j < _M; j ++){
        W2[j] = Wn1[_M * (i + 1) + j];
      }
      Rusanov(W, W2, norx, flux1);
      Rusanov(W1, W, norx, flux2);
    }

    // est ce que ça vaut 0 en dehors de la grille ?
    // j ai mis oui pour l instant mais cest a verifier
    else if(i == 0){
      for(j = 0; j < _M; j ++){
        W2[j] = Wn1[_M * (i + 1) + j];
      }
      Rusanov(Wi, W2, norx, flux1);
      Rusanov(0, Wi, norx, flux2);
    }

    else{
      for(j = 0; j < _M; j ++){
        W1[j] = Wn1[_M * (i - 1) + j];
      }
      Rusanov(Wi, 0, norx, flux1);
      Rusanov(W1, Wi, norx, flux2);
    }

    if(i < _NXTRANSBLOCK){
      Rusanov(W, W4, norx, flux3);
      Rusanov(0, W, norx, flux4);
    }
    else if(i > _NXTRANSBLOCK * (_NYTRANSBLOCK - 1)){
      Rusanov(W, 0, norx, flux3);
      Rusanov(W3, W, norx, flux4);
    }
    else{
      Rusanov(W, W4, norx, flux3);
      Rusanov(W3, W, norx, flux4);
    }

    for(j = 0; j < _M; j ++){
      Wns[i * _M + j] = (real)Wn1[i * _M + j] - (real)(dtt/dx)*((real)flux1[j] - (real)flux2[j]) - (real)(dtt/dy)((real)flux3[j] - (real)flux4[j]);
    }

  }

  for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; I ++){
    Wn1[i] = Wns[i];
  }

}
