for(int i=0; i<_NXTRANSBLOCK; i++){
  for(int j=0; j<_NYTRANSBLOCK; j++){
    //printf("i2 : %f; j2 : %f; x : %f; y : %f || ",i2,j2,x,y);
    //printf("i : %d; j : %d; i2 : %d; j2 : %d; x : %d; y : %d\n",i,j,i2,j2,x,y);
    for(int k=0;k<_M;k++){
      //printf("%f |", wtmp[k]);
      Wi[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
    }

    //attention si i = 0 ou i = _NXTRANSBLOCK*_NYTRANSBLOCK*_M - 1
    if(i > 0 && i < _NXTRANSBLOCK*_NYTRANSBLOCK*(_M - 1)){
      for(k = 0; k < _M; k ++){
        W1[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
      }
      for(k = 0; k < _M; k ++){
        W2[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
      }
      Rusanov(Wi, W2, norx, flux1);
      Rusanov(W1, Wi, norx, flux2);
    }

      // est ce que ça vaut 0 en dehors de la grille ?
      // j ai mis oui pour l instant mais cest a verifier
      else if(i == 0){
        for(k = 0; k < _M; k ++){
          W2[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
        }
        Rusanov(Wi, W2, norx, flux1);
        Rusanov(zero, Wi, norx, flux2);
      }

      else{
        for(k = 0; k < _M; k ++){
          W1[k] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i];
        }
        Rusanov(Wi, zero, norx, flux1);
        Rusanov(W1, Wi, norx, flux2);
      }
      for(k = 0; k < _M; k ++){
        Wns[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i ] = Wns[i + j] = Wn1[k*_NXTRANSBLOCK*_NYTRANSBLOCK+ j*_NXTRANSBLOCK + i ] - (*dtt/dx)*((real)flux1[j] - flux2[j]);;
      }
    }
    //printf("\n");
    //printf("\n");

    //printf("NX : %i; NY : %i\n", _NXTRANSBLOCK, _NYTRANSBLOCK);
  }
for(i = 0; i < _NXTRANSBLOCK*_NYTRANSBLOCK*_M; i ++){
  Wn1[i] = Wns[i];
}
