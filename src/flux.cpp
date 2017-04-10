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
