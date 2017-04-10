void conservatives(real* Y, real* W){
	
	real gam = _GAM;
	
	W = Y;
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
	
	Y = W;
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
