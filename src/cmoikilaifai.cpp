void conservatives(real* Y, real* W){
	W = Y;
	int i;
	for(i = 1; i <= 3; i ++){
		W[i] = W[0] * W[i];
	}
}

void primitives(real* Y, real* W){
	Y = W;
	int i;
	for(i = 1; i <= 3; i ++){
		Y[i] = Y[i] / Y[0];
	}
}
