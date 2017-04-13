#ifndef RESOLVE_H_
#define RESOLVE_H_

typedef float real;

void Wexact(real* x, real* y, real* W);
void copy(real* A, real* B);
void conservatives(real* Y, real* W);
void primitives(real* Y, real* W);
void flux(real* W, real* vn, real* flux);
void Rusanov(real* Wl, real* Wr, real* n, real* rus);
void TimesStepCPU1D(real *Wn1, real* dtt);
void TimesStepCPU2D(real *Wn1, real* dtt);

#endif /*RESOLVE_H_*/
