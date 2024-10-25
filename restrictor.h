#include "tridiagonal.h"

double dr;

int Nr = 1e5;

double rhoo = 997.075116517, muu = 0.000890438981615;
double Pch = -800., LenR = 8.6e-4;
double Cch = 1.08e-20, alphaV = 0.135e-15;

double *rL, *rR, *rS;

double *A, *B, *C, *D, *U;

double Qr = 0., Qf = 0.;

double *ptdata, *pVdot;

double intval (double val, double *tdat, double *Vdat) {
  double t0, t1, interp = 0.;
  for (int i = 0; i < 22; i++) {
    t0 = tdat[i] - tdat[0];
    t1 = tdat[i+1] - tdat[0];
    if (val == t0)
      interp = Vdat[i];
    else if (val > t0 && val < t1)
      interp = Vdat[i] + (val - t0)*(Vdat[i+1] - Vdat[i])/(t1 - t0);
  }
  return interp;
}

event restrictor (i++) {
  double dtt = dt*NozzR;
  
  B[0] = rhoo*rS[0]/2. + muu*rR[0]*dtt/dr;
  C[0] = -muu*rR[0]*dtt/dr;
  A[Nr-1] = -muu*rL[Nr-1]*dtt/dr;
  B[Nr-1] = rhoo*rS[Nr-1]/2. + muu*(rL[Nr-1] + 2.*rR[Nr-1])*dtt/dr;
  D[0] = rhoo*rS[0]*U[0]/2. + dtt*rS[0]*Pch/(2.*LenR);
  D[Nr-1] = rhoo*rS[Nr-1]*U[Nr-1]/2. + dtt*rS[Nr-1]*Pch/(2.*LenR);

  for (int i = 1; i < Nr-1; i++) {
    A[i] = -muu*rL[i]*dtt/dr;
    B[i] = rhoo*rS[i]/2. + muu*(rL[i] + rR[i])*dtt/dr;
    C[i] = -muu*rR[i]*dtt/dr;
    D[i] = rhoo*rS[i]*U[i]/2. + dtt*rS[i]*Pch/(2.*LenR);
  }

  tridiag(Nr,A,B,C,D,U);

  Qr = 0.; Qf = 0.;

  // Restrictor flow rate
  for (int i = 0; i < Nr; i++) {
    Qr += U[i]*(rL[i] + rR[i])/2;
  }
  Qr *= 2*pi*dr;

  // Feedthrough flow rate
  foreach(reduction(+:Qf)) {
    if (x == Delta/2.)
#if EMBED
      //      if(cs[] >= 1.)
      Qf += embed_face_value_x(point,u.x,0)*y*Delta/2.;
    //      Qf += fs.x[]*(u.x[] + u.x[-1])*y*Delta/2.;
#else
      Qf += (u.x[] + u.x[-1])*y*Delta/2.;
#endif
  }
  Qf *= 2*pi;

  // Pressure PDE
  Pch += dtt*(-alphaV*intval(NozzR*1.e6*t,ptdata,pVdot)*1e6 - Qf*pow(NozzR,2.) - Qr)/Cch;
}
