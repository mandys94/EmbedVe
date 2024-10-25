#define JACOBI 1

#define NozzR 8.e-6
#define NozzL 15.e-6
#define FeedThR 50.786e-6
#define FeedThL 100.e-6

#define sigmaa 71.709675e-3
//#define sigmaa 68.3e-3

#include "axi.h"
#include "navier-stokes/centered.h"
/* #if CONTACT */
/* #include "contact.h" */
/* #endif */
#include "two-phase.h"
//#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "restrictor.h"

/* #if CONTACT */
/* vector h[]; */
/* double theta0 = 90; */
/* h.t[bottom] = contact_angle (theta0*pi/180.); */
/* #endif */

#define Weber (rhoo*NozzR/sigmaa)
#define Reynolds (rhoo*NozzR/muu)

scalar kappa[];
int countLines (FILE * file) {
  int count = 0;
  char ch;
  if (file == NULL) {
    printf ("Error in opening file.\n");
    return 1;
  }
  else {
    printf ("Lines successfully counted.\n");
    while ((ch = fgetc(file)) != EOF)
      if (ch == '\n')
	count++;
  }
  fclose(file);
  return count;
}

#define LEVEL 12

u.n[left] = neumann(0.);
p[left] = dirichlet(Pch/rhoo);
pf[left] = dirichlet(Pch/rhoo);
f[left] = dirichlet(1.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
f[right] = neumann(0.);

u.t[top] = dirichlet(0.);
u.n[top] = dirichlet(0.);
p[top] = neumann(0.);
pf[top] = neumann(0.);
#if CONTACT
f[top] = neumann(0);
#else
f[top] = x < (FeedThL+NozzL)/NozzR ? 1. : 0.;
#endif
int main () {
  // Restrictor mesh
  dr = 4.84448e-5/Nr;
  rL = calloc(Nr, sizeof(double));
  rR = calloc(Nr, sizeof(double));
  rS = calloc(Nr, sizeof(double));
  for (int i = 0; i < Nr; i++) {
    *(rL + i) = i*dr;
    *(rR + i) = (i + 1)*dr;
    *(rS + i) = pow((i + 1)*dr,2.) - pow(i*dr,2.);
  }

  // Tridiagonal matrix algorithm data
  A = calloc(Nr, sizeof(double));
  B = calloc(Nr, sizeof(double));
  C = calloc(Nr, sizeof(double));
  D = calloc(Nr, sizeof(double));
  U = calloc(Nr, sizeof(double));

  L0 = 32.;
  init_grid (1 << (LEVEL - 4));

  f.sigma = 1./Weber;

  mu1 = 1./Reynolds;
  mu2 = 0.01*mu1;

  rho1 = 1.;
  rho2 = 1e-3;
  
  TOLERANCE = 1e-6;

  // Identifying the # of lines in the waveform
  FILE * fr;
  fr = fopen ("waveform.txt","r");
  int lineCount = countLines(fr);
  printf("Number of lines: %d\n", lineCount);

  // Reading the waveform
  char wtype[6];
  int wtypeI[lineCount];
  float wdur[lineCount], wvolt[lineCount];
  
  fr = fopen ("waveform.txt","r");
  if (fr == NULL) {
    printf ("Error in opening file.\n");
    return 1;
  }
  else {
    printf ("Waveform successfully read.\n");
    for (int i = 0 ; i < lineCount; i++) {
      fscanf(fr,"%s\t%f\t%f", wtype, &wdur[i], &wvolt[i]);
      wtypeI[i] = strcmp(wtype,"Const") == 0 ? 0 :
	strcmp(wtype,"Flank") == 0 ? 1 :
	none;
    }
  }
  fclose(fr);

  // Calculating \dot{V}(t)
  double tdata[2*lineCount], Vdot[2*lineCount];
  ptdata = tdata; pVdot = Vdot;
  for (int i = 0; i < lineCount; i++) {
   tdata[2*i] = !i ? 0. : tdata[2*i - 1] + 1e-6;
    tdata[2*i + 1] = !i ? wdur[0] : tdata[2*i - 1] + wdur[i];
    Vdot[2*i] = wtypeI[i] == 0 ? 0. : (wvolt[i+1] - wvolt[i-1])/wdur[i];
    Vdot[2*i + 1] = wtypeI[i] == 0 ? 0. : (wvolt[i+1] - wvolt[i-1])/wdur[i];
  }
  fr = fopen ("Vdot.txt","w");
  for (int i = 0; i < 2*lineCount; i++)
    fprintf (fr,"%.6f\t%.6f\n",tdata[i],Vdot[i]);
  fclose(fr);

/* #if CONTACT */
/*   f.height = h; */
/* #endif */
  
  run();
}

event init (t = 0) {
  double rad = 2.*sigmaa/(fabs(Pch)*NozzR); // need to be checked
  mask ((x > (FeedThL/NozzR) && (x < ((FeedThL + NozzL)/NozzR + 1.*Delta)) && y >= 1. ) || y > FeedThR/NozzR ? top : none);
  refine (fabs(x - (FeedThL + NozzL)/NozzR) < 0.2 && y < FeedThR/NozzR && level <= (LEVEL));

  fraction (f, difference((FeedThL + NozzL)/NozzR - x, sq(rad) - sq(y) - sq(x - ((FeedThL + NozzL)/NozzR + sqrt(sq(rad) - 1)))));
  
  foreach() {
    p[] = f[]*Pch/rhoo;
    foreach_dimension()
      u.x[] = 0;
  }
}

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-3,5.e-3,5.e-3}, LEVEL, 8);
}

event clean (i++) {
  remove_droplets (f);
}

/* event logs(i += 1){ */
/*   face vector ss = {{-1}}; */
/*   double fmax = 0.; */
/*   foreach(reduction(max:fmax)) */
/*     if (y < Delta) */
/*       if (f[] > 1e-6 && f[] < 1. - 1e-6) { */
/* 	coord n = facet_normal (point, f, ss); */
/* 	double alpha = plane_alpha (f[], n); */
/* 	coord segment[2];       */
/* 	if(facets (n, alpha, segment) == 2) */
/* 	  fmax = max(x + segment[1].x*Delta,x + segment[0].x*Delta); */
/*       } */
/*   fprintf(ferr,"%10.9f %10.9f %10.9f %10.9f %10.9f %10.9f\n",t,fmax,Pch,Qr,Qf,-alphaV*intval(NozzR*1.e+6*t,ptdata,pVdot)*1e6); */
/*   //  fprintf (stdout,"%g %g %g\n", t, intval(15.3*t,ptdata,pVdot), Pch);   */
/* } */
event logs(t += 0.01){
  face vector ss = {{-1}};
  double fmax = 0.;
  foreach(reduction(max:fmax))
    if (y < Delta)
      if (f[] > 1e-6 && f[] < 1. - 1e-6) {
	coord n = facet_normal (point, f, ss);
	double alpha = plane_alpha (f[], n);
	coord segment[2];      
	if(facets (n, alpha, segment) == 2)
	  fmax = max(x + segment[1].x*Delta,x + segment[0].x*Delta);
      }
  
  static FILE * fp;
  fp = fopen ("drops.dat", "a");

  scalar m[];
  foreach()
    m[] = (f[]) > 1e-3;
  int n = tag (m);
  double v[n];
  coord b[n];
  coord ub[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = ub[j].x = ub[j].y =0.;
  foreach (serial)
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension(){
	b[j].x += dv()*f[]*p.x;
	ub[j].x += dv()*f[]*u.x[];
      }
    }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ub, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  if(pid() ==0){
    for (int j = 0; j < n; j++)
      fprintf (fp, "%d %g %d %g %g %g %g %g\n", i, t,
	       j, v[j], b[j].x/v[j], b[j].y/v[j],ub[j].x/v[j], ub[j].y/v[j]);
    fclose(fp);
  }
  fprintf(ferr,"%10.9f %10.9f %10.9f %10.9f %10.9f %10.9f\n",t,fmax,Pch,Qr,Qf,-alphaV*intval(NozzR*1.e+6*t,ptdata,pVdot)*1e6);
  //  fprintf (stdout,"%g %g %g\n", t, intval(15.3*t,ptdata,pVdot), Pch);  
}

event snapshot (t += 0.01; t <= 2.5) {

  curvature (f, kappa, 1., add = false);
  
  char name[80];
  sprintf (name,"facets/facets-%g",t);
  FILE * fp = fopen (name,"w");
  output_facets (f,fp);
  fclose (fp);
  
  sprintf (name,"fields/fields-%g",t);
  FILE * fpp = fopen (name,"w");
  output_field ({u,p,kappa}, fpp, n = 2048, linear = true, box = {{X0,Y0},{X0 + L0,Y0 + FeedThR/NozzR}});
  fclose (fpp);
}
