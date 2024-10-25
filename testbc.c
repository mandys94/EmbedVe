#define NozzR 8.e-6
#define NozzL 15.e-6
#define FeedThR 50.78646481e-6
#define FeedThL 100.e-6

int LEVEL = 9;

#include "./embed.h"
#include "axi.h"
#include "./centered.h"
#include "./two-phase.h"
#include "log-conform.h"
#include "./tension.h"
#if CEMB
#include "./contact-embed.h"
#endif
#include "tag.h"
#include "./restrictor.h"

double muup=0.;
double relaxtime = 0.;

#define sigmaa 70.e-3
#define mus 0.001
#define Weber (rhoo*NozzR/sigmaa)
#define Reynolds (rhoo*NozzR/muu)
#define De (relaxtime*1.e-6/sqrt(rhoo*(NozzR*NozzR*NozzR)/sigmaa))
#define Beta (mus/muu)
scalar lambdav[], mupv[];
scalar toto[];

#define EPS 1e-10

scalar kappa[];
double nozzle (double X, double Y) {
  double Deltaa = L0/(1 << LEVEL);
  double
    x1 = round(FeedThL/NozzR/Deltaa)*Deltaa,
    x2 = round((FeedThL + 1.05*NozzL)/NozzR/Deltaa)*Deltaa,
    y1 = round(1./Deltaa)*Deltaa,
    y2 = round(FeedThR/NozzR/Deltaa)*Deltaa;
  
  return
    union(
	  difference(
		     intersection(
				  X - x1 + EPS,
				  -X + x2 + EPS
				  ),
		     -Y + y1 - EPS
		     ),
	  Y - y2 + EPS
	  );
}


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

u.n[left] = neumann(0.);
p[left] = dirichlet(Pch/rhoo);
pf[left] = dirichlet(Pch/rhoo);
f[left] = dirichlet(1.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
f[right] = neumann(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
uf.n[embed] = dirichlet(0.);
uf.t[embed] = dirichlet(0.);
p[embed] = neumann(0.);
pf[embed] = neumann(0.);
#if CEMB
f[embed] = neumann(0);
#else
f[embed] = x < (FeedThL + NozzL)/NozzR ? dirichlet(1.) : dirichlet(0.);
#endif
toto[embed] = neumann(0.);

event defaults (i=0){
  for (scalar s in (scalar *){tau_p}) {
    //    s.v.x.i = -1; // just a scalar, not the component of a vector
    //    foreach_dimension()
    s[embed] = neumann(0.);
    //      s[embed] = f[ghost] > 0 ? neumann(0) : dirichlet(0);
    //      s[embed] = x < (FeedThL + NozzL)/NozzR ? neumann(0) : dirichlet(0);
  }
#if AXI
  scalar s = tau_qq;
  s[embed] = neumann(0.);
  //  s[embed] = f[ghost] > 0 ? neumann(0) : dirichlet(0);
  //  s[embed] = x < (FeedThL + NozzL)/NozzR ? neumann(0) : dirichlet(0);
#endif  
}

int main (int argc, char * argv[]) {
  relaxtime = atof(argv[1]);
  double muratio = atof(argv[2]);
  LEVEL = atoi(argv[3]);
  muup = muu*(muratio);
  muu = (mus + muup);
  fprintf(ferr,"# Beta = %g De = %g Re = %g We = %g\n",Beta,De,Reynolds,Weber);

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
  init_grid (1 << (LEVEL /* - 4 */));

  f.sigma = 1./Weber;
  
  mu1 = 1./Reynolds*Beta;
  mu2 = 3.33e-3*mu1;
  
  rho1 = 1.;
  rho2 = 1e-3;
  
  TOLERANCE = 1e-6;
  
  // Identifying the # of lines in the waveform
  int lineCount;
  if (pid() == 0) {
    FILE * fr;
    fr = fopen ("waveform.txt","r");
    lineCount = countLines(fr);
    printf("Number of lines: %d\n", lineCount);
  }
  MPI_Bcast(&lineCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Reading the waveform
  double tdata[2*lineCount], Vdot[2*lineCount];
  if (pid() == 0) {
    FILE * fr;
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
    for (int i = 0; i < lineCount; i++) {
      tdata[2*i] = !i ? 0. : tdata[2*i - 1] + 1e-6;
      tdata[2*i + 1] = !i ? wdur[0] : tdata[2*i - 1] + wdur[i];
      Vdot[2*i] = wtypeI[i] == 0 ? 0. : (wvolt[i+1] - wvolt[i-1])/wdur[i];
      Vdot[2*i + 1] = wtypeI[i] == 0 ? 0. : (wvolt[i+1] - wvolt[i-1])/wdur[i];
    }
    fr = fopen ("Vdot.txt","a");
    for (int i = 0; i < 2*lineCount; i++)
      fprintf (fr,"%.6f\t%.6f\n",tdata[i],Vdot[i]);
    fclose(fr);
  }
  MPI_Bcast(&tdata[0], 2*lineCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Vdot[0], 2*lineCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  ptdata = tdata; pVdot = Vdot;

#if CEMB
  const scalar c[] = 90*pi/180.;
  contact_angle = c;
#endif

  mup = mupv;
  lambda = lambdav; 

  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) { 
    double rad = 2.*sigmaa/(fabs(Pch)*NozzR);
    
/* #if TREE */
/*     refine (fabs(nozzle(x,y)) < 0.3 && level < LEVEL); */
/*     refine (fabs(x - (FeedThL + NozzL)/NozzR) < 0.3 && y < 1. && level < LEVEL); */
/* #endif */
    solid (cs, fs, -nozzle(x,y));

    restriction ({cs, fs});
#if AXI
    cm_update (cm, cs, fs);
    fm_update (fm, cs, fs);
    restriction ({cm, fm, cs, fs});
#endif
    
    fraction (f, difference((FeedThL + NozzL)/NozzR - x, sq(rad) - sq(y) - sq(x - ((FeedThL + NozzL)/NozzR + sqrt(sq(rad) - 1)))));

    for (scalar s in {u})
      s.third = true;
    
    foreach() {
      p[] = f[]*Pch/rhoo;
      foreach_dimension()
	u.x[] = 0;
    }
    dump("init");
  }
}

#if TREE
event adapt (i++) {
  scalar cs_temp[];
  foreach()
    cs_temp[] = cs[];
  
  adapt_wavelet ({f, u, cs_temp}, (double[]){1e-3,5e-3,5e-3,1e-3}, LEVEL, LEVEL - 5);

  solid (cs, fs, -nozzle(x,y));

  restriction ({cs, fs});

#if AXI
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
  restriction ({cm, fm, cs, fs});
#endif
  
}
#endif

event properties (i++) {
  foreach() {
    mupv[] = (1. - Beta)*clamp(f[],0,1)/Reynolds;
    lambdav[] = De*clamp(f[],0,1)*sqrt(rho1/f.sigma);
    if(cs[] > 0.)
      toto[] = pid();
  }
  boundary({toto});
}

event clean (i++) {
  remove_droplets (f);
}

event snapshot (i+=10/* t=0.0025; t+=0.0025 ; t <= 2.5 */) {
  fprintf (stderr,"%g %g %g\n", t, intval(NozzR*1.e6*t,ptdata,pVdot), Pch);

  scalar omega[];
  vorticity(u, omega);

  curvature (f, kappa, 1., add = false);
  
  char name[80];
  sprintf(name, "snapshot-%d", i);
  dump(name, (scalar *){f,u,p,cs,omega,kappa,tau_p,tau_qq,toto});
  //  dump(name, (scalar *){f,u,p,cs,omega,kappa,tau_p,tau_qq,trA});
  
  sprintf(name, "interfaces%d", pid());
  FILE * fp = fopen(name,"w");
  output_facets(f,fp);
  fclose(fp);
  char command[80];
  sprintf(command, "LC_ALL=C cat interfa* > facets-%d",i);
  system(command);
  
  sprintf(name, "fields-%d",i);
  FILE * fpp = fopen(name,"w");
  if (t == 0)
    output_field({u,p,omega,cs,tau_p,tau_qq,toto}, fpp, n = 2048, linear = true, box = {{X0,Y0},{X0 + L0,Y0 + FeedThR/NozzR}});
  else
    output_field({u,p,omega,tau_p,tau_qq,toto}, fpp, n = 2048, linear = true, box = {{X0,Y0},{X0 + L0,Y0 + FeedThR/NozzR}});
  fclose(fpp);
}

event end(i=20){

}
