#define NozzR 8.e-6
#define NozzL 15.e-6
#define FeedThR 50.78646481e-6
#define FeedThL 100.e-6

int LEVEL = 9;

#include "./embed.h"
#include "axi.h"
#include "./centered.h"
#include "./two-phase.h"
#include "log-conformV1.h"
#include "./tension.h"
#include "tag.h"
#include "./restrictor.h"

char filename[150];
char name[80];

int main(int argc, char * argv[]){
  sprintf(filename, "%s", argv[1]);
  fprintf(ferr,"%s\n",filename);
  restore (file = filename);

  event ("metric");

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
  if(pid() ==0)
    for (int j = 0; j < n; j++)
      fprintf (fp, "%g %d %g %g %g %g %g %g\n", t,
	       j, fmax, v[j], b[j].x/v[j], b[j].y/v[j],ub[j].x/v[j], ub[j].y/v[j]);
    fclose(fp);

}
