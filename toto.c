#include "embed.h"
#include "navier-stokes/centered.h"
//#include "utils.h"
#include "log-conform.h"
//vector u[];

#undef center_gradient_x
#define center_gradient_x(a,i) (fs.x[0,i] && fs.x[1,i] ? (a[1,i] - a[-1,i])/(2.*Delta) : \
				fs.x[1,i] ? (a[1,i] - a[0,i])/Delta :	\
				fs.x[0,i]  ? (a[0,i] - a[-1,i])/Delta : 0)

#undef center_gradient_y
#define center_gradient_y(a,i) (fs.y[i] && fs.y[i,1] ? (a[i,1] - a[i,-1])/(2.*Delta) : \
			      fs.y[i,1] ? (a[i,1] - a[i])/Delta :        \
			      fs.y[i]  ? (a[i] - a[i,-1])/Delta : 0)

int main(){
  
  init_grid(32);
  
  solid (cs, fs, sq(x)+sq(y)-sq(1./2.));

  /* output_facets(cs,ferr); */
    foreach_face(){
      double shear1 = (tau_p.x.y[0,1] + tau_p.x.y[-1,1] -
		      tau_p.x.y[0,-1] - tau_p.x.y[-1,-1])/4.;
    }
    
    foreach_face(){
      double shear2 = (face_gradient_y(tau_p.x.y,1) + face_gradient_y(tau_p.x.y,-1))/2.;

    }
  /* dump("toto"); */
  /* foreach_face() */
  /*   fprintf(stdout,"%g %g %g %g\n",x,y,fs.x[],fs.y[]); */
  /* output_cells(ferr); */
}
