/* Implementation of Thomas' Algorithm as solver for tridiagonal */
/* linear system of equations. */

void tridiag(int n, double *a, double *b, double *c, double *d, double *T) {
  double w;
  for (int i = 1; i <= n-1; i++) {
    w = *(a+i)/(*(b+i-1));
    *(b+i) -= w*(*(c+i-1));
    *(d+i) -= w*(*(d+i-1));
  }
  
  *(T+n-1) = *(d+n-1)/(*(b+n-1));
  for (int i = n-2; i >= 0; i--) {
    *(T+i) = (*(d+i) - (*(c+i))*(*(T+i+1)))/(*(b+i));
  }
  
  return;
}
