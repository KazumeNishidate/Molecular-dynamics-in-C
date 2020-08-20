#include       "sk.h"

void hanning_window(double *x, int n)
{
  int i;
  double factor;

  factor = 2.0*PI/(double)n;

  for(i=0; i<n; i++) {
    *(x+i)=*(x+i)*(1.0-cos(factor*(double)i))/2.0; 
  }
}

void fft(double* x, double* y, int step, double flag)
{
  /* flag=-1.0 -> FFT  */
  /* flag= 1.0 -> RFT  */

  int i, i0, i1, i2, i3, j, ij, ir, n, n2;
  double s, c, sc, *xx, *yy, tx, ty, arg;

  n = step;  sc = 2.0 * PI / (double)n;

  xx = (double *)calloc(n, sizeof(double));
  yy = (double *)calloc(n, sizeof(double));

  n2 = n;

  for(j=0;j<n;j++){     /* main loop  */
    n2 /= 2;
    if(n2 != 0){
      for(i=0;i<n/2/n2;i++){
	ij = i * n2;    arg = flag * sc * (double)ij;
	c = cos(arg);   s = sin(arg);
	for(ir=0;ir<n2;ir++){
	  i0 = ir + ij;    i1 = i0 + ij;    i2 = i1 + n2;    i3 = i0 + n/2;
	  tx = c * x[i2] - s * y[i2];
	  ty = c * y[i2] + s * x[i2];
	  xx[i0] = x[i1] + tx;    yy[i0] = y[i1] + ty;
	  xx[i3] = x[i1] - tx;    yy[i3] = y[i1] - ty;
	}
      }
      for(i=0;i<n;i++){
	x[i] = xx[i];    y[i] = yy[i];
      }
    }
  }

  if(flag < 0.0)        /*  test fft or rft  */
    for(i=0;i<n;i++){
      x[i] /= (double)n;    y[i] /= (double)n;
    }
  free(xx);   free(yy);
}



