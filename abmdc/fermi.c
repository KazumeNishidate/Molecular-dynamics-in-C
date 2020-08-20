#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/***************************************************/
/* set occupation number by the fermi distribution */
/***************************************************/
void  set_fermi(void)
{
  int i;

  for(i=0;i<atom.nvband_all;i++) { 
    atom.fermi[i] = 1.0;  /* only for the valence band calculation yet */
  }
}

void  get_chgdns(void)
{
  int nb, ng, index;
  int i, im, pts;
  double factor;

  /*  printf(">>> get the charge density \n"); */

  for(i=0;i<ft.mesh;i++) { ft.rho[i]=0.0; }

  for(nb=0;nb<atom.nvband_all;nb++){ /* band by band reverse FFT3D */
    for(i=0;i<2*ft.mesh;i++) { wv.fftdat[i]=0.0; }
    for(ng=0;ng<wv.nplw;ng++) {
      index = wv.index[ng];
      wv.fftdat[index]   = wv.cg2r[nb][ng];
      wv.fftdat[index+1] = wv.cg2i[nb][ng];    /* + or - [CHECK HERE] */
    }

    fft3d(-1); /* (-1) : reverse FFT [CHECK HERE] */

    factor = (2.0/cell.vol)*atom.fermi[nb];

    for(im=0;im<ft.mesh;im++) {
      pts = 2*im;
      ft.rho[im] += factor*(wv.fftdat[pts]*wv.fftdat[pts] +
			    wv.fftdat[pts+1]*wv.fftdat[pts+1]);
    }

  }
}
