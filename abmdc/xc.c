#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

void  set_HVxc_hmat(int npw, int *nx, int *ny, int *nz)
{
  int nq1, nq2, cnt=0;
  int ndx, ndy, ndz, id;

  for(nq1=0;nq1<npw;nq1++) {
    for(nq2=0;nq2<=nq1;nq2++) {
      ndx = nx[nq1] - nx[nq2];
      ndy = ny[nq1] - ny[nq2];
      ndz = nz[nq1] - nz[nq2];
      
      if(ndx < 0) ndx += FFTPTSX;
      if(ndy < 0) ndy += FFTPTSY;
      if(ndz < 0) ndz += FFTPTSZ;
      
      id = FFTPTSX*FFTPTSY*ndz + FFTPTSX*ndy + ndx;
      /*-- fill up the Hamiltonian with Hartree and Vxc to diagonalize --*/
      hm.mat_re[cnt] += hxc.r[id];
      hm.mat_im[cnt] += hxc.i[id];
      cnt++;
      /*-----------------------------------------------------------------*/
    }
  }

}

void  set_Hartree(int npw)
{
  int i, cnt, id1, id2;
  double factor;

  /*  printf(">>> set Hartree energy \n"); */

  cnt=0;
  for(i=0;i<ft.mesh;i++) {
    wv.fftdat[cnt] = ft.rho[i];cnt++;
    wv.fftdat[cnt] = 0.0;cnt++;
  }

  fft3d(1); /* (1) : FFT [real -> G space] rho(r) -> rho(G) */

  wv.fftdat[0] = 0.0; /* cut the diverting part (G=0) */
  wv.fftdat[1] = 0.0; 

  for(i=1;i<npw;i++) { /* except the Gamma (G=0) */

    id2 = wv.index[i];
    id1 = id2/2;

    factor = 4.0*PI/wv.norm2[i];  /* get the VH(G) = (4Pi/|G|^2)*rho(G) */

    wv.fftdat[id2]   *= factor;
    wv.fftdat[id2+1] *= factor;

    hxc.r[id1] += wv.fftdat[id2];   /* add the Hartree in G space       */
    hxc.i[id1] += wv.fftdat[id2+1]; /* (exchange correlation) + Hartree */
  }

  fft3d(-1); /* (-1) : reverse FFT [ VH(G) -> VH(real space) ] */  

  for(i=0;i<ft.mesh;i++) { /* hxc.VH : for the Hartree calc. */
    hxc.VH[i] = wv.fftdat[i*2];
  }

}

void  set_Vxc(int npw)
{
  int i, cnt;
  int id1, id2;

  /*  printf(">>> get the exchange-correlation energy \n"); */

  for(i=0;i<ft.mesh;i++) {
    hxc.r[i] = 0.0;
    hxc.i[i] = 0.0;
  }

  calc_Vxc(); /* real-space */

  cnt=0;
  for(i=0;i<ft.mesh;i++) {
    wv.fftdat[cnt] = hxc.V[i];cnt++;
    wv.fftdat[cnt] = 0.0;cnt++;
  }

  fft3d(1); /* (1) : FFT [real space -> G space] */

  for(i=0;i<npw;i++) {  /* (hxc.r, hxc.i) = exchange correlation */
    id1 = wv.index[i];
    id2 = id1/2;
    hxc.r[id2] = wv.fftdat[id1];
    hxc.i[id2] = wv.fftdat[id1+1];
  }

}

/*----------------------------------------------------------------*/
/* exchange-correlation energy calculation                        */
/* - D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980) */
/* - J. Perdew and A. Zunger, Phys. Rev. B23, 5048 (1981)         */
/*----------------------------------------------------------------*/
void calc_Vxc()
{
  static const double  aa = 1.0,     ab = 1.0529,  acc = 0.3334;
  static const double  ad = 0.1423,  ae = 0.0480,   af = 0.0311;
  static const double  ag = 0.0116,  ah = 0.0020,   ai = 0.4582;
  double rs, ex, dex, ec, dec, t;

  int i;

  for(i=0;i<ft.mesh;i++) {

    /* x^y = exp(y*log(x))  [y = integer, x><0] or [y=float, x>0]  */
    rs = exp( (1.0/3.0)*log(3.0/(4.0*PI*ft.rho[i])) );
    ex = -ai/rs;
    dex = ai/(rs*rs);
    if(rs >= 1.0) {
      t   =  aa+ab*sqrt(rs)+acc*rs;
      ec  = -ad/t;
      dec =  ad*( ab/(2.0*sqrt(rs))+acc )/(t*t);
    } else {
      ec  = -ae+af*log(rs)-ag*rs+ah*rs*log(rs);
      dec = af/rs-ag+ah*(1.0+log(rs));
    }
    hxc.E[i] = ex+ec;
    hxc.V[i] = ex+ec-(dex+dec)*rs/3.0;
  }

}

void  set_VhxcC(int npw, int *nx, int *ny, int *nz)
{
  int ib, nq1, nq2;
  double ndx, ndy, ndz;
  int id;
  
  for(nq1=0;nq1<npw;nq1++) {
    
    for(ib=0;ib<atom.nvband_all;ib++) { /* band by band */
      wv.hvxr[ib][nq1] = 0.0;
      wv.hvxi[ib][nq1] = 0.0;
    }
    
    for(nq2=0;nq2<=npw;nq2++) {
      ndx = nx[nq1] - nx[nq2];
      ndy = ny[nq1] - ny[nq2];
      ndz = nz[nq1] - nz[nq2];
      
      if(ndx < 0) ndx += FFTPTSX;
      if(ndy < 0) ndy += FFTPTSY;
      if(ndz < 0) ndz += FFTPTSZ;
      
      id = FFTPTSX*FFTPTSY*ndz + FFTPTSX*ndy + ndx;
      
      for(ib=0;ib<atom.nvband_all;ib++) { /* band by band */
	wv.hvxr[ib][nq1] += hxc.r[id]*wv.cg2r[ib][nq2] 
	                  - hxc.i[id]*wv.cg2i[ib][nq2];
	wv.hvxi[ib][nq1] += hxc.r[id]*wv.cg2i[ib][nq2] 
	                  + hxc.i[id]*wv.cg2r[ib][nq2];
      }
    }
    
  }
}
