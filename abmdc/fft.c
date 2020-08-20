#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "headers/md.h"

/*********************************************************/
/* FFT3D routine                                         */
/* taken from: "Numerical Recipes in C",Japanese edition */
/*             GIJUTSU HYOURONSYA, (1992), p.440         */
/*********************************************************/
void fft3d(int isign)
{
#define DATA(i) wv.fftdat[(i)-1]
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

  int ndim, idim;
  unsigned long nn[4];
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit, k1, k2, n, nprev, nrem, ntot=1;
  double tempi, tempr;
  double theta, wi, wpi, wpr, wr, wtemp;
  int i;

  ndim = 3; /* 3D FFT */
  ntot = FFTPTSX*FFTPTSY*FFTPTSZ;
  nprev=1;
  nn[1] = FFTPTSX; nn[2] = FFTPTSY; nn[3] = FFTPTSZ;

  for(idim=ndim;idim>=1;idim--) { /* loop for the dimension */

    n = nn[idim];
    nrem=ntot/(n*nprev);
    ip1 = nprev << 1; /* 1 bit shift */
    ip2 = ip1*n;
    ip3 = ip2*nrem;
    i2rev = 1;
    for(i2=1;i2<=ip2;i2+=ip1) {   /* bit reverse */
      if(i2 < i2rev){
	for(i1=i2;i1<=i2+ip1-2;i1+=2){
	  for(i3=i1;i3<=ip3;i3+=ip2){
	    i3rev=i2rev+i3-i2;
	    SWAP(DATA(i3),DATA(i3rev));
	    SWAP(DATA(i3+1),DATA(i3rev+1));
	  }
	}
      }
      ibit=ip2 >> 1;
      while(ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }

    ifp1 = ip1;         /* Danielson-Lanczos part */
    while(ifp1<ip2){

      ifp2  = ifp1 << 1;
      theta = isign*PI2/(ifp2/ip1); /* initialization */
      wtemp = sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);

      wr = 1.0;
      wi = 0.0;

      for(i3=1;i3<=ifp1;i3+=ip1) {
	for(i1=i3;i1<=i3+ip1-2;i1+=2){
	  for(i2=i1;i2<=ip3;i2+=ifp2){
	    k1 = i2;
	    k2 = k1+ifp1;

	    tempr = wr*DATA(k2)   - wi*DATA(k2+1);
	    tempi = wr*DATA(k2+1) + wi*DATA(k2);
	    DATA(k2)   = DATA(k1) - tempr;
	    DATA(k2+1) = DATA(k1+1) - tempi;
	    DATA(k1)   += tempr;
	    DATA(k1+1) += tempi;
	  }
	}	
	wr = (wtemp=wr)*wpr - wi*wpi + wr;
	wi = wi*wpr + wtemp*wpi + wi;
      }
      ifp1 = ifp2;
    }
    nprev *= n;
  }

  if(isign > 0) { /* when in FFT [CHECK HERE] */
    for(i=0;i<2*ft.mesh;i++) {
      wv.fftdat[i]/=(double)ft.mesh; 

      /* CHECK HERE */
      /*      if((i+1)%2==0) wv.fftdat[i]*=-1.0; */

    }
  }

#undef SWAP
#undef DATA
}

void check_fft(void)
{
  int i;
  int cnt = 0;

  for(i=0;i<ft.mesh;i++) {
    printf(" %15.7e  %15.7e \n",wv.fftdat[cnt],wv.fftdat[cnt+1]);cnt+=2;
  }
}
