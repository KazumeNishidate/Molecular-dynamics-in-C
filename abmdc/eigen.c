#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/*------------------------------------------------------------------*/
/* eigen value and eigen vector solver functions                    */
/*                                                                  */
/* REF: - CAMM FORUM, CAMP-Atami / Car Parrinello method program    */
/*        JCPE, (1994) "diag.ff"                                    */
/*      - K. Murata, T. Oguni, Y. Karaki, "SUPERCOMPUTER",          */
/*        MARUZEN, (1985) [Japanese]                                */
/*------------------------------------------------------------------*/
void  get_eigen()
{

  printf(">>> plane wave set initialization\n");

  /* reduction to tridiagonal form by extended Housholder's method */
  househ();

  /* get eigen values by Bisection method */
  bisect();

  /* get eigen vectors by Inverse iteration method for symmetric band matrix */
  invitr();

  /* inverse transformation of the eigenvectors */
  invitr_vh();

  check_eigen_values_small();

  /*  check routines -------------*/
  /*  check_eigen_values_small(); */
  /*  check_ckg_small();          */
  /*  check_eigen_values_small(); */
  /*  check_hamiltonian_small();  */
  /*------------------------------*/

}

void  check_ckg_small()
{
#define ER(i,j) wv.cg2r[(j)-1][(i)-1] /* wfn : j-band i-wavevector */
#define EI(i,j) wv.cg2i[(j)-1][(i)-1] 

  int np, nb, i, j;
  np = wv.nplwin;        /* number of plane waves */
  nb = atom.nvband_all;  /* number of eigenvalues (band) */  

  for(j=1;j<=nb;j++) {    
    for(i=1;i<=np;i++) {
      printf("[band:%2d  npl:%2d] %10.7f %10.7f\n",j,i,ER(i,j),EI(i,j));
    }
  }

#undef ER
#undef EI
}

void  check_hamiltonian_small()
{
#define A(i) hm.mat_re[(i)-1]
#define B(i) hm.mat_im[(i)-1]
  int i, np;
  np = wv.nplwin;

  printf("wk.nplwin = %d\n",np);

  for(i=1;i<=(np)*(np+1)/2;i++) {

    printf("%9.6f \n",A(i));
  }
  printf("\n");

#undef A
#undef B
}

void  invitr_vh()
{
#define IWK(i) eg.iwk[(i)-1]
#define ER(i,j) wv.cg2r[(j)-1][(i)-1] /* wfn : j-band i-wavevector */
#define EI(i,j) wv.cg2i[(j)-1][(i)-1] 
#define VR(i) eg.d2[(i)-1]
#define VI(i) eg.d3[(i)-1]
#define A(i) hm.mat_re[(i)-1]
#define B(i) hm.mat_im[(i)-1]

  int i, j, np, nb, k;
  double tr, ti, sr, si;

  np = wv.nplwin;        /* number of plane waves */
  nb = atom.nvband_all;  /* number of eigenvalues (band) */  

  for(i=1;i<=np;i++) { /*=== [10] ===*/
    IWK(i)=i*(i-1)/2;
  }
  /* vr,vi are the diagonal unitary matrix */
  
  for(j=1;j<=nb;j++) { /*=== [60] ===*/
    for(i=1;i<=np;i++) { /*=== [20] ===*/
      EI(i,j) = -ER(i,j)*VI(i);
      ER(i,j) =  ER(i,j)*VR(i);
    }
    for(k=3;k<=np;k++) { /*=== [50] ===*/
      if( A(IWK(k)+k) <= 0.0) continue; /* [goto 50] CHECK HERE */
      tr = 0.0;
      ti = 0.0;
      for(i=1;i<=k-1;i++) { /*=== [30] ===*/
	tr+= A(IWK(k)+i)*ER(i,j)-B(IWK(k)+i)*EI(i,j);
	ti+= A(IWK(k)+i)*EI(i,j)+B(IWK(k)+i)*ER(i,j);
      }
      sr = tr/A(IWK(k)+k);
      si = ti/A(IWK(k)+k);
      for(i=1;i<=k-1;i++) { /*=== [40] ===*/
	ER(i,j) += -sr*A(IWK(k)+i)-si*B(IWK(k)+i);
	EI(i,j) += -si*A(IWK(k)+i)+sr*B(IWK(k)+i);
      } /*=== end [40] ===*/
    } /*=== end [50] ===*/
  } /*=== end [60] ===*/

#undef IWK
#undef ER
#undef EI
#undef VR
#undef VI
#undef A
#undef B
  
}

/*----------------------------------------------------*/
/* Inverse iteration method:                          */
/* originally written by K. Murata (HITACHI Ltd.)     */
/*                                                    */
/* REF:                                               */
/* K. Murata, T. Oguni, Y. Karaki, "SUPERCOMPUTER",   */
/* MARUZEN, (1985) [Japanese], p275                   */
/*----------------------------------------------------*/
void  invitr()
{
#define D(i) eg.d[(i)-1]
#define D1(i) eg.d1[(i)-1]
#define A(i, j) a[(i)-1][(j)-1]
#define E(i) hm.egv[(i)-1]
#define EPS 1.0e-15
#define DMAX(x,y) ((x) > (y) ? (x) : (y))
#define DMIN(x,y) ((x) < (y) ? (x) : (y))
#define V(i,j) wv.cg2r[(j)-1][(i)-1] /* wavefunction : j-band i-wavevector */
#define IWK(i) eg.iwk[(i)-1]
#define WK(i,j) eg.wk2[(i)-1][(j)-1]
#define IFLG(i) eg.iflg[(i)-1]

  double a[2][MAX_WAVE_VECTOR1];
  int i, n, nld, m, k, js, mmk, j, mpk;
  /* CHECK HERE: non-used variable n2 and m1 */
  /*  int n2, m1; */
  int ie, ngrp, ipnt, jj, nv, iter, nscnt, nsrow;
  int mpi, mmi, je, kmax, mk, kk, mult, nband;
  double anorm, eps, ueps, s, fn, eps1, eps2, eps3, eps4, eps5;
  double ejj, t, amax, piv, sum, vnorm;

  n = wv.nplwin;          /* number of plane waves */
  nv = atom.nvband_all; /* number of total valence bands of the system */

  nld = 1;  
  /*   n2 = 2;    m1 = 4; */
  anorm = 0.0;  m = nld+1;  eps = EPS;  ueps = eps;

  A(2,1) = D(1);
  for(i=2;i<=n;i++) {
    A(1,i)=D1(i-1);
    A(2,i)=D(i);
  }

  for(k=1;k<=n;k++) { /*=== [50] ===*/
    js=DMAX(1,k-nld);
    mmk=m-k;
    s=0.0;
    for(j=js;j<=k;j++) { /*=== [35] ===*/
      s += fabs(A(mmk+j,k));
    }
    mpk=m+k;
    ie=DMIN(k+nld,n);
    if((k+1) > ie) goto gt45;
    for(i=k+1;i<=ie;i++) { /*=== [40] ===*/
      s += fabs(A(mpk-i,i));
    }
  gt45:
    if(s > anorm) anorm=s;
  }
  fn=(double)n;
  if(anorm==0.0) anorm = 1.0;
  eps1=DMAX(anorm*ueps,eps);
  eps2=anorm*1.0e-3;
  eps3=fn*eps1;
  eps4=(10.0)*eps1;
  eps5=fn*eps4;
  ngrp=0;
  ipnt=1;

  for(jj=1;jj<=nv;jj++) { /*=== [420] ===*/
    iter=1;
    nscnt=0;
    nsrow=0;
    ejj=E(jj);
    if(jj==1) goto gt80;
    t=fabs(ejj-E(jj-1));
    if(t > eps1) goto gt60;
    ngrp+=1;
    ejj+=((double)ngrp)*eps1;
    goto gt80;
  gt60:
    if(t<eps2) goto gt70;
    ngrp=0;
    goto gt80;
  gt70:
    ngrp+=1;
  gt80:  /* initialize the vector */
    for(i=1;i<=n;i++) { /*=== [90] ===*/
      V(i,jj)=eps4;
    }
    /* triangular decomposition with interchanges */
    for(i=1;i<=n;i++) { /*=== [120] ===*/
      js=DMAX(1,i-nld);
      mpi = m+i;
      mmi = m-i;
      for(j=js;j<=i;j++) { /*=== [100] ===*/
	WK(mmi+j,i)=A(mmi+j,i);
	WK(mpi-j,j)=A(mmi+j,i);
      }
      WK(m,i)=A(m,i)-ejj;
      js=i+m;
      je=DMIN((i+2*nld),n);
      if(js > je) continue; /* [goto 120] CHECK HERE */
      for(j=js;j<=je;j++) {
	WK(mmi+j,i)=0.0;
      }
    } /*=== end [120] ===*/
    for(k=1;k<=n;k++) { /*=== [200] ===*/
      amax=0.0;
      kmax=k;
      ie=DMIN((k+nld),n);
      mpk=m+k;
      for(i=k;i<=ie;i++) { /*=== [130] ===*/
	if(fabs(WK(mpk-i,i)) <= amax ) continue; /* [goto 130] CHECK HERE */
	amax = fabs(WK(mpk-i,i));
	kmax=i;
      } /*=== [130] ===*/
      IWK(k)=kmax;
      amax=WK(mpk-kmax,kmax);
      if(amax == 0.0) amax = eps1;
      WK(mpk-kmax,kmax)=1.0/amax;
      if( fabs(amax) > eps3 ) goto gt140;
      if( nscnt == ngrp ) nsrow=k;
      nscnt+=1;
    gt140:
      if(k == n) continue; /* [goto 200] CHECK HERE */
      je=DMIN((k+2*nld),n);
      mk=m-kmax;
      mmk=m-k;
      if(kmax == k) goto gt160;
      for(j=k;j<=je;j++) { /*=== [150] ===*/
	t=WK(mk+j,kmax);
	WK(mk+j,kmax)=WK(mmk+j,k);
	WK(mmk+j,k)=t;
      } /*=== end [150] ===*/
    gt160:
      piv=1.0/amax;
      for(i=k;i<=ie;i++) { /*=== [170] ===*/
	WK(mpk-i,i)=piv*WK(mpk-i,i);
      }
      WK(m,k)=piv;
      if(nld == 0) continue; /* [goto 200] CHECK HERE */
      for(i=k+1;i<=ie;i++) { /*=== [190] ===*/
	t=-WK(mpk-i,i);
	mmi=m-i;
	for(j=k+1;j<=je;j++) { /*=== [180] ===*/
	  WK(mmi+j,i)+=t*WK(mmk+j,k);
	}
      } /*=== end [190] ===*/
    } /*=== end [200] ===*/
    if(nsrow != 0) V(nsrow,jj)+=eps4;
  gt210:
    for(kk=1;kk<=n;kk++) { /*=== [240] ===*/
      k=n+1-kk;
      sum=-V(k,jj);
      je=DMIN((k+2*nld),n);
      if( (k+1) > je ) goto gt230;
      mmk=m-k;
      for(j=k+1;j<=je;j++) { /*=== [220] ===*/
	sum+=WK(mmk+j,k)*V(j,jj);
      }
    gt230:
      V(k,jj)=-sum*WK(m,k);
    } /*=== end [240] ===*/
    if(ngrp==0) goto gt280;
    
    /* orthogonalize with respect to the previous member of the group */
    js=jj-ngrp;
    je=jj-1;
    for(j=js;j<=je;j++) { /*=== [270] ===*/
      sum = 0.0;
      for(i=1;i<=n;i++) { /*=== [250] ===*/
	sum += V(i,jj)*V(i,j);
      }
      if(sum == 0.0) continue; /* [goto 270] CHECK HERE */
      s = -sum;
      for(i=1;i<=n;i++) { /*=== [260] ===*/
	V(i,jj)+=s*V(i,j);
      }
    } /*=== end [270] ===*/

    /* test for convergence */
  gt280:
    vnorm=0.0;
    for(i=1;i<=n;i++) { /*=== [290] ===*/
      vnorm+=fabs(V(i,jj));
    }
    if(vnorm >= 1.0) goto gt380;
    if(iter < 10) goto gt300;
    /* no convergence */
    IFLG(jj)=0;
    goto gt390;
    /* change the initial vector */
  gt300:
    if(vnorm != 0.0) goto gt310;
    V(ipnt,jj)=eps3;
    ipnt+=1;
    if(ipnt > n) ipnt=1;
    goto gt330;
    /* normalize the initial vector */
  gt310:
    s=eps5/vnorm;
    if(iter >= 5 ) s*=100.0;
    for(i=1;i<=n;i++) { /*=== [320] ===*/
      V(i,jj)*=s;
    }
    /* forward elimination */
  gt330:
    if(nld == 0) goto gt370;
    for(k=1;k<=n;k++) { /*=== [360] ===*/
      kmax=IWK(k);
      sum=-V(kmax,jj);
      V(kmax,jj)=V(k,jj);
      js=DMAX(1,(k-nld));
      if(js > (k-1)) goto gt350;
      mmk=m-k;
      for(j=js;j<=k-1;j++) { /*=== [340] ===*/
	sum += WK(mmk+j,k)*V(j,jj);
      }
    gt350: V(k,jj)=-sum;
    } /*=== end [360] ===*/
  gt370:
    iter +=1;
    goto gt210;
    /* normalize as norm=1 */
  gt380:
    IFLG(jj)=1;
  gt390:
    sum=0.0;
    for(i=1;i<=n;i++) { /*=== [400] ===*/
      sum += V(i,jj)*V(i,jj);
    }
    if(sum==0) continue; /* [goto 420] CHECK HERE */
    s=1.0/sqrt(sum);
    for(i=1;i<=n;i++) { /*=== [410] ===*/
      V(i,jj) *= s;
    }
  } /*=== end [420] ===*/

  mult=1; /*--------------------------------------------------*/
  nband = atom.nvband_all; /* number of eigenvalues */  
  for(i=1;i<=nband;i++) {
    mult *=IFLG(i);
  }
  if(mult != 1) {
    printf("failed in invitr (mult = %d) : exiting\n",mult);
    exit(0);
  } /*---------------------------00---------------------------*/

#undef D
#undef D1
#undef A
#undef E
#undef EPS
#undef DMAX
#undef DMIN
#undef V
#undef IWK
#undef WK

}


/*----------------------------------------------------*/
/* Householder's reduction:                           */
/* originally written by A. Nagahori (HITACHI Ltd.)   */
/*                                                    */
/* REF:                                               */
/* K. Murata, T. Oguni, Y. Karaki, "SUPERCOMPUTER",   */
/* MARUZEN, (1985) [Japanese], p268                   */
/*                                                    */
/* modified by A. Fukumoto for CAMP-Atami : eghhlh()  */
/*----------------------------------------------------*/
void  househ()
{
#define A(i) hm.mat_re[(i)-1]
#define B(i) hm.mat_im[(i)-1]
#define D(i) eg.d[(i)-1]
#define D1(i) eg.d1[(i)-1]
#define IWK(i) eg.iwk[(i)-1]
#define D2(i) eg.d2[(i)-1]
#define D3(i) eg.d3[(i)-1]
#define IMOD(x,y) (int)(x)-((int)(x)/(int)(y))*(int)(y)

  int i, k, n, kp, j, kmod;
  int is, ids, ip, ip1;
  double sc, t, hk=0.0, s2, alm, tau;
  double as, d1r, d1i, s1, r, sum, dd;

  n = wv.nplwin; /* number of plane waves */

  if(n==2) goto gt155;
  for(i=1;i<=n;i++) { /*=== [5] ===*/
    IWK(i)=i*(i-1)/2;
  }
  /*--- set U-vector ---*/
  for(k=n;k>=3;k--) {  /*=== str [150] ===============================*/
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
    kp = IWK(k);
    D(k)=A(kp+k);
    sc = 0.0;
    for(j=1;j<=k-1;j++) { /*== [10] ==*/
      sc += fabs(A(kp+j))+fabs(B(kp+j));
    } 
    if(sc == 0.0) {
      D1(k-1) = 0.0; D3(k-1) = 0.0;
      goto gt140;
    }
    t = 1.0/sc;
    s2=0.0;

    for(j=1;j<=k-1;j++) { /*== [20] ==*/
      A(kp+j) *= t;
      B(kp+j) *= t;
      s2 += A(kp+j)*A(kp+j)+B(kp+j)*B(kp+j);
    }

    alm = A(kp+k-1)*A(kp+k-1)+B(kp+k-1)*B(kp+k-1);

    if(alm>0.0) {
      tau=sqrt(alm*s2);
      hk=s2+tau;
      as=1.0/alm;
      d1r=A(kp+k-1)*as*tau;
      d1i=B(kp+k-1)*as*tau;
      A(kp+k-1) *= (1.0+s2/tau);
      B(kp+k-1) *= (1.0+s2/tau);
    } else {
      hk=s2;
      s1=sqrt(s2);
      A(kp+k-1)=s1;
      B(kp+k-1)=0.0;
      d1r=s1;
      d1i=0.0;
    }
    r = 1.0/hk;
  /*--- set p-vector ---*/
    kmod=IMOD(k-1,2);

    if(kmod==0) {
      is=1;
      ids=2;
    } else {
      is=2;
      ids=1;
      D1(1)=A(kp+1)*A(1);
      D3(1)=B(kp+1)*A(1);
    }
    for(i=is;i<=k-2;i+=2) { /*=== [30] ===*/
      D1(i) = A(kp+i)*A(IWK(i+1)) + A(kp+i+1)*A(IWK(i+1)+i)
	- B(kp+i+1)*B(IWK(i+1)+i);
      D3(i) = B(kp+i)*A(IWK(i+1)) + A(kp+i+1)*B(IWK(i+1)+i)
	+ B(kp+i+1)*A(IWK(i+1)+i);
      D1(i+1) = A(kp+i)*A(IWK(i+1)+i) + A(kp+i+1)*A(IWK(i+1)+i+1)
	+ B(kp+i)*B(IWK(i+1)+i);
      D3(i+1) = -A(kp+i)*B(IWK(i+1)+i) + B(kp+i+1)*A(IWK(i+1)+i+1)
	+ B(kp+i)*A(IWK(i+1)+i);
    }

    for(i=1;i<=k-1;i++) { /*=== [40] ===*/
      D(i)=0.0;
      D2(i)=0.0;
    }

    for(i=is;i<=k-2;i+=2) { /*=== [60] ===*/
      ip=IWK(i);
      ip1=IWK(i+1);
      
      for(j=1;j<=i-1;j++) { /*=== [50] ===*/
	D(i)    +=  A(kp+j)*A(ip+j)  + B(kp+j)*B(ip+j);
	D2(i)   += -A(kp+j)*B(ip+j)  + B(kp+j)*A(ip+j);
	D(i+1)  +=  A(kp+j)*A(ip1+j) + B(kp+j)*B(ip1+j);
	D2(i+1) += -A(kp+j)*B(ip1+j) + B(kp+j)*A(ip1+j);

	D1(j)   += A(kp+i)*A(ip+j)   + A(kp+i+1)*A(ip1+j)
	  -B(kp+i)*B(ip+j) - B(kp+i+1)*B(ip1+j);

	D3(j)   += A(kp+i)*B(ip+j)   + A(kp+i+1)*B(ip1+j)
	  +B(kp+i)*A(ip+j) + B(kp+i+1)*A(ip1+j);
      }
    }

    for(i=is;i<=k-1;i++) { /*=== [70] ===*/
      D1(i) += D(i);
      D3(i) += D2(i);
    }

    /*--- set p-vector ---*/
    sum=0.0;
    for(i=1;i<=k-1;i++) { /*=== [80] ===*/
      D1(i) *= r;
      D3(i) *= r;
      sum += A(kp+i)*D1(i)+B(kp+i)*D3(i);
    }
    t = 0.5*sum*r;

    for(i=1;i<=k-1;i++) { /*=== [90] ===*/
      D1(i) += -t*A(kp+i);
      D3(i) += -t*B(kp+i);
    }
    /*---------------- set the mateix for the next reduction ----------*/
    for(i=is;i<=k-2;i+=2) { /*=== [110] ===*/
      ip = IWK(i);
      ip1 = IWK(i+1);
      for(j=1;j<=i;j++) { /*=== [100] ===*/
	A(ip+j) += -A(kp+i)*D1(j) - D1(i)*A(kp+j) 
                  - B(kp+i)*D3(j) - D3(i)*B(kp+j);
	B(ip+j) += -A(kp+i)*D3(j) + D3(i)*A(kp+j) 
                  + B(kp+i)*D1(j) - D1(i)*B(kp+j);
	A(ip1+j) += -A(kp+i+1)*D1(j) - D1(i+1)*A(kp+j) 
                   - B(kp+i+1)*D3(j) - D3(i+1)*B(kp+j);
	B(ip1+j) += -A(kp+i+1)*D3(j) + D3(i+1)*A(kp+j) 
                   + B(kp+i+1)*D1(j) - D1(i+1)*B(kp+j);
      }
    }

    for(i=ids;i<=k-1;i+=2) { /*=== [120] ===*/
      A(IWK(i)+i) += -2.0*( A(kp+i)*D1(i) + B(kp+i)*D3(i) );
    }

    /*------------------------------------------------------------------*/
    for(j=1;j<=k-1;j++) { /*=== [130] ===*/
      A(kp+j) *= sc;
      B(kp+j) *= sc;
    }
    D1(k-1) = -d1r*sc;
    D3(k-1) = -d1i*sc;

  gt140: 
    A(kp+k)=sc*sc*hk;

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  }  /*=== end [150] ===================================================*/

  gt155: 
  D(1) = A(1);
  D1(1) = A(2);
  D3(1) = B(2);
  D(2) = A(3);

  /*---------------------------------------*/  
  /* D1 : real subdiagonal element         */
  /* D2, D3 : diagonal unitary matrix      */
  /*---------------------------------------*/
  D2(n) = 1.0;
  D3(n) = 0.0;

  for(i=n-1;i>=1;i--) { /*=== [200] ===*/
    dd = D1(i)*D1(i)+D3(i)*D3(i);
    if(dd <= 0.0) {
      D2(i) = 1.0;
      D3(i) = 0.0;
    } else {
      r = sqrt(dd);
      D2(i) = ( D1(i)*D2(i+1)-D3(i)*D3(i+1) )/r;
      D3(i) = ( D1(i)*D3(i+1)+D3(i)*D2(i+1) )/r;
      D1(i) = r;
    }
  } /*=== end [200] ===*/

#undef A
#undef B
#undef D
#undef D1
#undef IWK
#undef D2
#undef D3
#undef IMOD
}

/*----------------------------------------------------*/
/* Bisection method:                                  */
/* originally written by K. Murata (HITACHI Ltd.)     */
/*                                                    */
/* REF:                                               */
/* K. Murata, T. Oguni, Y. Karaki, "SUPERCOMPUTER",   */
/* MARUZEN, (1985) [Japanese], p274                   */
/*----------------------------------------------------*/
void  bisect()
{
#define A(i) hm.mat_re[(i)-1]
#define B(i) hm.mat_im[(i)-1]
#define D(i) eg.d[(i)-1]
#define D1(i) eg.d1[(i)-1]
#define IWK(i) eg.iwk[(i)-1]
#define D2(i) eg.d2[(i)-1]
#define D3(i) eg.d3[(i)-1]
#define EPS 1.0e-15
#define DMAX(x,y) ((x) > (y) ? (x) : (y))
#define DMIN(x,y) ((x) < (y) ? (x) : (y))
#define W(i) eg.w[(i)-1]
#define E(i) hm.egv[(i)-1]

  double ueps, xh, xl, t, r, f, uepsr, dd, tfe, eps, q, qq;
  int n, i, ne, k, j, m;

  n = wv.nplwin; /* number of plane waves */
  ne = atom.nvband_all; /* number of eigenvalues */

  eps  = EPS;
  ueps = EPS;
  xh = DMAX(D(1)+fabs(D1(1)), D(n)+fabs(D1(n-1)));
  xl = DMIN(D(1)-fabs(D1(1)), D(n)-fabs(D1(n-1)));

  for(i=2;i<=n-1;i++) { /*=== [20] ===*/
    t = fabs(D1(i-1))+fabs(D1(i));
    if((D(i)+t) > xh ) xh = D(i)+t;
    if((D(i)-t) < xl ) xl = D(i)-t;
  }

  r = DMAX(fabs(xh), fabs(xl));
  t = r*ueps;
  xh += t;
  xl -= t;
  for(i=1;i<=n-1;i++) { /*=== [40] ===*/
    W(i) = D1(i)*D1(i);
  }

  f = xl; /* get the eigenvalues in small order */
  t = xh;

  for(i=1;i<=ne;i++) { /*=== [50] ===*/
    E(i) = t;
  }
  uepsr=1.0/ueps;

  for(k=1;k<=ne;k++) { /*=== [120] ===*/
    dd = E(k);
  gt60: 
    t = 0.5*(dd+f);
    tfe = 2.0*(fabs(dd)+fabs(f))*ueps+eps;

    if(fabs(dd-f) < tfe) goto gt115;
    j=0;
    q=D(1)-t;
    if(q>=0.0) j+=1;

    /***************************/
    for(i=2;i<=n;i++) { /*=== [90] ===*/
      if(q != 0.0) goto gt70;
      qq=fabs(D1(i-1))*uepsr;
      goto gt80;
    gt70:
      qq=W(i-1)/q;
    gt80:
      q=D(i)-t-qq;
      if(q>=0.0) j+=1;
    } /* end of [90] */
    /***************************/
    j=n-j; /* get the eigenvalues in small order */
    if(j >= k) goto gt100;
    f = t;
    goto gt60;

  gt100:
    dd=t;
    m=DMIN(j,ne);
    for(i=k;i<=m;i++) { /*=== [110] ===*/
      E(i) = t;
    }
    goto gt60;

  gt115:
    E(k)=t;
  } /*=== end of [120] ===*/ 

#undef A
#undef B
#undef D
#undef D1
#undef IWK
#undef D2
#undef D3
#undef EPS
#undef DMAX
#undef DMIN
#undef W
#undef E
}

void check_eigen_values_small()
{
#define E(i) hm.egv[(i)-1]

  int i;

  printf("----------------------------\n");
  printf("[band]   eigen value\n");
  printf("----------------------------\n");
  for(i=1;i<=atom.nvband_all;i++) printf("[%4d] %20.16f\n",i,E(i));

#undef E
}

