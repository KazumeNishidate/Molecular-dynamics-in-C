#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

void jacobi_trans(void) {
  int ii, jj, kk, mm, atom;
  int k, m, i=0, j=0, iteration_count=0;
  double akm, amax, theta, co, si, *a;
  double w, w_fx, w_fy, w_fz; /* working valiables */

  static double accuracy;
  static int NN, IMAX, cnt=0;

  /* --- taka ( Wednesday, June 7, 2000 ) --- */
  double fxx=0.0, fyy=0.0, fzz=0.0;
  double pxx=0.0, pyy=0.0, pzz=0.0;
  double w_px, w_py, w_pz;
  /* ---------------------------------------- */

  a = htb.mat1;  


/*======================== JACOBI rotation =========================*/
 /*\_\_\_\_\_\_\_\_\_\_ initialization \_\_\_\_\_\_\_\_\_\_\_*/
  if(cnt==0){
    accuracy = 1.0e-5;
    NN = sys.N*4;
    IMAX = 2*sys.N*sys.N*sys.N;
   /*--- prepare a Unit matrix which will be a Unitary ---*/
    for(k=0; k<NN; k++){
      for(m=0; m<NN; m++)
	if(k != m){ jacobi.Unitary[NN*k+m]=0.0; }
      jacobi.Unitary[NN*k+k]=1.0;
    }
   /*-----------------------------------------------------*/  
    cnt=1;
  }
 /*\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_*/

 /*--------------------- Ut_tb = Ut*tb ----------------------*/
  for(ii=0;ii<NN;ii++){
    for(jj=0;jj<NN;jj++){
      w = 0.0;
      for(kk=0;kk<NN;kk++){
	w += jacobi.Unitary[kk*NN+ii]*a[kk*NN+jj];
      }
      jacobi.Ut_tb[ii*NN+jj] = w;
    }
  }
 /*---------------- tb = Ut_tb*U = Ut*tb*U ------------------*/
  for(ii=0;ii<NN;ii++){
    for(jj=0;jj<NN;jj++){
      w = 0.0;
      for(kk=0;kk<NN;kk++){
	w += jacobi.Ut_tb[ii*NN+kk]*jacobi.Unitary[kk*NN+jj];
      }
      a[ii*NN+jj] = w;
    }
  }
 /*----------------------------------------------------------*/  

  while(IMAX > iteration_count){
    iteration_count++;
 /*----------- search of maximum -----------*/
    amax = 0.0;
    for(k=0; k<NN-1; k++){
      for(m=k+1; m<NN; m++){
	akm = fabs(a[NN*k+m]);
	if(akm > amax){
	  amax = akm;
	  i = k;
	  j = m;
	}
      }
    }
 /*-----------------------------------------*/
    if(amax <= accuracy)
      break;
    else{

      if(a[NN*i+i]==a[NN*j+j])
	theta = PI/4.0;
      else
	theta = 0.5*atan(2.0*a[NN*j+i]/(a[NN*i+i]-a[NN*j+j]));
      
      co = cos(theta);
      si = sin(theta);

      for(k=0; k<NN; k++){
	w = a[NN*i+k];
	a[NN*i+k] =  w*co + a[NN*j+k]*si;
	a[NN*j+k] = -w*si + a[NN*j+k]*co;
	if(k!=i && k!=j){
	  a[NN*k+i] = a[NN*i+k];
	  a[NN*k+j] = a[NN*j+k];
	}
      }

      a[NN*i+i] =  a[NN*i+i]*co + a[NN*i+j]*si;
      a[NN*j+j] = -a[NN*j+i]*si + a[NN*j+j]*co;
      a[NN*i+j] = 0.0;
      a[NN*j+i] = 0.0;
 /*--------- calculate the eigen-vector => Unitary matrix ---------*/
      for(k=0;k<NN;k++){
	w = jacobi.Unitary[NN*k+i];
	jacobi.Unitary[NN*k+i] =  w*co + jacobi.Unitary[NN*k+j]*si;
	jacobi.Unitary[NN*k+j] = -w*si + jacobi.Unitary[NN*k+j]*co;
      }
 /*----------------------------------------------------------------*/
    }
  }
/*====================================== end of JACOBI rotation ====*/

  sort_eigenvalue();  /* in "real.c" */

/*=================== electric FORCE evaluation ====================*/
  for(atom=0;atom<sys.N;atom++){ /* atom = target atom */

    fxx = fyy = fzz = 0.0;
    pxx = pyy = pzz = 0.0;
    
    for(ii=0;ii<NN;ii++){ 
      for(jj=0;jj<NN;jj++){
	jacobi.fx[NN*ii+jj] = 0.0;
	jacobi.fy[NN*ii+jj] = 0.0;
	jacobi.fz[NN*ii+jj] = 0.0;

	press.matX[NN*ii+jj] = 0.0;
	press.matY[NN*ii+jj] = 0.0;
	press.matZ[NN*ii+jj] = 0.0;
      }
    }
    
    for(ii=0;ii<sys.N-1;ii++){
      for(jj=ii+1;jj<sys.N;jj++){
	
	if( (ii==atom)||(jj==atom) ){
	  for(kk=0;kk<4;kk++){
	    for(mm=0;mm<4;mm++){
	      jacobi.fx[16*sys.N*ii+4*jj+NN*mm+kk] =
		htb.mat2x[16*sys.N*ii+4*jj+NN*mm+kk]; 
	      jacobi.fy[16*sys.N*ii+4*jj+NN*mm+kk] = 
		htb.mat2y[16*sys.N*ii+4*jj+NN*mm+kk];
	      jacobi.fz[16*sys.N*ii+4*jj+NN*mm+kk] = 
		htb.mat2z[16*sys.N*ii+4*jj+NN*mm+kk];

              /* ----- taka (Monday, June 5, 2000) ----- */
              press.matX[16*sys.N*ii+4*jj+NN*mm+kk] =
                htb.mat2x[16*sys.N*ii+4*jj+NN*mm+kk]*press.dx[ii*sys.N+jj];
              press.matY[16*sys.N*ii+4*jj+NN*mm+kk] =
                htb.mat2y[16*sys.N*ii+4*jj+NN*mm+kk]*press.dy[ii*sys.N+jj];
              press.matZ[16*sys.N*ii+4*jj+NN*mm+kk] =
                htb.mat2z[16*sys.N*ii+4*jj+NN*mm+kk]*press.dz[ii*sys.N+jj];
              /*-----------------------------------------*/
	      
	      if(ii!=atom) {
		jacobi.fx[16*sys.N*ii+4*jj+4*sys.N*mm+kk] *= -1.0;
		jacobi.fy[16*sys.N*ii+4*jj+4*sys.N*mm+kk] *= -1.0;
		jacobi.fz[16*sys.N*ii+4*jj+4*sys.N*mm+kk] *= -1.0;
	      }
	      
	    }
	  }
	}
	
      }
    }
    
    for(ii=0;ii<NN;ii++){
      for(jj=ii+1;jj<NN;jj++){
	jacobi.fx[NN*jj+ii] = jacobi.fx[NN*ii+jj];
	jacobi.fy[NN*jj+ii] = jacobi.fy[NN*ii+jj];
	jacobi.fz[NN*jj+ii] = jacobi.fz[NN*ii+jj];

        /* ----- taka (Monday, June 5, 2000) ----- */
        press.matX[NN*jj+ii] = press.matX[NN*ii+jj];
        press.matY[NN*jj+ii] = press.matY[NN*ii+jj];
        press.matZ[NN*jj+ii] = press.matZ[NN*ii+jj];
        /*-----------------------------------------*/
      }
    }
/*============================ end of electric FORCE evaluation ====*/

/*=================== unitary transformation =======================*/
    /* taka ( Wednesday, June 7, 2000 ) */
    for(ii=0;ii<NN;ii++){
      if( htb.occupied[ii] == 0 ) continue; /* <--- unoccupied AO */
      for(jj=0;jj<NN;jj++){
	w_fx = 0.0;
	w_fy = 0.0;
	w_fz = 0.0;

        /* ----- taka (Monday, June 5, 2000) ----- */
        w_px = 0.0;
        w_py = 0.0;
        w_pz = 0.0;
        /*-----------------------------------------*/

	for(kk=0;kk<NN;kk++){
	  w_fx += jacobi.Unitary[kk*NN+ii]*jacobi.fx[kk*NN+jj];
	  w_fy += jacobi.Unitary[kk*NN+ii]*jacobi.fy[kk*NN+jj];
	  w_fz += jacobi.Unitary[kk*NN+ii]*jacobi.fz[kk*NN+jj];

          /* ----- taka (Monday, June 5, 2000) ----- */
          w_px += jacobi.Unitary[kk*NN+ii]*press.matX[kk*NN+jj];
          w_py += jacobi.Unitary[kk*NN+ii]*press.matY[kk*NN+jj];
          w_pz += jacobi.Unitary[kk*NN+ii]*press.matZ[kk*NN+jj];
          /*-----------------------------------------*/
	}
        /* ===== Take a trace only for occupied AOs. ===== */
	fxx += w_fx*jacobi.Unitary[jj*NN+ii];
	fyy += w_fy*jacobi.Unitary[jj*NN+ii];
	fzz += w_fz*jacobi.Unitary[jj*NN+ii];

        /* ----- taka (Monday, June 5, 2000) ----- */
        pxx += w_px*jacobi.Unitary[jj*NN+ii];
        pyy += w_py*jacobi.Unitary[jj*NN+ii];
        pzz += w_pz*jacobi.Unitary[jj*NN+ii];
        /*-----------------------------------------*/
      }
    }
/*=============================== end of unitary transformation ====*/

      sys.fx[atom] -= 2.0*fxx;
      sys.fy[atom] -= 2.0*fyy;
      sys.fz[atom] -= 2.0*fzz;

      /* ----- taka (Monday, June 5, 2000) ----- */
      press.virX -= 2.0*pxx;
      press.virY -= 2.0*pyy;
      press.virZ -= 2.0*pzz;
      /*-----------------------------------------*/
  } /* end of atom loop */

}
