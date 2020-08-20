#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  void	hamiltonian(int target_atom)
*****/
void	hamiltonian(int target_atom)
{
  int h, i, j, k, cnt=0;
  int ii, jj;
  double dr=0.0;
  double dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;

  double s_R=0.0, s_dash_R, x, y=0.0; /* disposable variable */
  double eru, emu, enu, eru2, emu2, enu2; 
  double v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r;

  /*---------------------  LOCAL  ------------------------------*/
  int local_N; /* total number of particles in LOCAL region */

  local_N = local.num_of_nn_atoms[target_atom];
  /*------------------------------------------------------------*/

  /* make up the vector index [atom ID] of Htb in LOCAL region */
  cnt = 0;
  local.atoms[cnt++] = target_atom;
  for(k=local.address[target_atom];k<local.address[target_atom+1];k++) {
    local.atoms[cnt++] = local.lookup[k];
  }

/*
  for(k=0;k<local_N;k++) {
    printf("ID = %d\n",local.atoms[k]);
  }
*/

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  /* fill up the diagonal elements of local Htb */
  for(k=0;k<local_N;k++) {  
    diagonal(local_N, k);
  }  

  /* fill up the off diagonal elements [right-up part] */
  for(h=0;h<local_N;h++) {
    i = local.atoms[h];

    rxi = sys.rx[i];
    ryi = sys.ry[i];
    rzi = sys.rz[i];

  for(k=h+1;k<local_N;k++) { /* k-j loop */

    j = local.atoms[k];  /* neighbor atom of the target atom */

    dx = rxi - sys.rx[j];
    dy = ryi - sys.ry[j];
    dz = rzi - sys.rz[j];

    /* cyclic boundary condition */
    if(dx > rdx) dx -= sys.Lx;
    if(dy > rdy) dy -= sys.Ly;
    if(dz > rdz) dz -= sys.Lz;
    if(dx < -rdx) dx += sys.Lx;
    if(dy < -rdy) dy += sys.Ly;
    if(dz < -rdz) dz += sys.Lz;

    dr = dx*dx + dy*dy + dz*dz; /* dr' = dr^2 */

    if(dr >= htb.s_rm_2) { /* out-of-range */
      for(ii=0;ii<4;ii++){
	for(jj=0;jj<4;jj++){
	  htb.mat1[16*local_N*h+4*k+4*local_N*ii+jj] = 0.0;
	  htb.mat2x[16*local_N*h+4*k+4*local_N*ii+jj] = 0.0;
	  htb.mat2y[16*local_N*h+4*k+4*local_N*ii+jj] = 0.0;
	  htb.mat2z[16*local_N*h+4*k+4*local_N*ii+jj] = 0.0;
	}
      }
      continue;  /* check here [jump to NEXT the j-loop] */
    }

    else {  /* the range is I or II */

      dr = sqrt(dr);  /* dr' = dr^2 */
      eru = dx/dr;
      emu = dy/dr;
      enu = dz/dr;
      eru2 = eru*eru;
      emu2 = emu*emu;
      enu2 = enu*enu;

      if( dr <= htb.s_r1) { /* range I */

	x = exp(htb.s_nc*log(dr/htb.s_rc)); /* = (r0/rc)^nc */
	y = exp(htb.s_n*log(htb.s_r0/dr));  /* = (r0/dr)^n  */
	s_R = y * exp(htb.s_n*(-x + htb.r0_rc_nc ));
	s_dash_R = -(htb.s_n/dr)*(1.0+htb.s_nc * x)*s_R;
	v_sssig_r = htb.v_sssig * s_R;
	v_spsig_r = htb.v_spsig * s_R;
	v_ppsig_r = htb.v_ppsig * s_R;
	v_pppi_r  = htb.v_pppi * s_R;

	fill_htb1(local_N, h, k, dr, s_R, eru, emu, enu, eru2, 
		  emu2, enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);

	/* non-zero part of htb.mat2[h,k] is {h==0} */
	if(h==0) fill_htb2(local_N, h, k, dr, s_R, s_dash_R, eru, emu, enu, 
	  eru2, emu2, enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r, i, j);
      }

      else { /* range II */

	x = (dr - htb.s_r1);
	s_R = htb.s_c0 + htb.s_c1*x + htb.s_c2*x*x + htb.s_c3*x*x*x;
	s_dash_R = htb.s_c1 + 2.0*htb.s_c2*x + 3.0*htb.s_c3*x*x;
	v_sssig_r = htb.v_sssig * s_R;
	v_spsig_r = htb.v_spsig * s_R;
	v_ppsig_r = htb.v_ppsig * s_R;
	v_pppi_r  = htb.v_pppi * s_R;

	fill_htb1(local_N, h, k, dr, s_R, eru, emu, enu, eru2, 
		  emu2, enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);

	/* non-zero part of htb.mat2[h,k] is {h==0} */
	if(h==0) fill_htb2(local_N, h, k, dr, s_R, s_dash_R, eru, emu, enu, 
	  eru2, emu2, enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r, i, j);
      }

    }
  } /* end of "k-j" loop for neighbor atoms in LOCAL region */

  } /* end of "h-i" loop for neighbor atoms in LOCAL region */

  /* fill up the Htb symmetric matrix [left-down part] */
  for(ii=0;ii<local_N*4;ii++){
    for(jj=ii+1;jj<local_N*4;jj++){
      htb.mat1[4*local_N*jj+ii] = htb.mat1[4*local_N*ii+jj];

      if(ii<4&&jj<4) continue;  /* zero elements */
      if(ii>3&&jj>3) continue;  /* zero elements */
      htb.mat2x[4*local_N*jj+ii] = htb.mat2x[4*local_N*ii+jj];
      htb.mat2y[4*local_N*jj+ii] = htb.mat2y[4*local_N*ii+jj];
      htb.mat2z[4*local_N*jj+ii] = htb.mat2z[4*local_N*ii+jj];
      /*
      press.matX[4*local_N*jj+ii] = press.matX[4*local_N*ii+jj];
      press.matY[4*local_N*jj+ii] = press.matY[4*local_N*ii+jj];
      press.matZ[4*local_N*jj+ii] = press.matZ[4*local_N*ii+jj];
      */
    }
  }

/*********************************************************************/
  /*
  for(ii=0;ii<local_N*4;ii++){
    for(jj=0;jj<local_N*4;jj++){
      printf("%4.3f ",press.matZ[4*local_N*ii+jj]);
    }
    printf("\n");
  }
  exit(0);
  */
/*********************************************************************/


}

