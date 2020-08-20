#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/* tb.c */
void   fill_htb1(int i, int j, double dr, double s_R, 
		 double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r);
void   fill_htb2(int i, int j, double dr, double s_R, double s_dash_R,
		 double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r);
void   diagonal(int i);


/*****
  void	real_space(void)
*****/
void	real_space(void)
{
  int i, j;
  int ii, jj;

  int ij; /* ----- taka (Monday, June 5, 2000) ----- */

  double dr=0.0;
  double dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;

  double s_R=0.0, s_dash_R, x, y=0.0; /* disposable variable */
  double eru, emu, enu, eru2, emu2, enu2; 
  double v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r;

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  sys.pot = 0.0;

  for(i=0; i<sys.N; i++) {  /* "i" loop */
    rxi = sys.rx[i];
    ryi = sys.ry[i];
    rzi = sys.rz[i];

    phi.sum_phi[i] = 0.0; /* initialization for the repulsive term */
    phi.sum_phi_dash_x[i] = 0.0;
    phi.sum_phi_dash_y[i] = 0.0;
    phi.sum_phi_dash_z[i] = 0.0;
    phi.f_dash_sum_phi[i] = 0.0;

    for(j=0; j<sys.N; j++) {  /* "j" loop */

      phi.phi_dash_x[i*sys.N+j] = 0.0;
      phi.phi_dash_y[i*sys.N+j] = 0.0;
      phi.phi_dash_z[i*sys.N+j] = 0.0;

      if(i==j) {  /* fill up the diagonal Htb matrix elements */
	diagonal(i);
	continue;  /* jump to next "j" when i==j */
      }

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

      /* ----- taka (Monday, June 5, 2000) ----- */
      ij = i*sys.N+j;
      press.dx[ij] = dx;
      press.dy[ij] = dy;
      press.dz[ij] = dz;
      /*
      printf( "distance[%d %d]=[%f,%f,%f]\n",
          i, j, press.dx[ij],press.dy[ij],press.dz[ij]);
      */

      if((dr >= htb.s_rm_2)&&(j>i)) { /* out-of-range */
	for(ii=0;ii<4;ii++){
	  for(jj=0;jj<4;jj++){
	    htb.mat1[16*sys.N*i+4*j+4*sys.N*ii+jj] = 0.0;
	    htb.mat2x[16*sys.N*i+4*j+4*sys.N*ii+jj] = 0.0;
	    htb.mat2y[16*sys.N*i+4*j+4*sys.N*ii+jj] = 0.0;
	    htb.mat2z[16*sys.N*i+4*j+4*sys.N*ii+jj] = 0.0;
	  }
	}
	continue;  /* check here */
      }

      else {  /* the range is I or II or III */

	dr = sqrt(dr);  /* dr' = dr^2 */
	eru = dx/dr;
	emu = dy/dr;
	enu = dz/dr;
	eru2 = eru*eru;
	emu2 = emu*emu;
	enu2 = enu*enu;

	if( dr <= htb.s_r1){ /* range I */

	  if(j>i) { /* make up the Htb [right-up part] */
	    x = exp(htb.s_nc*log(dr/htb.s_rc)); /* = (r0/rc)^nc */
	    y = exp(htb.s_n*log(htb.s_r0/dr));  /* = (r0/dr)^n  */
	    s_R = y * exp(htb.s_n*(-x + htb.r0_rc_nc ));
	    s_dash_R = -(htb.s_n/dr)*(1.0+htb.s_nc * x)*s_R;
	    v_sssig_r = htb.v_sssig * s_R;
	    v_spsig_r = htb.v_spsig * s_R;
	    v_ppsig_r = htb.v_ppsig * s_R;
	    v_pppi_r  = htb.v_pppi * s_R;

	    /* printf("dr = %f  s_R = %f  s_dash_R = %f\n",dr,s_R,s_dash_R); */

	    fill_htb1(i, j, dr, s_R, eru, emu, enu, eru2, emu2, enu2, 
		      v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);
	    fill_htb2(i, j, dr, s_R, s_dash_R, eru, emu, enu, eru2, emu2, 
		      enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);
	  }

	  x = phi.phi * exp(phi.m*log(phi.d0/dr)) * 
	    exp(phi.m*( -exp(phi.mc*log(dr/phi.dc)) + phi.d0_dc_mc ));
	  phi.sum_phi[i] += x;
	  y = x * (-phi.m/dr)*(1.0+phi.mc * exp(phi.mc*log(dr/phi.dc)));

          /* printf("dr = %f  x = %f  y = %f\n",dr,x,y); */

	  phi.phi_dash_x[i*sys.N+j] = y*eru;
	  phi.phi_dash_y[i*sys.N+j] = y*emu;
	  phi.phi_dash_z[i*sys.N+j] = y*enu;
	  phi.sum_phi_dash_x[i] += y*eru;
	  phi.sum_phi_dash_y[i] += y*emu;
	  phi.sum_phi_dash_z[i] += y*enu;
	}

	else if((dr > htb.s_r1)&&(dr <= phi.d1)) { /* range II */

	  if(j>i){ /* make up the Htb [right-up part] */
	    x = (dr - htb.s_r1);
	    s_R = htb.s_c0 + htb.s_c1*x + htb.s_c2*x*x + htb.s_c3*x*x*x;
	    s_dash_R = htb.s_c1 + 2.0*htb.s_c2*x + 3.0*htb.s_c3*x*x;
	    v_sssig_r = htb.v_sssig * s_R;
	    v_spsig_r = htb.v_spsig * s_R;
	    v_ppsig_r = htb.v_ppsig * s_R;
	    v_pppi_r  = htb.v_pppi * s_R;

	    /* printf("dr = %f  s_R = %f  s_dash_R = %f\n",dr,s_R,s_dash_R); */

	    fill_htb1(i, j, dr, s_R, eru, emu, enu, eru2, emu2, enu2, 
		      v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);
	    fill_htb2(i, j, dr, s_R, s_dash_R, eru, emu, enu, eru2, emu2, 
		      enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);
	  }

	  x = phi.phi * exp(phi.m*log(phi.d0/dr)) * 
	    exp(phi.m*( -exp(phi.mc*log(dr/phi.dc)) + phi.d0_dc_mc ));
	  phi.sum_phi[i] += x;
	  y = x * (-phi.m/dr)*(1.0+phi.mc* exp(phi.mc*log(dr/phi.dc)));

          /* printf("dr = %f  x = %f  y = %f\n",dr,x,y); */

	  phi.phi_dash_x[i*sys.N+j] = y*eru;
	  phi.phi_dash_y[i*sys.N+j] = y*emu;
	  phi.phi_dash_z[i*sys.N+j] = y*enu;
	  phi.sum_phi_dash_x[i] += y*eru;
	  phi.sum_phi_dash_y[i] += y*emu;
	  phi.sum_phi_dash_z[i] += y*enu;
	}

	else if((dr > phi.d1)&&(dr < htb.s_rm)) { /* range III */
	  
	  if(j>i){ /* make up the Htb [right-up part] */
	    x = (dr - htb.s_r1);
	    s_R = htb.s_c0 + htb.s_c1*x + htb.s_c2*x*x + htb.s_c3*x*x*x;
	    s_dash_R = htb.s_c1 + 2.0*htb.s_c2*x + 3.0*htb.s_c3*x*x;

	    v_sssig_r = htb.v_sssig * s_R;
	    v_spsig_r = htb.v_spsig * s_R;
	    v_ppsig_r = htb.v_ppsig * s_R;
	    v_pppi_r  = htb.v_pppi * s_R;

	    /* printf("dr = %f  s_R = %f  s_dash_R = %f\n",dr,s_R,s_dash_R); */

	    fill_htb1(i, j, dr, s_R, eru, emu, enu, eru2, emu2, enu2, 
		      v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);
	    fill_htb2(i, j, dr, s_R, s_dash_R, eru, emu, enu, eru2, emu2, 
		      enu2, v_sssig_r, v_spsig_r, v_ppsig_r, v_pppi_r);
	  }

	  y = (dr - phi.d1);
	  x = phi.c0p + phi.c1p*y + phi.c2p*y*y + phi.c3p*y*y*y;
	  phi.sum_phi[i] += x;

          /* printf("dr = %f  x = %f  ",dr,x); */

	  x = phi.c1p + 2.0*phi.c2p*y + 3.0*phi.c3p*y*y;

          /* printf("y = %f\n",y); */

	  phi.phi_dash_x[i*sys.N+j] = x*eru;
	  phi.phi_dash_y[i*sys.N+j] = x*emu;
	  phi.phi_dash_z[i*sys.N+j] = x*enu;
	  phi.sum_phi_dash_x[i] += x*eru;
	  phi.sum_phi_dash_y[i] += x*emu;
	  phi.sum_phi_dash_z[i] += x*enu;
	}
      }
      /* getchar(); */
    } /* end of "j" loop */

    /*--- repulsive force evaluation [1st term] ---*/

    y = phi.sum_phi[i];

/*-------------------------------------------------------------*/
/*    printf("y = %8.4f\n",y);                                 */
/*    getchar();                                               */
/*-------------------------------------------------------------*/

    x = phi.c1f + 2.0*phi.c2f*y + 3.0*phi.c3f*y*y + 4.0*phi.c4f*y*y*y;
    phi.f_dash_sum_phi[i] = x;

    sys.fx[i] -= x * phi.sum_phi_dash_x[i];
    sys.fy[i] -= x * phi.sum_phi_dash_y[i];
    sys.fz[i] -= x * phi.sum_phi_dash_z[i];

    sys.pot += 
      phi.c0f + phi.c1f*y + phi.c2f*y*y + phi.c3f*y*y*y + phi.c4f*y*y*y*y;

  } /* end of "i" loop */

/*--------------------------------------------------------------*/
  /*printf("repulsive energy = %8.4f\n",sys.pot);               */
  /*getchar();                                                  */
  /*
  for(ii=0;ii<sys.N;ii++){
    printf("f[%d] = %.10f %.10f %.10f\n",ii,sys.fx[ii],sys.fy[ii],sys.fz[ii]);
  }
  printf( "\n" );
  */
/*--------------------------------------------------------------*/

  x = 0.0; y = 0.0; /* ----- taka (Monday, June 5, 2000) ----- */
  press.virX = press.virY = press.virZ = 0.0;

  /*--- repulsive force evaluation [2nd term] ---*/
  for(i=0; i<sys.N; i++) {   /* i = alpha */

    y = phi.f_dash_sum_phi[i]; /* ----- taka (Monday, June 5, 2000) ----- */

    for(j=0; j<sys.N; j++) { /* j = beta  */

      if(j==i) continue; /* jump to next "j" when i==j */

      ij = i*sys.N+j; /* ----- taka (Monday, June 5, 2000) ----- */

      x = phi.f_dash_sum_phi[j];
      sys.fx[i] -= x * phi.phi_dash_x[ij];
      sys.fy[i] -= x * phi.phi_dash_y[ij];
      sys.fz[i] -= x * phi.phi_dash_z[ij];

      /* ----- taka (Monday, June 5, 2000) ----- */
      press.virX -= (x+y)*phi.phi_dash_x[ij]*press.dx[ij];
      press.virY -= (x+y)*phi.phi_dash_y[ij]*press.dy[ij];
      press.virZ -= (x+y)*phi.phi_dash_z[ij]*press.dz[ij];

      /* -----
      printf( "repulsive force [%d,%d]\n",i, j );
      printf( "x component = %.8f,y component = %.8f,z component = %.8f\n",
          -(x+y)*phi.phi_dash_x[ij],
          -(x+y)*phi.phi_dash_y[ij],
          -(x+y)*phi.phi_dash_z[ij] );
      ----- */
    }
    /* -----
    printf("f[%d] = %.10f %.10f %.10f\n",i,
        press.virX, press.virY, press.virZ );
    press.virX = press.virY = press.virZ = 0.0;
    ----- */
  }
  /* getchar(); */
/*********************************************************************/
  /*
  for(ii=0;ii<sys.N;ii++){
    printf("f[%d] = %.10f %.10f %.10f\n",ii,sys.fx[ii],sys.fy[ii],sys.fz[ii]);
  }
  getchar();
  */
/*********************************************************************/

  /* fill up the Htb symmetric matrix [left-down part] */
  for(ii=0;ii<sys.N*4;ii++){
    for(jj=ii+1;jj<sys.N*4;jj++){
      htb.mat1[4*sys.N*jj+ii] = htb.mat1[4*sys.N*ii+jj];
      htb.mat2x[4*sys.N*jj+ii] = htb.mat2x[4*sys.N*ii+jj];
      htb.mat2y[4*sys.N*jj+ii] = htb.mat2y[4*sys.N*ii+jj];
      htb.mat2z[4*sys.N*jj+ii] = htb.mat2z[4*sys.N*ii+jj];
    }
  }

/*********************************************************************/
/*
  for(ii=0;ii<sys.N*4;ii++){
    for(jj=0;jj<sys.N*4;jj++){
      printf("%4.3f ",htb.mat2z[4*sys.N*ii+jj]);
    }
    printf("\n");
  }
exit(0);
*/
/*********************************************************************/

}

/*-------------------------------------------------------*/
/* determine the occupied atomic orbit                   */
/*-------------------------------------------------------*/
void sort_eigenvalue(void)
{
  int i, j, NN, address;
  double temp;

  NN=sys.N*4;

  for(i=0;i<NN;i++) {
    wv.eigen_values_address[i] = i;
    wv.eigen_values[i] = htb.mat1[NN*i+i];
    htb.occupied[i] = 0; /* taka ( Wednesday, June 7, 2000 ) */
  }

  for(i=0;i<NN-1;i++) {
      for(j=i+1;j<NN;j++) {
          if(wv.eigen_values[i] > wv.eigen_values[j])
            {
              temp = wv.eigen_values[i];
	      address = wv.eigen_values_address[i];

              wv.eigen_values[i] = wv.eigen_values[j];
	      wv.eigen_values_address[i] = wv.eigen_values_address[j];

              wv.eigen_values[j] = temp;
	      wv.eigen_values_address[j] = address;
            }
        }
    }
  for(i=0;i<sys.N*2;i++) { /* taka ( Wednesday, June 7, 2000 ) */
    htb.occupied[wv.eigen_values_address[i]] = 1;
  }

}

/* electronic potential evaluation for occupied AOs */
void calc_elec_E(void)
{
  int i;
  double tb_potential=0.0;

  for(i=0;i<sys.N*4;i++) {
    if( htb.occupied[i] == 0 ) continue;
    tb_potential += htb.mat1[i*sys.N*4+i];
  }
  sys.pot += 2.0 * tb_potential;
/*
printf(" 2.0 * tb_potential = %f\n",2.0 * tb_potential);
getchar();
*/
/*--------------------------------------------------------------------------*/
/*
  printf("total potential / number of atoms  = %f\n",sys.pot/(double)sys.N);
  getchar();
*/
/*--------------------------------------------------------------------------*/

/*********************************************************************/
/*
  printf("--- Eigenvalue --- \n");
  for(i=0;i<sys.N*2;i++) printf("%i : %10.7f \n",i+1,htb.mat1[htb.occupied[i]]);
  printf("Tight-Binding energy(sys.pot) = %f\n",sys.pot);
  getchar();
*/
/*********************************************************************/
}


