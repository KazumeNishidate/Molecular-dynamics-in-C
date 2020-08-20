#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/* tb.c */
void   fill_htb1(int i, int j,
		 double eru,  double emu,  double enu,
		 double eru2, double emu2, double enu2, 
		 double v_A_r, double v_B_r, double v_C_r, double v_D_r);
void   fill_htb2(int i, int j, double dr,
		 double ds_R_A, double ds_R_B, double ds_R_C, double ds_R_D,
		 double eru,  double emu,  double enu,
		 double eru2, double emu2, double enu2, 
		 double v_A_r, double v_B_r, double v_C_r, double v_D_r);
void   diagonal(int i);


void	real_space(void)
{
  int i, j;
  int ii, jj;
  double dr=0.0;
  double dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;
  double eru, emu, enu, eru2, emu2, enu2;      /* l, m, n */

  double x, y;  /* disposable variable */
  double phi_R, phi_dash_R;

  double s_R;
  double ds_R_A, v_A_r;  /* ss sigma --> A */
  double ds_R_B, v_B_r;  /* sp sigma --> B */
  double ds_R_C, v_C_r;  /* pp sigma --> C */
  double ds_R_D, v_D_r;  /* pp pi    --> D */
 
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
      
/*-- fill up the diagonal Htb matrix elements --*/
      if(i==j) {diagonal(i); continue;}
/*----------------------------------------------*/      

      dx = rxi - sys.rx[j];
      dy = ryi - sys.ry[j];
      dz = rzi - sys.rz[j];
      
      /* cyclic boundary condition */
      if(cluster.type==0){
	if(dx > rdx) dx -= sys.Lx;
	if(dy > rdy) dy -= sys.Ly;
	if(dz > rdz) dz -= sys.Lz;
	if(dx < -rdx) dx += sys.Lx;
	if(dy < -rdy) dy += sys.Ly;
	if(dz < -rdz) dz += sys.Lz;
      }
      dr = dx*dx + dy*dy + dz*dz; /* dr' = dr^2 */
      
      if( (dr >= htb.s_rm_2)&&(j>i) ) { /* out-of-range */
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
      
      else {  /* the range is I or II */
	dr = sqrt(dr);  /* dr' = dr^2 */
	eru = dx/dr;
	emu = dy/dr;
	enu = dz/dr;
	eru2 = eru*eru;
	emu2 = emu*emu;
	enu2 = enu*enu;
	
	if(dr <= htb.s_r1){ /* range I */
	  
	  if(j>i) { /* make up the Htb [right-up part] */

	    y = exp(htb.s_n*log(htb.s_r0/dr));  /* = (r0/r)^n */

	    /*------------------ ss sigma ------------------*/
	    x = exp(htb.s_nc_A*log(dr/htb.s_rc_A));  /* = (r/rc)^nc */
	    s_R = y*exp(htb.s_n*(-x + htb.r0_rc_nc_A));
	    ds_R_A = -(htb.s_n/dr)*(1.0 + htb.s_nc_A*x)*s_R;
	    v_A_r = htb.v_A*s_R;
	    /*------------------ sp sigma ------------------*/	    
	    x = exp(htb.s_nc_B*log(dr/htb.s_rc_B));
	    s_R = y*exp(htb.s_n*(-x + htb.r0_rc_nc_B));
	    ds_R_B = -(htb.s_n/dr)*(1.0 + htb.s_nc_B*x)*s_R;
	    v_B_r = htb.v_B*s_R;
	    /*------------------ pp sigma ------------------*/	    
	    x = exp(htb.s_nc_C*log(dr/htb.s_rc_C));
	    s_R = y*exp(htb.s_n*(-x + htb.r0_rc_nc_C));
	    ds_R_C = -(htb.s_n/dr)*(1.0 + htb.s_nc_C*x)*s_R;
	    v_C_r = htb.v_C*s_R;
	    /*------------------ pp pi ---------------------*/
	    x = exp(htb.s_nc_D*log(dr/htb.s_rc_D));
	    s_R = y*exp(htb.s_n*(-x + htb.r0_rc_nc_D));
	    ds_R_D = -(htb.s_n/dr)*(1.0 + htb.s_nc_D*x)*s_R;
	    v_D_r = htb.v_D*s_R;
	    /*----------------------------------------------*/

/*--------------------------------------------------*/
/*	    printf("dr = %f  s_R = %f \n",dr,s_R);  */
/*--------------------------------------------------*/

	    fill_htb1(i, j, eru, emu, enu, eru2, emu2, enu2, 
		      v_A_r, v_B_r, v_C_r, v_D_r);
	    fill_htb2(i, j, dr, ds_R_A, ds_R_B, ds_R_C, ds_R_D,
		      eru, emu, enu, eru2, emu2, enu2, 
		      v_A_r, v_B_r, v_C_r, v_D_r);
	  }
	  y = exp(phi.m*log(phi.r0/dr));   /* = (r0/r)^m  */
	  x = exp(phi.mc*log(dr/phi.dc));  /* = (r/dc)^mc */
	  phi_R = y*exp(phi.m*(-x + phi.r0_dc_mc));
	  phi.sum_phi[i] += phi_R;
	  phi_dash_R = phi_R*(-phi.m/dr)*(1.0 + phi.mc*x);
	  
	  phi.phi_dash_x[i*sys.N+j] = phi_dash_R*eru;
	  phi.phi_dash_y[i*sys.N+j] = phi_dash_R*emu;
	  phi.phi_dash_z[i*sys.N+j] = phi_dash_R*enu;
	  phi.sum_phi_dash_x[i] += phi_dash_R*eru;
	  phi.sum_phi_dash_y[i] += phi_dash_R*emu;
	  phi.sum_phi_dash_z[i] += phi_dash_R*enu;
	}
	
	else if((dr > htb.s_r1)&&(dr < htb.s_rm)) { /* range II */
	  
	  if(j>i){ /* make up the Htb [right-up part] */

	    x = (dr - htb.s_r1);

	    /*-------------------------- ss sigma --------------------------*/
	    s_R = htb.s_c0A + htb.s_c1A*x + htb.s_c2A*x*x + htb.s_c3A*x*x*x;
	    ds_R_A = htb.s_c1A + 2.0*htb.s_c2A*x + 3.0*htb.s_c3A*x*x;
	    v_A_r = htb.v_A*s_R;
	    /*-------------------------- sp sigma --------------------------*/
	    s_R = htb.s_c0B + htb.s_c1B*x + htb.s_c2B*x*x + htb.s_c3B*x*x*x;
	    ds_R_B = htb.s_c1B + 2.0*htb.s_c2B*x + 3.0*htb.s_c3B*x*x;
	    v_B_r = htb.v_B*s_R;
	    /*-------------------------- pp sigma --------------------------*/
	    s_R = htb.s_c0C + htb.s_c1C*x + htb.s_c2C*x*x + htb.s_c3C*x*x*x;
	    ds_R_C = htb.s_c1C + 2.0*htb.s_c2C*x + 3.0*htb.s_c3C*x*x;
	    v_C_r = htb.v_C*s_R;
	    /*-------------------------- pp pi -----------------------------*/
	    s_R = htb.s_c0D + htb.s_c1D*x + htb.s_c2D*x*x + htb.s_c3D*x*x*x;
	    ds_R_D = htb.s_c1D + 2.0*htb.s_c2D*x + 3.0*htb.s_c3D*x*x;
	    v_D_r = htb.v_D*s_R;
	    /*--------------------------------------------------------------*/

/*--------------------------------------------------*/
/*	    printf("dr = %f  s_R = %f \n",dr,s_R);  */
/*--------------------------------------------------*/

	    fill_htb1(i, j, eru, emu, enu, eru2, emu2, enu2, 
		      v_A_r, v_B_r, v_C_r, v_D_r);
	    fill_htb2(i, j, dr, ds_R_A, ds_R_B, ds_R_C, ds_R_D,
		      eru, emu, enu, eru2, emu2, enu2, 
		      v_A_r, v_B_r, v_C_r, v_D_r);
	  }
	  y = (dr - phi.r1);
	  phi_R = phi.c0p + phi.c1p*y + phi.c2p*y*y + phi.c3p*y*y*y;
	  phi.sum_phi[i] += phi_R;
	  phi_dash_R = phi.c1p + 2.0*phi.c2p*y + 3.0*phi.c3p*y*y;
	  
	  phi.phi_dash_x[i*sys.N+j] = phi_dash_R*eru;
	  phi.phi_dash_y[i*sys.N+j] = phi_dash_R*emu;
	  phi.phi_dash_z[i*sys.N+j] = phi_dash_R*enu;
	  phi.sum_phi_dash_x[i] += phi_dash_R*eru;
	  phi.sum_phi_dash_y[i] += phi_dash_R*emu;
	  phi.sum_phi_dash_z[i] += phi_dash_R*enu;
	}
      }
    } /* end of "j" loop */
    
/*================= repulsive force evaluation [1st term] =================*/
    y = phi.sum_phi[i];
/*------------------------------*/
/*    printf("y = %8.4f\n",y);  */
/*    getchar();                */
/*------------------------------*/
    x = phi.c1f + 2.0*phi.c2f*y + 3.0*phi.c3f*y*y + 4.0*phi.c4f*y*y*y;
    phi.f_dash_sum_phi[i] = x;
    
    sys.fx[i] -= x*phi.sum_phi_dash_x[i];
    sys.fy[i] -= x*phi.sum_phi_dash_y[i];
    sys.fz[i] -= x*phi.sum_phi_dash_z[i];
    
    sys.pot += phi.c1f*y + phi.c2f*y*y + phi.c3f*y*y*y + phi.c4f*y*y*y*y;
/*=========================================================================*/
    
  } /* end of "i" loop */
/*-------------------------------------------------*/
/*  printf("repulsive energy = %8.4f\n",sys.pot);  */
/*  getchar();                                     */
/*-------------------------------------------------*/
/*================= repulsive force evaluation [2nd term] =================*/
  for(i=0; i<sys.N; i++) {   /* i = alpha */

    y = phi.f_dash_sum_phi[i];

    for(j=0; j<sys.N; j++) { /* j = beta  */
      
      if(j==i){continue;} /* jump to next "j" when i==j */
      
      x = phi.f_dash_sum_phi[j];
      sys.fx[i] -= x * phi.phi_dash_x[i*sys.N+j];
      sys.fy[i] -= x * phi.phi_dash_y[i*sys.N+j];
      sys.fz[i] -= x * phi.phi_dash_z[i*sys.N+j];

      press.virX -= (x+y)*phi.phi_dash_x[i*sys.N+j]*press.dx[i*sys.N+j];
      press.virY -= (x+y)*phi.phi_dash_y[i*sys.N+j]*press.dy[i*sys.N+j];
      press.virZ -= (x+y)*phi.phi_dash_z[i*sys.N+j]*press.dz[i*sys.N+j];
    }
  }
/*=========================================================================*/

/*-------------------------------------------------------------------*/
/*  
  for(ii=0;ii<sys.N;ii++){
    printf("f[%d] = %f %f %f\n",ii,sys.fx[ii],sys.fy[ii],sys.fz[ii]);
  }
  getchar();                                                         
*/  
/*-------------------------------------------------------------------*/

  /* fill up the Htb symmetric matrix [left-down part] */
  for(ii=0;ii<sys.N*4;ii++){
    for(jj=ii+1;jj<sys.N*4;jj++){
      htb.mat1[4*sys.N*jj+ii] = htb.mat1[4*sys.N*ii+jj];
      htb.mat2x[4*sys.N*jj+ii] = htb.mat2x[4*sys.N*ii+jj];
      htb.mat2y[4*sys.N*jj+ii] = htb.mat2y[4*sys.N*ii+jj];
      htb.mat2z[4*sys.N*jj+ii] = htb.mat2z[4*sys.N*ii+jj];
    }
  }
  
/*-----------------------------------------------*/
/* 
   for(ii=0;ii<sys.N*4;ii++){
    for(jj=0;jj<sys.N*4;jj++){
      printf("%4.3f ",htb.mat2z[4*sys.N*ii+jj]);
    }
    printf("\n");
  }
  exit(0);                                       
*/
/*-----------------------------------------------*/

/*  
  for(i=0; i<sys.N*4*sys.N*4; i++){
    if(i%32==0){printf("\n");}
    if(htb.mat1[i]>=0){printf(" ");}
    printf("%f",htb.mat1[i]);
  }
  getchar();
*/
}


/****************************************************************************/
/* determine the occupied atomic orbit                                      */
/****************************************************************************/
void sort_eigenvalue(void)
{
  int i, j, NN, address;
  double temp;

  NN = sys.N*4;

  for(i=0;i<NN;i++) {
    wv.eigen_values_address[i] = i;
    wv.eigen_values[i] = htb.mat1[NN*i+i];
    htb.occupied[i] = 0;
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
  for(i=0;i<sys.N*2;i++) {
    htb.occupied[wv.eigen_values_address[i]] = 1;
  }
}


/****************************************************************************/
/* electronic potential evaluation for occupied AOs                         */
/****************************************************************************/
void calc_elec_E(void)
{
  int i;
  double tb_potential=0.0;

  for(i=0;i<sys.N*2;i++) {
    if( htb.occupied[i] == 0 ) continue;
    tb_potential += htb.mat1[i*sys.N*4+i];
  }
  sys.pot += 2.0 * tb_potential;
  sys.pot += htb.e_0*sys.N;  

/*------------------------------------------------------------*/
/*  printf(" 2.0 * tb_potential = %f\n",2.0 * tb_potential);  */
/*  getchar();                                                */
/*------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*  printf("total potential/number of atoms = %f\n",sys.pot/(double)sys.N);  */
/*  getchar();                                                               */
/*---------------------------------------------------------------------------*/

/*--------------------------------------------------------*/
/*  
  printf("--- Eigenvalue --- \n");
  for(i=0;i<sys.N*2;i++){
    printf("%i : %10.7f \n",i+1,htb.mat1[htb.occupied[i]]);
  }    
  printf("Tight-Binding energy(sys.pot) = %f\n",sys.pot);
  getchar();                                              
*/    
/*--------------------------------------------------------*/
}
