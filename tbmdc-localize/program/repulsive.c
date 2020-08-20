#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  void	repulsive(void)
*****/
void	repulsive(void)
{
  int i, j, k;
  double dr=0.0;
  double dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;

  double x, y=0.0; /* disposable variable */
  double eru, emu, enu, eru2, emu2, enu2; 

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  sys.pot = 0.0;

  for(i=0; i<sys.N; i++) {  /* "i" loop */
    rxi = sys.rx[i];
    ryi = sys.ry[i];
    rzi = sys.rz[i];

    phi.sum_phi[i] = 0.0; /* initialization for the repulsive term */
    phi.f_dash_sum_phi[i] = 0.0;

    for(k=local.address[i];k<local.address[i+1];k++) { /* j loop */

      j = local.lookup[k];  /* neighbor atom of the target "i" atom */

      phi.phi_dash_x[i*sys.N+j] = 0.0;
      phi.phi_dash_y[i*sys.N+j] = 0.0;
      phi.phi_dash_z[i*sys.N+j] = 0.0;

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

      /* --- taka (Wednesday, June 7, 2000) --- */
      press.dx[i*sys.N+j] = dx;
      press.dy[i*sys.N+j] = dy;
      press.dz[i*sys.N+j] = dz;

      /*
      printf("distance[%d %d] %f %f %f\n", i, j,
              press.dx[i*sys.N+j],press.dy[i*sys.N+j],press.dz[i*sys.N+j]);
      getchar();
      */

      if(dr >= htb.s_rm_2) { /* out-of-range [htb.s_rm = phi_rm]*/
	continue;  /* jump to NEXT j-loop */
      }

      else {  /* the range is I or II */

	dr = sqrt(dr);  /* dr' = dr^2 */
	eru = dx/dr;
	emu = dy/dr;
	enu = dz/dr;
	eru2 = eru*eru;
	emu2 = emu*emu;
	enu2 = enu*enu;

	if( dr <= phi.d1 ){ /* range I */

	  x = phi.phi * exp(phi.m*log(phi.d0/dr)) * 
	    exp(phi.m*( -exp(phi.mc*log(dr/phi.dc)) + phi.d0_dc_mc ));
	  phi.sum_phi[i] += x;
	  y = x * (-phi.m/dr)*(1.0+phi.mc* exp(phi.mc*log(dr/phi.dc)));

	  phi.phi_dash_x[i*sys.N+j] = y*eru;
	  phi.phi_dash_y[i*sys.N+j] = y*emu;
	  phi.phi_dash_z[i*sys.N+j] = y*enu;
	}

	else { /* range II */
	  
	  y = (dr - phi.d1);
	  x = phi.c0p + phi.c1p*y + phi.c2p*y*y + phi.c3p*y*y*y;
	  phi.sum_phi[i] += x;

	  x = phi.c1p + 2.0*phi.c2p*y + 3.0*phi.c3p*y*y;

	  phi.phi_dash_x[i*sys.N+j] = x*eru;
	  phi.phi_dash_y[i*sys.N+j] = x*emu;
	  phi.phi_dash_z[i*sys.N+j] = x*enu;
	}
      }
    } /* end of "j" loop */

    /*--- repulsive force evaluation [1st term] ---*/

    y = phi.sum_phi[i];

/*-------------------------------------------------------------*/
/*    printf("y = %8.4f\n",y);                                 */
/*    getchar();                                               */
/*-------------------------------------------------------------*/

    x = phi.c1f + 2.0*phi.c2f*y + 3.0*phi.c3f*y*y + 4.0*phi.c4f*y*y*y;
    phi.f_dash_sum_phi[i] = x;

    sys.pot += 
      phi.c0f + phi.c1f*y + phi.c2f*y*y + phi.c3f*y*y*y + phi.c4f*y*y*y*y;

  } /* end of "i" loop */

/*--------------------------------------------------------------*/
/*  printf("repulsive energy = %8.4f\n",sys.pot);               */
/*  getchar();                                                  */
/*--------------------------------------------------------------*/

  press.virX = press.virY = press.virZ = 0.0;

  /*--- repulsive force evaluation [2nd term] ---*/
  for(i=0; i<sys.N-1; i++) {   /* i = alpha */

    y = phi.f_dash_sum_phi[i];
    for(j=i+1; j<sys.N; j++) { /* j = beta  */

      x = phi.f_dash_sum_phi[j];
      sys.fx[i] -= (x+y) * phi.phi_dash_x[i*sys.N+j];
      sys.fy[i] -= (x+y) * phi.phi_dash_y[i*sys.N+j];
      sys.fz[i] -= (x+y) * phi.phi_dash_z[i*sys.N+j];

      /* --- DELETE ---
      sys.fx[j] -= (x+y) * phi.phi_dash_x[j*sys.N+i];
      sys.fy[j] -= (x+y) * phi.phi_dash_y[j*sys.N+i];
      sys.fz[j] -= (x+y) * phi.phi_dash_z[j*sys.N+i];
      ----------------- */

      sys.fx[j] += (x+y) * phi.phi_dash_x[i*sys.N+j];
      sys.fy[j] += (x+y) * phi.phi_dash_y[i*sys.N+j];
      sys.fz[j] += (x+y) * phi.phi_dash_z[i*sys.N+j];
 
      /* --- taka (Wednesday, June 7, 2000) --- */
      press.virX -= 2.0*(x+y)*phi.phi_dash_x[i*sys.N+j]*press.dx[i*sys.N+j];
      press.virY -= 2.0*(x+y)*phi.phi_dash_y[i*sys.N+j]*press.dy[i*sys.N+j];
      press.virZ -= 2.0*(x+y)*phi.phi_dash_z[i*sys.N+j]*press.dz[i*sys.N+j];

      /* --- DELETE ---
      press.virX -= (x+y)*phi.phi_dash_x[i*sys.N+j]*press.dx[i*sys.N+j];
      press.virY -= (x+y)*phi.phi_dash_y[i*sys.N+j]*press.dy[i*sys.N+j];
      press.virZ -= (x+y)*phi.phi_dash_z[i*sys.N+j]*press.dz[i*sys.N+j];

      press.virX -= (x+y)*phi.phi_dash_x[j*sys.N+i]*press.dx[j*sys.N+i];
      press.virY -= (x+y)*phi.phi_dash_y[j*sys.N+i]*press.dy[j*sys.N+i];
      press.virZ -= (x+y)*phi.phi_dash_z[j*sys.N+i]*press.dz[j*sys.N+i];
      ----------------- */
    }
  }

/*********************************************************************/
  /*
  for(i=0;i<sys.N;i++){
    printf("f[%d] = %f %f %f\n",i,sys.fx[i],sys.fy[i],sys.fz[i]);
  }
  getchar();
  printf("virial(rep) = %f %f %f\n",
          press.virX,press.virY,press.virZ);
  getchar();
  */
/*********************************************************************/

}

