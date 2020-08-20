#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  void	real_space(void)

  The function for potential and force calculation of Si diamond structure
  with Stillinger Weber potentail set.

*****/
void	real_space(void)
{
  int i, j, k; 
  double rdx, rdy, rdz;
  double x_ij, y_ij, z_ij;
  double x_jk, y_jk, z_jk;
  double x_ki, y_ki, z_ki;
  double r_ij, r_jk, r_ki;

  /*----- sigma unit ------*/
  double X_ij, Y_ij, Z_ij;
  double X_jk, Y_jk, Z_jk;
  double X_ki, Y_ki, Z_ki;
  double Rij, Rjk, Rki, Rij4;

  double pot_2, pot_3;
  double for_2, f2_1, f2_2;
  double cs_i, cs_j, cs_k;
  double ex_i, ex_j, ex_k;
  double h_i, h_j, h_k;
  double t_i, t_j, t_k;
  double r_i, r_j, r_k;

  /*----- 3Body Force ----*/
  double F1x, F1y, F1z;
  double F2x, F2y, F2z;
  double F3x, F3y, F3z;
  double F4x, F4y, F4z;
  double F5x, F5y, F5z;
  double F6x, F6y, F6z;
  double F7x, F7y, F7z;
  double F8x, F8y, F8z;
  double F9x, F9y, F9z;

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  /* initial set */
  sys.pot = 0.0;
  sys.virX = 0.0;
  sys.virY = 0.0;
  sys.virZ = 0.0;

  pot_3 = 0.0;

  for(i=0; i<sys.N; i++) {
    for(j=i+1; j<sys.N; j++) {
      x_ij = sys.rx[i] - sys.rx[j];
      y_ij = sys.ry[i] - sys.ry[j];
      z_ij = sys.rz[i] - sys.rz[j];
      
      /* cyclic boundary condition */
      if(x_ij > rdx) x_ij -= sys.Lx;
      if(y_ij > rdy) y_ij -= sys.Ly;
      if(z_ij > rdz) z_ij -= sys.Lz;
      if(x_ij < -rdx) x_ij += sys.Lx;
      if(y_ij < -rdy) y_ij += sys.Ly;
      if(z_ij < -rdz) z_ij += sys.Lz;

      if(x_ij/sw.sig>sys.radius || x_ij/sw.sig<-sys.radius) continue;
      if(y_ij/sw.sig>sys.radius || y_ij/sw.sig<-sys.radius) continue;
      if(z_ij/sw.sig>sys.radius || z_ij/sw.sig<-sys.radius) continue;
      
      /*  ion_i - ion_j distance  */
      r_ij = sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
      
      /*  Rij --> sigma UNIT */
      Rij = r_ij/sw.sig;
      
      /* sys.radius = cut-off radius */
      if(Rij > sys.radius) continue;

      /*  Rij4 = Rij^4 */
      Rij4 = Rij * Rij * Rij * Rij;
      
      /* ----------------  SET 2-BODY POTENTIAL ---------------------*/
      pot_2 = sw.ep*sw.A*(sw.B/Rij4-1.0)*exp(1.0/(Rij-sw.a));
      
      /* ----------------  SET 2-BODY FORCE ------------------------*/
      f2_1 = -sw.ep*sw.A*exp(1.0/(Rij-sw.a))/Rij;
      f2_2 = -(4.0*sw.B/(Rij4*Rij)) + (1.0-sw.B/Rij4)/((Rij-sw.a)*(Rij-sw.a));
      
      for_2 = f2_1 * f2_2;

      /*  ion_i - ion_j distance -> sigma UNIT */
      X_ij = x_ij/sw.sig;
      Y_ij = y_ij/sw.sig;
      Z_ij = z_ij/sw.sig;
      
      /* ------------------  SET 3-BODY ------------------------------*/
      
      for(k=j+1;k<sys.N;k++){
	x_jk = sys.rx[j] - sys.rx[k];  
	y_jk = sys.ry[j] - sys.ry[k];  
	z_jk = sys.rz[j] - sys.rz[k];  
	
	/* cyclic boundary condition */
	if(x_jk > rdx) x_jk -= sys.Lx;
	if(y_jk > rdy) y_jk -= sys.Ly;
	if(z_jk > rdz) z_jk -= sys.Lz;
	if(x_jk < -rdx) x_jk += sys.Lx;
	if(y_jk < -rdy) y_jk += sys.Ly;
	if(z_jk < -rdz) z_jk += sys.Lz;
	
	x_ki = sys.rx[k] - sys.rx[i];  
	y_ki = sys.ry[k] - sys.ry[i];  
	z_ki = sys.rz[k] - sys.rz[i];  
	
	/* cyclic boundary condition */
	if(x_ki > rdx) x_ki -= sys.Lx;
	if(y_ki > rdy) y_ki -= sys.Ly;
	if(z_ki > rdz) z_ki -= sys.Lz;
	if(x_ki < -rdx) x_ki += sys.Lx;
	if(y_ki < -rdy) y_ki += sys.Ly;
	if(z_ki < -rdz) z_ki += sys.Lz;

	if(x_jk/sw.sig>sys.radius || x_jk/sw.sig<-sys.radius) continue;
	if(y_jk/sw.sig>sys.radius || y_jk/sw.sig<-sys.radius) continue;
	if(z_jk/sw.sig>sys.radius || z_jk/sw.sig<-sys.radius) continue;
	if(x_ki/sw.sig>sys.radius || x_ki/sw.sig<-sys.radius) continue;
	if(y_ki/sw.sig>sys.radius || y_ki/sw.sig<-sys.radius) continue;
	if(z_ki/sw.sig>sys.radius || z_ki/sw.sig<-sys.radius) continue;

	/* ion_i - ion_j - ion_k distance */
	r_jk = sqrt(x_jk*x_jk + y_jk*y_jk + z_jk*z_jk);
	r_ki = sqrt(x_ki*x_ki + y_ki*y_ki + z_ki*z_ki);
 	
	/* sigma UNIT */
	Rjk = r_jk/sw.sig;
	Rki = r_ki/sw.sig;
	
	/* sys.radius = cut-off radius */
	if(Rjk > sys.radius || Rki > sys.radius) continue;
	
	/* ------- sigma UNIT ---------*/
	
	X_jk = x_jk/sw.sig;
	Y_jk = y_jk/sw.sig;
	Z_jk = z_jk/sw.sig;
	
	X_ki = x_ki/sw.sig;
	Y_ki = y_ki/sw.sig;
	Z_ki = z_ki/sw.sig;

	/* ----------------------------*/

	cs_i = -(x_ij*x_ki + y_ij*y_ki + z_ij*z_ki)/(r_ij * r_ki);
	cs_j = -(x_ij*x_jk + y_ij*y_jk + z_ij*z_jk)/(r_ij * r_jk);
	cs_k = -(x_ki*x_jk + y_ki*y_jk + z_ki*z_jk)/(r_ki * r_jk);
	
	ex_i = exp(sw.gam/(Rij-sw.a)+sw.gam/(Rki-sw.a));
	ex_j = exp(sw.gam/(Rij-sw.a)+sw.gam/(Rjk-sw.a));
	ex_k = exp(sw.gam/(Rki-sw.a)+sw.gam/(Rjk-sw.a));

	r_i = -sw.gam/((Rij-sw.a)*(Rij-sw.a));
	r_j = -sw.gam/((Rjk-sw.a)*(Rjk-sw.a));
	r_k = -sw.gam/((Rki-sw.a)*(Rki-sw.a));

	/*------------------- Note ------------------------------*/
	/*     sw.beta = -2.0*sw.ep*sw.lam -> constant           */
	/*-------------------------------------------------------*/

	t_i = sw.beta*(cs_i+1.0/3.0)*ex_i/(Rij*Rki);
	t_j = sw.beta*(cs_j+1.0/3.0)*ex_j/(Rij*Rjk);
	t_k = sw.beta*(cs_k+1.0/3.0)*ex_k/(Rki*Rjk);

	/* cut-off condition [ion_j_i_k  -> ion i center position] */
	h_i = sw.lam * ex_i * (cs_i+1.0/3.0)*(cs_i+1.0/3.0);
	h_j = sw.lam * ex_j * (cs_j+1.0/3.0)*(cs_j+1.0/3.0);
	h_k = sw.lam * ex_k * (cs_k+1.0/3.0)*(cs_k+1.0/3.0);

	/* ----------------  SET 3-BODY POTENTIAL ---------------------*/
	pot_3 += sw.ep*(h_i + h_j + h_k);

	/* ----------------  SET 3-BODY FORCE -------------------------*/

	F1x = -(X_ij*sw.ep/Rij)*r_i*(h_i + h_j);
	F1y = -(Y_ij*sw.ep/Rij)*r_i*(h_i + h_j);
	F1z = -(Z_ij*sw.ep/Rij)*r_i*(h_i + h_j);

	F2x = (X_ki*sw.ep/Rki)*r_k*(h_i + h_k);
	F2y = (Y_ki*sw.ep/Rki)*r_k*(h_i + h_k);
	F2z = (Z_ki*sw.ep/Rki)*r_k*(h_i + h_k);

	F3x = (-X_jk-Rjk*X_ij*cs_j/Rij)*t_j;
	F3y = (-Y_jk-Rjk*Y_ij*cs_j/Rij)*t_j;
	F3z = (-Z_jk-Rjk*Z_ij*cs_j/Rij)*t_j;
	
	F4x = (X_jk+Rjk*X_ki*cs_k/Rki)*t_k;
	F4y = (Y_jk+Rjk*Y_ki*cs_k/Rki)*t_k;
	F4z = (Z_jk+Rjk*Z_ki*cs_k/Rki)*t_k;
	
	F5x = -(X_jk*sw.ep/Rjk)*r_j*(h_j + h_k);
	F5y = -(Y_jk*sw.ep/Rjk)*r_j*(h_j + h_k);
	F5z = -(Z_jk*sw.ep/Rjk)*r_j*(h_j + h_k);

	F6x = (X_ki+Rki*X_ij*cs_i/Rij)*t_i;
	F6y = (Y_ki+Rki*Y_ij*cs_i/Rij)*t_i;
	F6z = (Z_ki+Rki*Z_ij*cs_i/Rij)*t_i;
	
	F7x = (-X_ki-Rki*X_jk*cs_k/Rjk)*t_k;
	F7y = (-Y_ki-Rki*Y_jk*cs_k/Rjk)*t_k;
	F7z = (-Z_ki-Rki*Z_jk*cs_k/Rjk)*t_k;

	F8x = (-X_ij-Rij*X_ki*cs_i/Rki)*t_i;
	F8y = (-Y_ij-Rij*Y_ki*cs_i/Rki)*t_i;
	F8z = (-Z_ij-Rij*Z_ki*cs_i/Rki)*t_i;

	F9x = (X_ij+Rij*X_jk*cs_j/Rjk)*t_j;
	F9y = (Y_ij+Rij*Y_jk*cs_j/Rjk)*t_j;
	F9z = (Z_ij+Rij*Z_jk*cs_j/Rjk)*t_j;

	sys.fx[i]+=(F1x + F2x + F3x + F4x -F6x -F8x)/sw.sig;
	sys.fy[i]+=(F1y + F2y + F3y + F4y -F6y -F8y)/sw.sig;
	sys.fz[i]+=(F1z + F2z + F3z + F4z -F6z -F8z)/sw.sig;

	sys.fx[j]+=(-F1x + F5x + F6x -F3x -F9x + F7x)/sw.sig; 
	sys.fy[j]+=(-F1y + F5y + F6y -F3y -F9y + F7y)/sw.sig; 
	sys.fz[j]+=(-F1z + F5z + F6z -F3z -F9z + F7z)/sw.sig; 

	sys.fx[k]+=(-F2x -F5x + F8x + F9x -F4x -F7x)/sw.sig;
	sys.fy[k]+=(-F2y -F5y + F8y + F9y -F4y -F7y)/sw.sig;
	sys.fz[k]+=(-F2z -F5z + F8z + F9z -F4z -F7z)/sw.sig;

      }

      sys.pot += pot_2 + pot_3;

      sys.fx[i] += (for_2 * X_ij)/sw.sig;
      sys.fy[i] += (for_2 * Y_ij)/sw.sig;
      sys.fz[i] += (for_2 * Z_ij)/sw.sig;

      sys.fx[j] -= (for_2 * X_ij)/sw.sig;
      sys.fy[j] -= (for_2 * Y_ij)/sw.sig;
      sys.fz[j] -= (for_2 * Z_ij)/sw.sig;

      /* 2-body part for VIRIAL : 3-body contribution is not considered */
      sys.virX += ((for_2 * X_ij)/sw.sig) * x_ij; 
      sys.virY += ((for_2 * Y_ij)/sw.sig) * y_ij; 
      sys.virZ += ((for_2 * Z_ij)/sw.sig) * z_ij; 
    }
  }
}














