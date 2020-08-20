#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

FILE *chgdns, *positions, *energy, *force, *velocity, *eigen, *temp, *for_mesa, *pres, *fmsd;

void  open_files()
{
  chgdns    = fopen("./files/chgdns","w");
  positions = fopen("./files/positions","w");
  energy    = fopen("./files/energy","w");
  force     = fopen("./files/force","w");
  velocity  = fopen("./files/velocity","w");
  eigen     = fopen("./files/eigen","w");
  temp      = fopen("./files/temperature","w");
  pres      = fopen("./files/pressure","w");
  for_mesa  = fopen("./files/for_mesa","w");
  fmsd       = fopen("./files/msd","w");
}

void  close_files()
{
  fclose(chgdns);
  fclose(positions);
  fclose(energy);
  fclose(force);
  fclose(velocity);
  fclose(eigen);
  fclose(temp);
  fclose(pres);
  fclose(for_mesa);
  fclose(fmsd);
}

void dump_temp()
{
  static int cnt=0;
  if(cnt==0) {
    fprintf(temp,"# [MD step] [ionic temp.] [electronic temp.] [fictitious temp.]\n");
    cnt=1;
  }
  fprintf(temp," %d  %e  %e  %e\n", cpa.step, E.cTemp, E.eTemp, E.efTemp);
  fflush(temp);  
}

void dump_pres()
{
  static int cnt=0;
  if(cnt==0) {
    fprintf(pres,"# MDstep  Pxx   Pyy   Pzz   P   [Pascal] \n");
    cnt=1;
  }
  fprintf(pres,"%d   %e  %e  %e  %e\n", cpa.step, ctl.Pxx, ctl.Pyy, ctl.Pzz, ctl.pres);
  fflush(pres);  
}

void dump_eigen(void)
{
  int ib;

  fprintf(eigen,"%4d  ",cpa.step);
  for(ib=0;ib<atom.nvband_all;ib++) { 
    fprintf(eigen,"%.15f ",cp.egvl[ib]);
  }
  fprintf(eigen,"\n");
  fflush(eigen);  
}

void dump_chgdns(void)
{
  int nkx, nky, nkz;
  int id;

  printf(">>> charge density dump\n");
  
  for(nkx=0;nkx<FFTPTSX;nkx++)
    for(nky=0;nky<FFTPTSY;nky++)
      for(nkz=0;nkz<FFTPTSZ;nkz++){
	id = FFTPTSX*FFTPTSY*nkz + FFTPTSX*nky + nkx;
	fprintf(chgdns,"%.15f  ",ft.rho[id]);
	if((nkz+1)%4==0){ fprintf(chgdns,"\n"); }
      }
  fflush(chgdns);
}

void dump_position(void)
{
  int i;

  fprintf(positions,"%4d  ",cpa.step);
  for(i=0; i<atom.N; i++){
    fprintf(positions,"%.15f %.15f %.15f ",
	    cell.rx[i]/au.r, cell.ry[i]/au.r, cell.rz[i]/au.r);
  }
  fprintf(positions,"\n");
  fflush(positions);  
}

void dump_position_mesa(void)
{
  int i;

  for(i=0; i<atom.N; i++){
    fprintf(for_mesa," %d  %f %f %f\n",
            atom.type[i], cell.rx[i]/au.r, cell.ry[i]/au.r, cell.rz[i]/au.r);
  }
  fprintf(for_mesa,"\n");
  fflush(for_mesa);  
}

void dump_msd(void)
{
  int i;

  fprintf(fmsd,"%10d  %10.4e ",cpa.step, cpa.dt*(double)cpa.step);

  for(i=0; i<atom.ntypes; i++){
    fprintf(fmsd," % 8.4e",msd.value[i]);
  } 
  fprintf(fmsd,"\n");
  fflush(fmsd);  
}

/*---------------------------------------------------------------------------*/
/* LDA calculation -> occupied electron = [2]*n                              */
/*                                                                           */
/* E.ekin = [2]*SUM[ n*<p| -(1/2) D^2 |p> ] : electronic kinetic energy      */
/* E.fkin = myu*SUM[ <p'|p'> ]            : fictitious kinetic energy [Real] */
/* E.VpsR = [2]*SUM[ n*<p| Vps |p> ]        : Vps energy [Real]              */
/* E.hart = (1/2) Integral[ rho(r)rho(r')/(r-r') ] : Hartree energy (bug)    */
/* E.xc   = Integral[ Exc(r)rho(r) ]      : exchange + correlation energy    */
/*                                                                           */
/* E.EW   = P(I)+P(II)+P(III)             : EWALD potential                  */
/* E.ckin = (1/2)*m*v^2                   : core (ionic) kinetic energy      */
/*                                                                           */
/* E.eTotal = E.ekin + E.VpsR + E.hart + E.xc                                */
/* E.cTotal = E.EW + E.ckin                                                  */
/* E.TotalF = E.eTotal + E.cTotal + E.fkin                                   */
/* E.Total  = E.eTotal + E.cTotal                                            */
/*                                                                           */
/* E.cTemp  = ionic temperature per atom                                     */
/* E.eTemp  = electronic temperature per electron                            */
/* E.efTemp = fictitious temperature per electron                            */
/*---------------------------------------------------------------------------*/
void dump_energy(void)
{
  static int cnt=0;
  if(cnt==0) {
    fprintf(energy,"# (1)cpa.step (2)cpq.step (3)E.Total (4)E.ekin (5)E.VpsR (6)E.hart (7)E.xc (8)E.EW (9)E.ckin (10)E.fkin (11)E.TotalF (12)E.cTotal (13)E.eTotal (14)ionicTemp (15)eTemp (16)fictitious eTemp  (17)cp.residual_max\n");
    cnt=1;
  }
  fprintf(energy," %d %d % .8e % .8e % .8e  % .8e % .8e % .8e % .8e % .8e  % .8e % .8e % .8e % .8e % .8e % .8e % .8e\n",
	  cpa.step, cpq.step,
	  E.Total, E.ekin, E.VpsR, E.hart, E.xc, E.EW, E.ckin, E.fkin,
	  E.TotalF, E.cTotal, E.eTotal, E.cTemp, E.eTemp, E.efTemp,
	  cp.residual_max);
  fflush(energy);  
}

void dump_force(void)
{
  int i;

  fprintf(force," %4d  ",cpa.step);
  for(i=0; i<atom.N; i++){
    fprintf(force,"%.15e %.15e %.15e ",
	    wk.fx[i], wk.fy[i], wk.fz[i]);
  }
  fprintf(force,"\n");
  fflush(force);  
}

void dump_velocity(void)
{
  int i;

  fprintf(velocity," %4d  ",cpa.step);
  for(i=0; i<atom.N; i++){
    fprintf(velocity,"%.15e %.15e %.15e ",
	    cell.vx[i], cell.vy[i], cell.vz[i]);
  }
  fprintf(velocity,"\n");
  fflush(velocity);  
}

void display_md_start(void)
{
    printf("**********************************************************\n");
    printf("***                      MD start                      ***\n");  
    printf("**********************************************************\n");
}

void display_bo_start(void)
{
    printf("==========================================================\n");
    printf(">>>        approaching to the BO potential surface        \n");
    printf("==========================================================\n");
}

void display_bo_end(void)
{
    printf("==========================================================\n");
    printf(">>>        END of approaching to the BO surface           \n");
    printf("==========================================================\n");
    cpq.step=0;
}

void  display(void)
{
  int i;

  printf(" STEP [%5d:%5d] ( residual = %e )\n",
	 cpa.step, cpq.step, cp.residual_max);
  printf("   E.Total      E.cTotal     E.eTotal     E.EW         E.ckin\n");
  printf("  % 8.4e  % 8.4e  % 8.4e  % 8.4e  % 8.4e\n", 
	 E.Total,E.cTotal,E.eTotal,E.EW,E.ckin);
  printf("   E.ekin       E.fkin       E.VpsR       E.hart       E.xc\n");
  printf("  % 8.4e  % 8.4e  % 8.4e  % 8.4e  % 8.4e\n",
	 E.ekin,E.fkin,E.VpsR,E.hart,E.xc);

  printf("   MSD (rb)^2:");
  for(i=0; i<atom.ntypes; i++){
    printf(" % 8.4e ",msd.value[i]);
  } 
  printf("\n");

  /*  printf("   Pxx          Pyy          Pzz          P  [Pa]\n"); */
  /*  printf("  % 8.4e  % 8.4e  % 8.4e  % 8.4e\n", */
  /*	 ctl.Pxx, ctl.Pyy, ctl.Pzz, ctl.pres); */

  /*  printf("   TEMP [K]: % 9.4f(ion)  % 9.4f(e)  % 9.4f(ef)\n",E.cTemp,E.eTemp,E.efTemp); */

  printf("-------------------------------------------------------------------\n");

}
