#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

void  calc_press(void)  
{ 

  /* note : this function is wrong. and must be fixed to include */
  /* "stress" calculated by first-principle manner.              */

  int i;
  double vir_xx=0.0, vir_yy=0.0, vir_zz=0.0;  /* Virial term */

  for(i=0;i<atom.N;i++) { 
    vir_xx += cell.rx[i] * wk.fx[i];
    vir_yy += cell.ry[i] * wk.fy[i];
    vir_zz += cell.rz[i] * wk.fz[i];
  }

  /* Pxx * V = (m*v^2 + Sum[ri*Fi] ) */
  ctl.Pxx = (2.0*E.ckinx + vir_xx)/cell.vol;
  ctl.Pyy = (2.0*E.ckiny + vir_yy)/cell.vol;
  ctl.Pzz = (2.0*E.ckinz + vir_zz)/cell.vol;
  ctl.Pxx *= au.toPa;  /* convert to the [Pa] unit */
  ctl.Pyy *= au.toPa;
  ctl.Pzz *= au.toPa;

  ctl.pres = (ctl.Pxx + ctl.Pyy + ctl.Pzz)/3.0; /* current pressure */

}

void  ctl_temp(int control_step, double Tset)
{
  int i;
  double scale, Eaverage, T;
  static double Esum=0.0;
  static int cnt=0;

  cnt++;
  Esum += E.ckin;

  if(cpa.step % control_step == 0){ 
    Eaverage = Esum/(double)cnt;  /* evaluate the averaged Ekinetic */
    T = Eaverage*au.H2K;         /* convert [Hartree]->[K]         */
    scale = sqrt(Tset/T);         /* evaluate the scaling factor    */

    for(i=0; i<atom.N; i++) {
      cell.vx[i] *= scale;
      cell.vy[i] *= scale;
      cell.vz[i] *= scale;
    }

    Esum = 0.0;
    cnt = 0;
  }

}

