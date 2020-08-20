#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "headers/md.h"

/*****
 *  function: void init_atomic_unit(void)
 *  purpose: atomic_unit setting
 *  [Hartree Unit] is <mainly> used throughout this program.
 *****/
void   init_atomic_unit(void)  
{
  /* Atomic Unit
   * [Rydberg Unit] ------------------------------------------------------
   * definition     : m = 1/2,   e^2 = 2,   hbar = 1
   * Mass           = 1 [mRy] = 1/2 [emass] = (1/2)*0.91093897x10^(-30) [Kg]
   * Bohr radius    = 1 [aR] = hbar^2/(m*e^2) = 0.529177249 [A]
   * Energy         = 1 [Ry] = epsiron_1s = m*e^4/(2*hbar^2) = 13.6058 [eV]
   * Time           = 1 [tR] = 2*hbar^3/(m*e^4) = 4.837x10^(-17) [sec]
   *
   * [Hartree Unit] -------------------------------------------------------
   * definition     : m = e^2 = hbar = 1
   * Mass           = 1 [mH] = 1 [emass] = 0.91093897x10^(-30) [Kg]
   * Bohr radius    = 1 [aH] = hbar^2/(m*e^2) = 0.529177249 [A]
   * Energy         = 1 [H]  = 2*epsiron_1s = m*e^4/(hbar^2) = 27.2116 [eV]
   * Time           = 1 [tH] = hbar^3/(m*e^4) = 2.418x10^(-17) [sec]
   * 
   * [unit conversion] ----------------------------------------------------
   * Energy         : 1 [H] = 2 [Ry]
   * Time           : 1 [tH] = 2 [tR]
   * ----------------------------------------------------------------------
   */

  au.hbar = 1.05457266;   /* 1.05457266 x 10^(-34) [J s] = 1 [a.u.]   */
  au.mH   = 9.109389754;  /* 9.109389754 x 10^(-31) [Kg] = 1 [mass H] */
  au.mRy  = au.mH/2.0;    /* 1 [mass H] = 1/2 [mass Ry]               */
  au.e2H  = 2.5669722;    /* x 10^(-19) [C^2]                         */
  au.e2Ry = 2.0*au.e2H;   /* x 10^(-19) [C^2]                         */

  au.r = 1.0/(RBOHR);  /* 0.529177249 [A] = 1 [Bohr radius a0] */
                       /* Aungstrome to A.U.                   */

  au.timeH  = 2.418884337e-17; /* 1 [Hartree time] = 2.418e-17       [sec] */
  au.timeRy = 2.0 * au.timeH;  /* 1 [Rydberg time] = 2 x [Hartree time]    */
  au.EH2J   = 4.3597482e-18;   /* 1 [EH] = 4.3597482e-18 [J]               */
  au.ERy2J  = au.EH2J/2.0;     /* 1/2 [Ry] =  1 [EH]                       */

  /* Ek = (1/2)SUM[mv^2] = (3/2)N Kb T,  H2K : [Hartree]->[K] */
  au.H2K  = EPOT*EZ * (2.0/(3.0*KB));      

  /* [au length] = 0.529177249 (=RBOR) * 10^(-10) [m]          */
  /* [au mass]   = 9.109389754 x 10^(-31) [Kg]                 */
  /* [au time]   = 2.418x10^(-17) [sec]                        */
  /* 1 [Pa] = [N/m^2] = [m^(-1)*kg*sec^(-2)] : [Giga] = 10^(9) */
  au.toPa = (1.0/RBOHR) * au.mH * (1.0/(2.418884337*2.418884337)) * (1.0e13);

}

void  init_velocity(void)
{
  int i;

  for(i=0;i<atom.N;i++){
    cell.vx[i] = 0.0;
    cell.vy[i] = 0.0;
    cell.vz[i] = 0.0;
  }
}
