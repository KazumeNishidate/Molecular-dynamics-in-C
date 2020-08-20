#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"
#include  <time.h>

#define RANDOMIZE() srand(time(NULL))      /*  set seed random number  */
#define RANDOM(x) (rand()%(x))             /*  0 =< RANDOM(x) =< x-1   */

/*********************************************************************/
/*  AgI soft-core potential setting                                  */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void   get_control_param(void)
{
  ctl.calc_max = 10000;      /* maximum MD time step                     */
  ctl.delta_time_fs = 1.0;   /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 100.0;          /* Temperature setting [K] */
  ctl.t_control_step = 10;   /* scale at every ctl.t_control_step steps  */

  ctl.set_press_GPa_X = 0.0001;     /* Pressure setting  [GPa] */
  ctl.set_press_GPa_Y = 0.0001; 
  ctl.set_press_GPa_Z = 0.0001; 
  ctl.p_control_step = 20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = 5.08;   /* initial lattice constant [AgI]      */
  sys.Ay = 5.08;
  sys.Az = 5.08;
  sys.nx = 3;      /* number of unit cells of x-direction */
  sys.ny = 3;
  sys.nz = 3;

  ctl.natoms_in_unit_cell = 4;    /* number of atoms in unit cell            */
  ctl.natoms_in_mol_unit  = 2;    /* number of atoms in primitive Mol unit   */
  ctl.kinds_of_ions       = 2;    /* Ag or I */

  sys.a1   = 0.1;    /* alpha setting for EWALD calculation */
  sys.hm   = 23;     /* = |n^2| : cutoff for the reciprocal lattice vector */
  sys.radius = 7.5;  /* cutoff radious [A] in real-space                   */
                     /* sys.radius < Min[MD-basic-cell dimension]/2 [A]    */
}

/***************************************************************************/
/*                     AgI specific parameters setting                     */
/***************************************************************************/

/*****
  void   set_loc(void)
  [AgI] Set-up the Ag and I particle locations on AgI lattice.
*****/

void  identify_ion(void){ /* do nothing here */
    /* see "sys.ion[cnt]" in set_loc() */  }   

void  set_loc(void)
{
  int ii, jj, kk, ran_Ag_men, cnt=0;
  double shiftx, shifty, shiftz; 

  /*  I  */  
  for(ii=0; ii<sys.nx; ii++){        
    for(jj=0; jj<sys.ny; jj++){
      for(kk=0; kk<sys.nz; kk++){

        shiftx = (double)ii;
        shifty = (double)jj;
        shiftz = (double)kk;

        sys.rx[cnt] = (0.0 + shiftx)*sys.Ax;
        sys.ry[cnt] = (0.0 + shifty)*sys.Ay;
        sys.rz[cnt] = (0.0 + shiftz)*sys.Az;
        sys.ion[cnt] = 0;
        cnt++;

        sys.rx[cnt] = (0.5 + shiftx)*sys.Ax;
        sys.ry[cnt] = (0.5 + shifty)*sys.Ay;
        sys.rz[cnt] = (0.5 + shiftz)*sys.Az;
        sys.ion[cnt] = 0;
        cnt++;
      }
    }
  }
  
  /*  Ag  */
  for(ii=0; ii<sys.nx; ii++){         
    for(jj=0; jj<sys.ny; jj++){
      for(kk=0; kk<sys.nz; kk++){

        sys.ion[cnt] = 1;
        sys.ion[cnt+1] = 1;

        shiftx = (double)ii;
        shifty = (double)jj;
        shiftz = (double)kk;

        ran_Ag_men = RANDOM(3);           /* {0, 1, 2} */

        if(ran_Ag_men == 1 || ran_Ag_men == 2){
          switch(RANDOM(4)){
          case 0:
            sys.rx[cnt] = (0.5  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.0  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.25 + shiftz)*sys.Az;
            cnt++;
            break;
          case 1:           
            sys.rx[cnt] = (0.75 + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.0  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.5  + shiftz)*sys.Az;
            cnt++;
            break;
          case 2:
            sys.rx[cnt] = (0.25 + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.0  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.5  + shiftz)*sys.Az;
            cnt++;
            break;
          case 3:
            sys.rx[cnt] = (0.5  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.0  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.75 + shiftz)*sys.Az;
            cnt++;
            break;
          }
        }

        if(ran_Ag_men == 0 || ran_Ag_men == 2){
          switch(RANDOM(4)){
          case 0:
            sys.rx[cnt] = (0.25 + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.5  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.0  + shiftz)*sys.Az;
            cnt++;
            break;
          case 1:
            sys.rx[cnt] = (0.5  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.25 + shifty)*sys.Ay;
            sys.rz[cnt] = (0.0  + shiftz)*sys.Az;
            cnt++;
            break;
          case 2:
            sys.rx[cnt] = (0.5  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.75 + shifty)*sys.Ay;
            sys.rz[cnt] = (0.0  + shiftz)*sys.Az;
            cnt++;
            break;
          case 3:
            sys.rx[cnt] = (0.75 + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.5  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.0  + shiftz)*sys.Az;
            cnt++;
            break;
          }
        }

        if(ran_Ag_men == 1 || ran_Ag_men == 0){
          switch(RANDOM(4)){
          case 0:
            sys.rx[cnt] = (0.0  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.25 + shifty)*sys.Ay;
            sys.rz[cnt] = (0.5  + shiftz)*sys.Az;
            cnt++;
            break;
          case 1:
            sys.rx[cnt] = (0.0  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.5  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.25 + shiftz)*sys.Az;
            cnt++;
            break;
          case 2:
            sys.rx[cnt] = (0.0  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.5  + shifty)*sys.Ay;
            sys.rz[cnt] = (0.75 + shiftz)*sys.Az;
            cnt++;
            break;
          case 3:
            sys.rx[cnt] = (0.0  + shiftx)*sys.Ax;
            sys.ry[cnt] = (0.75 + shifty)*sys.Ay;
            sys.rz[cnt] = (0.5  + shiftz)*sys.Az;
            cnt++;
            break;
          }
        }
      }
    }
  }
}
    
void   set_potential(void)
{
  double m_Ag, m_I, z_Ag, z_I;
  double *soft_core_sigma_mem;
  
  /*--- memory allocation ---*/
  soft_core_sigma_mem =
    (double *)calloc(ctl.kinds_of_interaction,sizeof(double));
  soft_core.sigma = soft_core_sigma_mem;
  /*-------------------------*/

  z_Ag =  1.0*sqrt(sys.kk)*0.6;
  z_I  = -1.0*sqrt(sys.kk)*0.6;
  m_Ag = 107.87/(6.0221367*1.6605402/10.0);   /* [m_u] */
  m_I  = 126.9 /(6.0221367*1.6605402/10.0);   /* [m_u] */

  ion.Z[0] = z_I;
  ion.m[0] = m_I;
  ion.Z[1] = z_Ag;
  ion.m[1] = m_Ag;

  /* Soft-Core-Ion potential parameter set  */
  soft_core.sigma[0] = 2.2;  /* [A] */
  soft_core.sigma[1] = 0.63; /* [A] */
}

void   mk_table(void)  /* make a look up table for Soft Core potential [AgI] */
{
  int ion_i, ion_j, ddr;
  int ij = 0;
  double dr, adr, erfc_adr_per_dr, ZZ;
  double p_e1, p_r1, f_e1, f_r1;
  double p_3, p_7, S_ij;
  double epsilon;

/*--------------------------------------------------------------*/
/*                    (sigma_i + sigma_j)^n   Z_i* Z_j*(fe)^2   */
/*  u_ij(r) = epsilon --------------------- + ---------------   */
/*                             r^n                   r          */
/*--------------------------------------------------------------*/

  epsilon = 0.0851*1.60217733*0.1;  /* 0.0851 [eV] => [energy] */

  for(ion_i=0; ion_i<ctl.kinds_of_ions; ion_i++) {
    for(ion_j=ion_i; ion_j<ctl.kinds_of_ions; ion_j++) {

      ZZ = ion.Z[ion_i]*ion.Z[ion_j];

      S_ij = soft_core.sigma[ion_i] + soft_core.sigma[ion_j];

      for(ddr=0; ddr<sys.table_division+1; ddr++) { /* division = 0.001 [A] */

        dr = ((double)ddr)/1000.0 + 0.5; /* 0.5 [A] to sys.radius (cut-off) */
        adr = sys.a1*dr;  
        erfc_adr_per_dr = erfcc(adr)/dr;

        /* Potential */
        p_3 = (S_ij/dr)*(S_ij/dr)*(S_ij/dr);
        p_7 = p_3*p_3*(S_ij/dr);
        p_e1 = ZZ*erfc_adr_per_dr;
        p_r1 = epsilon*p_7;        /* soft-core-ion specific part */

        /* Force */
        f_e1 = ZZ*(2.0*sys.a3*exp(-adr*adr) + erfc_adr_per_dr)/(dr*dr);
        f_r1 = ((7.0*p_r1)/dr)/dr;           /* soft-core-ion specific part */

        sys.pe1r1[ij*(sys.table_division+1)+ddr] = p_e1 + p_r1;
        sys.fe1r1[ij*(sys.table_division+1)+ddr] = f_e1 + f_r1;
      }

      /* fill the lookup table index */
      sys.lookup[ctl.kinds_of_ions*ion_i + ion_j] = ij;
      sys.lookup[ctl.kinds_of_ions*ion_j + ion_i] = ij;
      ij++;  /* shift the look-up table column */
    }
  }
}

