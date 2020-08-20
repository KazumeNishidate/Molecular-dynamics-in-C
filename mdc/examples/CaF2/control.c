#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/***
*  CaF2 soft-core potential setting
***/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void get_control_param(void)
{
  ctl.calc_max = 10000;      /* maximum MD time step                     */
  ctl.delta_time_fs = 2.34;  /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 100.0;         /* Temperature setting [K] */
  ctl.t_control_step = 10;   /* scale at every ctl.t_control_step steps  */

  ctl.set_press_GPa_X = 0.0001;  /* Pressure setting  [GPa] */
  ctl.set_press_GPa_Y = 0.0001; 
  ctl.set_press_GPa_Z = 0.0001; 
  ctl.p_control_step = 20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = 5.9;      /* initial lattice constant [CaF2]              */
  sys.Ay = 5.9;
  sys.Az = 5.9;
  sys.nx = 3;        /* number of unit cells in x-direction */
  sys.ny = 3;
  sys.nz = 3;

  ctl.natoms_in_unit_cell = 12;    /* number of atoms in unit cell          */
  ctl.natoms_in_mol_unit  = 3;     /* number of atoms in primitive Mol unit */
  ctl.kinds_of_ions       = 2;     /* Ca or F */

  sys.a1   = 0.2;    /* alpha setting for EWALD calculation */
  sys.hm   = 23;     /* = |n^2| : cutoff for the reciprocal lattice vector */
  sys.radius = 7.5;  /* cutoff radious [A] in real-space                   */
                     /* sys.radius < Min[MD-basic-cell dimension]/2 [A]    */
}

/***************************************************************************/
/*                   CaF2 specific parameters setting                      */
/***************************************************************************/
void	identify_ion(void)
{
  short x, y, z, ion_cnt, i=0;

  for(z=0; z<sys.nz*2; z++)
    for(y=0; y<sys.ny*2; y++)
      for(x=0; x<sys.nx*2; x++){
	ion_cnt = (x + y + z)%2;
	if(ion_cnt==0){
	  sys.ion[i++] = 0;    /*  Ca --> red  */
	}
	else{
	  sys.ion[i++] = 1;    /*  F  --> blue */
	}
      }
}

/*****
  void set_loc(void)
  [CaF2] Set-up the Ca and F particle locations on CaF2 lattice.
*****/
void	set_loc(void)
{
  short	x, y, z, i=0;
  short nx, ny, nz;
  double nnx, nny, nnz;
  double nx3, ny3, nz3;
  double shiftx, shifty, shiftz;

  nx = sys.nx * 2;
  ny = sys.ny * 2;
  nz = sys.nz * 2;
  shiftx = (sys.Ax)/2.0; /*   X display shift */
  shifty = (sys.Ay)/2.0;
  shiftz = (sys.Az)/2.0;

  /****  Ca ion ****/
  for(z=0; z<nz; z++)
    for(y=0; y<ny; y++)
      for(x=0; x<nx; x++) {
	  if((x+y+z)%2==0){
	    nnx = (double)(x)/2.0;
	    nny = (double)(y)/2.0;
	    nnz = (double)(z)/2.0;
	    sys.rx[i] = sys.Ax * nnx +shiftx;
	    sys.ry[i] = sys.Ay * nny +shifty;
	    sys.rz[i] = sys.Az * nnz +shiftz;
	    sys.ion[i] = 0;
		    i++;
	  }
	}

  /****  F ion ****/
  for(z=0; z<nz; z++)
    for(y=0; y<ny; y++)
      for(x=0; x<nx; x++) {
	  nx3 = ((double)x + 0.5)/2.0;
	  ny3 = ((double)y + 0.5)/2.0;
	  nz3 = ((double)z + 0.5)/2.0;

	  sys.rx[i] = sys.Ax * nx3 +shiftx;
	  sys.ry[i] = sys.Ay * ny3 +shifty;
	  sys.rz[i] = sys.Az * nz3 +shiftz;
	  sys.ion[i] = 1;
	  i++;
	}
}

void set_potential(void)
{
  double m_Ca, m_F, z_Ca, z_F;
/***
*  function form:
*  
*  phy_ij(r) = epsilon((sigma_i + sigma_j)/r)^n + Z_i Z_j(f e)^2/r
*
***/

  z_Ca = 2.0 * sqrt(sys.kk);    /* Ca */
  z_F = -1.0 * sqrt(sys.kk);    /* F  */
  m_Ca = 40.08 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */      /* Ca */
  m_F = 19.00  /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */      /* F  */

  ion.Z[0] = z_Ca;    
  ion.m[0] = m_Ca;
  ion.Z[1] = z_F;
  ion.m[1] = m_F;

  /* soft.epsilon = 0.154 * 1.60217733 * 0.1; */
  soft.epsilon = 0.280 * 1.60217733 * 0.1;  /* 0.28 [eV] */
  soft.sigma_ca = 1.28;
  soft.sigma_f = 1.28;
}

void mk_table(void)   /* make a look up table for [CaF2] potential */
{
  int ion_i, ion_j, ddr;
  int ij = 0;
  double dr, adr, erfc_adr_per_dr, ZZ;
  double p_e1, p_r1, f_e1, f_r1;
  double x = 0;

  for(ion_i=0;ion_i<ctl.kinds_of_ions;ion_i++) {
    for(ion_j=ion_i;ion_j<ctl.kinds_of_ions;ion_j++) {

      ZZ = ion.Z[ion_i] * ion.Z[ion_j];

      for(ddr=0;ddr<sys.table_division+1;ddr++) {  /* division = 0.001 [A]  */

	dr = ((double)ddr)/1000.0 + 0.5; /* 0.5 [A] to sys.radius (cut-off) */
	adr = sys.a1 * dr;  
	erfc_adr_per_dr = erfcc(adr)/dr;

	/* Potential */
	p_e1 = ZZ * erfc_adr_per_dr;
	if(ion_i==0&&ion_j==0){
	  x = (soft.sigma_ca + soft.sigma_ca)/dr;
	}
	if(ion_i==0&&ion_j==1){
	  x = (soft.sigma_ca + soft.sigma_f)/dr;
	}
	if(ion_i==1&&ion_j==1){
	  x = (soft.sigma_f + soft.sigma_f)/dr;
	}

	/*	p_r1 = soft.epsilon *x*x*x*x*x*x*x*x*x*x*x*x;   */
	p_r1 = soft.epsilon *x*x*x*x*x*x*x;

	/* Force */
	f_e1 = ZZ* (2.0 * sys.a3 * exp(- adr * adr) +
		    erfc_adr_per_dr )/(dr*dr);
	f_r1 = 7.0 * soft.epsilon *x*x*x*x*x*x*x /(dr*dr);

	sys.pe1r1[ij*(sys.table_division+1)+ddr] = p_e1 + p_r1;
	sys.fe1r1[ij*(sys.table_division+1)+ddr] = f_e1 + f_r1;
      }

      /* fill the lookup table index */
      sys.lookup[ctl.kinds_of_ions * ion_i + ion_j] = ij;
      sys.lookup[ctl.kinds_of_ions * ion_j + ion_i] = ij;
      ij++;  /* shift the look-up table column */
    }
  }
}
