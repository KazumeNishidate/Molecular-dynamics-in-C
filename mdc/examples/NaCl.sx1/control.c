#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/*  NaCl SX1 potential setting                                       */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void get_control_param(void)
{
  ctl.calc_max = 10000;       /* maximum MD time step                       */
  ctl.delta_time_fs = 3.0;    /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s]   */

  ctl.temp = 300.0;           /* Temperature setting [K]                    */
  ctl.t_control_step = 5;     /* scale at every ctl.t_control_step steps    */

  ctl.set_press_GPa_X = 1.0;  /* Pressure setting  [GPa]                    */
  ctl.set_press_GPa_Y = 1.0; 
  ctl.set_press_GPa_Z = 1.0; 
  ctl.p_control_step  =  20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = 5.63;              /* lattice constant of unit cell [NaCl]       */
  sys.Ay = 5.63;
  sys.Az = 5.63;
  sys.nx = 3;	              /* number of unit cells in x-direction        */
  sys.ny = 3;
  sys.nz = 3;

  ctl.natoms_in_unit_cell = 8;     /* number of atoms in unit cell           */
  ctl.natoms_in_mol_unit  = 2;     /* number of atoms in primitive Mol unit  */
  ctl.kinds_of_ions       = 2;     /* Na or Cl */

  sys.a1     = 0.3; /* alpha setting for EWALD calculation               */
  sys.hm     = 20;      /* = |n^2| : cutoff for the reciprocal lattice vector*/
  sys.radius = 8.0;     /* cutoff radious [A] in real-space                  */
                        /* sys.radius < Min[MD-basic-cell dimension]/2 [A]   */
}

void	identify_ion(void)
{
  short	x, y, z, i=0;

  for(z=0; z<sys.nz*2; z++)
    for(y=0; y<sys.ny*2; y++)
      for(x=0; x<sys.nx*2; x++)
	sys.ion[i++] = (x%2 + y%2 + z%2)%2;  /* [Na=0 or Cl=1] */
}

/*****
  void set_loc(void)
  [NaCl] Set-up the Na and Cl particle locations on NaCl lattice.
*****/
void	set_loc(void)
{
  short	x, y, z, i=0;
  short nx, ny, nz;
  double nnx, nny, nnz;
  double shiftx, shifty, shiftz;

  nx = sys.nx * 2;
  ny = sys.ny * 2;
  nz = sys.nz * 2;

  nnx = (double)nx;
  nny = (double)ny; 
  nnz = (double)nz; 

  shiftx = (sys.Lx/nnx)/2.0;
  shifty = (sys.Ly/nny)/2.0;
  shiftz = (sys.Lz/nnz)/2.0;

  for(z=0; z<nz; z++)
    for(y=0; y<ny; y++)
      for(x=0; x<nx; x++) {
	sys.rx[i] = sys.Lx*(double)(x/nnx)+shiftx;
	sys.ry[i] = sys.Ly*(double)(y/nny)+shifty;
	sys.rz[i] = sys.Lz*(double)(z/nnz)+shiftz;
	i++;
      }
}

void set_potential(void)
{
  double m_na, m_cl, z_na, z_cl;
  double *sx1_alpha_mem, *sx1_beta_mem;

  /*--- memory allocation ---*/
  sx1_alpha_mem = (double *)calloc(ctl.kinds_of_interaction,
				   sizeof(double));
  sx1_beta_mem = (double *)calloc(ctl.kinds_of_interaction,
				  sizeof(double));
  sx1.alpha = sx1_alpha_mem;
  sx1.beta  = sx1_beta_mem;
  /*-------------------------*/

  /*------------------------------------------------------------------------*/
  /* << mass unit conversion >>                                             */ 
  /*                                                                        */
  /*                  [g]                          10^{-3} [kg]             */
  /*  m_na = 22.9898 ------ = 22.9898 --------------------------------------*/
  /*                 [mol]            [mol] x 6.0221367 x 10^{23} [mol^{-1}]*/
  /*           22.9898 x 10^{-3}                    1                       */ 
  /*       = ---------------------- [kg] x ---------------------- [m_u/kg]  */
  /*          6.0221367 x 10^{23}           1.6605402 x 10^{-27}            */
  /*                                     1                                  */
  /*       = 22.9898 x --------------------------------------- [m_u]        */
  /*                     6.0221367 x 1.6605402 x 10^{3+23-27}               */
  /*                                  1                                     */
  /*       = 22.9898 x  ---------------------------------- [m_u]            */
  /*                     6.0221367 x 1.6605402 x 10^{-1}                    */
  /*                                  1                                     */
  /*  m_cl = 35.4528 x  ---------------------------------- [m_u]            */
  /*                     6.0221367 x 1.6605402 x 10^{-1}                    */
  /*                                                                        */
  /*  # see init_sys() in "init.c" for the definition of sys.kk.            */
  /*                                                                        */
  /*------------------------------------------------------------------------*/
  /* ---------------------------------------------------------- */
  /* SX-1 potential parameter set: MEMO                         */
  /*                                                            */
  /* function form:                                             */
  /* u_(i,j) = Z_(i)*Z_(j)*e^(2)/r_(i,j) +                      */
  /*    f(b_(i)+b_(j))*Exp[(a_(i)+a_(j)-r_(i,j))/(b_(i)+b_(j))  */
  /*                                                            */
  /* f0=1 [kcal/(A*mol)]                                        */ 
  /* 1 [kcal] = 4.1868*10^3  [J=N*m]                            */
  /* 1 [A] = 10^(-10) [m]                                       */
  /*                                                            */
  /* for NaCl crystal system:                                   */
  /* a_Na = 1.260 [A]                                           */ 
  /* a_Cl = 1.950 [A]                                           */
  /* b_Na = 0.080 [A] or 0.085 [A]                              */
  /* b_Cl = 0.090 [A]                                           */
  /* ---------------------------------------------------------- */

  z_na = 1.0 * sqrt(sys.kk);
  z_cl = -1.0 * sqrt(sys.kk);
  m_na = 22.9898 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */
  m_cl = 35.4528 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.Z[0] = z_na;
  ion.m[0] = m_na;
  ion.Z[1] = z_cl;
  ion.m[1] = m_cl;

  sx1.alpha[0] = 1.260;
  sx1.alpha[1] = 1.950;
  sx1.beta[0]  = 0.080;
  sx1.beta[1]  = 0.090;

  /*------------------------------------------------------------------------*/
  /* << f0 unit conversion >>                                               */
  /*     =>  1 [kcal][mol^{-1}] = 4.18605 x 10^{3} [J][mol^{-1}]            */
  /*     =>	 avogadro constant: 6.0221367e23                            */
  /*                                                                        */
  /*             [kcal]                               [J]                   */
  /*  f0 = 1.0 ---------- = 1.0 x 4.18605 x 10^{3} ----------               */
  /*            [A][mol]                            [A][mol]                */
  /*                               [J] x 10^{18} [energy]/[J]      1        */
  /*     = 1.0 x 4.18605 x 10^{3} --------------------------- x -------     */
  /*                                      [A]                    [mol]      */
  /*                                                                        */
  /*     => [mol^{-1}] => [ (one atom)^{-1}] => per one atom                */
  /*                                                                        */
  /*                                 10^{18} [energy]          1            */
  /*     => 1.0 x 4.18605 x 10^{3} ------------------ x ----------------    */
  /*                                      [A]            (avogadro cnst)    */ 
  /*                                   10^{18} [energy]                     */
  /*     = 1.0 x 4.18605 x 10^{3} --------------------------                */
  /*                                6.0221367 x 10^{23}  [A]                */
  /*     = 1.0 x (4.18605/6.0221367) x 10^{-2}  [energy][A^{-1}]            */
  /*------------------------------------------------------------------------*/
  sx1.f0 = 1.0 * 4.18605/(6.0221367 * 100.0);
}

void mk_table(void)   /* make a look up table for SX1 [NaCl] potential */
{
  int ion_i, ion_j, ddr;
  int ij = 0;
  double dr, adr, erfc_adr_per_dr, ZZ, a, b;
  double p_e1, p_r1, f_e1, f_r1;

  for(ion_i=0;ion_i<ctl.kinds_of_ions;ion_i++) {
    for(ion_j=ion_i;ion_j<ctl.kinds_of_ions;ion_j++) {

      a = sx1.alpha[ion_i] + sx1.alpha[ion_j];
      b = sx1.beta[ion_i] + sx1.beta[ion_j];
      ZZ = ion.Z[ion_i] * ion.Z[ion_j];

      for(ddr=0;ddr<sys.table_division+1;ddr++) {  /* division = 0.001 [A] */

	dr = ((double)ddr)/1000.0+0.5;  /* 0.5 [A] to sys.radius (cut-off) */
	adr = sys.a1 * dr;  
	erfc_adr_per_dr = erfcc(adr)/dr;

	/* Potential */
	p_e1 = ZZ * erfc_adr_per_dr;
	p_r1 = sx1.f0 * b * exp((a - dr) / b);  /* SX1 specific part */

	/* Force */
	f_e1 = ZZ* (2.0 * sys.a3 * exp(- adr * adr) +
		    erfc_adr_per_dr )/(dr*dr);
	f_r1 = sx1.f0*exp((a -dr) / b)/dr;      /* SX1 specific part */

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



