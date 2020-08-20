#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/*  NaCl Born-Mayer-Huggins potential setting                        */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*
*  [note] cutoff length of the real space force is taken as the radius
*  of the sphere in MD super cell box  in "set_cutoff()".
*****/
void get_control_param(void)
{
  ctl.skip_ewald      = 0;     /* 0 => calc Ewald, 1 => skip Ewald         */

  ctl.calc_max        = 100000; /* maximum MD time step                     */
  ctl.delta_time_fs   = 2.0;   /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  pt.temp_set_initial = 600.0; /* initial temperature setting [K]          */
  pt.temp_set         = 600.0; /* temperature setting [K]                  */
  pt.t_control_step   = 10;     /* scale at every ctl.t_control_step steps  */

  pt.pres_set[0][0]   = 10.0;  /* Pxx : Pressure setting [GPa]             */
  pt.pres_set[1][1]   = 10.0;  /* Pyy : the unit will be chenged to [GPa'] */
  pt.pres_set[2][2]   = 10.0;  /* Pzz : "unit_converter(void)" of "init.c" */
  pt.p_control_step   =  20;
  pt.pres_scalling_factor = 0.5;

  cell.hu[0][0] = 5.63;          /* unit cell matrix H (3x3) = {a,b,c} [A] */
  cell.hu[1][1] = 5.63;              /* initial setting = cubic unit cell  */
  cell.hu[2][2] = 5.63;
  /*
  cell.hu[0][0] = 4.4;
  cell.hu[1][1] = 4.4;
  cell.hu[2][2] = 4.4;
  */
  cell.nx = 3;	                 /* number of unit cells in x-direction    */
  cell.ny = 3;
  cell.nz = 3;

  ctl.natoms_in_unit_cell = 8;   /* number of atoms in an unit cell        */
  ctl.natoms_in_mol_unit  = 2;   /* number of atoms in primitive Mol unit  */
  ctl.kinds_of_ions       = 2;   /* Na or Cl */

  sys.hm     = 27;      /* = |n^2| : cutoff for the reciprocal lattice vector*/
  sys.radius = 20.0;    /* cutoff radious [A] for the potential lookup table */
}

void	identify_ion(void)
{
  int	x, y, z, i=0;

  for(z=0; z<cell.nz*2; z++)
    for(y=0; y<cell.ny*2; y++)
      for(x=0; x<cell.nx*2; x++)
	sys.ion[i++] = (x%2 + y%2 + z%2)%2;  /* [Na=0 or Cl=1] */
}

/*****
  void set_loc(void)
  [NaCl] Set-up the Na and Cl particle locations on NaCl lattice.
*****/
void	set_loc(void)
{
  int	x, y, z, i=0;
  int  nx, ny, nz;
  double nnx, nny, nnz;
  double shiftx, shifty, shiftz;

  nx  = cell.nx * 2;  ny  = cell.ny * 2;   nz  = cell.nz * 2;
  nnx = (double)nx;   nny = (double)ny;    nnz = (double)nz; 

  shiftx = (cell.hm[0][0]/nnx)/2.0;
  shifty = (cell.hm[1][1]/nny)/2.0;
  shiftz = (cell.hm[2][2]/nnz)/2.0;

  for(z=0; z<nz; z++)
    for(y=0; y<ny; y++)
      for(x=0; x<nx; x++) {
	cell.rx[i] = cell.hm[0][0]*(double)(x/nnx)+shiftx;
	cell.ry[i] = cell.hm[1][1]*(double)(y/nny)+shifty;
	cell.rz[i] = cell.hm[2][2]*(double)(z/nnz)+shiftz;
	i++;
      }
}

/*----------------------- Born-Mayer-Huggins potential [NaCl] ------------- */
void set_potential(void)
{
  double m_na, m_cl, z_na, z_cl;

  /* --------------------------------------------------------------- */
  /*   Huggins-Mayer potential with Tosi-Fumi parameter set          */
  /*   see references for more detail                                */
  /* --------------------------------------------------------------- */
  /*                                                                 */
  /*  Function Form :                                                */
  /* u_(i,j) = Z_(i)*Z_(j)*e^(2)/r_(i,j) +                           */
  /*     c_(i,j)*b*Exp((sigma_(i) + sigma_(j) - r_(i,j))/rho)        */
  /*         - ( C_(i,j)/(r_(i,j)^(6)) ) - ( D_(i,j)/(r_(i,j)^(8)) ) */
  /*                                                                 */
  /*  Parameter Table :                                              */
  /*    b = 3.38e-13 [erg]                                           */
  /*                                                                 */
  /*    c_na_na = 1.25                                               */
  /*    c_na_cl = 1                                                  */
  /*    c_cl_cl = 0.75                                               */
  /*                                                                 */
  /*    sigma_na = 1.170 [A]                                         */
  /*    sigma_cl = 1.585 [A]                                         */
  /*    rho      = 0.317 [A]                                         */
  /*                                                                 */
  /*    C_na_na  = 1.68e-60  [erg][cm^6]                             */
  /*    C_na_cl  = 11.2e-60  [erg][cm^6]                             */
  /*    C_cl_cl  = 116.0e-60 [erg][cm^6]                             */
  /*                                                                 */
  /*    D_na_na  = 0.8e-76   [erg][cm^8]                             */
  /*    D_na_cl  = 13.9e-76  [erg][cm^8]                             */
  /*    D_cl_cl  = 233.0e-76 [erg][cm^8]                             */
  /*                                                                 */
  /*                 |                                               */
  /*                 V                                               */
  /*        UNIT  CONVERTED  [erg]->[energy],[cm]->[A]               */
  /*                 |                                               */
  /*                 V                                               */
  /*                                                                 */
  /*----------- UNIT CONVERSION  ------------------------------------*/
  /*  1[erg] = 1.0e-7 [J]                                            */
  /*                                                                 */
  /*  b = 3.38e-13 [erg]                                             */
  /*                                                                 */
  /*    =  3.383e-13 * 1.0e-7 [J]                                    */
  /*                                                                 */
  /*       3.38e-13 * 1.0e-07 [J]              [energy]              */
  /*    = ------------------------  *   1.0e18 ----------            */
  /*                1                           [ J ]                */
  /*    = 3.38e-2 [energy]                                           */
  /*                                                                 */
  /*-----------------------------------------------------------------*/
  /*                                                                 */
  /*    b = 3.38e-2 [energy]                                         */
  /*                                                                 */
  /*    c_na_na = 1.25                                               */
  /*    c_na_cl = 1                                                  */
  /*    c_cl_cl = 0.75                                               */
  /*                                                                 */
  /*    sigma_na = 1.170  [A]                                        */
  /*    sigma_cl = 1.585  [A]                                        */
  /*    rho      = 0.317  [A]                                        */
  /*                                                                 */
  /*    C_na_na  = 1.68e-01  [energy][A]                             */
  /*    C_na_cl  = 11.2e-01  [energy][A]                             */
  /*    C_cl_cl  = 116.0e-01 [energy][A]                             */
  /*                                                                 */
  /*    D_na_na  = 0.8e-01   [energy][A]                             */
  /*    D_na_cl  = 13.9e-01  [energy][A]                             */
  /*    D_cl_cl  = 233.0e-01 [energy][A]                             */
  /*                                                                 */
  /* ----------------------------------------------------------------*/

  z_na =  1.0 * sqrt(sys.kk);
  z_cl = -1.0 * sqrt(sys.kk);
  m_na = 22.9898 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */
  m_cl = 35.4528 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.Z[0] = z_na;
  ion.m[0] = m_na;
  ion.Z[1] = z_cl;
  ion.m[1] = m_cl;
  
  hm.b = 3.38/100.0;
  hm.c_na_na = 1.25;
  hm.c_na_cl = 1.0;
  hm.c_cl_cl = 0.75;
  hm.sigma_na = 1.170;
  hm.sigma_cl = 1.585;
  hm.rho      = 0.317;
  hm.C_na_na  = 1.68/10.0;
  hm.C_na_cl  = 11.2/10.0;
  hm.C_cl_cl  = 116.0/10.0;
  hm.D_na_na  = 0.8/10.0;
  hm.D_na_cl  = 13.9/10.0;
  hm.D_cl_cl  = 233.0/10.0; 
}

void mk_table(void)   /* make a look up table for BMH potential */
{
  int ion_i, ion_j, ddr;
  int ij = 0, iijj;
  double dr, dr6, dr8;
  double hm_c_ij, hm_sigma_ij, hm_C_ij, hm_D_ij;
  double p_r1, f_r1;

  for(ion_i=0;ion_i<ctl.kinds_of_ions;ion_i++) {
    for(ion_j=ion_i;ion_j<ctl.kinds_of_ions;ion_j++) {

      iijj = ion_i + ion_j;

      switch (iijj) {
      case 0:
        hm_c_ij = hm.c_na_na;
        hm_sigma_ij = hm.sigma_na + hm.sigma_na;
        hm_C_ij = hm.C_na_na;
        hm_D_ij = hm.D_na_na;
        break;
      case 1:
        hm_c_ij = hm.c_na_cl;
        hm_sigma_ij = hm.sigma_na + hm.sigma_cl;
        hm_C_ij = hm.C_na_cl;
        hm_D_ij = hm.D_na_cl;
        break;
      default:
        hm_c_ij = hm.c_cl_cl;
        hm_sigma_ij = hm.sigma_cl + hm.sigma_cl;
        hm_C_ij = hm.C_cl_cl;
        hm_D_ij = hm.D_cl_cl;
        break;
      }

      for(ddr=0;ddr<sys.table_division+1;ddr++) {  /* division = 0.001 [A] */

	dr = ((double)ddr)/1000.0 + 0.5; /* 0.5 [A] to sys.radius (cut-off) */
	dr6 = dr*dr*dr*dr*dr*dr;
	dr8 = dr6*dr*dr;

	/* Potential */
	p_r1 = hm_c_ij*hm.b*exp((hm_sigma_ij - dr)/hm.rho) -
	  hm_C_ij/dr6 - hm_D_ij/dr8;

	/* Force */
	f_r1 = (hm_c_ij*hm.b*exp((hm_sigma_ij - dr)/hm.rho)/(hm.rho) -
		6.0*hm_C_ij/(dr6*dr) - 8.0*hm_D_ij/(dr8*dr))/dr;

	sys.pe1r1[ij*(sys.table_division+1)+ddr] = p_r1;
	sys.fe1r1[ij*(sys.table_division+1)+ddr] = f_r1;
      }

      /* fill the lookup table index */
      sys.lookup[ctl.kinds_of_ions * ion_i + ion_j] = ij;
      sys.lookup[ctl.kinds_of_ions * ion_j + ion_i] = ij;
      ij++;  /* shift the look-up table column */
    }
  }
}




