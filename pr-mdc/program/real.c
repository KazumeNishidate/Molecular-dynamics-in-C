#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  The function for potential and force calculation of <repulsion> term 
  using the look-up table created by mk_table() in control.c.
*****/
void	real_space(void)
{
  int i, j, ij, ddr;
  double dr;
  double pij, fij;
  double dx, dy, dz;
  double rxi, ryi, rzi;

  sys.pot = 0.0;
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {  
      pt.virial[i][j] = 0.0;
    }
  }

  for(i=0; i<sys.N; i++) {
    rxi = cell.rx[i];
    ryi = cell.ry[i];
    rzi = cell.rz[i];

    for(j=i+1; j<sys.N; j++) {
      dx = rxi - cell.rx[j];
      dy = ryi - cell.ry[j];
      dz = rzi - cell.rz[j];

      cyclic_bc_dr(&dx, &dy, &dz);

      dr = dx*dx + dy*dy + dz*dz; /* dr' = dr^2 */

      /* sys.cutoff2 = (cut-off radius)^2 */
      if(dr > sys.cutoff2) continue; 

      /* shift the look-up table column */
      ij = sys.lookup[sys.ion[i]*(ctl.kinds_of_ions) + sys.ion[j]];

      /* dr' -> sqrt(dr) : 0.001 [A] division */
      ddr = (int)((sqrt(dr)-0.5)*1000.0);

      /* look-up the table */
      pij = *(sys.pe1r1+ij*(sys.table_division+1)+ddr);
      fij = *(sys.fe1r1+ij*(sys.table_division+1)+ddr);

      sys.pot += pij;

      cell.fx[i] += fij * dx;
      cell.fy[i] += fij * dy;
      cell.fz[i] += fij * dz;
      cell.fx[j] -= fij * dx;
      cell.fy[j] -= fij * dy;
      cell.fz[j] -= fij * dz;
      /* (note): fij' => -Fij/dr, fijx => -fij'*dx => -Fij*dx/dr  */

      /* VIRIAL : repulsion */
      pt.virial[0][0] += fij * dx * dx; 
      pt.virial[0][1] += fij * dx * dy; 
      pt.virial[0][2] += fij * dx * dz; 

      pt.virial[1][1] += fij * dy * dy; 
      pt.virial[1][2] += fij * dy * dz; 

      pt.virial[2][2] += fij * dz * dz; 
    }
  }
  pt.virial[1][0] = pt.virial[0][1];
  pt.virial[2][0] = pt.virial[0][2];
  pt.virial[2][1] = pt.virial[1][2];
}

