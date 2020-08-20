#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  void	real_space(void)

  The function for potential and force calculation of
  <Ewald-first-term> + <repulsion> terms using the look-up table
  created by mk_table() in control.c.

*****/
void	real_space(void)
{
  int i, j, ij, ddr;
  double dr;
  double pij, fij;
  double dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  sys.pot = 0.0;
  sys.virX = 0.0;
  sys.virY = 0.0;
  sys.virZ = 0.0;

  for(i=0; i<sys.N; i++) {
    rxi = sys.rx[i];
    ryi = sys.ry[i];
    rzi = sys.rz[i];

    for(j=i+1; j<sys.N; j++) {
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

      /* dr' = dr^2 */
      dr = dx*dx + dy*dy + dz*dz;

      /* sys.radius = cut-off radius */
      /* sys.radius2 = (cut-off radius)^2 */
      if(dr > sys.radius2) continue; 

      /* shift the look-up table column */
      ij = sys.lookup[sys.ion[i]*(ctl.kinds_of_ions) + sys.ion[j]];

      /* dr' -> sqrt(dr) */
      ddr = (int)((sqrt(dr)-0.5)*1000.0);

      /* look-up the table */
      pij = *(sys.pe1r1+ij*(sys.table_division+1)+ddr);
      fij = *(sys.fe1r1+ij*(sys.table_division+1)+ddr);

      sys.pot += pij;

      sys.fx[i] += fij * dx;
      sys.fy[i] += fij * dy;
      sys.fz[i] += fij * dz;
      sys.fx[j] -= fij * dx;
      sys.fy[j] -= fij * dy;
      sys.fz[j] -= fij * dz;
      /* (note): fij' => -Fij/dr, fijx => -fij'*dx => -Fij*dx/dr  */

      /* VIRIAL : (repulsion) + (EWALD 1st term) */
      sys.virX += fij * dx * dx; 
      sys.virY += fij * dy * dy; 
      sys.virZ += fij * dz * dz; 
    }	     
  }
}
