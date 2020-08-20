#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

void  get_energy(void)
{
  int i, ib, nq;
  double xr;

  E.ekin = 0.0;
  E.VpsR = 0.0;
  for(ib=0;ib<atom.nvband_all;ib++) {  /* band by band */
    xr = 0.0;  
    for(nq=0;nq<wv.nplw;nq++) {
      /*---------------------- kinetic ----------------------*/
      E.ekin += atom.fermi[ib]*wv.norm2[nq]
	*(wv.cg2r[ib][nq]*wv.cg2r[ib][nq] + wv.cg2i[ib][nq]*wv.cg2i[ib][nq]);
      
      /*------------------ pseudopotential ------------------*/
      E.VpsR += 2.0*atom.fermi[ib]
	*(wv.cg2r[ib][nq]*wv.hcgr[ib][nq] + wv.cg2i[ib][nq]*wv.hcgi[ib][nq]);
      
      /*----------- Fictitious kinetic energy ---------------*/
      xr +=(wv.cg2r[ib][nq]-wv.cg1r[ib][nq])*(wv.cg2r[ib][nq]-wv.cg1r[ib][nq])
	  +(wv.cg2i[ib][nq]-wv.cg1i[ib][nq])*(wv.cg2i[ib][nq]-wv.cg1i[ib][nq]);
    }
    E.fkin = xr/cpq.dt2myu;
  }
  
  /*-------------- hartree and exchnge-correlation ----------------*/
  E.hart = 0.0;  
  E.xc   = 0.0;
  for(i=0; i<ft.mesh; i++) {
    E.hart += 0.5*ft.rho[i]*hxc.VH[i];
    E.xc   +=     ft.rho[i]*hxc.E[i];
  }
  E.hart *= (cell.vol/(double)ft.mesh);
  E.xc   *= (cell.vol/(double)ft.mesh);
  /*---------------------------------------------------------------*/

  calc_kinetic(); /* calculate the core kinetic energy */

  E.eTotal = E.ekin + E.VpsR + E.hart + E.xc;
  E.cTotal = E.ckin + E.EW;
  
  E.TotalF = E.eTotal + E.cTotal + E.fkin;
  E.Total  = E.eTotal + E.cTotal;

  /*------- temperature conversion --------------------------------*/
  /* Ek = (1/2)SUM[mv^2] = (3/2)N Kb T */
  /* T  = Ek*(2.0/(3.0*N*KB))          */

  /* ionic temperature per atom */
  E.cTemp  = E.ckin*au.H2K/((double)atom.N);

  /* electronic temperature per electron */
  E.eTemp  = ME*E.ekin*au.H2K/((double)atom.nvband_all); 

  /* fictitious temperature per electron */
  E.efTemp = E.fkin*au.H2K/((double)atom.nvband_all);
  /*---------------------------------------------------------------*/

}

void  calc_kinetic(void)
{
  int i;
  double px,py,pz,m2;

  E.ckin = 0.0;
  px=py=pz=0.0;

  for(i=0; i<atom.N; i++){
    m2 = atom.m[atom.type[i]]/2.0;
    px += m2*cell.vx[i]*cell.vx[i];
    py += m2*cell.vy[i]*cell.vy[i];
    pz += m2*cell.vz[i]*cell.vz[i];
  }

  E.ckinx = px;
  E.ckiny = py;
  E.ckinz = pz;

  E.ckin  = E.ckinx + E.ckiny + E.ckinz;

}

