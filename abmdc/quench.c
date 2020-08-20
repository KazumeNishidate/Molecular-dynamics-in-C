#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

void Born_Oppenheimer(void)  /* approach to the BO potential surface */
{
  int i, j;

  /* [Cg] 2 -> 1 :[1->2->3] : 2 -> 1 */ 

  for(i=0;i<atom.nvband_all;i++)  /* C(t-dt) = C(t) : in quench loop */
    for(j=0;j<wv.nplw;j++) {
      wv.cg1r[i][j] = wv.cg2r[i][j];
      wv.cg1i[i][j] = wv.cg2i[i][j];
    }

  display_bo_start();

  /*================= approach to the BO surface ========================*/
  for(cpq.step=1; cpq.step<=cpq.max_qstep; cpq.step++){

    for(i=0;i<atom.nvband_all;i++) { /* store the current wave functions */
      for(j=0;j<wv.nplw;j++) {
	wv.cg3r[i][j] = wv.cg2r[i][j];
	wv.cg3i[i][j] = wv.cg2i[i][j];
      }
    }

    /*-------------------------------------------------------------------*/
    set_Vps(wv.nplw, wv.kgx, wv.kgy, wv.kgz);     /* H*C(t) with Vps     */
    get_chgdns();                                 /* calc rho            */
    set_Vxc(wv.nplw);                             /* H*C(t) with Vxc     */
    set_Hartree(wv.nplw);                         /* H*C(t) with Hartree */
    set_VhxcC(wv.nplw, wv.ngx, wv.ngy, wv.ngz);   /* Vhxc*C(t)           */

    mk_HxC();         /* making up the H*C                               */
                      /* -> wv.rhcr wv.rhci                              */

    set_Legendre();   /* setting up the Legendre Multiplier              */ 
                      /* wv.rhcr wv.cg2r wv.cg2i -> cp.LI cp.LR          */

    quench_wv();      /* scaling the wave function                       */
                      /* calculate the residual                          */

    sort_eigen();     /* sort the eigenvalues in reverse order           */
                      /* preparation for the Gram_Schmidt                */

    Gram_Schmidt();   /* Gram_Schmidt orthogonalization                  */
    /*-------------------------------------------------------------------*/

    /* orthogonality_check(); */
    /* show_eigen();          */

    for(i=0;i<atom.nvband_all;i++) { /* recover the wave functions */
      for(j=0;j<wv.nplw;j++) {
	wv.cg1r[i][j] = wv.cg3r[i][j];
	wv.cg1i[i][j] = wv.cg3i[i][j];
      }
    }

    get_energy();    
    dump_energy();
    display();  

    if(cp.residual_max < cpq.residual) break;  /* exit this loop */
  }
  /*=====================================================================*/

  for(i=0;i<atom.nvband_all;i++) {    /* initialize the wave function */ 
    for(j=0;j<wv.nplw;j++) {
      wv.cg1r[i][j] = wv.cg2r[i][j];
      wv.cg1i[i][j] = wv.cg2i[i][j];
    }
  }
  display_bo_end();
}

void  quench_wv(void)
{
  int i, j;
  double quench, zr, zi;

  quench = 1.0; /* setting up the quench rate */
  if(cpq.step%cpq.qstep_interval==0){ quench = cpq.quench; }

  for(i=0;i<atom.nvband_all;i++) {
    cp.residual[i] = 0.0;
    for(j=0;j<wv.nplw;j++) {
      zr = -wv.rhcr[i][j] + cp.Lr[i]*wv.cg2r[i][j] - cp.Li[i]*wv.cg2i[i][j];
      zi = -wv.rhci[i][j] + cp.Lr[i]*wv.cg2i[i][j] + cp.Li[i]*wv.cg2r[i][j];

      wv.cg2r[i][j] += quench*(wv.cg2r[i][j]-wv.cg1r[i][j]) + cpq.dt2myu*zr;
      wv.cg2i[i][j] += quench*(wv.cg2i[i][j]-wv.cg1i[i][j]) + cpq.dt2myu*zi;

      cp.residual[i] += zr*zr + zi*zi;
    }
  }

  cp.residual_max = cp.residual[0];
  for(i=1;i<atom.nvband_all;i++) {
    if(cp.residual[i] > cp.residual_max){ cp.residual_max = cp.residual[i]; }
  }

}

void  set_Legendre()  /* setting up the Legendre Multiplier */
{                     /* for the orthogonalized basis set   */
  int i, j;

  for(i=0;i<atom.nvband_all;i++) {
    cp.Lr[i] = 0.0;
    cp.Li[i] = 0.0;
    for(j=(wv.nplw-1);j>=0;j--) { /* fill up only the diagonal elements */
      cp.Lr[i] += wv.cg2r[i][j]*wv.rhcr[i][j] + wv.cg2i[i][j]*wv.rhci[i][j];
      cp.Li[i] += wv.cg2r[i][j]*wv.rhci[i][j] - wv.cg2i[i][j]*wv.rhcr[i][j];
    }
    cp.egvl[i] = cp.Lr[i]; /* store the eigen values */
  }
  
}

void mk_HxC(void)
{
  int i, j;

  for(i=0;i<atom.nvband_all;i++)
    for(j=0;j<wv.nplw;j++) {
      wv.rhcr[i][j] = 0.5*wv.norm2[j]*wv.cg2r[i][j] 
	+ wv.hcgr[i][j] + wv.hvxr[i][j];
      wv.rhci[i][j] = 0.5*wv.norm2[j]*wv.cg2i[i][j] 
	+ wv.hcgi[i][j] + wv.hvxi[i][j];
    }
}

void sort_eigen()
{
  int i, tmp, cnt=0;

  for(i=0; i<atom.nvband_all; i++){ cp.oder[i] = i; }

  while(cnt==0){
    cnt=1;
    for(i=0; i<(atom.nvband_all-1); i++){
      if( cp.egvl[cp.oder[i]] > cp.egvl[cp.oder[i+1]] ){
	tmp          = cp.oder[i+1];
	cp.oder[i+1] = cp.oder[i];
	cp.oder[i]   = tmp;
	cnt          = 0;
      }
    }
  }

}

void show_eigen()
{
  int i, ne;
  ne = atom.nvband_all; /* number of eigenvalues */

  printf("eigen values\n");
  for(i=0; i<ne; i++){
    printf("[%4d] %20.16f\n",i,cp.egvl[cp.oder[i]]);
  }
}

void Gram_Schmidt(void) /* Gram-Schmidt orthogonalization */
{
  int i, j, k;
  double denominator, ar, ai;
  double vr[MAX_WAVE_VECTOR2], vi[MAX_WAVE_VECTOR2];

  denominator = 1.0/sqrt(inner_productR(cp.oder[0], cp.oder[0]));
  for(i=0;i<wv.nplw;i++) {
    wv.cg2r[cp.oder[0]][i] *= denominator;
    wv.cg2i[cp.oder[0]][i] *= denominator;
  }
  
  for(i=1;i<atom.nvband_all;i++) {
    
    for(k=0;k<wv.nplw;k++) {
      vr[k] = 0.0;
      vi[k] = 0.0;
    }
    
    for(j=0; j<i; j++) {
      ar = inner_productR(cp.oder[j], cp.oder[i]);
      ai = inner_productI(cp.oder[j], cp.oder[i]);
      for(k=0; k<wv.nplw; k++) {
	vr[k] += ar*wv.cg2r[cp.oder[j]][k] - ai*wv.cg2i[cp.oder[j]][k];
	vi[k] += ar*wv.cg2i[cp.oder[j]][k] + ai*wv.cg2r[cp.oder[j]][k];
      }
    }
    
    for(k=0;k<wv.nplw;k++) {
      wv.cg2r[cp.oder[i]][k] -= vr[k];
      wv.cg2i[cp.oder[i]][k] -= vi[k];
    }
    
    denominator = 1.0/sqrt(inner_productR(cp.oder[i], cp.oder[i]));
    for(k=0;k<wv.nplw;k++) {
      wv.cg2r[cp.oder[i]][k] *= denominator;
      wv.cg2i[cp.oder[i]][k] *= denominator;
    }

  }
  
}

double  inner_productR(int ib, int jb)
{
  int i;
  double z;
  
  z = 0.0;
  for(i=0; i<wv.nplw; i++){
    z += wv.cg2r[ib][i]*wv.cg2r[jb][i] - (-wv.cg2i[ib][i])*wv.cg2i[jb][i];
  }
  return z;
}

double  inner_productI(int ib, int jb)
{
  int i;
  double z;
  
  z = 0.0;
  for(i=0; i<wv.nplw; i++){
    z += wv.cg2r[ib][i]*wv.cg2i[jb][i] + (-wv.cg2i[ib][i])*wv.cg2r[jb][i];
  }
  return z;
}
