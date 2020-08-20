#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "headers/md.h"

/*----------------------------------------------------*/
/*               Verlet's velocity form               */
/*----------------------------------------------------*/
void  next_r(void)  
{
  int i;
  double dt2=cpa.dt*cpa.dt;
  double drx, dry, drz;
  double m2;

  for(i=0; i<atom.N; i++){
    m2 = 2.0*atom.m[atom.type[i]];
    drx = cpa.dt*cell.vx[i] + (dt2/m2)*wk.fx[i];
    dry = cpa.dt*cell.vy[i] + (dt2/m2)*wk.fy[i];
    drz = cpa.dt*cell.vz[i] + (dt2/m2)*wk.fz[i];

    cell.rx[i] += drx;
    cell.ry[i] += dry;
    cell.rz[i] += drz;

    /* periodic boundary condition */
    if(cell.rx[i]>cell.Lx) cell.rx[i] -= cell.Lx;
    if(cell.ry[i]>cell.Ly) cell.ry[i] -= cell.Ly;
    if(cell.rz[i]>cell.Lz) cell.rz[i] -= cell.Lz;
    if(cell.rx[i]<0.0) cell.rx[i] += cell.Lx;
    if(cell.ry[i]<0.0) cell.ry[i] += cell.Ly;
    if(cell.rz[i]<0.0) cell.rz[i] += cell.Lz;
  }
}

void  next_v(void)
{
  int i;
  double dvx, dvy, dvz;
  double m2, quench;
  
  quench = 1.0; /* setting up the quench rate */
  if(cpa.step%cpa.qstep_interval==0){ quench = cpa.quench; }

  for(i=0; i<atom.N; i++){
    m2 = 2.0*atom.m[atom.type[i]];
    dvx = (cpa.dt/m2)*(wk.fx[i] + wk.fx0[i]);
    dvy = (cpa.dt/m2)*(wk.fy[i] + wk.fy0[i]);
    dvz = (cpa.dt/m2)*(wk.fz[i] + wk.fz0[i]);

    cell.vx[i] += dvx;
    cell.vy[i] += dvy;
    cell.vz[i] += dvz;

    cell.vx[i] *= quench;
    cell.vy[i] *= quench;
    cell.vz[i] *= quench;
  }
}
/*----------------------------------------------------*/

void  next_wf(void)
{
  int i, j;

  for(i=0;i<atom.nvband_all;i++){   /* store the current wave functions */
    for(j=0;j<wv.nplw;j++) {
      wv.cg3r[i][j] = wv.cg2r[i][j];
      wv.cg3i[i][j] = wv.cg2i[i][j];
    }
  }

  mk_HxC();       /* making up the H*C                                    */
  Ryckaert();     /* orthogonalization with Ryckaert-method               */
  get_next_wf();  /* scaling the wave function and calculate the residual */

  /* orthogonality_check(); */

  for(i=0;i<atom.nvband_all;i++){   /* recover the wave functions */
    for(j=0;j<wv.nplw;j++) {
      wv.cg1r[i][j] = wv.cg3r[i][j];
      wv.cg1i[i][j] = wv.cg3i[i][j];
    }
  }
  
  for(i=0;i<atom.nvband_all;i++){   /* eigen value */
    cp.egvl[i] = cp.X0r[i][i]/cpa.dt2myu;
  }

}

void  get_next_wf(void)
{
  int i, j, k;
  double xcr, xci;
  
  for(i=0;i<atom.nvband_all;i++){
    cp.residual[i] = 0.0;
    for(k=0;k<wv.nplw;k++) {
      xcr = 0.0;
      xci = 0.0;
      for(j=0;j<atom.nvband_all;j++){
	xcr += cp.X0r[i][j]*wv.cg2r[j][k] + cp.X0i[i][j]*wv.cg2i[j][k];
	xci += cp.X0r[i][j]*wv.cg2i[j][k] - cp.X0i[i][j]*wv.cg2r[j][k];
      }
      wv.cg2r[i][k] = wv.ncgr[i][k] + xcr;
      wv.cg2i[i][k] = wv.ncgi[i][k] + xci;
      
      cp.residual[i] += 
	  (-wv.rhcr[i][k] + xcr/cpa.dt2myu)*(-wv.rhcr[i][k] + xcr/cpa.dt2myu)
	+ (-wv.rhci[i][k] + xci/cpa.dt2myu)*(-wv.rhci[i][k] + xci/cpa.dt2myu);
    }
  }  
  
  cp.residual_max = cp.residual[0];
  for(i=1;i<atom.nvband_all;i++){
    if(cp.residual[i] > cp.residual_max){ cp.residual_max = cp.residual[i]; }
  }
  
}

void  Ryckaert(void)
{
  int i, j, k, iterate;
  double residue;
  double Xnr[MAX_BANDS][MAX_BANDS], Xni[MAX_BANDS][MAX_BANDS];
  double  Ir[MAX_BANDS][MAX_BANDS],  Ii[MAX_BANDS][MAX_BANDS];
  double  Ar[MAX_BANDS][MAX_BANDS],  Ai[MAX_BANDS][MAX_BANDS];
  double  Br[MAX_BANDS][MAX_BANDS],  Bi[MAX_BANDS][MAX_BANDS];
  
  /*------------ make delta (i=j -> 1, i!=j -> 0) ------------*/
  for(j=0;j<atom.nvband_all;j++) {
    for(i=0;i<atom.nvband_all;i++) {  
      Ir[i][j] = 0.0;
      Ii[i][j] = 0.0;
    }
    Ir[j][j] = 1.0;
  }
  
  /*--------------- make non orthogonal C~(t+dt) -------------*/
  for(i=0;i<atom.nvband_all;i++) {  
    for(k=0;k<wv.nplw;k++) {
      wv.ncgr[i][k] = 2.0*wv.cg2r[i][k]-wv.cg1r[i][k]-cpa.dt2myu*wv.rhcr[i][k];
      wv.ncgi[i][k] = 2.0*wv.cg2i[i][k]-wv.cg1i[i][k]-cpa.dt2myu*wv.rhci[i][k];
    }
  }
  
  /*---------------------- make A and B ----------------------*/
  for(i=0;i<atom.nvband_all;i++) {  
    for(j=0;j<atom.nvband_all;j++) {     
      Ar[i][j] = 0.0;  Ai[i][j] = 0.0;
      Br[i][j] = 0.0;  Bi[i][j] = 0.0;
      for(k=(wv.nplw-1);k>=0;k--) {
	Ar[i][j] += wv.ncgr[i][k]*wv.ncgr[j][k] + wv.ncgi[i][k]*wv.ncgi[j][k];
	Ai[i][j] += wv.ncgr[i][k]*wv.ncgi[j][k] - wv.ncgi[i][k]*wv.ncgr[j][k];
	Br[i][j] += wv.ncgr[i][k]*wv.cg2r[j][k] + wv.ncgi[i][k]*wv.cg2i[j][k];
	Bi[i][j] += wv.ncgr[i][k]*wv.cg2i[j][k] - wv.ncgi[i][k]*wv.cg2r[j][k];
      }
    }
  }
    
  /*------------- make initial X  -> X[0]=(1-A)/2 ------------*/
  for(i=0;i<atom.nvband_all;i++) {  
    for(j=0;j<atom.nvband_all;j++) {  
      cp.X0r[i][j] = (Ir[i][j] - Ar[i][j])/2.0;
      cp.X0i[i][j] = (Ii[i][j] - Ai[i][j])/2.0;
    }
  }
  
  /*======================= Loop of self-consistent =======================*/
  for(iterate=1; iterate<=100; iterate++){

    for(i=0;i<atom.nvband_all;i++){
      for(j=0;j<atom.nvband_all;j++) {     
	Xnr[i][j] = Ir[i][j] - Ar[i][j];
	Xni[i][j] = Ii[i][j] - Ai[i][j];
	for(k=0;k<atom.nvband_all;k++) {     
	  Xnr[i][j] += (Ir[i][k] - Br[i][k])*cp.X0r[j][k]
	             + (Ii[i][k] - Bi[i][k])*cp.X0i[j][k]
	             + cp.X0r[i][k]*(Ir[j][k] - Br[j][k])
	             + cp.X0i[i][k]*(Ii[j][k] - Bi[j][k])
	             - cp.X0r[i][k]*cp.X0r[j][k]
	             - cp.X0i[i][k]*cp.X0i[j][k];
	  Xni[i][j] +=-(Ir[i][k] - Br[i][k])*cp.X0i[j][k]
	             + (Ii[i][k] - Bi[i][k])*cp.X0r[j][k]
	             - cp.X0r[i][k]*(Ii[j][k] - Bi[j][k])
	             + cp.X0i[i][k]*(Ir[j][k] - Br[j][k])
	             + cp.X0r[i][k]*cp.X0i[j][k]
	             - cp.X0i[i][k]*cp.X0r[j][k];
	}
	Xnr[i][j] *= 0.5;
	Xni[i][j] *= 0.5;
      }
    }

    residue=0.0;
    for(j=0;j<atom.nvband_all;j++){
      for(i=0;i<atom.nvband_all;i++){
	residue += (Xnr[i][j] - cp.X0r[i][j])*(Xnr[i][j] - cp.X0r[i][j])
	         + (Xni[i][j] - cp.X0i[i][j])*(Xni[i][j] - cp.X0i[i][j]);
      }
    }
    residue = sqrt(residue/((double)(atom.nvband_all*atom.nvband_all)));

    for(j=0;j<atom.nvband_all;j++){
      for(i=0;i<atom.nvband_all;i++){
	cp.X0r[i][j] = Xnr[i][j];  /* X(n) -> X(n-1) */
	cp.X0i[i][j] = Xni[i][j]; 
      }
    }
    
    if(residue <= 1.0e-12) break;
  }
  /*=======================================================================*/

  if(iterate>=100){
    printf("self-consistent ERROR (residue = %e)\n",residue);
    exit(0);
  }
  
}

void set_force(void)
{
  int i;

  get_HF_force(wv.nplw, wv.kgx, wv.kgy, wv.kgz);

  for(i=0;i<atom.N;i++){
    wk.fx0[i] = wk.fx[i];
    wk.fy0[i] = wk.fy[i];
    wk.fz0[i] = wk.fz[i];
    
    wk.fx[i] = wk.fccx[i] + wk.fecx[i];
    wk.fy[i] = wk.fccy[i] + wk.fecy[i];
    wk.fz[i] = wk.fccz[i] + wk.fecz[i];
  }

}
