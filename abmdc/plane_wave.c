#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "headers/md.h"

void  orthogonality_check(void)
{
  int i, j, g;
  double re, im;

  for(i=0;i<atom.nvband_all;i++){  
    for(j=0;j<atom.nvband_all;j++){  
      re = 0.0;
      im = 0.0;
      for(g=0;g<wv.nplw;g++) {
	re += wv.cg2r[i][g]*wv.cg2r[j][g] + wv.cg2i[i][g]*wv.cg2i[j][g];
	im += wv.cg2r[i][g]*wv.cg2i[j][g] - wv.cg2i[i][g]*wv.cg2r[j][g];
      }
      printf("[%d %d]  re = %.15e  im = %.15e\n", i, j, re, im);
      getchar();
    }
  }

}

void  mkindex_fft3d(void)
{
  int i, j, index, cnt=1;
  int nkx, nky, nkz;
  double kx, ky, kz;
  
  for(i=0;i<wv.nplw;i++) {
    nkx = wv.ngx[i];
    nky = wv.ngy[i];
    nkz = wv.ngz[i];

    if(nkx < 0) nkx += FFTPTSX;
    if(nky < 0) nky += FFTPTSY;
    if(nkz < 0) nkz += FFTPTSZ;
    index = 2*(FFTPTSX*FFTPTSY*nkz + FFTPTSX*nky + nkx);
    wv.index[i] = index;

    kx = wv.kgx[i];
    ky = wv.kgy[i];
    kz = wv.kgz[i];
    wv.norm2[i] = kx*kx+ky*ky+kz*kz;

    /* printf("(%3d, %3d, %3d) %5d\n",wv.ngx[i],wv.ngy[i],wv.ngz[i],index); */
    /* printf("%f %f %f %f\n",kx,ky,kz,wv.norm2[i]); */
  }

  /*---------- for mk_Bessel() ----------*/
  for(i=0;i<wv.nplw;i++){
    bhs.id[i]=0;  /* initialize */
  }
  for(i=1;i<wv.nplw;i++){
    if(bhs.id[i]!=0) continue;
    for(j=1;j<wv.nplw;j++){
      if(wv.norm2[j]==wv.norm2[i]) bhs.id[j]=cnt;
    }
    cnt++;
    if(cnt==BESSEL_MAX){
      printf("\nBESSEL_MAX = %d : size.h\n\n",BESSEL_MAX);
      exit(0);
    }
  }
  
  /* check */
  /*
  for(i=0;i<wv.nplw;i++){
    printf("%3d  %10f [%3d]\n",i,wv.norm2[i],bhs.id[i]);
  }exit(0);
  */
  /*-------------------------------------*/

}

void   show_vectors(int npl)
{
  int i;

  printf("generated reciprocal lattice vectors = %d\n",npl);
  printf("[num]    h   k   l  :     Gx         Gy         Gz\n");
  for(i=0;i<npl;i++) {
    printf("[%3d]  %3d %3d %3d   %10.6f %10.6f %10.6f\n",
	   i+1,wv.ngx[i],wv.ngy[i],wv.ngz[i],
	   wv.kgx[i],wv.kgy[i],wv.kgz[i]);
  }
}

void   init_plane_wave_coefficients(void)
{
  /*  int i; */

  /* randomize(); */
  /* normalization */
  /*
  for(i=0;i<atom.nvband_all;i++){
    normalize(&wv.cg2r[i][0], &wv.cg2i[i][0]);
  }
  */
}

void  normalize(double *re, double *im)
{
  int i;
  double norm=0.0, norm_sqrt;
  
  for(i=0;i<wv.nplw;i++) {
    norm += re[i]*re[i]+im[i]*im[i];
  }
  
  norm_sqrt = sqrt(norm);
  if(norm_sqrt==0.0) return;
  
  for(i=0;i<wv.nplw;i++) {
    re[i] /= norm_sqrt;
    im[i] /= norm_sqrt;
  }

}

void  randomize()
{
  int i, j;
  unsigned int seed;

  seed = (unsigned int)time( 0 );
  srand( seed );

  for(i=0;i<atom.nvband_all;i++) { /* initialization : random setting [re] */
    for(j=0;j<wv.nplw;j++) {
      wv.cg1r[i][j] = (rand()/(RAND_MAX + 1.0)) /200.0;
      wv.cg1i[i][j] = (rand()/(RAND_MAX + 1.0)) /200.0;
    }
  }
}

/*****
 * reciprocal lattice vector gererator for the plane wave basis set
 *
 *   condition:   (1/2)|k+g|^2  < E_cut
 *
 *   E_cut = cutoff energy [Hartree]
 *   reciprocal lattice vector "g" = 2Pi {nx/Lx, ny/Ly, nz/Lz}
 *   wave vector "k" = only the gamma point {0,0,0} is considered
 *
 *   total number of plane waves to represent one wave function
 *   (one atomic orbital) is "wv.nplw" (cnt).
 *
 *   # restriction: only for a cubic MD super cell system.
 *
 *   note: "#define MAX_WAVE_VECTORS 3000" must be increased for the 
 *         large scale simulation using many plane waves.
 *****/
void   reciprocal_lattice_vector_generator(void)
{
  int i;
  int cnt = 0, cnti = 1; /* check here */
  int gx, gy, gz;
  int gx_max, gy_max, gz_max;
  double ggx, ggy, ggz, gg2;

  /* temporary array to store the generated plane waves */
  double kgx[MAX_WAVE_VECTORS],  kgxi[MAX_WAVE_VECTORS];
  double kgy[MAX_WAVE_VECTORS],  kgyi[MAX_WAVE_VECTORS];
  double kgz[MAX_WAVE_VECTORS],  kgzi[MAX_WAVE_VECTORS];

  int ngx[MAX_WAVE_VECTORS],  ngxi[MAX_WAVE_VECTORS];
  int ngy[MAX_WAVE_VECTORS],  ngyi[MAX_WAVE_VECTORS];
  int ngz[MAX_WAVE_VECTORS],  ngzi[MAX_WAVE_VECTORS];

  /* maximum radius of the reciprocal lattice sphere space */
  gx_max = (int)(sqrt(2.0*wv.E_cut)*cell.Lx/PI2);
  gy_max = (int)(sqrt(2.0*wv.E_cut)*cell.Ly/PI2);
  gz_max = (int)(sqrt(2.0*wv.E_cut)*cell.Lz/PI2);

  for(gx=-gx_max; gx<=gx_max; gx++) {  /* reciprocal lattice vector "G" */
    for(gy=-gy_max; gy<=gy_max; gy++) {  
      for(gz=-gz_max; gz<=gz_max; gz++) {  

	check_array_size(cnt, cnti);

	ggx = PI2*((double)gx)/cell.Lx;
	ggy = PI2*((double)gy)/cell.Ly;
	ggz = PI2*((double)gz)/cell.Lz;
	gg2 = ggx*ggx + ggy*ggy + ggz*ggz;

	if(gx==0&&gy==0&&gz==0) { /* Gamma point = G[0] */
	  ggx = 0.0;
	  ggy = 0.0;
	  ggz = 0.0;
	  gg2 = ggx*ggx + ggy*ggy + ggz*ggz;

	  kgxi[0] = ggx;   ngxi[0] = 0;
	  kgyi[0] = ggy;   ngyi[0] = 0;
	  kgzi[0] = ggz;   ngzi[0] = 0;
	  continue;
	}

	if( gg2 <= wv.E_cuti*2.0 ) {
	  kgxi[cnti] = ggx;   ngxi[cnti] = gx;
	  kgyi[cnti] = ggy;   ngyi[cnti] = gy;
	  kgzi[cnti] = ggz;   ngzi[cnti] = gz;
	  cnti++;
	}

	if( gg2 > wv.E_cuti*2.0 && gg2 <= wv.E_cut*2.0 ) {
	  kgx[cnt] = ggx;   ngx[cnt] = gx;
	  kgy[cnt] = ggy;   ngy[cnt] = gy;
	  kgz[cnt] = ggz;   ngz[cnt] = gz;
	  cnt++;
	}

      }
    }
  }

  check_array_size(cnt, cnti);

  for(i=0;i<cnti;i++) {
    wv.kgx[i] = kgxi[i];   wv.ngx[i] = ngxi[i];
    wv.kgy[i] = kgyi[i];   wv.ngy[i] = ngyi[i];
    wv.kgz[i] = kgzi[i];   wv.ngz[i] = ngzi[i];
  }
  for(i=0;i<cnt;i++) {
    wv.kgx[cnti+i] = kgx[i];   wv.ngx[cnti+i] = ngx[i];
    wv.kgy[cnti+i] = kgy[i];   wv.ngy[cnti+i] = ngy[i];
    wv.kgz[cnti+i] = kgz[i];   wv.ngz[cnti+i] = ngz[i];
  }

  wv.nplw  = cnt+cnti;    /* total number of reciprocal lattice vectors "G" */
  wv.nplwin = cnti;       /* --- for initial plane wave set to diagonalize  */

  printf("Ecut_0 = %f  plane waves at initial = %d\n",wv.E_cuti,wv.nplwin);
  printf("Ecut   = %f  number of plane waves  = %d\n",wv.E_cut,wv.nplw);
}

void check_array_size(int cnt, int cnti)
{
  if(cnt >= MAX_WAVE_VECTORS) { /* error handling */
    printf("MAX_WAVE_VECTORS %d exceeded : ",MAX_WAVE_VECTORS);
    printf("size.h\n");
    exit(0);
  }
  if(cnti >= MAX_WAVE_VECTOR1) { /* error handling */
    printf("MAX_WAVE_VECTOR1 %d exceeded : ",MAX_WAVE_VECTOR1);
    printf("size.h\n");
    exit(0);
  }
  if(cnt+cnti >= MAX_WAVE_VECTOR2) { /* error handling */
    printf("MAX_WAVE_VECTOR2 %d exceeded : ",MAX_WAVE_VECTOR2);
    printf("size.h\n");
    exit(0);
  }
}
