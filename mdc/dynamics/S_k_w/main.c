#include "sk.h"
#include "prototypes.h"

int main(void)
{
  init_param();  /* set the calculational control parameters [control.c] */
  init_mem();    /* dynamic memory allocation [control.c] */ 
  openfiles();   /* file manupulation [files.c] */
  calc_Skw();    /* dynamical structure factor S(k,w) calculation [main.c] */
  record_data(); /* write a calculated S(k,w) to the file "out" [files.c] */
  close_files(); /* file closing operation [files.c] */
  return 0;      /* do nothing */
}

void calc_Skw(void)     /* [Na] -> positive ion  [Cl] -> negative ion */
{
  int i, j, ion_kind, Tstep;
  float  x, y, z;     
  double omega_zero, g2w, gw, g0; 
  double g7_4w, g3_2w, g5_4w, g3_4w, g1_2w, g1_4w;
  double sum=0;

  for(Tstep=0; Tstep<total_time_step; Tstep++){

    if(Tstep%10==0) printf("reading time-step = %d\n",Tstep);

    for(i=0;i<num_of_particles;i++){
      fscanf(fp_positions,"%d %f %f %f", &ion_kind, &x, &y, &z);  
      switch(ion_kind){
      case 0:
	rho_k_c_Na[Tstep] += cos((double)x*hh+(double)y*kk+(double)z*ll);
	rho_k_s_Na[Tstep] += sin((double)x*hh+(double)y*kk+(double)z*ll);
	break;
      case 1:
	rho_k_c_Cl[Tstep] += cos((double)x*hh+(double)y*kk+(double)z*ll);
	rho_k_s_Cl[Tstep] += sin((double)x*hh+(double)y*kk+(double)z*ll);
	break;
      default:
	printf("new ION number detected\n");
	exit(0);
      }
    }
  }

  hanning_window(rho_k_c_Na, total_time_step);
  hanning_window(rho_k_s_Na, total_time_step);
  hanning_window(rho_k_c_Cl, total_time_step);
  hanning_window(rho_k_s_Cl, total_time_step);

  /*-------   F F T  ------------------------------------*/
  fft(rho_k_c_Na, rho_k_c_Na_Im, total_time_step, -1.0);
  fft(rho_k_s_Na, rho_k_s_Na_Im, total_time_step, -1.0);
  fft(rho_k_c_Cl, rho_k_c_Cl_Im, total_time_step, -1.0);
  fft(rho_k_s_Cl, rho_k_s_Cl_Im, total_time_step, -1.0);
  /*-----------------------------------------------------*/
  
  for(Tstep=0; Tstep<total_time_step; Tstep++){

    /*** S(k,w)++ [Na-Na] ***/
    skw_PP[Tstep] = (rho_k_c_Na[Tstep]-rho_k_s_Na_Im[Tstep])*
      (rho_k_c_Na[Tstep]-rho_k_s_Na_Im[Tstep]) + 
	(rho_k_c_Na_Im[Tstep]+rho_k_s_Na[Tstep])*
	  (rho_k_c_Na_Im[Tstep]+rho_k_s_Na[Tstep]);  

    /*** S(k,w)-- [Cl-Cl] ***/
    skw_NN[Tstep] = (rho_k_c_Cl[Tstep]-rho_k_s_Cl_Im[Tstep])*
      (rho_k_c_Cl[Tstep]-rho_k_s_Cl_Im[Tstep]) + 
	(rho_k_c_Cl_Im[Tstep]+rho_k_s_Cl[Tstep])*
	  (rho_k_c_Cl_Im[Tstep]+rho_k_s_Cl[Tstep]);

    /*** S(k,w)+- [Na-Cl] ***/
    skw_PN[Tstep] = (rho_k_c_Na[Tstep]-rho_k_s_Na_Im[Tstep])*
      (rho_k_c_Cl[Tstep]-rho_k_s_Cl_Im[Tstep]) + 
	(rho_k_s_Na[Tstep]+rho_k_c_Na_Im[Tstep])*
	  (rho_k_s_Cl[Tstep]+rho_k_c_Cl_Im[Tstep]);

    /*===  S(q,w) evaluation  ============================================*/
    /* S(q,w) = Sz(Q,w) [charge weighting] */
    skw[Tstep] = (Z_ion0*Z_ion0*num_of_ion0*skw_PP[Tstep] + 
		  Z_ion1*Z_ion1*num_of_ion1*skw_NN[Tstep] + 
		  2.0*Z_ion0*Z_ion1*sqrt(num_of_ion0*num_of_ion1)*
		  skw_PN[Tstep])/num_of_particles; 

    /* S(q,w) = Sm(Q,w) [mass weighting] */
/*
    skw[Tstep] = (M_ion0*M_ion0*num_of_ion0*skw_PP[Tstep] + 
		  M_ion1*M_ion1*num_of_ion1*skw_NN[Tstep] + 
		  2.0*M_ion0*M_ion1*sqrt(num_of_ion0*num_of_ion1)*
		  skw_PN[Tstep])/num_of_particles; 
*/

    /* S(q,w) = Sn(Q,w) [nertron scattering length weighting] */
/*
       skw[Tstep] = B_ion0*B_ion0*skw_PP[Tstep] + 
       B_ion1*B_ion1*skw_NN[Tstep] + 2.0*B_ion0*B_ion1*skw_PN[Tstep];       
*/

    /*====================================================================*/
  }

  /*---  Gaussian window filtering for noise reduction -------------------*/
  omega_zero = 3.0*delta_w;  /* omega-0 for Gauss window  */
  g2w    = (2.0/omega_zero)*0.3095493881*exp(-1.204119983*4.0);
  g7_4w  = (2.0/omega_zero)*0.3095493881*exp(-1.204119983*49.0/16.0);
  g3_2w  = (2.0/omega_zero)*0.3095493881*exp(-1.204119983*9.0/4.0);
  g5_4w  = (2.0/omega_zero)*0.3095493881*exp(-1.204119983*25.0/16.0);
  gw     = (2.0/omega_zero)*0.3095493881*exp(-1.204119983);
  g3_4w  = (2.0/omega_zero)*0.3095493881*exp(-1.204119983*9.0/16.0);
  g1_2w  = (2.0/omega_zero)*0.3095493881*exp(-1.204119983/4.0);
  g1_4w  = (2.0/omega_zero)*0.3095493881*exp(-1.204119983/16.0);
  g0     = (2.0/omega_zero)*0.3095493881;

  printf("skw[0] = %f\n",skw[0]);
  printf("skw[1] = %f\n",skw[1]);
  printf("skw[2] = %f\n",skw[2]);
  printf("skw[total_time_step-1] = %f\n",skw[total_time_step-1]);
  printf("skw[total_time_step-2] = %f\n",skw[total_time_step-2]);
  printf("skw[total_time_step-3] = %f\n",skw[total_time_step-3]);

  for(i=0;i<total_time_step;i++){

    j = i;
    if(i<8) j = (total_time_step-1)-i; /* cyclic BC in FFT */
    if(i>total_time_step-9) j = i-(total_time_step-9);

    gskw[j] = skw[j-8]*g2w+skw[j-7]*g7_4w+
      skw[j-6]*g3_2w+skw[j-5]*g5_4w+skw[j-4]*gw+
	skw[j-3]*g3_4w+skw[j-2]*g1_2w+skw[j-1]*g1_4w+
	  skw[j]*g0+skw[j+1]*g1_4w+skw[j+2]*g1_2w+
	    skw[j+3]*g3_4w+skw[j+4]*gw+
	      skw[j+5]*g5_4w+skw[j+6]*g3_2w+
		skw[j+7]*g7_4w+skw[j+8]*g2w;
    if(j!=0) sum +=gskw[j];  /* sum up the normalize factor except 0 point */
  }

  for(i=0;i<total_time_step;i++){ /* normalize */
    gskw[i] /= sum;
  }  
}

