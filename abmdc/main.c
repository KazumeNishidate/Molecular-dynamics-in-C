#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/*---------------------------------------------------*/
/* ab-initio Car-Parrinello Molecular Dynamics code  */
/*---------------------------------------------------*/
int  main(void)
{
  open_files();

  init_atomic_unit();
  set_atoms();
  set_param();

  dump_position_mesa(); /* dump initial position */

  init_bhs_pseudo_potential();
  /* check_bhs_dV_in_real_space(0, 0); */
  
  reciprocal_lattice_vector_generator();
  mkindex_fft3d();

  /* show_vectors(wv.nplw); */

  /*--- set Hamiltonian matrix with kinetic and pseudopotential ---*/
  clear_hmat_small();
  mk_Bessel(wv.nplw, wv.kgx, wv.kgy, wv.kgz);
  set_Vps_hmat(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);
  set_Ek(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);

  get_eigen();

  set_fermi();  /*--- set occupation number [under development] ---*/

  /*--- add hartree and correlation and exchange potential ---*/
  clear_hmat_small();
  set_Vps_hmat(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);
  set_Ek(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);

  get_chgdns();
  dump_chgdns();
  set_Vxc(wv.nplwin);
  set_Hartree(wv.nplwin);
  set_HVxc_hmat(wv.nplwin, wv.ngx, wv.ngy, wv.ngz);

  get_eigen();
  /*-------------------------------------------------------------*/
  
  init_plane_wave_coefficients();  

  init_ewald();
  ewald();

  Born_Oppenheimer();  /* approaching to the BO potential surface        */
  set_force();         /* store the F(t=0)                               */

  /*================================================ MD start ===*/
  for(cpa.step=1; cpa.step<=cpa.max_step; cpa.step++){

    ewald();
    set_force();  /* F(t+dt), store F(t) for next_v()  */
    calc_press();

    set_Vps(wv.nplw, wv.kgx, wv.kgy, wv.kgz);     /* H*C(t) with Vps     */
    get_chgdns();                                 /* calc rho            */
    set_Vxc(wv.nplw);                             /* H*C(t) with Vxc     */
    set_Hartree(wv.nplw);                         /* H*C(t) with Hartree */
    set_VhxcC(wv.nplw, wv.ngx, wv.ngy, wv.ngz);   /* Vhxc*C(t)           */

    get_energy();
    display();
    calc_msd();   /* MSD calculation for ions          */

    next_r();     /* Verlet's velocity form -> r(t+dt) */
    next_v();     /* Verlet's velocity form -> v(t+dt) */
    next_wf();    /* get the C(t+dt) */


    /* control temperature with velocity scaling */
    /* ctl_temp(ctl.tstep, ctl.temp); */

    if(cp.residual_max > cpa.residual) {
      Born_Oppenheimer();   /* approaching to the BO potential surface */
      set_force();
    }

    /*--- dump data ---*/
    dump_eigen();
    dump_position();
    dump_energy();
    dump_force();
    dump_velocity();
    dump_temp();
    dump_pres();
    dump_position_mesa();
    dump_msd();
    /*-----------------*/
    
  } /*======================================== End of MD Loop ===*/  

  close_files();  
  return 0;
}
