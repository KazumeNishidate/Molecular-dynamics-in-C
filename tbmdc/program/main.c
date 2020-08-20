#include  <stdlib.h>
#include  <stdio.h>
#include  "md.h"
#include  "prototypes.h"

void   xd(void);


void  newton(void)
{
/* see "debug.c" ---------------------------------------------*/
/*   calc_R_vs_E();                                           */
/*------------------------------------------------------------*/

/* initialize the system -------------------------------------*/
  init();
/*------------------------------------------------------------*/

/* momentum correction for the initial velocities ------------*/
  moment_correction();
/*------------------------------------------------------------*/

/*================================================ MD calculation ===*/
  for(sys.step=1; sys.step<=ctl.calc_max; sys.step++) {

/* time integration [select one of the following 2 methods] --*/
    next_rv_verlet();                                       
/*    next_rv_gear();                                         */
/*------------------------------------------------------------*/

    calc_kin();   
    calc_press(); 
    calc_msd();  

/* record data -----------------------------------------------*/
/*    r_position();                                           */
/*    r_velocities();                                         */
/*    print_to_file();                                        */
/*------------------------------------------------------------*/

/* terminal output -------------------------------------------*/
    display4();                                             
/*------------------------------------------------------------*/

/* pressure and temprature control ---------------------------*/
/*    control_press(ctl.p_control_step);                      */
/*    control_temp(ctl.t_control_step, ctl.temp);             */
/*------------------------------------------------------------*/

/* momentum correction for every 1000 MD steps ---------------*/
    if(sys.step % 1000 == 0){moment_correction();};         
/*------------------------------------------------------------*/

/* X interface -----------------------------------------------*/
#ifdef XD
    if(sys.step % 2 == 0){xd();}
#endif
/*------------------------------------------------------------*/
  }
/*========================================= end of MD calculation ===*/
}


int  main(void)
{
  open_files();
  newton();
  close_files();
  return 0;
}

