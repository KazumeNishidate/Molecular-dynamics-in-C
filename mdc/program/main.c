#include  <stdlib.h>
#include  <stdio.h>
#include  "md.h"
#include  "prototypes.h"

int  main(void)
{
  open_files();
  newton();
  close_files();
  return 0;
}

void  newton(void)
{

  init();       /* initialize the system            */

  /* for the estimation of alpha value [optional]     */
  /* calc_alpha(); */
  
  mk_table();   /* make up the table for force and potential */

  /*********************************************************************/
  /*================================================ MD calculation ===*/
  for(sys.step=1; sys.step<=ctl.calc_max; sys.step++) {

    /* time integration [select one of the following 2 methods] */
    next_rv_verlet(); 
    /* next_rv_gear(); */

    calc_kin();  
    calc_press();

    /* record position data in [A] unit "files/positions"         */
    /* r_position(); */
    /* record velocity data in [A]/[fsec] unit "files/velocities" */
    /* r_velocities(); */
    /* record position data in MACMOL format "files/MACMOL" */
    /* if(sys.step==1){r_macmol();}  */
    /* if(sys.step%1000==0){r_macmol();}  */

    /* TTY [terminal] output  ====================================*/
    /* display1(); */
    display2();
    /*    no_display(); */
    /*============================================================*/

    /* MSD. calculation [optional] */
    calc_msd();  

    /* file output              */
    print_to_file();

    /* pressure and temprature control */
    /*    control_press(ctl.p_control_step); */
    if(sys.step<10000){ control_temp(ctl.t_control_step, ctl.temp); };

    /* momentum correction for every 1000 MD steps */
    if(sys.step % 1000 == 0){moment_correction();}; 

    /*-------------- X interface -------------------*/
#ifdef XD 
    if(sys.step % 2 == 0){xd();}; 
#endif
    /*----------------------------------------------*/
  }
  /*========================================= end of MD calculation ===*/
  /*********************************************************************/
}

