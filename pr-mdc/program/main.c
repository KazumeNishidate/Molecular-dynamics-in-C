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

  /* for the estimation of alpha value [optional]   */
  /*  calc_alpha(); */
  
  mk_table();  /* make up the table for force and potential */

  /*********************************************************************/
  /*================================================ MD calculation ===*/
  /*                                                                   */
  for(sys.step=1; sys.step<=ctl.calc_max; sys.step++) {

    if(sys.step == 10000) {
      pt.pres_set[0][0]   = 20.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 20.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 20.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }

    if(sys.step == 20000) {
      pt.pres_set[0][0]   = 30.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 30.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 30.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }

    if(sys.step == 30000) {
      pt.pres_set[0][0]   = 40.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 40.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 40.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }

    if(sys.step == 40000) {
      pt.pres_set[0][0]   = 50.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 50.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 50.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }

    if(sys.step == 50000) {
      pt.pres_set[0][0]   = 60.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 60.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 60.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }

    if(sys.step == 60000) {
      pt.pres_set[0][0]   = 70.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 70.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 70.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }

    if(sys.step == 70000) {
      pt.pres_set[0][0]   = 80.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[1][1]   = 80.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
      pt.pres_set[2][2]   = 80.0 *1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );
    }



    /* time integration [select one of the following 2 methods] */
    next_rv_verlet(); 
    /* next_rv_gear(); */

    calc_msd();   /* MSD. calculation */

    /* Pressure and Temperature control --------------------------*/

    control_press_PR();

    /* control_press(pt.p_control_step); */

    control_temp(pt.t_control_step, pt.temp_set);

    /*------------------------------------------------------------*/

    calc_kin();  
    calc_press();

    /* TTY [terminal] output  ------------------------------------*/
    display1();
    /* display2(); */
    /*    no_display(); */
    /*------------------------------------------------------------*/

    /* record MD calculation history -----------------------------*/
    print_to_file();
    r_cell();
    /* r_position(); */
    /* r_velocities(); */
    /*------------------------------------------------------------*/

    /* momentum correction for every 1000 MD steps [optional] */
    /*    if(sys.step % 1000 == 0){moment_correction();};  */

#ifdef XD 
    if(sys.step % 2 == 0){xd();};  /*--- X interface ---*/
#endif

  }
  /*                                                                   */
  /*========================================= end of MD calculation ===*/
  /*********************************************************************/
}


