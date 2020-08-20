#include "sk.h"
#include "prototypes.h"

/*-------------------   Open Files  ----------------------*/
void openfiles(void)
{
  if((fp_positions = fopen("../../files/positions","r"))==NULL){
    printf("Cannot open 'positions'. Abort\n");
    exit(EXIT_FAILURE);
  }
  if((fpout = fopen("out","w"))==NULL){
    printf("Can't open out. Abort\n");
    exit(EXIT_FAILURE);
  }
}

/*--------------------  Close Files  ---------------------*/
void close_files(void)
{
  fclose(fp_positions);
  fclose(fpout);
}

/*--------------------  Record Data ----------------------*/
void record_data(void)
{
  int i;  /* X-axis = omega  [rad]/[ps] */

  for(i=0;i<total_time_step;i++){
    fprintf(fpout,"%f   %f\n",
	    (float)delta_w*((float)(i)), (float)gskw[i]);
  }
}





