#include  <stdlib.h>
#include  <stdio.h>
#include  "md.h"
#include  "prototypes.h"

/*------------ open files ------------*/
void   open_files(void)
{
  if((fpout = fopen("../files/out","w"))==NULL){
    printf("cannot open out. Abort\n");
    exit(EXIT_FAILURE);
  }
  if((fpout_positions = fopen("../files/positions","w"))==NULL){
    printf("cannot open fpout_positions. Abort\n");
    exit(EXIT_FAILURE);
  }
  if((fpout_velocities = fopen("../files/velocities","w"))==NULL){
    printf("cannot open fpout_velocities. Abort\n");
    exit(EXIT_FAILURE);
  }

  /* to WRITE a blank " " at first time            */
  /* especially for DEC (digital-UNIX) system ;-). */

  fprintf(fpout," ");  
  fprintf(fpout_positions," ");
  fprintf(fpout_velocities," ");
}

/*------------ close files ------------*/
void   close_files(void)
{
  fclose(fpout);
  fclose(fpout_positions);
  fclose(fpout_velocities);
}

/*------------ print data to the file(s) ------------*/
void   print_to_file(void)
{
  int atom_kind;
  
  fprintf(fpout,"%d  %f %f %f  %.2f %.2f  %.3f %.3f %.3f ",
	  sys.step,
	  sys.kin/(double)sys.N/htb.eV2E,
	  sys.pot/(double)sys.N/htb.eV2E,
	  (sys.kin+sys.pot)/(double)sys.N/htb.eV2E,
	  sys.kin*sys.e2t,
	  sys.pres*sys.pp2gpa,
	  sys.Lx/sys.nx,sys.Ly/sys.ny, sys.Lz/sys.nz);
  
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    fprintf(fpout," %f ",msd.value[atom_kind]);
  } 
  fprintf(fpout,"\n");
  if(sys.step % 10 == 1){fflush(fpout);}
}

/*------------ print data to display FORM4 ------------*/
void   display4(void)
{
  int atom_kind;

  if(sys.step % 10 == 1){
    printf("  -----------------------------------------------------------------\
-------\n");
    printf("  step:  Ek  : Up[eV] : Et[eV] : T[K] : P[Gpa]  ( Ax   Ay   Az )  MSD \
[A^2]\n");
    printf("  -----------------------------------------------------------------\
-------\n");
  };

  printf("%6u %6.2f  %6.2f  %6.2f   %6.2f  %6.3f   %3.2f %3.2f %3.2f   ",
         sys.step,
         sys.kin/(double)sys.N/htb.eV2E,
         sys.pot/(double)sys.N/htb.eV2E,
         (sys.kin+sys.pot)/(double)sys.N/htb.eV2E,
         sys.kin*sys.e2t,
         sys.pres*sys.pp2gpa,
         sys.Lx/sys.nx,sys.Ly/sys.ny, sys.Lz/sys.nz);
 
  /* MSD output */
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    printf("%5.3f ",msd.value[atom_kind]);
  }
  printf("\n");
}



/*--- record position data to file in [A] unit ---*/
void   r_position(void)
{
  int i;
  int cnt=0; 

  for(i=0;i<sys.N;i++)
    {
      fprintf(fpout_positions,"%d %5.3f %5.3f %5.3f ", 
	      sys.ion[i],sys.rx[i],sys.ry[i],sys.rz[i]);
      cnt++;
      if(cnt == 4){
	fprintf(fpout_positions," \n");
	cnt = 0;
      }
    }
  fprintf(fpout_positions,"\n");
  fflush(fpout_positions);
}

/*--- record velocity data to file in [A]/[fsec] unit ---*/
void   r_velocities(void)
{
  int i;
  int cnt=0; 

  for(i=0;i<sys.N;i++)
    {
      fprintf(fpout_velocities,"%d %5.3f %5.3f %5.3f ", 
	      sys.ion[i],sys.vx[i],sys.vy[i],sys.vz[i]);
      cnt++;
      if(cnt == 4){
	fprintf(fpout_velocities," \n");
	cnt = 0;
      }
    }
  fprintf(fpout_velocities,"\n");
}

