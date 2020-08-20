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
  if((fpmac = fopen("../files/MACMOL","w"))==NULL){
    printf("cannot open MACMOL. Abort\n");
    exit(EXIT_FAILURE);
  }

  /* to WRITE a blank " " at first time ! */

  fprintf(fpout," ");  
  fprintf(fpout_positions," ");
  fprintf(fpout_velocities," ");
  fprintf(fpmac," ");  
}

/*------------ close files ------------*/
void   close_files(void)
{
  fclose(fpout);
  fclose(fpout_positions);
  fclose(fpout_velocities);
  fclose(fpmac);
}

/*------------ print data to the file(s) ------------*/
void   print_to_file(void)
{
  int atom_kind;

  fprintf(fpout,"%d  %.2f %.2f %.2f  %.2f %.2f  %.3f %.3f %.3f ",
	  sys.step,
	  sys.kin*sys.perMol,
	  sys.pot*sys.perMol,
	  (sys.kin+sys.pot)*sys.perMol,
	  sys.kin*sys.e2t,
	  sys.pres*sys.pp2gpa,
	  sys.Lx/sys.nx,sys.Ly/sys.ny,sys.Lz/sys.nz);

  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    fprintf(fpout," %.3f ",msd.value[atom_kind]);
  } 
  fprintf(fpout,"\n");

  if(sys.step % 10 == 1){ 
    fflush(fpout);
  }

}

/*------------ print data to display FORM1 ------------*/
void   display1(void)
{
  int atom_kind;

  if(sys.step % 10 == 1){
    printf("  -------------------------------------------------------------------------\n");
    printf("  step: Ek  :  Up   :  Et   :  T[K] : P[Gpa]  ( Ax   Ay   Az )  MSD [A^2]\n");
    printf("  -------------------------------------------------------------------------\n");
  };

  printf("%6u % 5.2f % 5.2f % 5.2f % 7.2f % 6.3f  % 3.2f % 3.2f % 3.2f  ",
	 sys.step,
	 sys.kin*sys.perMol/1000.0,
	 sys.pot*sys.perMol/1000.0,
	 (sys.kin+sys.pot)*sys.perMol/1000.0,
	 sys.kin*sys.e2t,
	 sys.pres*sys.pp2gpa,
	 sys.Lx/sys.nx,sys.Ly/sys.ny, sys.Lz/sys.nz);

  /* MSD output */
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    printf("%5.3f ",msd.value[atom_kind]);
  }
  printf("\n");
}

/*------------ print data to display FORM2 ------------*/
void   display2(void)
{
  int atom_kind;

  if(sys.step % 10 == 1){
    printf("------------------------------------------------------------------------\n");
    printf("     step :   Ek  :    Up   :    Et   :  T[K] :P[GPa] : ( Ax   Ay   Az )  msd[A^2]\n");
    printf("------------------------------------------------------------------------\n");
  };
  printf("%10u %7.1f %9.1f %9.1f %7.1f %7.3f  (%4.1f %4.1f %4.1f)",
	 sys.step,
	 sys.kin*sys.perMol/1000.0,
	 sys.pot*sys.perMol/1000.0,
	 (sys.kin+sys.pot)*sys.perMol/1000.0,
	 sys.kin*sys.e2t,
	 sys.pres*sys.pp2gpa,
	 sys.Lx/sys.nx,sys.Ly/sys.ny, sys.Lz/sys.nz);

  /* MSD output */
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    printf(" %6.2f\n",msd.value[atom_kind]);
  }
}

/*------------ print data to display FORM3 ------------*/
void   no_display(void)  /* display only a numerical value list */
{
  int atom_kind;

  printf("%6u %5.2f %5.2f %5.2f  %6.2f %6.3f   %3.2f %3.2f %3.2f   ",
	 sys.step,
	 sys.kin*sys.perMol/1000.0,
	 sys.pot*sys.perMol/1000.0,
	 (sys.kin+sys.pot)*sys.perMol/1000.0,
	 sys.kin*sys.e2t,
	 sys.pres*sys.pp2gpa,
	 sys.Lx/sys.nx,sys.Ly/sys.ny, sys.Lz/sys.nz);

  /* MSD output */
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    printf("%4.2f ",msd.value[atom_kind]);
  }
  printf("\n");

}

/*--- record position data to file in [A] unit ---*/
void   r_position(void)
{
  int i;
  int cnt=0; 

  for(i=0;i<sys.N;i++) { /* [A] unit */
    fprintf(fpout_positions,"%d %.5f %.5f %.5f ", 
	    sys.ion[i],sys.rx[i],sys.ry[i],sys.rz[i]);
    cnt++;
    if(cnt == 4){
      fprintf(fpout_positions," \n");
      cnt = 0;
    }
  }
  fprintf(fpout_positions,"\n");
}

/*--- record velocity data to file in [A]/[fsec] unit ---*/
void   r_velocities(void)
{
  int i;
  int cnt=0; 
  double vx, vy, vz;

  for(i=0;i<sys.N;i++) {
    vx = sys.vx[i]/SECD_2_FS;  /* velocity :   [A/fsec'] -> [A/fsec]     */
    vy = sys.vy[i]/SECD_2_FS;  /* unit     :   1 [fsec] = 10^(-15) [sec] */
    vz = sys.vz[i]/SECD_2_FS;

    fprintf(fpout_velocities,"%d %.8e %.8e %.8e ", 
	    sys.ion[i], vx, vy, vz); /* [A/fsec] */
    cnt++;
    if(cnt == 4){
      fprintf(fpout_velocities," \n");
      cnt = 0;
    }
  }
  fprintf(fpout_velocities,"\n");
}

/*--- record position data to file in MACMOL format ---*/
void   r_macmol(void)
{
  int i;
  int cnt=0; 

  fprintf(fpmac,"MD-STEP=%d \n",sys.step);
  fprintf(fpmac,"%d \n",sys.N);
  for(i=0;i<sys.N;i++) { /* [A] unit */
    fprintf(fpmac,"ION[%d]    %5d    %3.6f   %3.6f   %3.6f \n", 
	    sys.ion[i],i+1, sys.rx[i],sys.ry[i],sys.rz[i]);
  }
  fprintf(fpmac,"\n");
}

