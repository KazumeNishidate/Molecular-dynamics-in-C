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
  if((fpout_cell = fopen("../files/cell","w"))==NULL){
    printf("cannot open fpout_cell. Abort\n");
    exit(EXIT_FAILURE);
  }

  /* to WRITE a blank " " at first time ! */
  fprintf(fpout," ");  
  fprintf(fpout_positions," ");
  fprintf(fpout_velocities," ");
  fprintf(fpout_cell," ");

}

/*------------ close files ------------*/
void   close_files(void)
{
  fclose(fpout);
  fclose(fpout_positions);
  fclose(fpout_velocities);
  fclose(fpout_cell);
}

/*------------ print data to the file(s) ------------*/
void   print_to_file(void)
{
  int atom_kind;
  double enthalpy;

  fprintf(fpout,"%d  %.2f %.2f %.2f  %.2f %.2f  %.3f %.3f %.3f ",
	  sys.step,
	  sys.kin*sys.perMol,
	  sys.pot*sys.perMol,
	  (sys.kin+sys.pot)*sys.perMol,
	  sys.kin*sys.e2t,
	  pt.pres_trace*sys.pp2gpa,
	  cell.Lx/(double)cell.nx,
	  cell.Ly/(double)cell.ny,
	  cell.Lz/(double)cell.nz);

  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    fprintf(fpout," %.3f ",msd.value[atom_kind]);
  } 

  /* enthalpy = E + PV */
  enthalpy = (sys.kin + sys.pot) + pt.pres_trace*cell.vol;
  fprintf(fpout,"%.8e\n", enthalpy*sys.perMol);

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
    printf("  -----------------------------------------------------------------------------\n");
    printf("  step:  Ek  :  Up   :E[kJ/Mol]:  T[K] : P[Gpa]  ( Ax   Ay   Az )  MSD [A^2]\n");
    printf("  -----------------------------------------------------------------------------\n");
  };

  printf("%6u % 6.2f % 5.2f % 8.2f % 7.2f % 7.3f  % 3.2f % 3.2f % 3.2f  ",
	 sys.step,
	 sys.kin*sys.perMol/1000.0,
	 sys.pot*sys.perMol/1000.0,
	 (sys.kin+sys.pot)*sys.perMol/1000.0,
	 sys.kin*sys.e2t,
	 pt.pres_trace*sys.pp2gpa,
	 cell.Lx/(double)cell.nx,
	 cell.Ly/(double)cell.ny,
	 cell.Lz/(double)cell.nz);

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
    printf("  ------------------------------------------------------------------------\n");
    printf("  step:  Ek  :  Up   :  Et   :  T[K] : P[Gpa] ( Ax   Ay   Az )  MSD [A^2]\n");
    printf("  ------------------------------------------------------------------------\n");
  };
  printf("%6u % 6.2f % 5.2f % 5.2f % 7.2f % 6.3f  % 3.2f % 3.2f % 3.2f\n",
	 sys.step,
	 sys.kin*sys.perMol/1000.0,
	 sys.pot*sys.perMol/1000.0,
	 (sys.kin+sys.pot)*sys.perMol/1000.0,
	 sys.kin*sys.e2t,
	 pt.pres_trace*sys.pp2gpa,
	 cell.Lx/(double)cell.nx,
	 cell.Ly/(double)cell.ny,
	 cell.Lz/(double)cell.nz);

  /* MSD output */
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    printf(" %10.5f\n",msd.value[atom_kind]);
  }
}

/*------------ print data to display FORM3 ------------*/
void   no_display(void)  /* display only a numerical value list */
{
  int atom_kind;

  printf("%6u %6.2f %5.2f %5.2f  %6.2f %6.3f   %3.2f %3.2f %3.2f   ",
	 sys.step,
	 sys.kin*sys.perMol/1000.0,
	 sys.pot*sys.perMol/1000.0,
	 (sys.kin+sys.pot)*sys.perMol/1000.0,
	 sys.kin*sys.e2t,
	 pt.pres_trace*sys.pp2gpa,
	 cell.Lx/(double)cell.nx,
	 cell.Ly/(double)cell.ny,
	 cell.Lz/(double)cell.nz);

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
    fprintf(fpout_positions,"%d %.5f %.5f %.5f", 
	    sys.ion[i],cell.rx[i],cell.ry[i],cell.rz[i]);
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
    vx = cell.vx[i]/SECD_2_FS;  /* velocity :   [A/fsec'] -> [A/fsec]     */
    vy = cell.vy[i]/SECD_2_FS;  /* unit     :   1 [fsec] = 10^(-15) [sec] */
    vz = cell.vz[i]/SECD_2_FS;

    fprintf(fpout_velocities,"%d %.8e %.8e %.8e", 
	    sys.ion[i], vx, vy, vz); /* [A/fsec] */
    cnt++;
    if(cnt == 4){
      fprintf(fpout_velocities," \n");
      cnt = 0;
    }
  }
  fprintf(fpout_velocities,"\n");
}

/*--- record cell data to file in [A] unit ---*/
void   r_cell(void)
{

  fprintf(fpout_cell,"%6u %.3f %.3f %.3f %.3f %.3f %.3f   %.3f\n",
	  sys.step,
	  cell.Lx, cell.Ly, cell.Lz,
	  cell.angle_ab, cell.angle_bc, cell.angle_ca, cell.vol);

  if(sys.step%10==1) fflush(fpout_cell);
}
