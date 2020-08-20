### Compile command #########################################
###  for FreeBSD or LINUX (GCC)
CC= gcc
CFLAGS= -Wall -fstrength-reduce -O3
###  for DEC: degital UNIX
# CC= cc
# CFLAGS= -O4 -w0 -fast 
###  for HPUX   << "-Aa": compile under ANSI mode [HPUX] >>
# CC=cc
# CFLAGS= -Aa +O3 
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0)
# CC=cc
# CFLAGS= -Aa -O
##############################################################

### Include files and library ##############################
LIBS= -lm
############################################################

SOURCES= main.c ext.c pt.c control.c init.c\
	 files.c rv.c msd.c tb.c chebyshev.c. fermi.c\
	 chemical.c tb_init.c lookup.c hamiltonian.c repulsive.c

OBJECTS= main.o ext.o pt.o control.o init.o\
	 files.o rv.o msd.o tb.o chebyshev.o fermi.o\
	 chemical.o tb_init.o lookup.o hamiltonian.o repulsive.o

TARGET= md 

.c.o	:$(SOURCES)
	$(CC) -c $(CFLAGS) $(INC) $< -o $@ 

$(TARGET)	:$(OBJECTS)
	$(CC) $(CFLAGS) $(RUN_OPTIONS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean	:
	rm -f $(OBJECTS) *core *~ "#"* *.ln ./xsrc/*~ ./xsrc/"#"* ./xsrc/*.ln

fclean	:
	rm -f ../files/positions ../files/out ../files/velocities

lint	:
	lint -cbx $(SOURCES)

cflow   :
	cflow -ix $(SOURCES)
