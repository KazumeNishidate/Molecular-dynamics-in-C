### Compile command #########################################
###  for FreeBSD, Linux
CC= gcc
CFLAGS= -Wall -O2 -finline-functions -fomit-frame-pointer -funroll-loops
# CFLAGS= -fstrength-reduce -Wall -O4
###  for Linux/Alpha with Compaq C Compiler
# CC= ccc
# CFLAGS= -arch ev56 -O3 -g3
###  for HPUX   << "-Aa": compile under ANSI mode [HPUX] >>
# CC=cc
# CFLAGS= -Aa +O3 
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0)
# CC=cc
# CFLAGS= -Aa -O
##############################################################

### Include files and library ##############################
###  for FreeBSD, Linux
LIBS= -lm
###  for Linux/Alpha with libcpml.a or libffm.a
# LIBS= -lcpml
# LIBS= -lffm -lm
############################################################

SOURCES= main.c control.c init.c ext.c bhs_init.c\
	 plane_wave.c pseudo.c eigen.c fermi.c\
	 fft.c ewald.c xc.c energy.c quench.c next.c monitor.c pt.c\
	 msd.c

OBJECTS= $(SOURCES:.c=.o)

HEADERS=./headers/physical_constants.h ./headers/function_macros.h\
	./headers/prototypes.h ./headers/size.h	./headers/md.h

TARGET= md 

$(TARGET):$(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

$(OBJECTS):$(HEADERS)

clean:
	rm -f $(TARGET) *.o *~ "#"* *.ln *core ./headers/*~

fclean:
	rm -f ./files/[a-z]*

lint:
	lint -ix $(SOURCES)

cflow:
	cflow -ix $(SOURCES)

tag:
	etags $(SOURCES)
