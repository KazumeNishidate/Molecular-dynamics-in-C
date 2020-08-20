### Compile command ##########################################
###  for FreeBSD
CC= gcc
CFLAGS= -Wall -fstrength-reduce -O 
###  for DEC: degital UNIX
# CC= cc
# CFLAGS= -O3 -w0
###  for HPUX   << "-Aa": compile under ANSI mode [HPUX] >>
# CC=cc 
# CFLAGS= -Aa +O3 
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0)
# CC=cc
# CFLAGS= -Aa -O
##############################################################

### Include files and library ################################
###  where is the X include file "Xlib.h" and "Xutil.h"?
INC=  -I/usr/include/X11R5 -I/usr/include/X11 -I/usr/X11R6/include
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0: OpenWindows)
# INC=-I/usr/openwin/include

###  where is the X library "libX11.a"?
LIBS= -lm -L/usr/X11R6/lib -lX11
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0: OpenWindows)
# LIBS= -lm -L/usr/openwin/lib -lX11
##############################################################

###  Run options  ############################################
RUN_OPTIONS= -DXD
## to display in X Window
# -DXD
#############################################################

SOURCES= main.c ext.c pt.c control.c init.c\
	files.c rv.c msd.c tb.c chebyshev.c. fermi.c\
	chemical.c tb_init.c lookup.c hamiltonian.c repulsive.c\
	./xsrc/xdmain.c ./xsrc/xdinit.c ./xsrc/xdicon.c\
	./xsrc/xhist.c ./xsrc/xmsd.c ./xsrc/xnet.c\
	./xsrc/network.c ./xsrc/xunit.c ./xsrc/xmycolor.c 
OBJECTS= main.o ext.o pt.o control.o init.o\
	 files.o rv.o msd.o tb.o chebyshev.o fermi.o\
	 chemical.o tb_init.o lookup.o hamiltonian.o repulsive.o\
	 ./xsrc/xdmain.o ./xsrc/xdinit.o ./xsrc/xdicon.o\
	 ./xsrc/xhist.o ./xsrc/xmsd.o ./xsrc/xnet.o\
	 ./xsrc/network.o ./xsrc/xunit.o ./xsrc/xmycolor.o 

TARGET= md 

.c.o	:$(SOURCES)
	$(CC) -c $(CFLAGS) $(INC) $(RUN_OPTIONS) $< -o $@ 

$(TARGET)	:$(OBJECTS)
	$(CC) $(CFLAGS) $(RUN_OPTIONS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean	:
	rm -f $(OBJECTS) *~ "#"* *.ln ./xsrc/*~ ./xsrc/"#"* ./xsrc/*.ln

fclean	:
	rm -f ../files/positions ../files/out ../files/velocities

lint	:
	lint -cbx $(SOURCES)

cflow   :
	cflow -ix $(SOURCES)

