### Compile command ##########################################
###  for FreeBSD, Linux, Linux/Alpha
CC= gcc
CFLAGS= -Wall -O2 -finline-functions -fomit-frame-pointer -funroll-loops
###  for Linux/Alpha with Compaq C Compiler
# CC= ccc
# CFLAGS= -arch ev56 -O3 -g3
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

### Math Library ################################
###  for FreeBSD, Linux, Linux/Alpha
ML= -lm
###  for Linux/Alpha with libcpml.a or libffm.a
# ML= -lcpml
# ML= -lffm -lm
#################################################

### Include files and library ################################
###  where is the X include file "Xlib.h" and "Xutil.h"?
# INC=  -I/usr/include/X11R5 -I/usr/include/X11 -I/usr/X11R6/include
INC= -I/usr/X11R6/include
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0: OpenWindows)
# INC=-I/usr/openwin/include

###  where is the X library "libX11.a"?
# LIBS= $(ML) -L/usr/lib -L/usr/X11R6/lib -L/usr/lib/X11R5 -lX11
LIBS= $(ML) -L/usr/X11R6/lib -lX11
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0: OpenWindows)
# LIBS= $(ML) -L/usr/openwin/lib -lX11
##############################################################

###  Run options  ############################################
RUN_OPTIONS= -DXD
##############################################################

SOURCES= main.c ext.c pt.c control.c init.c\
	 files.c real.c rv.c msd.c tb.c jacobi.c debug.c\
	 ./xsrc/xdmain.c ./xsrc/xdinit.c ./xsrc/xdicon.c\
	 ./xsrc/xhist.c ./xsrc/xmsd.c ./xsrc/xnet.c\
	 ./xsrc/network.c ./xsrc/xunit.c ./xsrc/xmycolor.c 

OBJECTS= $(SOURCES:.c=.o)

HEADERS= md.h potential.h prototypes.h ./xsrc/xr.h

TARGET= md 

.c.o:
	$(CC) $(CFLAGS) -c $(INC) $(RUN_OPTIONS) $< -o $@ 

$(TARGET):$(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(TARGET)

$(OBJECTS):$(HEADERS)

clean:
	rm -f *.o *~ "#"* *.ln ./xsrc/*.o ./xsrc/*~ ./xsrc/"#"* ./xsrc/*.ln $(TARGET)

lint:
	lint -cbx $(SOURCES)

cflow:
	cflow -ix $(SOURCES)

tag:
	etags $(SOURCES)
