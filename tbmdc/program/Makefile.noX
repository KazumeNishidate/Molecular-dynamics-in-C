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

### Include files and library ##############################
###  for FreeBSD, Linux, Linux/Alpha
LIBS= -lm
###  for Linux/Alpha with libcpml.a or libffm.a
# LIBS= -lcpml
# LIBS= -lffm -lm
############################################################

SOURCES= main.c ext.c pt.c control.c init.c\
	 files.c real.c rv.c msd.c tb.c jacobi.c debug.c

OBJECTS= $(SOURCES:.c=.o)

HEADERS= md.h potential.h prototypes.h

TARGET= md 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@ 

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
