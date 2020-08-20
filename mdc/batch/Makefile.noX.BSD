### Compile command #########################################
###  for FreeBSD of Linux with gcc
CC= gcc
CFLAGS= -Wall -O2 -finline-functions -fomit-frame-pointer -funroll-loops
##############################################################

### Include files and library ##############################
LIBS= -lm
############################################################

SOURCES= main.c ext.c pt.c control.c init.c\
	 files.c real.c rcprcl.c rv.c msd.c

OBJECTS= main.o ext.o pt.o control.o init.o\
	 files.o real.o rcprcl.o rv.o msd.o

TARGET= md 

.c.o	:$(SOURCES)
	$(CC) -c $(CFLAGS) $(INC) $< -o $@ 

$(TARGET)	:$(OBJECTS)
	$(CC) $(CFLAGS) $(RUN_OPTIONS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean	:
	rm -f $(TARGET) $(OBJECTS) *~ "#"* *.ln ./xsrc/*~ ./xsrc/"#"* ./xsrc/*.ln

lint	:
	lint $(SOURCES)

cflow   :
	cflow $(SOURCES)

fclean  :
	rm -f ../files/[a-z]*

