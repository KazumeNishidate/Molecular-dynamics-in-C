### Compile command ##########################################
###  for FreeBSD
CC= gcc
CFLAGS=  -Wall -fstrength-reduce -O 
######################################
###  for HPUX   << "-Aa": compile under ANSI mode [HPUX] >>
# CC=cc 
# CFLAGS= -Aa +O3 
######################################
###  for Sun Ultra1 Model200E (SPARCompiler C 3.0)
# CC=cc
# CFLAGS= -Aa -O
######################################
###  for DEC: degital UNIX
###  note: to make program, type "make run" in "dynamics" directory.
###
# CC= cc
# CFLAGS= -O3 -w0
##############################################################
LIBS= -lm

SOURCES= fft1.c control.c files.c main.c 
OBJECTS= fft1.o control.o files.o main.o 

TARGET= run

.c.o	:$(HEADERS) $(SOURCES)
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

$(TARGET)	:$(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean	:
	rm -f $(OBJECTS) *~ "#"* *.ln out
lint	:
	lint -cbx $(SOURCES)







