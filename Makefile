CC      = g++ 
CFLAGS  = -O2 -Wall -ferror-limit=3
EXEC    = a.out

INCDIR += -I$(GSL_DIR)/include 
LIBDIR += -L$(GSL_DIR)/lib -lgsl

all: $(EXEC) 

$(EXEC): main.o dump.o SFR.o misc.o tridiag.o cosmo_progs.o ps.o reionization.o utilities.o
	$(CC) $(CFLAGS) $(INCDIR) -o $@ $^ $(LIBDIR)

%.o: %.cpp constants.h
	$(CC) $(CFLAGS) $(INCDIR) -c -o $@ $< 	

.PHONY: clean

clean:
	@rm -vf *.o
	@rm -vf $(EXEC)
