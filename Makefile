CC      = g++ 
CFLAGS  = -O2 -Wall
EXEC    = a.out

INCDIR += 
LIBDIR += 

all: $(EXEC) 

$(EXEC): main.o utilities.o
	$(CC) $(CFLAGS) $(INCDIR) -o $@ $^ $(LIBDIR)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c -o $@ $< 	

.PHONY: clean

clean:
	@rm -vf *.o
	@rm -vf $(EXEC)
