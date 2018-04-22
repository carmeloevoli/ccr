# C compiler and flags
CPPFLAGS = -I${GSL_DIR}/include 
LDFLAGS  = -L${GSL_DIR}/lib -lgsl -lgslcblas 
CC       = /usr/local/Cellar/gcc/5.2.0/bin/gcc-5 -fopenmp

# non-project dependencies
COSMO_FILES = ps.o cosmo_progs.o misc.o 

main: main.c 
	${CC} ${CPPFLAGS} -o a.out ${LDFLAGS} main.c

integral: main_integral.c
	${CC} ${CPPFLAGS} -o a.out ${LDFLAGS} main_integral.c

clean:
	rm -vf *.o *~ a.out
