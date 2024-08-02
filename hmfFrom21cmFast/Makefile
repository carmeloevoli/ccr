# C compiler and flags
GSL_DIR  = /usr/local/Cellar/gsl/2.8
CPPFLAGS = -I${GSL_DIR}/include 
LDFLAGS  = -L${GSL_DIR}/lib -lgsl -lgslcblas 
CC       = /usr/local/Cellar/gcc/14.1.0_2/bin/gcc-14 -fopenmp

# non-project dependencies
COSMO_FILES = ps.o cosmo_progs.o misc.o 

main: main.c 
	${CC} ${CPPFLAGS} -o a.out ${LDFLAGS} main.c

integral: main_integral.c
	${CC} ${CPPFLAGS} -o a.out ${LDFLAGS} main_integral.c

clean:
	rm -vf *.o *~ a.out
