CC = g++
CFLAGS = -Wall -O2 -std=c++11 -fopenmp
MKLINC = -I${MKLROOT}/include
CPLUSINC =
LFGSL = -lgsl -lgslcblas
LFMKL = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -liomp5
LFARP = -larpack
LFLAGS = -lpthread -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas

SOURCES = main.cpp \
		  basis.cpp \
		  matrix.cpp \
		  operator.cpp \
		  diag.cpp

OBJS = ${SOURCES:.cpp=.o}


all: $(OBJS)
	$(CC) $(OBJS) -o vortexUR0 $(LFGSL) $(LFMKL) $(LFARP)
	rm -f *.o


$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) $(MKLINC) $(CPLUSINC) $< -o $@

clean:
	rm -f *.o vortexUR*

