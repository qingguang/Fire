####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

#CC = mpicc 
CC = scorep --user mpicc
#CC= mpicc
CFLAGS = -Wall -g -O3 
LIBS = -lm

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

CFLAGS += $(METIS_INC)
LIBS   += $(METIS_LIB)
CFLAGS += $(PAPI_INC)
LIBS += $(PAPI_LIB)
#LIBS += -lpapi
#LIBPOS  += libmetis.a
#LIBS += -lmetis
SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c test_functions.c
#SRCS = initialization.c util_read_files.c util_write_files.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS) 
#test_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg *.out *.vtk $(LIBPOS)
