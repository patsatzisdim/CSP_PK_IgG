SHELL = /bin/sh
ARCH = macx
TARGET1= main
TARGETS= $(TARGET1)
SRCS1= SUB_main.f  \
dlsode.f  \
fdjac1.f \
rg.f  \
balanc.f  \
balbak.f  \
cdiv.f  \
d1mach.f  \
elmhes.f  \
eltran.f  \
hqr.f  \
hqr2.f  \
SUB_eigen.f  \
SUB_sinve.f \
SUB_smult.f \
SUB_unitary.f \
InitCond.f \
SUB_ProbDef.f \
SUB_sort.f \
SUB_diagnostics.f \
CSP_kern.f  
CSRCS1= 
OBJS1= $(SRCS1:.f=.o) $(CSRCS1:.c=.o)
DEPSRC = $(SRCS1) 
DEPCSRC = $(CSRCS1)  
LIBDIR =
LIBS1=
DEFS =
########### the following flags are for debugging #############
#FFLAGS =   -O -fno-automatic -Wunused -Wsurprising $(DEFS) -fexternal-blas -fblas-matmul-limit=50
########### the following flags are for optimized running #############
FFLAGS =  -O2 $(DEFS) -ffixed-line-length-132 -cpp
LDFLAGS = 
INCL = 
F77 = gfortran	
FC = $(F77)	
CC = cc	
CFLAGS =   $(DEFS) -O	
CPP = gcc	
MKDEP = 	

all: $(TARGETS)	

$(TARGET1): $(OBJS1)	
	$(F77) $(LDFLAGS) $(FFLAGS) $(INCL) -o $(TARGET1) $(OBJS1) $(LIBDIR) $(LIBS1)

clean:	
	rm -f $(OBJS1)  Makefile.bak

depend:	
	@(sed '/^# DO NOT DELETE THIS LINE/q' Makefile && \
	 $(MKDEP) CPP=$(CPP) $(DEPCSRC) $(DEPSRC) \
	) > Makefile.new
	@mv Makefile Makefile.bak
	@mv Makefile.new Makefile

.f.o:	
	$(F77) $(FFLAGS) -c $*.f

#	Dependencies of F77 object files
SUB_main.o:		paramet.i
#SUB_diagnostics.o:        ptMerge.i

# DO NOT DELETE THIS LINE	


