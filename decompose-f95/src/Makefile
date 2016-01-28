.SUFFIXES:
.SUFFIXES: .f .f90 .f95 .o .mod .F

# Location of ARPACK
#ARPACKHOME = $(HOME)/local/src/ARPACK
#ARPACKLIBDIR = $(ARPACKHOME)/lib
#ARPACKLIB = arpack_OPTERON

# Location of XMLF90
XMLF90HOME = $(HOME)/local/src/xmlf90-1.2g/macros/
XMLF90LIB = $(XMLF90HOME)/lib
XMLF90INC = $(XMLF90HOME)/modules

# Location of NetCDF
NETCDFHOME = /scratch/groups/pande/amber9/src/netcdf/
NETCDFLIB = $(NETCDFHOME)/lib
NETCDFINC = $(NETCDFHOME)/include

# For opteron with Intel Fortran Compiler 9.1 and OpenMP
#CC = icc
#FC = ifort
#FC90 = ifort
#FC95 = ifort
#OPTFLAGS = -O3 -static -no-prec-div -openmp -g -check bounds -fpstkchk -inline-debug-info -debug extended -traceback -threads
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include -I$(HOME)/local/src/acml3.5.0/gfortran64/include -I$(XMLF90INC)
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -L$(HOME)/local/src/acml3.5.0/gfortran64/lib -lacml -lgfortran -lnetcdf -g -L$(XMLF90LIB) -lflib -openmp -static -inline-debug-info -debug extended -traceback -threads

# For opteron with Intel Fortran Compiler 9.1
#CC = gcc
#FC = ifort
#FC90 = ifort
#FC95 = ifort
#OPTFLAGS = -O3 -static -no-prec-div -ipo -g
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -lnetcdf -static -g

# For opteron with gfortran
#CC = gcc
#FC = gfortran
#FC90 = gfortran
#FC95 = gfortran
#OPTFLAGS = -march=opteron -O3 -ffast-math
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -lnetcdf

# For opteron with g95
#CC = gcc
#FC = g95
#FC90 = g95
#FC95 = g95
#OPTFLAGS = -march=opteron -O3 -ffast-math
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -lnetcdf

# For PowerPC with gcc and g95
#CC = gcc
#FC = g95
#FC90 = g95
#FC95 = g95
#OPTFLAGS = -mpowerpc -O3 -ffast-math -fno-second-underscore # -fno-underscoring -fbounds-check
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include -I$(HOME)/local/src/xmlf90-1.2g/macros/modules
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -lnetcdf -L$(HOME)/local/src/xmlf90-1.2g/macros/lib -lflib -framework vecLib

# For PowerPC with gcc and gfortran
#CC = gcc
#FC = gfortran
#FC90 = gfortran
#FC95 = gfortran
#OPTFLAGS = -mpowerpc -O3 -ffast-math -fno-second-underscore # -fno-underscoring -fbounds-check
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include -I$(HOME)/local/src/xmlf90-1.2g/macros/modules
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -lnetcdf -L$(HOME)/local/src/xmlf90-1.2g/macros/lib -lflib -framework vecLib

# For Intel Mac with gfortran
#CC = gcc
#FC = gfortran
#FC90 = gfortran
#FC95 = gfortran
#OPTFLAGS = -O3 -ffast-math -fno-second-underscore # -fno-underscoring -fbounds-check -m64
#INCFLAGS = -I$(HOME)/local/src/netcdf-3.6.1/include -I$(HOME)/local/src/xmlf90-1.2g/macros/modules
#LDOPTS = -L$(HOME)/local/src/netcdf-3.6.1/lib -lnetcdf -L$(HOME)/local/src/xmlf90-1.2g/macros/lib -lflib -framework vecLib

# For Intel Mac with Intel Fortran Compiler 9.1 and OpenMP
#CC = icc
#FC = ifort
#FC90 = ifort
#FC95 = ifort
#OPTFLAGS = -O3 -no-prec-div -openmp -g -check all -inline-debug-info -debug extended -traceback -threads
#INCFLAGS = -I$(HOME)/local/include -I$(XMLF90INC)
#LDOPTS = -L$(XMLF90LIB) -lflib -framework vecLib -openmp -inline-debug-info -debug extended -traceback -threads -L$(HOME)/local/lib -lnetcdff -lnetcdf

# For Xeon cluster node with Intel Fortran compiler.
CC = icc
FC = ifort
FC90 = ifort
FC95 = ifort
OPTFLAGS = -O3 -no-prec-div -openmp -g -check all,noarg_temp_created -inline-debug-info -debug extended -traceback -threads
INCFLAGS = -I$(HOME)/local/include -I$(XMLF90INC) -I$(NETCDFINC)
LDOPTS = -L$(NETCDFLIB) -L$(XMLF90LIB) -lflib -openmp -inline-debug-info -debug extended -traceback -threads -L$(HOME)/local/lib -lnetcdff -lnetcdf -llapack

# Set flags
FC90FLAGS = $(OPTFLAGS)
FC95FLAGS = $(OPTFLAGS)
FFLAGS = $(OPTFLAGS)
CCOPTS = $(OPTFLAGS)

LIBPATH =
SYSLIBS = -lm

# Object code (like F90SRCS but with .o suffix):
#F90OBJS = numeric_kinds.o constants.o timer.o eigs_arpack.o eigs_lapack.o eigs_trlan.o eispack.o kabsch_rmsd.o theobald_rmsd.o netcdfio.o pdbio.o utilities.o transitionmatrix.o decompose.o
F90OBJS = numeric_kinds.o constants.o timer.o eigs_lapack.o kabsch_rmsd.o theobald_rmsd.o netcdfio.o pdbio.o utilities.o transitionmatrix.o decompose.o

# Default targets:
default: decompose kmeans

# Automatic state decomposition driver
decompose: $(F90OBJS) decomposedriver.o 
	$(FC90) -o decompose decomposedriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# K-means clustering driver
kmeans: $(F90OBJS) kmeansdriver.o 
	$(FC90) -o kmeans kmeansdriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# Compute timescales for existing state decomposition
compute-timescales: $(F90OBJS) timescalesdriver.o 
	$(FC90) -o compute-timescales timescalesdriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# ?
compute-distances: $(F90OBJS) distancedriver.o 
	$(FC90) -o compute-distances distancedriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# Compute metastability metric for given decomposition
compute-metastability: $(F90OBJS) metastabilitydriver.o 
	$(FC90) -o compute-metastability metastabilitydriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# Compute transition matrix for given decomposition
compute-transitionmatrix: $(F90OBJS) transitionmatrixdriver.o 
	$(FC90) -o compute-transitionmatrix transitionmatrixdriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# Generate bootstrap sampled transition matrices for given decomposition
compute-transitionmatrixbootstrap: $(F90OBJS) transitionmatrixbootstrapdriver.o 
	$(FC90) -o compute-transitionmatrixbootstrap transitionmatrixbootstrapdriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

# Test functions in decomposition class
testdriver: $(F90OBJS) testdriver.o
	$(FC90) -o testdriver testdriver.o $(F90OBJS) $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

clean:
	rm *.o *.mod

# Fixed source, if any:
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<

.F.o:
#       $(FC) -c $(FFLAGS) $<
	$(FC) -c $(CCOPTS) $<

.f90.mod:
	$(FC90) -c $(FC90FLAGS) $(INCFLAGS) $<

.f95.mod:
	$(FC95) -c $(FC95FLAGS) $(INCFLAGS) $<

# Free form source:
.f90.o:
	$(FC90) -c $(FC90FLAGS) $(INCFLAGS) $<

.f95.o:
	$(FC95) -c $(FC95FLAGS) $(INCFLAGS) $<

# Dependencies based on "use" statements:
#mod_samp.o:     booly.o
