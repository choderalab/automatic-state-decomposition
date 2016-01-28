.SUFFIXES:
.SUFFIXES: .c .o

# Location of gromacs header and library files.
GROMACSINC = $(HOME)/local/gromacs-3.3.1/include
GROMACSLIB = $(HOME)/local/gromacs-3.3.1/lib


# Location of netCDF header and library files.
NETCDFINC = $(HOME)/local/src/netcdf-3.6.1/include
NETCDFLIB = $(HOME)/local/src/netcdf-3.6.1/lib

# For opteron with Intel C compiler
#CC = icc
#OPTFLAGS = -O -g

# For Mac OS X with GCC
CC = gcc
OPTFLAGS = -O -g

INCFLAGS = -I$(GROMACSINC) -I$(NETCDFINC)
LDOPTS = -L$(GROMACSLIB) -lgmx -L$(NETCDFLIB) -lnetcdf -g

# Set flags
CCOPTS = $(OPTFLAGS)

LIBPATH =
SYSLIBS = -lm

# Default target: mod_samp executable:
xtc2netcdf: xtc2netcdf.o
	$(CC) -o xtc2netcdf xtc2netcdf.o $(OPTFLAGS) $(LDOPTS) $(LIBPATH) $(SYSLIBS)

clean:
	rm *.o *.mod

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $<
