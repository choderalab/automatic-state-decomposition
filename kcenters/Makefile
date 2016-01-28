all: kcenters

kcenters.o: kcenters.cpp
	../../../mpich-1.2.7p1/bin/mpicxx -c kcenters.cpp -O3 -ffast-math -mfpmath=sse -funroll-loops -I/usr/local/gromacs/include/gromacs

kcenters: kcenters.o theobald_rmsd.o
	../../../mpich-1.2.7p1/bin/mpicxx -o kcenters kcenters.o theobald_rmsd.o -L/usr/local/gromacs/lib -lgmx -lm

theobald_rmsd.o: theobald_rmsd.cpp
	g++ -c theobald_rmsd.cpp -O3 -ffast-math -mfpmath=sse -funroll-loops -I/usr/local/gromacs/include/gromacs
