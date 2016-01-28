#!/usr/bin/python

#=============================================================================================
# Process trpzip2 data from IBM Almaden to prepare for automatic state decomposition.
#
# This code generates a number of AMBER format netCDF trajectories containing 10 ps trajectory segments.
# The first and last snapshot in each trajectory segment is written to the netcdf file, withot duplication of first/last configurations.
# Data is written starting from the first complete trajectory segment.
#=============================================================================================

#=============================================================================================
# REQUIREMENTS
#
# This code requires the 'pynetcdf' package, containing the Scientific.IO.NetCDF package built for numpy.
#
# http://pypi.python.org/pypi/pynetcdf/
# http://sourceforge.net/project/showfiles.php?group_id=1315&package_id=185504
#=============================================================================================

#=============================================================================================
# TODO
#=============================================================================================

#=============================================================================================
# CHAGELOG
#=============================================================================================

#=============================================================================================
# VERSION CONTROL INFORMATION
# * 2008-04-01 JDC
# Created file.
#=============================================================================================

__version__ = "$Revision: $ $Date: $"
# $Date: $
# $Revision:  $
# $LastChangedBy: $
# $HeadURL: $
# $Id: $

#=============================================================================================
# IMPORTS
#=============================================================================================

from numpy import * # numerical objects
#import Scientific.IO.NetCDF # netCDF file access
#from Scientific.IO import NetCDF
from pynetcdf import NetCDF
import os
import os.path
import gzip

#=============================================================================================
# PARAMETERS
#=============================================================================================

# system information 
natoms = 219

# list of parallel tempering directories to crawl and extract replica trajectories
repex_directories = ['rep6a', 'rep6b', 'rep6c', 'rep6d', 'rep6e', 'rep6f', 'rep6g', 'rep6h', 'rep6i']

# directory to store netcdf trajectories
output_directory = 'rep6-extracted'

# number of frames per iteration
nframes_per_iteration = 40

# period for velocity randomizations
nframes_per_trajectory = 10

#=============================================================================================
# SUBROUTINES
#=============================================================================================
def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line)
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def read_file(filename):
    """Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines

def initialize_netcdf(netcdf_file, title, natoms, is_periodic, has_velocities):
  """Initialize the given NetCDF file according to the AMBER NetCDF Convention Version 1.0, Revision B.

  ARGUMENTS
    netcdf_file (NetCDFFile object) - the file to initialize global attributes, dimensions, and variables for
    title (string) - the title for the netCDF file
    natoms (integer) - the number of atoms in the trajectories to be written
    is_periodic (boolean) - if True, box coordinates will also be stored
    has_velocities (boolean) - if True, the velocity trajectory variables will also be created

  NOTES
    The AMBER NetCDF convention is defined here:

    http://amber.scripps.edu/netcdf/nctraj.html    

  """

  # Create dimensions.
  netcdf_file.createDimension('frame', 0)        # unlimited number of frames in trajectory
  netcdf_file.createDimension('spatial', 3)      # number of spatial coordinates
  netcdf_file.createDimension('atom', natoms)    # number of atoms in the trajectory
  netcdf_file.createDimension('label', 5)        # label lengths for cell dimensions
  netcdf_file.createDimension('cell_spatial', 3) # naming conventions for cell spatial dimensions
  netcdf_file.createDimension('cell_angular', 3) # naming conventions for cell angular dimensions
  
  # Set attributes.
  setattr(netcdf_file, 'title', title)
  setattr(netcdf_file, 'application', 'AMBER')
  setattr(netcdf_file, 'program', 'sander')
  setattr(netcdf_file, 'programVersion', '8')
  setattr(netcdf_file, 'Conventions', 'AMBER')
  setattr(netcdf_file, 'ConventionVersion', '1.0')

  # Define variables to store unit cell data, if specified.
  if is_periodic:
    cell_spatial = netcdf_file.createVariable('cell_spatial', 'c', ('cell_spatial',))
    cell_angular = netcdf_file.createVariable('cell_angular', 'c', ('cell_spatial', 'label'))
    cell_lengths = netcdf_file.createVariable('cell_lengths', 'd', ('frame', 'cell_spatial'))
    setattr(cell_lengths, 'units', 'angstrom')
    cell_angles = netcdf_file.createVariable('cell_angles', 'd', ('frame', 'cell_angular'))
    setattr(cell_angles, 'units', 'degree')  

    netcdf_file.variables['cell_spatial'][0] = 'x'
    netcdf_file.variables['cell_spatial'][1] = 'y'
    netcdf_file.variables['cell_spatial'][2] = 'z'

    netcdf_file.variables['cell_angular'][0] = 'alpha'
    netcdf_file.variables['cell_angular'][1] = 'beta '
    netcdf_file.variables['cell_angular'][2] = 'gamma'

  # Define variables to store velocity data, if specified.
  if has_velocities:
    velocities = netcdf_file.createVariable('velocities', 'd', ('frame', 'atom', 'spatial'))    
    setattr(velocities, 'units', 'angstrom/picosecond')
    setattr(velocities, 'scale_factor', 20.455)  

  # Define coordinates and snapshot times.
  frame_times = netcdf_file.createVariable('time', 'f', ('frame',))
  setattr(frame_times, 'units', 'picosecond')
  frame_coordinates = netcdf_file.createVariable('coordinates', 'f', ('frame', 'atom', 'spatial'))
  setattr(frame_coordinates, 'units', 'angstrom')

  # Define optional data not specified in the AMBER NetCDF Convention that we will make use of.
  frame_energies = netcdf_file.createVariable('total_energy', 'f', ('frame',))
  setattr(frame_energies, 'units', 'kilocalorie/mole')
  frame_energies = netcdf_file.createVariable('potential_energy', 'f', ('frame',))
  setattr(frame_energies, 'units', 'kilocalorie/mole')  
  
  return

def write_netcdf_frame(netcdf_file, frame_index, time = None, coordinates = None, cell_lengths = None, cell_angles = None, total_energy = None, potential_energy = None):
  """Write a NetCDF frame.

  ARGUMENTS
    netcdf_file (NetCDFFile) - the file to write a frame to
    frame_index (integer) - the frame to be written

  OPTIONAL ARGUMENTS
    time (float) - time of frame (in picoseconds)
    coordinates (natom x nspatial NumPy array) - atomic coordinates (in Angstroms)
    cell_lengths (nspatial NumPy array) - cell lengths (Angstroms)
    cell_angles (nspatial NumPy array) - cell angles (degrees)
    total_energy (float) - total energy (kcal/mol)
    potential_energy (float) - potential energy (kcal/mol)

  """
  if time != None: netcdf_file.variables['time'][frame_index] = time      
  if coordinates != None: netcdf_file.variables['coordinates'][frame_index,:,:] = coordinates
  if cell_lengths != None: netcdf_file.variables['cell_lengths'][frame_index,:] = cell_lengths
  if cell_angles != None: netcdf_file.variables['cell_angles'][frame_index,:] = cell_angles
  if total_energy != None: netcdf_file.variables['total_energy'][frame_index] = total_energy
  if potential_energy != None: netcdf_file.variables['total_energy'][frame_index] = potential_energy
  
  return

def read_amber_energy_frame(infile):
    """Read a frame of energy components from the AMBER energy file.

    ARGUMENTS
      infile (Python file handle) - the file to read from

    RETURNS
      energies (Python dict) -- energies[keyword] contains the energy for the corresponding keyword
    """

    # number of lines per .ene block
    ene_lines_per_block = 10
    
    # energy keys
    energy_keys = [
        'Nsteps', 'time', 'Etot', 'EKinetic', # L0
        'Temp', 'T_solute', 'T_solv', 'Pres_scal_solu', # L1
        'Pres_scal_solv', 'BoxX', 'BoxY', 'BoxZ', # L2
        'volume', 'pres_X', 'pres_Y', 'pres_Z',
        'Pressure', 'EKCoM_x', 'EKCoM_y', 'EKCoM_z',
        'EKComTot', 'VIRIAL_x', 'VIRIAL_y', 'VIRIAL_z',
        'VIRIAL_tot', 'E_pot', 'E_vdw', 'E_el',
        'E_hbon', 'E_bon', 'E_angle', 'E_dih',
        'E_14vdw', 'E_14el', 'E_const', 'E_pol',
        'AV_permMoment', 'AV_indMoment', 'AV_totMoment', 'Density', 'dV/dlambda'
        ]

    # Read energy block.
    energies = dict()
    key_index = 0
    for line_counter in range(ene_lines_per_block):
        line = infile.readline() # read the line
        elements = line.split() # split into elements
        elements.pop(0) # drop the 'L#' initial element
        for element in elements:
            key = energy_keys[key_index] # get the key
            energies[key] = float(element) # store the energy
            key_index += 1 # increment index

    return energies

def read_amber_frame(amber_trj_file, natoms, is_periodic):
  """Read a coordinate frame from an AMBER text-formatted trajectory file.

  ARGUMENTS
    amber_trj_file (Python file handle) - the file to read from
    natoms (integer) - the number of atoms per frame
    is_periodic (boolean) - if True, box coordinates are read too

  RETURNS
    coordinates (numpy array) - or None if exception is encountered
    boxsize (numpy array) - or None if exception is encountered

  TODO
    Change return values to key-value pairs in a dict() to allow handling of cell lengths and cell angles.
  """

  # determine number of lines per frame
  nlines = int(ceil(natoms * 3 / 10.0))
  coordinates = zeros([natoms, 3], float32)

  # parse coordinates
  atom_index = 0
  dimension_index = 0    
  for line_index in range(nlines):
    line = amber_trj_file.readline()
    elements = line.split()
    for element in elements:
      coordinates[atom_index, dimension_index] = float(element)
      dimension_index += 1
      # wrap dimension
      if dimension_index == 3:
        dimension_index = 0
        atom_index += 1
    
  # skip box
  if is_periodic:
    # get box line
    line = amber_trj_file.readline()
    # determine how many numbers are needed to describe the box size
    elements = line.split()
    nbox = len(elements)
    boxsize = zeros([nbox], float32)
    index = 0
    for element in elements:
      boxsize[index] = float(element)
      index += 1
    # return coordinates and box size
    return (coordinates, boxsize)

  return coordinates
  
#=============================================================================================
# MAIN
#=============================================================================================

# TODO: Create netcdf directory to store output if it doesn't exist.

# Crawl replica directories
for repex_directory in repex_directories:
  # create subdirectory to store data
  netcdf_directory = os.path.join(output_directory, repex_directory, 'netcdf')
  os.makedirs(netcdf_directory)
  trajectory_energy_filename = os.path.join(output_directory, repex_directory, 'trajectory-energies')  

  # Open file to store all trajectory segment energies.
  trajectory_energy_file = open(trajectory_energy_filename, 'w')

  # Read replica-for-temperature files.
  temp_for_replica_lines = read_file(os.path.join(repex_directory, 'temp_for_replica.txt'))
  replica_for_temp_lines = read_file(os.path.join(repex_directory, 'replica_at_temp.txt'))

  # Read list of temperatures.
  lines = read_file(os.path.join(repex_directory, '0', 'temp_list'))
  temperatures = list()
  lookup_temperature_index = dict() # lookup_temperature_index[temperature_string] is numerical index (starting from 0) of temperature
  ntemperatures = 0 # number of temperatures
  for line in lines:
    temperature = line.strip()
    temperatures.append(temperature)
    lookup_temperature_index[temperature] = ntemperatures
    ntemperatures += 1 
  print temperatures

  # Determine number of replicas.
  elements = temp_for_replica_lines[0].split()
  nreplicas = len(elements)
  replicas = range(nreplicas)

  # Open all temperatures for reading.
  amber_temperature_file = dict()
  for temperature in temperatures:
    # Open the trajectory.    
    filename = os.path.join(repex_directory, '%(temperature)s.trj.gz' % vars())
    print "Opening AMBER trajectory file %(filename)s for reading..." % vars()
    infile = gzip.open(filename, 'r')
    # eat header line    
    line = infile.readline() 
    # Store the file handle.
    amber_temperature_file[temperature] = infile

  # Open all energy files for reading.
  amber_energy_file = dict()
  for temperature in temperatures:
    # Open the energy file.
    filename = os.path.join(repex_directory, '%(temperature)s.ene.gz' % vars())
    print "Opening AMBER energy file %(filename)s for reading..." % vars()
    infile = gzip.open(filename, 'r')
    # eat header line    
    line = infile.readline() 
    # Store the file handle.
    amber_energy_file[temperature] = infile

  # Open all replicas for writing.
  nframes = zeros([nreplicas], int32)
  ntrajectories = zeros([nreplicas], int32) # ntrajectories[replica_index] is the number of trajectories that have been written to replica 'replica_index'
  ntrajectories -= 1 # set all to -1
  netcdf_replica_file = dict()
  for replica in replicas:
    # Open the NetCDF trajectory file.
    filename = os.path.join(netcdf_directory, "%(repex_directory)s-%(replica)d.netcdf" % vars())
    print "Opening NetCDF trajectory file %(filename)s for writing..." % vars()
    outfile = NetCDF.NetCDFFile(filename, 'w')
    
    # Initialize the NetCDF file.
    title = "%(repex_directory)s replica %(replica)d" % vars()
    is_periodic = True
    has_velocities = False
    initialize_netcdf(outfile, title, natoms, is_periodic, has_velocities)

    # Store file handle.
    netcdf_replica_file[replica] = outfile

  # Process iterations.
  niterations = len(temp_for_replica_lines)
  for iteration in range(niterations):
    print "Iteration %d..." % iteration

    # Extract permuation information.
    line = temp_for_replica_lines[iteration]
    temperature_for_replica = line.split()
    for replica in replicas:
      # get the name of the temperature for this replica
      temperature = temperature_for_replica[replica]

      # read from this temperature
      total_energies = zeros([nframes_per_trajectory], float64) # storage for total energies
      trajectory_frame_index = 0
      for counter in range(nframes_per_iteration):
        # Read AMBER snapshot.
        (coordinates, cell_lengths) = read_amber_frame(amber_temperature_file[temperature], natoms, is_periodic)
        # Read AMBER energy snapshot.
        energies = read_amber_energy_frame(amber_energy_file[temperature])        
        # Compute mean and standard error of total energy for this trajectory segment.
        total_energies[trajectory_frame_index] = energies['Etot'] # store total energy
        trajectory_frame_index += 1
        # Write only last frame in each trajectory, along with trajectory energy
        if trajectory_frame_index == nframes_per_trajectory:
            # Write NetCDF frame.
            write_netcdf_frame(netcdf_replica_file[replica], nframes[replica], coordinates = coordinates, cell_lengths = cell_lengths, cell_angles = array([90., 90., 90.]), total_energy = energies['Etot'], potential_energy = energies['E_pot'])
            nframes[replica] += 1

            # Compute mean total energy and std error of mean over each trajectory.
            Etot_mean = total_energies.mean()
            Etot_std = total_energies.std()
            Etot_err = Etot_std / sqrt(float(nframes_per_trajectory))

            # DEBUG
            #print "replica %5d    temperature %8s K    trajectory_index %5d    Etot = %8.1f +- %.1f kcal/mol" % (replica, temperature, ntrajectories[replica], Etot_mean, Etot_err)

            # Write mean total energy over trajectory, provided this is not the first trajectory.
            if (ntrajectories[replica] > -1):
                trajectory_energy_file.write("%5d %5d %8d %16.10e %6.3f\n" % (replica, lookup_temperature_index[temperature], ntrajectories[replica], Etot_mean, Etot_err))

            # increment trajectory counter
            ntrajectories[replica] += 1

            # reset frame counter within trajectory
            trajectory_frame_index = 0
            
      # Sync replica netcdf file to disk.
      netcdf_replica_file[replica].sync()

    # Sync trajectory energy file to disk.
    trajectory_energy_file.flush()                                       
      
  # Close all files.
  print "Closing files..."  

  for replica in replicas:
    # Close the NetCDF trajectory file.
    netcdf_replica_file[replica].close()
    
  for temperature in temperatures:    
    amber_temperature_file[temperature].close()

# All done -- close all remaining files
print "Closing all remaining files..."
trajectory_energy_file.close()

