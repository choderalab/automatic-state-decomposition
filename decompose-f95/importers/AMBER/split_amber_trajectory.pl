#!/usr/bin/perl -w

# split_amber_trajectory.pl
# ==============================================================================================
# Split a concatenated AMBER trajectory into individual netCDF trajectory segments, and write
# a master trajectory index file with weights for use with decompose-f95.
#
# Written by John Chodera, Dill lab, UCSF 2006.
# ==============================================================================================#
# Copyright (C) 2006  John D. Chodera
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# ==============================================================================================
# COMMAND-LINE ARGUMENTS:
# ./split_amber_trajectory.pl
# ==============================================================================================
# REQUIREMENTS:
# - AMBER9 'ptraj'

# ==============================================================================================
# USE STATEMENTS
# ==============================================================================================

use strict;

# ==============================================================================================
# PARAMETERS
# ==============================================================================================

#my $source_trajectory = "400K-replica.trj";  # AMBER source trajectory
my $source_trajectory = "snapshots.netcdf";  # AMBER source trajectory
my $nsnapshots_per_trajectory = 200;   # number of snapshots per trajectory segment
my $ntrajectories = 1000;  # number of trajectory segments

my $destination_directory = "trajectories"; # where final trajectories should go
my $trajectory_index_filename = "trajectory-index"; # file to contain index of trajectories

# ==============================================================================================
# MAIN
# ==============================================================================================

# make directory to receive trajectories if it doesn't exist
if(! -e $destination_directory) {
  system("mkdir $destination_directory");
}

# open trajectory index file
open INDEX_FILE, ">", $trajectory_index_filename || die("Could not open $trajectory_index_filename ...");

# loop over trajectory
for(my $trajectory_index = 1; $trajectory_index <= $ntrajectories; $trajectory_index++) {
  # construct trajectory filename
  my $trajectory_filename = sprintf("%s/trajectory-%06d", $destination_directory, $trajectory_index);
  print "$trajectory_filename\n";
  
  # determine start and stop
  my $start = $nsnapshots_per_trajectory * ($trajectory_index-1) + 1;
  my $stop = $nsnapshots_per_trajectory * $trajectory_index;
  
  # write ptraj input file
  my $ptraj_input = <<END;
 trajin $source_trajectory $start $stop
 trajout $trajectory_filename netcdf
END

  # save to file
  open PTRAJ_INFILE, ">", "ptraj.in";
  print PTRAJ_INFILE "$ptraj_input\n";
  close PTRAJ_INFILE;

  # run ptraj
  system("ptraj a1e.prmtop_nowat < ptraj.in >& /dev/null");

  # add to trajectory index file
  # format: <weight> <directory>
  printf INDEX_FILE "%16.8f %s\n", 1.0, $trajectory_filename;
}

close INDEX_FILE;



