#!/usr/bin/perl -w

# compile-trajectory-index.pl
# ==============================================================================================
# Compile a master trajectory index file and convert trajectories to AMBER netCDF format for
# use with decompose-f95.
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

my $destination_directory = "netcdf-trajectories"; # where final trajectories should go
my $trajectory_index_filename = "trajectory-index"; # file to contain index of trajectories
my $prmtop = "prmtop.box";
my $nskip = 10;   # interval for writing data

# ==============================================================================================
# MAIN
# ==============================================================================================

# make directory to receive trajectories if it doesn't exist
if(! -e $destination_directory) {
  system("mkdir $destination_directory");
}

# open trajectory index file
open INDEX_FILE, ">", $trajectory_index_filename || die("Could not open $trajectory_index_filename ...");

# construct list of source trajectories
my $command = "ls | grep '^[0-9]\\{3\\}'";
#print "$command\n";
my @directories = `$command`;
my @trajectory_filenames = ();

# Make a list of all source trajectories.
foreach my $directory ( @directories ) {
    chomp $directory;
    print ">$directory\n";
    # my $trajectory_filename = "$directory/mdcrd.1.gz";
    my $trajectory_filename = "$directory/425_${directory}_cat.trj.gz";
    if( -e $trajectory_filename ) {
	push @trajectory_filenames, $trajectory_filename;
    }
}

# DEBUG
print "@trajectory_filenames\n";

# Process source trajectories.
my $trajectory_index = 1;
foreach my $source_trajectory ( @trajectory_filenames ) {
  # construct destination trajectory filename
  my $trajectory_filename = sprintf("%s/trajectory-%06d", $destination_directory, $trajectory_index);
  print "$trajectory_filename\n";
  
  # write ptraj input file
  my $ptraj_input = <<END;
 trajin $source_trajectory 1 10000000 $nskip
 trajout $trajectory_filename netcdf
END

  # save to file
  open PTRAJ_INFILE, ">", "ptraj.in";
  print PTRAJ_INFILE "$ptraj_input\n";
  close PTRAJ_INFILE;

  # run ptraj
  #system("ptraj $prmtop < ptraj.in >& /dev/null");
  system("ptraj $prmtop < ptraj.in");

  # add to trajectory index file
  # format: <weight> <directory>
  printf INDEX_FILE "%16.8f %s\n", 1.0, $trajectory_filename;

  # increment trajectory index
  $trajectory_index++;
}

close INDEX_FILE;



