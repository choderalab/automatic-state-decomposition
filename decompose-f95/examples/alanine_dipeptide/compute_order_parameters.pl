#!/usr/bin/perl

# Compute various order parameters for the alanine dipeptide using ptraj.

#===============================================================================
# USE DECLARATIONS
#===============================================================================

use strict;
use IO::File;

#===============================================================================
# PARAMETERS
#===============================================================================

my $trajdir = "trajectories";
my $outdir = "orderparameters";

#===============================================================================
# MAIN
#===============================================================================

# generate list of trajectory files
my @trjlist = `ls $trajdir/trajectory-*`;

# analyze each trajectory
foreach my $trjfile ( @trjlist )
{
    # extract index
    $trjfile =~ /trajectory-(\d+)/;
    my $trajectory_index = $1;

    print "$trajectory_index\n";

    # write ptraj input file
    my $ptraj_input = <<END;
 trajin $trjfile
 dihedral omega1 :1\@CH3 :1\@C :2\@N :2\@CA out $outdir/omega1-$trajectory_index
 dihedral phi :1\@C :2\@N :2\@CA :2\@C out $outdir/phi-$trajectory_index
 dihedral psi :2\@N :2\@CA :2\@C :3\@N out $outdir/psi-$trajectory_index
 dihedral omega2 :2\@CA :2\@C :3\@N :3\@CH3 out $outdir/omega2-$trajectory_index
END

    # DEBUG
    print $ptraj_input;

  # save to file
  open PTRAJ_INFILE, ">", "ptraj.in";
  print PTRAJ_INFILE "$ptraj_input\n";
  close PTRAJ_INFILE;

  # run ptraj
  system("ptraj a1e.prmtop_nowat < ptraj.in >& /dev/null");
}

