/*
 * $Id: xtc2netcdf.c,v 1.2 2006/09/18 07:20:05 jchodera Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "fatal.h"
#include "xtcio.h"
#include "enxio.h"
#include "smalloc.h"
#include "names.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "txtdump.h"

#include <netcdf.h>

/* Storage for netCDF trajectory variable IDs */
typedef struct _netcdfTrajectoryInfo {
  int ncid;
  int currentFrame;
  int frameDID;
  int spatialDID;
  int atomDID;
  int cell_spatialDID;
  int cell_angularDID;
  int labelDID;
  int spatialVID;
  int timeVID;
  int coordinateVID;
  int cell_spatialVID;
  int cell_angularVID;
  int cellLengthVID;
  int cellAngleVID;
  int velocityVID;
  double velocityScale;
  char *timeUnits;
  char *coordinateUnits;
  char *cellLengthUnits;
  char *cellAngleUnits;
  char *velocityUnits;
  char *Conventions;
  char *ConventionVersion;
  float *R;

} netcdfTrajectoryInfo;

#define INITIALIZE_netcdfTrajectoryInfo( _p_ ) \
  _p_->ncid = -1; \
  _p_->currentFrame = 0; \
  _p_->frameDID = 0; \
  _p_->spatialDID = 0; \
  _p_->atomDID = 0; \
  _p_->cell_spatialDID = 0; \
  _p_->cell_angularDID = 0; \
  _p_->labelDID = 0; \
  _p_->spatialVID = 0; \
  _p_->timeVID = 0; \
  _p_->coordinateVID = 0; \
  _p_->cell_spatialVID = 0; \
  _p_->cell_angularVID = 0; \
  _p_->cellLengthVID = 0; \
  _p_->cellAngleVID = 0; \
  _p_->velocityVID = 0; \
  _p_->timeUnits = NULL; \
  _p_->coordinateUnits = NULL; \
  _p_->cellLengthUnits = NULL; \
  _p_->cellAngleUnits = NULL; \
  _p_->velocityUnits = NULL; \
  _p_->velocityScale = 0.0; \
  _p_->Conventions = NULL; \
  _p_->ConventionVersion = NULL; \
  _p_->R = NULL

#define FREE_netcdfTrajectoryInfo( _p_ ) \
  safe_free(_p_->timeUnits); \
  safe_free(_p_->coordinateUnits); \
  safe_free(_p_->cellLengthUnits); \
  safe_free(_p_->cellAngleUnits); \
  safe_free(_p_->velocityUnits); \
  safe_free(_p_->Conventions); \
  safe_free(_p_->ConventionVersion); \
  safe_free(_p_->R); \
  safe_free(_p_)

#define AMBER_NETCDF_FRAME "frame"
#define AMBER_NETCDF_SPATIAL "spatial"
#define AMBER_NETCDF_ATOM "atom"
#define AMBER_NETCDF_CELL_SPATIAL "cell_spatial"
#define AMBER_NETCDF_CELL_ANGULAR "cell_angular"
#define AMBER_NETCDF_COORDS "coordinates"
#define AMBER_NETCDF_TIME "time"
#define AMBER_NETCDF_LABEL "label"
#define AMBER_NETCDF_LABELLEN 5

void netcdfPutAttributeText( int ncid, int vid, char *attribute, char *text )
{
  int err;

  err = nc_put_att_text(ncid, vid, attribute, strlen(text), text);
  if (err != NC_NOERR) {
    fprintf(stdout, "netcdfPutAttributeText: Error putting attribute (%s): %s\n", 
            attribute, nc_strerror(err));
  }
}

/*
 * Open a netCDF file for writing, returning the ncid.
 *
 * if bAppend is TRUE, we will append if the file exists.
 */
netcdfTrajectoryInfo * open_netcdf_file(char * filename, bool bAppend, int natoms, char * title, bool bPeriodic)
{
  int status; /* status flag returned by netCDF operations */
  netcdfTrajectoryInfo * NCInfo; /* storage for netCDF trajectory info */
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[3], count[3];
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
		   'b', 'e', 't', 'a', ' ', 
		   'g', 'a', 'm', 'm', 'a' };
  int i, j, err;
  float time;


  /* Allocate storage for netCDF file info. */
  NCInfo = (netcdfTrajectoryInfo *) malloc(sizeof(netcdfTrajectoryInfo));
  INITIALIZE_netcdfTrajectoryInfo(NCInfo);
  
  if(bAppend) {
    /* Attempt to open an existing dataset for writing. */
    status = nc_open(filename, NC_WRITE, &NCInfo->ncid);
  }
  if(!bAppend || (status != NC_NOERR)) {
    /* Attempt to open new file for writing. */
    status = nc_create(filename, NC_64BIT_OFFSET, &NCInfo->ncid);

    /* Die if there is an error. */
    if (status != NC_NOERR) {
      fprintf(stdout, "Error opening netCDF output coordinate file (%s): %s\n", filename, nc_strerror(status));
      exit(1);
    }
  }

  /* Create dimensions. */
  status = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_FRAME, NC_UNLIMITED, &NCInfo->frameDID);
  if (status != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF frame dimension ID: %s\n",
	    nc_strerror(status));
  
  status = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_SPATIAL, 3, &NCInfo->spatialDID);
  if (status != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF spatial dimension ID: %s\n",
	    nc_strerror(status));
  
  status = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_ATOM, natoms, &NCInfo->atomDID);
  if (status != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF atom dimension ID: %s\n",
	    nc_strerror(status));
  
  status = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_LABEL, AMBER_NETCDF_LABELLEN, &NCInfo->labelDID);
  if (status != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF label dimension ID: %s\n",
	    nc_strerror(status));
  
  status = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, 3, &NCInfo->cell_spatialDID);
  if (status != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF cell_spatial dimension ID: %s\n",
	    nc_strerror(status));
  
  status = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, 3, &NCInfo->cell_angularDID);
  if (status != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF cell_angular dimension ID: %s\n",
	    nc_strerror(status));

  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "title", title);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "application", "gromacs");
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "program", "xtc2netcdf");
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "programVersion", "1.0");
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "Conventions", "AMBER");
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "ConventionVersion", "1.0");

  dimensionID[0] = NCInfo->spatialDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->spatialVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining spatial VID: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->frameDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_TIME, NC_FLOAT, 1, dimensionID, &NCInfo->timeVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining time VID: %s\n", nc_strerror(err));
  netcdfPutAttributeText(NCInfo->ncid, NCInfo->timeVID, "units", "picosecond");
  if (err != NC_NOERR) fprintf(stdout, "NetCDF putting time VID units: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->frameDID;
  dimensionID[1] = NCInfo->atomDID;
  dimensionID[2] = NCInfo->spatialDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_COORDS, NC_FLOAT, 3, dimensionID, &NCInfo->coordinateVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining coordinate VID: %s\n", nc_strerror(err));
  netcdfPutAttributeText(NCInfo->ncid, NCInfo->coordinateVID, "units", "angstrom");
  if (err != NC_NOERR) fprintf(stdout, "NetCDF putting coordinate VID units: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->cell_spatialDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->cell_spatialVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell spatial VID: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->cell_angularDID;
  dimensionID[1] = NCInfo->labelDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &NCInfo->cell_angularVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell spatial VID: %s\n", nc_strerror(err));
  
  if(bPeriodic) {    
    dimensionID[0] = NCInfo->frameDID;
    dimensionID[1] = NCInfo->cell_spatialDID;
    err = nc_def_var(NCInfo->ncid, "cell_lengths", NC_DOUBLE, 2, dimensionID, &NCInfo->cellLengthVID);
    if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell length VID: %s\n", nc_strerror(err));
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellLengthVID, "units", "angstrom");
    if (err != NC_NOERR) fprintf(stdout, "NetCDF putting cell length VID units: %s\n", nc_strerror(err));
    
    dimensionID[1] = NCInfo->cell_angularDID;
    err = nc_def_var(NCInfo->ncid, "cell_angles", NC_DOUBLE, 2, dimensionID, &NCInfo->cellAngleVID);
    if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell angle VID: %s\n", nc_strerror(err));
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellAngleVID, "units", "degree");
    
    if (err != NC_NOERR) fprintf(stdout, "NetCDF putting cell angle VID units: %s\n", nc_strerror(err));
  }

  err = nc_set_fill(NCInfo->ncid, NC_NOFILL, &i);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF setting fill value: %s\n", nc_strerror(err));
  
  err = nc_enddef(NCInfo->ncid);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF error on ending definitions: %s\n", nc_strerror(err));
  
  start[0] = 0;
  count[0] = 3;
  
  xyz[0] = 'x';
  xyz[1] = 'y';
  xyz[2] = 'z';
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->spatialVID, start, count, xyz);
  if (err != NC_NOERR) 
    fprintf(stdout, "Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n", 
	    nc_strerror(err));
  
  xyz[0] = 'a';
  xyz[1] = 'b';
  xyz[2] = 'c';
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_spatialVID, start, count, xyz);
  if (err != NC_NOERR) 
    fprintf(stdout, "Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n", 
	    nc_strerror(err));
  
  
  start[0] = 0;
  start[1] = 0;
  count[0] = 3;
  count[1] = 5;
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_angularVID, start, count, abc);
  if (err != NC_NOERR) 
    fprintf(stdout, "Error on NetCDF output of angular VID 'alpha', 'beta ' and 'gamma': %s\n", 
	    nc_strerror(err));

  /* Return the netCDF file info. */
  return NCInfo;
}


/*
 * Close a netCDF file.
 */
void close_netcdf_file(netcdfTrajectoryInfo * NCInfo)
{
  int status; /* status flag returned by netCDF operations */

  /* Close the netCDF file. */
  status = nc_close(NCInfo->ncid);

  /* Die if error encountered. */
  if(status != NC_NOERR) {
    fprintf(stdout, "Error closing NetCDF file (nc_close), error: %s\n", nc_strerror(status));
    exit(1);
  }

  /* Free memory. */
  free(NCInfo);
}

void netcdf_write_frame(netcdfTrajectoryInfo * NCInfo, int natoms, rvec * x)
{
  static size_t start[3], count[3];
  int status;
  float * array;
  int i, d, index;

  /* Write atom coordinates */
  start[0] = 0;
  start[1] = 0;
  start[2] = NCInfo->currentFrame;
  
  count[0] = 3;
  count[1] = natoms;
  count[2] = 1;

  array = (float *)malloc(natoms * DIM * sizeof(float));

  index = 0;
  for(i=0; i<natoms; i++)
    for(d=0; d<DIM; d++)
      array[index++] = x[i][d];
  
  status = nc_put_vara_float(NCInfo->ncid, NCInfo->coordinateVID, start, count, array);
  
  free(array);
  
  NCInfo->currentFrame++;
}

void transcode_xtc(char *fn, char * netcdf_filename, bool bAppend)
{
  int    xd,indent;
  char   buf[256];
  rvec   *x;
  matrix box;
  int    nframe,natoms,step;
  real   prec,time;
  bool   bOK;
  int ncid; /* netCDF file id */
  netcdfTrajectoryInfo * NCInfo;
  bool   bPeriodic;

  xd = open_xtc(fn,"r");
  read_first_xtc(xd,&natoms,&step,&time,box,&x,&prec,&bOK);

  /* Open a netCDF file for writing. */
  NCInfo = open_netcdf_file(netcdf_filename, bAppend, natoms, "Transcoded from gromacs xtc trajectory.", TRUE);

  /* Read and write subsequent frames. */
  nframe=0;
  do {
    /*
    if (bXVG) {
      int i,d;
      
      fprintf(stdout,"%g",time);
      for(i=0; i<natoms; i++)
	for(d=0; d<DIM; d++)
	  fprintf(stdout," %g",x[i][d]);
      fprintf(stdout,"\n");
    } else {
    */
    sprintf(buf,"%s frame %d",fn,nframe);
    indent=0;
    indent=pr_title(stdout,indent,buf);
    pr_indent(stdout,indent);
    fprintf(stdout,"natoms=%10d  step=%10d  time=%10g  prec=%10g\n",
	    natoms,step,time,prec);

    /* TODO: Worry about precision. */
    /* status = nc_put_vara_float(NCInfo->ncid, NCInfo->cellLengthVID, start, count, box); */
    
    netcdf_write_frame(NCInfo, natoms, x);

    /* 
    pr_rvecs(stdout,indent,"box",box,DIM);
    pr_rvecs(stdout,indent,"x",x,natoms);
    */

    nframe++;
  } while (read_next_xtc(xd,natoms,&step,&time,box,x,&prec,&bOK) && nframe<2);
  if (!bOK)
    fprintf(stderr,"\nWARNING: Incomplete frame at time %g\n",time);
  close_xtc(xd);

  /* Close netCDF file. */
  close_netcdf_file(NCInfo);
}

void transcode_trx(char *fn, char * netcdf_filename, bool bAppend)
{
  int ftp; /* gromacs trajectory file id */
  
  ftp = fn2ftp(fn);
  if (ftp == efXTC)
    transcode_xtc(fn, netcdf_filename, bAppend);
  else
    fprintf(stderr,"File %s not supported.\n", fn);
}

/*
 * Write an AMBER format netCDF trajectory.
 */

  /*
   ! Create netCDF file
   cmode = nf90_64bit_offset
   exists = .false.
   if (owrite == 'N') cmode = ior(cmode, nf90_noclobber)
   if (owrite == 'U' .and. facc == 'A') cmode = ior(cmode, nf90_noclobber)
   status = nf90_create(path=filename,cmode=cmode,ncid=ncid)
   frame = 1
   if (status == nf90_eexist) then
      exists = .true.
      status = nf90_open(path = filename, mode = nf90_write, ncid = ncid)
      if (status == nf90_noerr) then !else fall through to error handling below and exit
         call checkerror(nf90_inq_dimid(ncid, "frame",FrameDimID), "query Frame ID")
         call checkerror(nf90_Inquire_Dimension(ncid,FrameDimID,name,frame),"query current frame")
         frame=frame + 1
         call checkerror(nf90_inq_varid(ncid, "time", TimeVarID), &
               "query time ID")         
         return 
      end if
   end if
   call checkerror(status, "create file")
   if (status /= nf90_noerr) then
#ifndef DUMB
      write (0,*) 'Error on opening ', filename
#endif
      call mexit(6,1)
   end if

   ! Define dimensions
   call checkerror(nf90_def_dim(ncid,"frame", nf90_unlimited, FrameDimID))
   call checkerror(nf90_def_dim(ncid,"spatial", 3, SpatialDimID))
   call checkerror(nf90_def_dim(ncid,"atom", natom, AtomDimID))
   call checkerror(nf90_def_dim(ncid,"label", NCLABELLEN, LabelDimID))


   ! Set global attributes
   call checkerror(nf90_put_att(ncid,nf90_global, "title", &
         title), "define title")
   call checkerror(nf90_put_att(ncid,nf90_global, "application", &
         'AMBER'), "define application")
   call checkerror(nf90_put_att(ncid,nf90_global, "program", &
         'sander'), "define program")
   call checkerror(nf90_put_att(ncid,nf90_global, "programVersion", &
         '9.0'), "define programVersion")
   call checkerror(nf90_put_att(ncid,nf90_global, "Conventions", &
         'AMBER'), "define Convention")
   call checkerror(nf90_put_att(ncid,nf90_global, "ConventionVersion", &
         '1.0'), "define ConventionVersion")

   ! Define non-optional variables
   call checkerror(nf90_def_var(ncid, "spatial", nf90_char, &
         (/  SpatialDimID /), SpatialVarID))
   
   call checkerror(nf90_def_var(ncid, "time", nf90_float, &
         (/  FrameDimID /), TimeVarID))
   call checkerror(nf90_put_att(ncid, TimeVarID, "units", &
         "picosecond"), "define time units")
  */

  /*   
   ! Define coordinate variables
   if (ntwx > 0 ) then
      call open_nc_file(crd_ncid,mdcrd,owrite,facc,atomCnt,title, &
            crd_file_exists,crd_frame,crd_TimeVarID)
      if (crd_file_exists) then
         call checkerror(nf90_inq_varid(crd_ncid, "coordinates", CoordVarID), &
               "query coordinates")
      else
         call checkerror(nf90_def_var(crd_ncid, "coordinates", nf90_float, &
               (/ SpatialDimID, AtomDimID, FrameDimID /), CoordVarID), "define coordinates")
         call checkerror(nf90_put_att(crd_ncid, CoordVarID, "units", &
               "angstrom"), "define coordinates units")
      end if
      ! Define unit cell data for periodic simulations
      if (ntb > 0 ) then
         if (crd_file_exists) then
            call checkerror(nf90_inq_varid(crd_ncid, "cell_lengths", Cell_lengthVarID), &
                  "query cell lengths")
            call checkerror(nf90_inq_varid(crd_ncid, "cell_angles",  Cell_angleVarID), &
                  "query cell angles")
         else
            ! Dimensions
            call checkerror(nf90_def_dim(crd_ncid,"cell_spatial", 3, Cell_spatialDimID), &
                  "define cell spatial dim")
            call checkerror(nf90_def_dim(crd_ncid,"cell_angular", 3, Cell_angularDimID), &
                  "define cell angular dim")
            ! Label variables
            call checkerror(nf90_def_var(crd_ncid, "cell_spatial", nf90_char, &
                  (/  Cell_spatialDimID /), Cell_spatialVarID),"define cell spatial var" )
            call checkerror(nf90_def_var(crd_ncid, "cell_angular", nf90_char, &
                  (/ LabelDimID, Cell_angularDimID /), Cell_angularVarID), "define cell angular var")
            ! Data variables and attributes
            call checkerror(nf90_def_var(crd_ncid, "cell_lengths", nf90_double, &
                  (/ Cell_spatialDimID, FrameDimID /), Cell_lengthVarID), "define cell lengths")
            call checkerror(nf90_put_att(crd_ncid, Cell_lengthVarID, "units", &
                  "angstrom"), "define cell length units")
            call checkerror(nf90_def_var(crd_ncid, "cell_angles", nf90_double, &
                  (/ Cell_angularDimID, FrameDimID /), Cell_angleVarID), "define cell angles")
            call checkerror(nf90_put_att(crd_ncid, Cell_angleVarID, "units", &
                  "degree"), "define cell angle units")
         end if
      end if
   end if



   ! Define velocity variables in file pointed to by vel_ncid
   if (ntwv /= 0 ) then
      if (vel_file_exists) then
         call checkerror(nf90_inq_varid(vel_ncid, "velocities", VelocVarID), &
               "query velocities")
      else
         call checkerror(nf90_def_var(vel_ncid, "velocities", nf90_float, &
               (/ SpatialDimID, AtomDimID, FrameDimID /), VelocVarID), &
               "define velocities")
         call checkerror(nf90_put_att(vel_ncid, VelocVarID, "units", &
               "angstrom/picosecond"), "define velocity units")
         call checkerror(nf90_put_att(vel_ncid, VelocVarID, "scale_factor", &
               20.455), "define velocity scale_factor")
      end if
   end if


   ! Prepare files for data writing
   if (ntwx > 0 .and. .not. crd_file_exists) then 
      call end_nc_define(crd_ncid)
      ! Fill dimension label variables
      if (ntb > 0) then
         call checkerror(nf90_put_var(crd_ncid, Cell_spatialVarID, &
	 (/ 'a','b','c' /), start = (/ 1 /), count = (/ 3 /)), &
               "write spatial variable")
         call checkerror(nf90_put_var(crd_ncid, Cell_angularVarID, &
               (/ 'alpha','beta ','gamma' /), &
               start = (/ 1, 1 /), count = (/ NCLABELLEN, 3 /)), &
               "write spatial variable")
      end if
   end if
   if (ntwv > 0 .and. .not. vel_file_exists) call end_nc_define(vel_ncid)

   ! Write frame.			      
   select case (unit)
      case (MDCRD_UNIT)
         if (n == 3 .and. ntb > 0) then       ! Assume this is box (fails on one atom systems)
            call checkerror(nf90_put_var(crd_ncid,Cell_lengthVarID,(/ a,b,c/), &
                  start = (/ 1, crd_frame /), count = (/ 3, 1 /)), 'write cell lengths')
            call checkerror(nf90_put_var(crd_ncid,Cell_angleVarID,(/ alpha,beta,gamma /), &
                  start = (/ 1, crd_frame /), count = (/ 3, 1 /)), 'write cell angles')
         else
            call checkerror(nf90_put_var(crd_ncid,CoordVarID, r(istart:n), &
                  start = (/ 1, 1, crd_frame /),count = (/ 3, (n-istart+1)/3, 1 /)), &
                  'write atom coords')
         end if
      case (MDVEL_UNIT)
         call checkerror(nf90_put_var(vel_ncid,VelocVarID, r(istart:n), &
               start = (/ 1, 1, vel_frame /), count = (/ 3, (n-istart+1)/3, 1 /)), &
               'write velocities')
      case default
         write (6,*) 'Error: unhandled unit ',unit,' selected for output in bintraj'
   end select

   ! close
   if (ntwx > 0) call checkerror(nf90_close(crd_ncid),"close mdcrd")   
   if (ntwv > 0) call checkerror(nf90_close(vel_ncid),"close mdvel")
  */

  /*
  NCInfo = (netcdfTrajectoryInfo *) safe_malloc(sizeof(netcdfTrajectoryInfo));
  INITIALIZE_netcdfTrajectoryInfo(NCInfo);
  
  err = nc_create(outInfo->filename, NC_64BIT_OFFSET, &NCInfo->ncid);
  if (err != NC_NOERR) {
    fprintf(stdout, "WARNING in ptrajOutputCoordinates(): Error opening\n");
    fprintf(stdout, "NetCDF output coordinate file (%s): %s\n", 
	    outInfo->filename, nc_strerror(err));
    exit(1);
  }

  err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_FRAME, NC_UNLIMITED, &NCInfo->frameDID);
  if (err != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF frame dimension ID: %s\n",
	    nc_strerror(err));
  
  err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_SPATIAL, 3, &NCInfo->spatialDID);
  if (err != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF spatial dimension ID: %s\n",
	    nc_strerror(err));
  
  err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_ATOM, atoms, &NCInfo->atomDID);
  if (err != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF atom dimension ID: %s\n",
	    nc_strerror(err));
  
  err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_LABEL, AMBER_NETCDF_LABELLEN, &NCInfo->labelDID);
  if (err != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF label dimension ID: %s\n",
	    nc_strerror(err));
  
  err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, 3, &NCInfo->cell_spatialDID);
  if (err != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF cell_spatial dimension ID: %s\n",
	    nc_strerror(err));
  
  err = nc_def_dim(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, 3, &NCInfo->cell_angularDID);
  if (err != NC_NOERR)
    fprintf(stdout, "ptrajOutputCoordinates() defining NetCDF cell_angular dimension ID: %s\n",
	    nc_strerror(err));

  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "title", outInfo->title);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "application", outInfo->application);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "program", outInfo->program);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "programVersion", outInfo->version);
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "Conventions", "AMBER");
  netcdfPutAttributeText(NCInfo->ncid, NC_GLOBAL, "ConventionVersion", "1.0");
  
  dimensionID[0] = NCInfo->spatialDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->spatialVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining spatial VID: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->frameDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_TIME, NC_FLOAT, 1, dimensionID, &NCInfo->timeVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining time VID: %s\n", nc_strerror(err));
  netcdfPutAttributeText(NCInfo->ncid, NCInfo->timeVID, "units", "picosecond");
  if (err != NC_NOERR) fprintf(stdout, "NetCDF putting time VID units: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->frameDID;
  dimensionID[1] = NCInfo->atomDID;
  dimensionID[2] = NCInfo->spatialDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_COORDS, NC_FLOAT, 3, dimensionID, &NCInfo->coordinateVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining coordinate VID: %s\n", nc_strerror(err));
  netcdfPutAttributeText(NCInfo->ncid, NCInfo->coordinateVID, "units", "angstrom");
  if (err != NC_NOERR) fprintf(stdout, "NetCDF putting coordinate VID units: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->cell_spatialDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_CELL_SPATIAL, NC_CHAR, 1, dimensionID, &NCInfo->cell_spatialVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell spatial VID: %s\n", nc_strerror(err));
  
  dimensionID[0] = NCInfo->cell_angularDID;
  dimensionID[1] = NCInfo->labelDID;
  err = nc_def_var(NCInfo->ncid, AMBER_NETCDF_CELL_ANGULAR, NC_CHAR, 2, dimensionID, &NCInfo->cell_angularVID);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell spatial VID: %s\n", nc_strerror(err));
  
  if (outInfo->isBox) {    
    dimensionID[0] = NCInfo->frameDID;
    dimensionID[1] = NCInfo->cell_spatialDID;
    err = nc_def_var(NCInfo->ncid, "cell_lengths", NC_DOUBLE, 2, dimensionID, &NCInfo->cellLengthVID);
    if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell length VID: %s\n", nc_strerror(err));
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellLengthVID, "units", "angstrom");
    if (err != NC_NOERR) fprintf(stdout, "NetCDF putting cell length VID units: %s\n", nc_strerror(err));
    
    dimensionID[1] = NCInfo->cell_angularDID;
    err = nc_def_var(NCInfo->ncid, "cell_angles", NC_DOUBLE, 2, dimensionID, &NCInfo->cellAngleVID);
    if (err != NC_NOERR) fprintf(stdout, "NetCDF defining cell angle VID: %s\n", nc_strerror(err));
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->cellAngleVID, "units", "degree");
    
    if (err != NC_NOERR) fprintf(stdout, "NetCDF putting cell angle VID units: %s\n", nc_strerror(err));
  }

  err = nc_set_fill(NCInfo->ncid, NC_NOFILL, &i);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF setting fill value: %s\n", nc_strerror(err));
  
  err = nc_enddef(NCInfo->ncid);
  if (err != NC_NOERR) fprintf(stdout, "NetCDF error on ending definitions: %s\n", nc_strerror(err));
  
  start[0] = 0;
  count[0] = 3;
  
  xyz[0] = 'x';
  xyz[1] = 'y';
  xyz[2] = 'z';
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->spatialVID, start, count, xyz);
  if (err != NC_NOERR) 
    fprintf(stdout, "Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n", 
	    nc_strerror(err));
  
  xyz[0] = 'a';
  xyz[1] = 'b';
  xyz[2] = 'c';
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_spatialVID, start, count, xyz);
  if (err != NC_NOERR) 
    fprintf(stdout, "Error on NetCDF output of spatial VID 'x', 'y' and 'z': %s\n", 
	    nc_strerror(err));
  
  
  start[0] = 0;
  start[1] = 0;
  count[0] = 3;
  count[1] = 5;
  err = nc_put_vara_text(NCInfo->ncid, NCInfo->cell_angularVID, start, count, abc);
  if (err != NC_NOERR) 
    fprintf(stdout, "Error on NetCDF output of angular VID 'alpha', 'beta ' and 'gamma': %s\n", 
	    nc_strerror(err));
  
  outInfo->info = (void *) NCInfo;
*/

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gmx2netcdf reads a trajectory ([TT].trj[tt]/[TT].trr[tt]/[TT].xtc[tt]) and writes an AMBER format netCDF trajectory.[PAR]"
  };
  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTRJ, "-o", NULL, ffWRITE },
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bAppend=TRUE;
  t_pargs pa[] = {
    { "-O", FALSE, etBOOL, {&bAppend}, "Overwrite (default is to append)" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  /* Transcode the trajectory to netCDF. */
  transcode_trx(ftp2fn(efTRX,NFILE,fnm), ftp2fn(efTRJ,NFILE,fnm), bAppend);
    
  thanx(stderr);

  return 0;
}
