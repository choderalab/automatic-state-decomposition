!=============================================================================================
! A module for reading snapshot data from AMBER trajectory netCDF files.
!
! The AMBER trajectory netCDF convention is documented here:
! http://amber.scripps.edu/netcdf/nctraj.html
!
! The code for reading data from the AMBER format netCDF file is based upon the code appearing
!  in the AMBER9 bintraj.f file, written by John Mongan, though no code is copied directly.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
! - Make thread-safe by returning a type(netcdfdata) containing all necessary information,
!   instead of storing as module data?
!=============================================================================================
module netcdfio
  use numeric_kinds ! for precision
  use constants ! for global constants
  use netcdf ! for netCDF file access  

  implicit none

  private
  public :: netcdf_open, netcdf_close, netcdf_read
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  integer, parameter :: NCLABELLEN = 5
  integer, parameter :: MAX_ATTRIBUTE_LEN = 256 ! maximum netCDF character string attribute length

  logical, parameter :: debug = .TRUE.

  !=============================================================================================
  ! Data types.
  !=============================================================================================
  type, public :: netcdfio_t
     character(len=MAX_FILENAME_LENGTH) :: filename ! the filename opened
     
     integer :: ncid ! netCDF dataset ids
     integer :: FrameDimID, SpatialDimID, AtomDimID ! netCDF dimension ids
     integer :: CoordVarID, SpatialVarID ! netCDF variable ids
    
     integer :: nspatial
     integer :: natoms
     integer :: nframes
  end type netcdfio_t

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

contains

  !=============================================================================================
  ! Open a snapshot file and return the number of atoms and stored frames.
  !=============================================================================================    
  function netcdf_open(filename, natoms, nframes) result(nc)
    ! Parameters.
    character(*), intent(in) :: filename
      ! the filename of the netCDF file to open
    integer, intent(out), optional :: natoms
      ! returns the number of atoms, if desired
    integer, intent(out), optional :: nframes
      ! returns the number of frames, if desired

    ! Return value.
    type(netcdfio_t) :: nc
      ! netcdf file 'handle'

    ! Local variables.
    integer :: status 
      ! status result from netCDF library call
    character(len = MAX_ATTRIBUTE_LEN) :: attribute
      ! storage for attribute read from file    

    if(debug) write(*,*) 'Opening snapshotfile ', trim(filename), ' for reading...'

    ! Store filename.
    nc%filename = filename

    ! Open the netCDF file for read-only access and obtainin the netCDF dataset id.
    status = nf90_open(path = trim(nc%filename), mode = nf90_nowrite, ncid = nc%ncid)
    call checkerror(status, 'Attempting to open netCDF snapshot file ' // filename)

    ! Check global attributes.   
    status = nf90_get_att(nc%ncid, nf90_global, "Conventions", attribute)
    call checkerror(status, "Attempting to read attribute value for global attribute 'Conventions'")   
    ! TODO: Check to make sure "Conventions" includes the token "AMBER".
    status = nf90_get_att(nc%ncid, nf90_global, "ConventionVersion", attribute)
    call checkerror(status, "Attempting to read attribute value for global attribute 'ConventionVersion'")  
    ! TODO: Check to make sure this version is "1.0".

    ! Get the netCDF id of the coordinates variable.
    status = nf90_inq_varid(nc%ncid, "coordinates", nc%CoordVarID)
    call checkerror(status, "Attempting to obtain id for 'coordinates'")   

    ! Get the NetCDF dimension ids for the spatial, atom, and frame dimensions.
    status = nf90_inq_dimid(nc%ncid, "spatial", nc%SpatialDimID)
    call checkerror(status, "Attempting to obtain id for 'spatial' dimension")   
    status = nf90_inq_dimid(nc%ncid, "atom", nc%AtomDimID)
    call checkerror(status, "Attempting to obtain id for 'atom' dimension")   
    status = nf90_inq_dimid(nc%ncid, "frame", nc%FrameDimID)
    call checkerror(status, "Attempting to obtain id for 'frame' dimension")   
    
    ! Get the dimension sizes.
    status = nf90_inquire_dimension(nc%ncid, nc%SpatialDimID, len = nc%nspatial)
    call checkerror(status, "Attempting to determine number of 'spatial' dimensions present")   
    status = nf90_inquire_dimension(nc%ncid, nc%AtomDimID, len = nc%natoms)
    call checkerror(status, "Attempting to determine number of atoms present")   
    status = nf90_inquire_dimension(nc%ncid, nc%FrameDimID, len = nc%nframes)
    call checkerror(status, "Attempting to determine number of frames present")   

    if(debug) write(*,*) 'nframes = ', nc%nframes, ', natoms = ', nc%natoms, ', nspatial ', nc%nspatial
    
    ! Return number of atoms and frames.
    if(present(natoms)) natoms = nc%natoms
    if(present(nframes)) nframes = nc%nframes

  end function netcdf_open

  !=============================================================================================
  ! Read snapshots from snapshotfile.
  !=============================================================================================
  subroutine netcdf_read(nc, snapshots)
    ! Parameters.
    type(netcdfio_t), intent(in) :: nc
      ! netcdf file 'handle'
    real(sp), dimension(nc%nframes, nc%natoms, nc%nspatial), intent(inout) :: snapshots
      ! Storage for the snapshots.
      ! snapshots(i,n,k) is coordinate k of atom n of snapshot i.

    ! Local variables.
    integer :: status ! status result from netCDF library call    
    
    if(debug) write(*,*) 'Reading ', nc%nframes, ' snapshots...'
    if(debug) write(*,*) 'Dimensions of snapshots: ', nc%nframes, nc%natoms, nc%nspatial
    if(debug) write(*,*) 'Dimensions of snapshots: ', shape(snapshots)

    

    ! Read all coordinates to snapshot store.
    status = nf90_get_var(ncid=nc%ncid, varid=nc%CoordVarID, start=(/1,1,1/), &
         count=(/ nc%nspatial, nc%natoms, nc%nframes/), &
         values=snapshots, map =(/ nc%nframes*nc%natoms, nc%nframes, 1 /))
    call checkerror(status, "Attempting to read coordinate snapshots.")   

    if(debug) write(*,*) 'Done.'

  end subroutine netcdf_read

  !=============================================================================================
  ! Close the snapshotfile.
  !=============================================================================================  
  subroutine netcdf_close(nc)    
    ! Parameters.
    type(netcdfio_t), intent(in) :: nc
      ! netcdf file 'handle'

    ! Local variables.
    integer :: status ! status result from netCDF library call    

    if(debug) write(*,*) 'Closing snapshotfile...'
    
    ! Close snapshot file.
    status = nf90_close(ncid=nc%ncid)
    call checkerror(status, "Closing snapshotfile.")

    if(debug) write(*,*) 'Closed.'
    
  end subroutine netcdf_close

  !=============================================================================================
  ! Check whether netCDF library call was successful.  
  ! If not, throw an appropriate error message and terminate.
  !=============================================================================================
  subroutine checkerror(status, message)
    
    ! Parameters.
    integer, intent(in) :: status 
      ! status returned by call to netCDF.
    character(*), optional, intent(in) :: message 
      ! error message to report

    ! If status reports an error, report the netCDF error message and supplied message.
    if(status /= nf90_noerr) then
       ! Report netCDF error message.
       write (6,*) 'NetCDF error: ', trim(nf90_strerror(status))
       
       ! Display developer message, if provided.
       if (present(message)) then
          write (6,*) '  Encountered during: ', message
       end if
       
       ! Halt execution.
       stop
    end if
    
  end subroutine checkerror
  
end module netcdfio
