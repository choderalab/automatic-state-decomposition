!=============================================================================================
! A module for reading gromacs XTC trajectories.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
!=============================================================================================
module xtcio
  use numeric_kinds ! for precision
  use constants ! for global constants

  implicit none

  private
  public :: xtc_open
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================

  logical, parameter :: debug = .FALSE.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

  !=============================================================================================
  ! External C functions from gromacs library.
  !=============================================================================================
  integer, external :: open_xtc, read_first_xtc, read_next_xtc, close_xtc

contains

  !=============================================================================================
  ! Open a snapshot file and return the number of atoms and stored frames.
  !=============================================================================================    
  function snapshotfile_open(filename, natoms, nframes) result(handle)
    ! Parameters.
    character(*), intent(in) :: filename
      ! the filename of the netCDF file to open
    integer, intent(out), optional :: natoms
      ! returns the number of atoms, if desired
    integer, intent(out), optional :: nframes
      ! returns the number of frames, if desired

    ! Return values.
    integer :: handle
      ! an integer handle to the currently-open XTC file

    ! Local variables.
    integer :: current_frame
      ! index of current frame
    integer :: natoms_, step_
    real(sp) :: time_

    ! Open the XTC file (gromacs call).
    !
    ! gromacs call from kernel/gmxdump.c:
    !  xd = open_xtc(fn,"r");
    handle = open_xtc(filename, 'r')

    ! Initialize frame counter.
    current_frame = 0

    ! Read first frame from trajectory (gromacs call).
    !
    ! C declaration:
    !  int read_first_xtc(int fp,int *natoms,int *step,real *time, matrix box,rvec **x,real *prec,bool *bOK)
    ! gromacs call from kernel/gmxdump.c:
    !  read_first_xtc(xd,&natoms,&step,&time,box,&x,&prec,&bOK);
    read_first_xtc(handle, natoms, step, time, box, &x, prec, bOK);
    current_frame = current_frame + 1

    ! Retrieve frames.
    status = 1
    do while (status)
       ! Read next frame from trajectory (gromacs call):
       !
       ! gromacs call from kernel/gmxdump.c:
       !  read_next_xtc(xd,natoms,&step,&time,box,x,&prec,&bOK)
       status = read_next_xtc(handle, natoms, &step, &time, box, x, &prec, &bOK)
       current_frame = current_frame + 1
    end do

    ! Close the file and reopen it to return to beginning.
    status = close_xtc(handle)

    xtc_open = xd

    ! Return number of atoms and frames if desired
    if(present(natoms)) natoms = 1
    if(present(nframes)) nframes = current_frame

    
  end function snapshotfile_open

  !=============================================================================================
  ! Read snapshots from snapshotfile.
  !=============================================================================================
  subroutine snapshotfile_read(nc, snapshots)
    ! Parameters.
    type(netcdfio_t), intent(in) :: nc
      ! netcdf file 'handle'
    real(sp), dimension(nc%nframes, nc%natoms, nc%nspatial) :: snapshots
      ! Storage for the snapshots.
      ! snapshots(i,n,k) is coordinate k of atom n of snapshot i.

    ! Local variables.
    integer :: status ! status result from netCDF library call    
    
    if(debug) write(*,*) 'Reading ', nc%nframes, ' snapshots...'
    
    ! Read all coordinates to snapshot store.
    status = nf90_get_var(ncid=nc%ncid, varid=nc%CoordVarID, start=(/1,1,1/), &
         count=(/ nc%nspatial, nc%natoms, nc%nframes/), &
         values=snapshots, map =(/ nc%nframes*nc%natoms, nc%nframes, 1 /))
    call checkerror(status, "Attempting to read coordinate snapshots.")   

    if(debug) write(*,*) 'Done.'

  end subroutine snapshotfile_read

  !=============================================================================================
  ! Close the snapshotfile.
  !=============================================================================================  
  subroutine snapshotfile_close(nc)    
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
    
  end subroutine snapshotfile_close

  
end module xtcio

program x


end program x
