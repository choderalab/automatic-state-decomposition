!=============================================================================================
! A module for reading snapshot data from delimited ASCII files.
!
! Written by John D. Chodera <jchodera@gmail.com>, Pande lab, Stanford, 2007.

! Copyright (c) 2007 Stanford University?
! All Rights Reserved.
!=============================================================================================
! TODO:
!=============================================================================================
module delimitedio
  use numeric_kinds ! for precision
  use constants ! for global constants
  use pdbio, only : new_unit ! for automatic unit numbers

  implicit none

  private
  public :: delimited_open, delimited_read, delimited_close
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  logical, parameter :: debug = .FALSE. ! if .TRUE., extra debug information will be printed

  !=============================================================================================
  ! Data types.
  !=============================================================================================
  type, public :: delimited_t
     character(len=MAX_FILENAME_LENGTH) :: filename ! the filename opened
     integer :: unit ! unit number
     integer :: lines_read ! number of lines read
     integer :: natoms ! number of atoms per frame (total)
     integer :: nframes ! number of frames in files (total)
     integer :: nspatial ! number of spatial dimensions per atom
  end type delimited_t

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

contains

  !=============================================================================================
  ! Open a snapshot file and return the number of atoms and stored frames.
  !=============================================================================================    
  function delimited_open(filename, format_string, nspatial, natoms, nframes) result(file)
    ! Parameters.
    character(*), intent(in) :: filename
      ! the filename of the netCDF file to open
    character(format_string) :: format_string
      ! the format string describing the data format (e.g. '108F5.3')
    integer, intent(in) :: nspatial
      ! number of spatial dimensions per atom
    integer, intent(out), optional :: natoms
      ! returns the number of atoms, if desired
    integer, intent(out), optional :: nframes
      ! returns the number of frames, if desired

    ! Return value.
    type(delimitiedio_t) :: file
      ! metadata about the opened file

    ! Local variables.
    character(len=MAX_LINE_LENGTH) :: line
      ! an entire line from the PDB file
      ! Note that MAX_LINE_LENGTH is defined in the 'constants' module.
    real(sp), dimension(sp) :: parsed_line
    integer :: iostat

    if(debug) write(*,*) 'Opening delimited ASCII snapshotfile ', trim(filename), ' for reading...'

    ! Store filename.
    file%filename = trim(filename)

    ! Store number of spatial dimensions
    file%nspatial = nspatial

    ! Open file for reading, ensuring it exists.
    call new_unit(iunit)
    open(unit=iunit, file=filename, status='OLD', err=11, position='REWIND', &
         form='FORMATTED', action='READ')

    ! TODO: Determine natoms and nframes
    
    do
       ! Read a line.
       lineNumber = lineNumber + 1
       read(unit=iunit, fmt='(a)', iostat=iostat, end=20) line
       ! TODO: Terminate loop if EOF encountered?  Here or after parsing?
       ! if(iostat < 0) exit


    end do

    ! read(unit=iunit, fmt='(a)', iostat=iostat, end=20) line
    
    ! End of file has been reached.
20  continue
    
    ! Close file.
    close(iunit)

    if(debug) write(*,*) 'nframes = ', file%nframes, ', natoms = ', file%natoms, ', nspatial ', file%nspatial

    ! Return number of atoms and frames.
    if(present(natoms)) natoms = file%natoms
    if(present(nframes)) nframes = file%nframes

    ! Return control.
    return

    ! Report error if file could not be opened.
11  continue
    write(*,*) 'delimitedio: delimited_read: Error opening file ', filename
    stop

  end function delimited_open

  !=============================================================================================
  ! Read snapshots from snapshotfile.
  !=============================================================================================
  subroutine delimited_read(file, snapshots)
    ! Parameters.
    type(delimited_t), intent(in) :: file
      ! file metadata
    real(sp), dimension(file%nframes, file%natoms, file%nspatial) :: snapshots
      ! Storage for the snapshots.
      ! snapshots(i,n,k) is coordinate k of atom n of snapshot i.

    if(debug) write(*,*) 'Reading ', file%nframes, ' snapshots...'
    
    ! Read all coordinates to snapshot store.
    ! TODO
    ! read(unit=iunit, fmt='(a)', iostat=iostat, end=20) line

    if(debug) write(*,*) 'Done.'

  end subroutine delimited_read

  !=============================================================================================
  ! Close the snapshotfile.
  !=============================================================================================  
  subroutine delimited_close(file)    
    ! Parameters.
    type(delimited_t), intent(in) :: file
      ! file metadata

    ! Local variables.
    integer :: status ! status result from netCDF library call    

    if(debug) write(*,*) 'Closing snapshotfile...'
    
    ! Close snapshot file.
    close(file%unit)

    if(debug) write(*,*) 'Closed.'
    
  end subroutine delimited_close

end module delimitedio

