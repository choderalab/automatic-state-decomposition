!=============================================================================================
! Main driver program for state decomposition algorithm.
!=============================================================================================
program main
  use decompose ! for state space decomposition
  use utilities, only : setRandomSeed ! to set the random seed from the current time
  use constants ! global constants

  implicit none

  ! Constants.
  character(len=*), parameter :: version = '$Id: decomposedriver.f90,v 1.7 2006/11/04 08:32:01 jchodera Exp $' ! version string

  ! Local variables.
  integer :: narguments
  character(len=MAX_FILENAME_LENGTH) :: filename

  ! Write copyright notice.
  write(*,*) 'decompose-f90'
  write(*,*) version
  write(*,*) 'Copyright (C) The Regents of the University of California'
  write(*,*) 'Written by John D. Chodera'
  write(*,*) 'This software comes with ABSOLUTELY NO WARRANTY.'
  write(*,*) 'This is free software, and you are welcome to redistribute it'
  write(*,*) 'under certain conditions.'
  write(*,*) 'See the included GPL license for details.'
  
  ! Initialize random number generator.
  call setRandomSeed()

  ! Get control XML filename from command line.
  narguments = command_argument_count()
  if(narguments .ne. 1) then
     write(*,*) ''
     write(*,*) 'usage: ./decompose control.xml'
     write(*,*) ''
     stop
  end if
  call get_command_argument(1, value=filename)

  ! Initialize.
  call initialize(trim(filename))

  ! Run algorithm.
  call iterate()

  ! Clean up.
  call finalize()

  ! Terminate
  stop

end program main

