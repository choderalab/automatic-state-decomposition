!=============================================================================================
! Driver to compute timescales from existing state decomposition.
!=============================================================================================
program main
  use decompose ! for state space decomposition
  use utilities, only : setRandomSeed ! to set the random seed from the current time
  use constants ! global constants
  use utilities, only : unique

  implicit none

  ! Constants.
  character(len=*), parameter :: version = '$Id: timescalesdriver.f90,v 1.2 2007/05/29 08:47:23 jchodera Exp $' ! version string

  ! Local variables.
  integer :: narguments
  character(len=MAX_FILENAME_LENGTH) :: control_filename
  character(len=MAX_FILENAME_LENGTH) :: stateassignments_filename
  character(len=MAX_FILENAME_LENGTH) :: timescales_prefix
  character(len=MAX_FILENAME_LENGTH) :: string
!  integer, dimension(:), allocatable :: stateassignments
!  integer :: tau_step
!  integer :: max_tau

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
  if(narguments .ne. 5) then
     write(*,*) ''
     write(*,*) 'usage: ./compute-timescales control.xml stateassignments-file timescales-prefix max_tau tau_step'
     write(*,*) ''
     stop
  end if
  call get_command_argument(1, value=control_filename)
  call get_command_argument(2, value=stateassignments_filename)

  ! Initialize.
  call initialize(trim(control_filename))

  ! Read state assignments.
  call readStateAssignments(trim(stateassignments_filename), macrostate_assignments)
  
  ! Determine macrostates and compacted macrostates.
  macrostates => unique(macrostate_assignments, compacted=compacted_macrostate_assignments)
  nmacrostates = size(macrostates,1)
  if(debug) write(*,*) 'There are ', nmacrostates, ' macrostates:'
  if(debug) write(*,*) macrostates
  
  ! DEBUG
  call checkstates(macrostate_assignments)
  call checkstates(compacted_macrostate_assignments)

  ! Parse remainder of arguments.
  call get_command_argument(3, value=timescales_prefix)
  call get_command_argument(4, value=string)
  read(string,'(I8)') max_tau
  call get_command_argument(5, value=string)
  read(string,'(I8)') tau_step

  write(*,*) 'Computing timescales with max_tau = ', max_tau, ', tau_step = ', tau_step
  ! Write timescales.
  call writeTimescales(trim(timescales_prefix), nmacrostates, compacted_macrostate_assignments, max_tau, tau_step)

  ! Write eigenvectors.
  call writeEigenvectors(trim(timescales_prefix), nmacrostates, compacted_macrostate_assignments, tau, &
       min(max_ntimescales+1,nmacrostates))

  ! Clean up.
  call finalize()

  ! Terminate
  stop

end program main

