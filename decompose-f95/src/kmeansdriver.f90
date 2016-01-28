!=============================================================================================
! Main driver program for K-means clustering.
!=============================================================================================
program main
  use decompose ! for K-means clustering functions
  use utilities, only : setRandomSeed ! to set the random seed from the current time
  use constants ! global constants
  use utilities, only : find, unique
  use timer ! for timing info

  implicit none

  ! Constants.
  character(len=*), parameter :: version = '$Id: kmeansdriver.f90,v 1.1 2007/04/04 03:02:49 jchodera Exp $' ! version string

  ! Local variables.
  integer :: narguments
  character(len=MAX_FILENAME_LENGTH) :: filename
  character(len=MAX_FILENAME_LENGTH) :: string

  ! Local variables.
  integer :: macrostate, compacted_macrostate_index
    ! macrostate index
  integer, dimension(:), pointer :: indices
    ! snapshot indices of macrostate members
  integer :: nindices
    ! number of members of the macrostate
  integer :: nmembers
    ! Number of members of this state.
  integer :: nmaxmicrostates
    ! Maximum possible number of microstates.
  integer :: target_nmicrostates
    ! Target number of microstates to split to.
  integer :: status
  integer :: counter

  ! Write copyright notice.
  write(*,*) 'decompose-f90: K-means clustering driver'
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
  if(narguments .ne. 2) then
     write(*,*) ''
     write(*,*) 'usage: ./kmeans control.xml nstates'
     write(*,*) ''
     stop
  end if
  ! First argument is control filename.
  call get_command_argument(1, value=filename)
  ! Second argument is number of microstates to split to.
  call get_command_argument(2, value=string)
  read(string,'(I8)') target_nmicrostates
  
  ! Initialize.
  call initialize(trim(filename), nocheck = .TRUE.)
  write(*,*) macrostate_assignments

  !=============================================================================================
  ! Split to microstates.
  !=============================================================================================
  
  ! Determine the distance metric to use for this iteration.
  rmsd_mode = rmsd_mode_first_iteration
  if(debug) write(*,*) 'Using rmsd_mode = ', rmsd_mode

  ! Perform K-means clustering.
  call resetTimer()

  ! Get snapshot indices of all snapshots belonging to microstate 1.
  indices => find(macrostate_assignments == 1)

  ! Determine number of snapshots to be clustered.
  nindices = size(indices,1)
  write(*,*) '======================================================================'
  write(*,*) 'Will cluster ', nindices, ' snapshots into ', target_nmicrostates, ' states.'
  write(*,*) '======================================================================'

  ! Copy existing state assignments over so that states that will not be split remain unchanged.
  microstate_assignments = macrostate_assignments

  ! Set all state representatives to zero.
  microstates = 0

  ! Split the macrostate into microstates, assigning them unique ids by their generator.
  call splitState(nindices, indices, target_nmicrostates, microstate_assignments)

  ! Free list of macrostate members.
  deallocate(indices)

  ! Generate compacted microstate assignments, where indices go from 1...nmicrostates.
  microstates => unique(microstate_assignments, compacted=compacted_microstate_assignments)
  nmicrostates = size(microstates,1)    
  write(*,*) 'there are ', nmicrostates, ' unique microstates'
  write(*,*) 'checking microstate assignments...'
  call checkstates(compacted_microstate_assignments)
  if(debug) write(*,*) readTimer(), ' seconds elapsed'

  ! Write state assignments.
  write(filename,*) iteration
  filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.stateassignments-split'
  call writeStateassignments(filename, nmicrostates, microstate_assignments)
  write(filename,*) iteration
  filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.stateassignments-split-compacted'
  call writeStateassignments(filename, nmicrostates, compacted_microstate_assignments)

  ! Write list of exemplars.
  write(filename,*) iteration
  filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.microstates'
  call writeStateassignments(filename, nmicrostates, microstates)

  ! Write state PDBs.
  call resetTimer()
  write(filename,*) iteration
  filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.ensemble-split'
  call writeStatePDBs(filename, nmicrostates, microstates, microstate_assignments, split_pdb_nmodels)
  if(debug) write(*,*) readTimer(), ' seconds elapsed'       

  ! Clean up.
  call finalize()

  ! Terminate
  stop

end program main

