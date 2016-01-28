!=============================================================================================
! Driver to compute bootstrap sample of the transition matrix at a given lag time for a given decomposition.
!=============================================================================================
program main
  use decompose ! for state space decomposition
  use utilities, only : setRandomSeed ! to set the random seed from the current time
  use constants ! global constants
  use utilities, only : unique, trace, std, generateUniformIntegers, compute_confidence_bounds
  use pdbio, only : new_unit
  use transitionmatrix, only : computeTransitionMatrix, printTransitionMatrix
    
  implicit none

  ! Constants.
  character(len=*), parameter :: version = '$Id: transitionmatrixbootstrapdriver.f90,v 1.1 2006/11/07 08:00:19 jchodera Exp $' ! version string

  ! Local variables.
  integer :: narguments
  character(len=MAX_FILENAME_LENGTH) :: control_filename
  character(len=MAX_FILENAME_LENGTH) :: stateassignments_filename
  character(len=MAX_FILENAME_LENGTH) :: filename
  character(len=MAX_FILENAME_LENGTH) :: string
!  integer, dimension(:), allocatable :: stateassignments
!  integer :: tau_step
!  integer :: max_tau
  integer :: trial ! current bootstrap trial
  integer :: ntrials ! number of bootstrap trials to conduct / samples to generate
  !integer :: tau ! lag time to compute transition matrix for
  integer, dimension(:), allocatable :: trajectory_indices ! indices of trajectories in bootstrap sample
  real(sp), dimension(:,:), allocatable :: Tji ! transition matrix
  real(sp), dimension(:,:,:), allocatable :: Tjin ! transition matrix for bootstrap sample n
  integer :: iunit ! file unit for writing
  integer :: i, j

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
     write(*,*) 'usage: ./compute-transitionmatrix control.xml stateassignments-file transitionmatrix-file tau nsamples'
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
  call get_command_argument(3, value=filename)
  call get_command_argument(4, value=string)
  read(string,'(I8)') tau
  call get_command_argument(5, value=string)
  read(string,'(I8)') ntrials
  write(*,*) 'tau = ', tau, ', ntrials = ', ntrials

  ! allocate storage for bootstrap samples, taus, and metastabilities
  allocate(trajectory_indices(ntrajectories))
  allocate(Tji(nmacrostates,nmacrostates),Tjin(nmacrostates,nmacrostates,ntrials))

  ! Compute metastabilities and uncertainties by bootstrap.
  write(*,*) 'Performing ', ntrials, ' bootstrap trials...'
  do trial = 1, ntrials
     write(*,*) 'trial = ', trial, '/', ntrials
     
     ! Generate bootstrap sample of trajectory indices.
     call generateUniformIntegers(ntrajectories, trajectory_indices)
     
     ! Compute transition matrix.
     Tji = computeTransitionMatrix(ntrajectories, trajectories(trajectory_indices), nmacrostates, &
          compacted_macrostate_assignments, tau)

     !call printTransitionMatrix(Tji)

     ! Store bootstrapped sample.
     Tjin(:,:,trial) = Tji
       
  end do

  ! Write to file.
  if(debug) write(*,*) 'Writing to file ', trim(filename), '...'
  call new_unit(iunit)
  open(unit=iunit, file=trim(filename), status='REPLACE', err=10, position='REWIND', &
       form='FORMATTED', action='WRITE')    
  do trial = 1,ntrials
     do j = 1,nmacrostates
        do i = 1,nmacrostates
           ! Write line.
           write(iunit,'(I5,X,I4,X,I4,X,F8.6,X,F8.6,X,F8.6,X,F8.6)') trial, j, i, Tjin(j,i,trial)
        end do
     end do
  end do
  close(unit=iunit)
  if(debug) write(*,*) 'Done.'

  ! Clean up.
  call finalize()
  deallocate(trajectory_indices,Tji,Tjin)

  ! Terminate
  stop

  ! Report error if file could not be opened.
10 continue
  write(*,*) 'Error opening file ', trim(filename)
  stop

end program main

