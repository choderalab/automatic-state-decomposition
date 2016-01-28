!=============================================================================================
! Driver to compute metastability (with uncertainties) as a function of lag time for given state decomposition.
!=============================================================================================
program main
  use decompose ! for state space decomposition
  use utilities, only : setRandomSeed ! to set the random seed from the current time
  use constants ! global constants
  use utilities, only : unique, trace, std, generateUniformIntegers
  use pdbio, only : new_unit
  use transitionmatrix, only : computeTransitionMatrix, printTransitionMatrix
    
  implicit none

  ! Constants.
  character(len=*), parameter :: version = '$Id: metastabilitydriver.f90,v 1.1 2006/11/06 00:51:37 jchodera Exp $' ! version string

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
  integer :: ntrials ! number of bootstrap trials to conduct
  integer :: ntau ! number of tau values to evaluate metastability for
  integer :: tau_index
  integer, dimension(:), allocatable :: taus ! tau values to evaluate metastability at
  integer, dimension(:), allocatable :: trajectory_indices ! indices of trajectories in bootstrap sample
  real(dp), dimension(:,:), allocatable :: Tji ! transition matrix
  real(dp), dimension(:,:), allocatable :: M_trial_tau ! M_trial_tau(trial,tau_index) is the metastability from trial TRIAL at lag time index TAU
  real(dp), dimension(:), allocatable :: M_tau ! M_tau(tau_index) is the metastability at lag time index TAU
  real(dp), dimension(:), allocatable :: dM_tau ! dM_tau(tau_index) is the uncertainty of the metastability at lag time index TAU
  integer :: iunit ! file unit for writing

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
     write(*,*) 'usage: ./compute-metastability control.xml stateassignments-file metastability-file max_tau tau_step'
     write(*,*) ''
     write(*,*) 'The same number of bootstrap trials specified for the timescales in control.xml will be used here.'
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
  read(string,'(I8)') max_tau
  call get_command_argument(5, value=string)
  read(string,'(I8)') tau_step
  write(*,*) 'max_tau = ', max_tau, ', tau_step = ', tau_step

  ! use number of bootstrap samples specified for timescales
  ntrials = timescales_bootstrap_ntrials  
  
  ! Determine number of lag times.
  ntau = floor(real(max_tau,sp) / real(tau_step,sp))
  if(ntau == 0) stop

  ! allocate storage for bootstrap samples, taus, and metastabilities
  allocate(trajectory_indices(ntrajectories))
  allocate(taus(ntau))
  allocate(M_trial_tau(ntrials,ntau))
  allocate(Tji(nmacrostates,nmacrostates))

  ! Calculate lag times.
  do tau_index = 1,ntau
     taus(tau_index) = tau_index * tau_step
  end do

  ! Compute metastabilities and uncertainties by bootstrap.
  write(*,*) 'Performing ', ntrials, ' bootstrap trials...'
  do trial = 1, ntrials
     write(*,*) 'trial = ', trial, '/', ntrials
     
     ! Generate bootstrap sample of trajectory indices.
     call generateUniformIntegers(ntrajectories, trajectory_indices)
     
     ! Compute metastabilities at all lag times.
     !!$omp parallel do default(shared) private(tau, Tji)
     do tau_index = 1,ntau
        ! Get lag time to compute transition matrix at.
        tau = taus(tau_index)
        
        ! Compute transition matrix.
        Tji = computeTransitionMatrix(ntrajectories, trajectories(trajectory_indices), nmacrostates, &
             compacted_macrostate_assignments, tau)

        !call printTransitionMatrix(Tji)
        
        ! Compute metastabilites.
        M_trial_tau(trial,tau_index) = trace(Tji)
       
     end do
     !!$omp end parallel do
     
  end do

  ! Compute maximum-likelihood estimate, and uncertainties as the std dev of the estimate over bootstrap trials.
  write(*,*) 'Computing MLE and std devs...'
  allocate(M_tau(ntau),dM_tau(ntau))
  do tau_index = 1,ntau
     ! Get lag time to compute transition matrix at.
     tau = taus(tau_index)

     ! Compute transition matrix.
     Tji = computeTransitionMatrix(ntrajectories, trajectories, nmacrostates, &
          compacted_macrostate_assignments, tau)

     ! Compute maximum-likelihood estimate.
     M_tau(tau_index) = trace(Tji)

     ! Compute std dev.
     dM_tau(tau_index) = std(M_trial_tau(:,tau_index))
  end do

  ! Write to file.
  if(debug) write(*,*) 'Writing to file ', trim(filename), '...'
  call new_unit(iunit)
  open(unit=iunit, file=trim(filename), status='REPLACE', err=10, position='REWIND', &
       form='FORMATTED', action='WRITE')    
  do tau_index = 1,ntau
     ! Get lag time to compute transition matrix at.
     tau = taus(tau_index)

     ! Write line.
     write(iunit,'(F8.3,X,F8.3,X,F8.3)') tau*tau_unit, M_tau, dM_tau
  end do
  close(unit=iunit)
  if(debug) write(*,*) 'Done.'

  ! Clean up.
  call finalize()
  deallocate(trajectory_indices,taus,M_trial_tau,M_tau,dM_tau,Tji)

  ! Terminate
  stop

  ! Report error if file could not be opened.
10 continue
  write(*,*) 'Error opening file ', trim(filename)
  stop

end program main

