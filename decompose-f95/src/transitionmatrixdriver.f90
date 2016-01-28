!=============================================================================================
! Driver to compute transition matrix elements (with uncertainties) as a function of lag time for given state decomposition.
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
  character(len=*), parameter :: version = '$Id: transitionmatrixdriver.f90,v 1.2 2006/11/06 05:56:03 jchodera Exp $' ! version string

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
  real(sp), dimension(:,:), allocatable :: Tji ! transition matrix
  real(sp), dimension(:,:,:,:), allocatable :: Tjimn ! transition matrix as a function of lag time with bootstrap trials
  real(sp), dimension(:,:,:), allocatable :: Tjim ! MLE of transition matrix as a function of lag time
  real(sp), dimension(:,:,:), allocatable :: dTjim ! uncertainty in transition matrix as a function of lag time
  real(sp), dimension(:,:,:), allocatable :: Tjim_lower, Tjim_upper ! confidence bounds
  integer :: iunit ! file unit for writing
  integer :: i, j
  real(sp) :: confidence_interval

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
     write(*,*) 'usage: ./compute-transitionmatrix control.xml stateassignments-file transitionmatrix-file max_tau tau_step'
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
  confidence_interval = timescales_confidence_interval
  
  ! Determine number of lag times.
  ntau = floor(real(max_tau,sp) / real(tau_step,sp))
  if(ntau == 0) stop

  ! allocate storage for bootstrap samples, taus, and metastabilities
  allocate(trajectory_indices(ntrajectories))
  allocate(taus(ntau))
  allocate(Tji(nmacrostates,nmacrostates))
  allocate(Tjim(nmacrostates,nmacrostates,ntau),dTjim(nmacrostates,nmacrostates,ntau))
  allocate(Tjim_lower(nmacrostates,nmacrostates,ntau),Tjim_upper(nmacrostates,nmacrostates,ntau))
  allocate(Tjimn(nmacrostates,nmacrostates,ntau,ntrials))

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
             compacted_macrostate_assignments, tau, ignore_empty_states = .TRUE.)

        !call printTransitionMatrix(Tji)

        ! Store bootstrapped sample.
        Tjimn(:,:,tau_index,trial) = Tji
       
     end do
     !!$omp end parallel do
     
  end do

  ! Compute maximum-likelihood estimate, and uncertainties as the std dev of the estimate over bootstrap trials.
  write(*,*) 'Computing MLE and std devs...'
  do tau_index = 1,ntau
     ! Get lag time to compute transition matrix at.
     tau = taus(tau_index)

     ! Compute transition matrix.
     Tji = computeTransitionMatrix(ntrajectories, trajectories, nmacrostates, &
          compacted_macrostate_assignments, tau, ignore_empty_states = .TRUE.)

     !call printTransitionMatrix(real(Tji,dp))

     ! Store maximum-likelihood estimate.
     Tjim(:,:,tau_index) = Tji

     ! Compute uncertainties.
     do i = 1,nmacrostates
        do j = 1,nmacrostates
           ! std dev
           dTjim(j,i,tau_index) = std(Tjimn(j,i,tau_index,:))

           ! confidence bounds
           call compute_confidence_bounds(Tjimn(j,i,tau_index,:), confidence_interval, &
                Tjim_lower(j,i,tau_index), Tjim_upper(j,i,tau_index))
        end do
     end do
  end do

  ! Write to file.
  if(debug) write(*,*) 'Writing to file ', trim(filename), '...'
  call new_unit(iunit)
  open(unit=iunit, file=trim(filename), status='REPLACE', err=10, position='REWIND', &
       form='FORMATTED', action='WRITE')    
  do tau_index = 1,ntau
     ! Get lag time at which to compute transition matrix.
     tau = taus(tau_index)

     do j = 1,nmacrostates
        do i = 1,nmacrostates
           ! Write line.
           write(iunit,'(I5,X,I4,X,I4,X,F8.6,X,F8.6,X,F8.6,X,F8.6)') tau, j, i, Tjim(j,i,tau_index), dTjim(j,i,tau_index), &
                Tjim_lower(j,i,tau_index), Tjim_upper(j,i,tau_index)
        end do
     end do
  end do
  close(unit=iunit)
  if(debug) write(*,*) 'Done.'

  ! Clean up.
  call finalize()
  deallocate(trajectory_indices,taus,Tji,Tjim,dTjim,Tjimn,Tjim_lower,Tjim_upper)

  ! Terminate
  stop

  ! Report error if file could not be opened.
10 continue
  write(*,*) 'Error opening file ', trim(filename)
  stop

end program main

