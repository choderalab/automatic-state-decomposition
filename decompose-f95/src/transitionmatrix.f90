!=============================================================================================
! A module to compute and perform operations on transition matrices.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
! - Allow computeTimescalesBootstrap to automatically determine number of trials required
! to compute requested confidence bounds to 10% accuracy.
!=============================================================================================
module transitionmatrix
  use numeric_kinds ! for precision
  use constants ! global constants

  implicit none

  private
  public :: computeTransitionMatrix, computeTimescales, computeTimescalesBootstrap, &
       transitionMatrixEigs, computeTransitionCountMatrix, computeLumpedTransitionMatrix, &
       printTransitionMatrix
              
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  
  logical, parameter :: debug = .TRUE. ! Flag to enable printing of debugging information.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  ! Trajectory information
  type, public :: trajectory_t
     character(len=MAX_FILENAME_LENGTH) :: filename ! filename of trajectory
     real(dp) :: log_weight ! the (unnormalized) log weight of the trajectory
     integer :: start ! the index of the first snapshot in the trajectory
     integer :: length ! the length (in snapshots) of a trajectory
  end type trajectory_t
     
  !=============================================================================================
  ! Private module data.
  !=============================================================================================

contains

  !=============================================================================================
  ! Print the elements of a transition matrix.
  !=============================================================================================
  subroutine printTransitionMatrix(Tji)
    ! Parameters.
    real(dp), dimension(:,:) :: Tji

    ! Local variables.
    integer :: i
      ! Loop index.
    integer :: N
      ! Size.    
    character(len=128) :: format_string

    ! Get size.
    N = size(Tji,1)

    ! Construct format string
    write(format_string,*) N
    format_string = '('//adjustl(trim(format_string))//'F8.5)'

    write(*,*) 'Tji = '
    do i = 1,N
       write(*,format_string) Tji(i,:)
    end do
    
  end subroutine printTransitionMatrix

  !=============================================================================================
  ! Write the elements of a transition matrix to a file.
  !=============================================================================================
  subroutine writeTransitionMatrix(filename, Tji)
    use pdbio, only : new_unit

    ! Parameters.
    character(len=*), intent(in) :: filename
      ! File to write stateassignments to.
    real(dp), dimension(:,:) :: Tji
      ! Transition matrix to write.

    ! Local variables.
    integer :: i
      ! Loop index.
    integer :: N
      ! Size.    
    character(len=128) :: format_string
    integer :: iunit
      ! unit to use for file management

    write(*,*) 'Writing matrix to ', trim(filename)

    ! Get size.
    N = size(Tji,1)

    ! Construct format string
    write(format_string,*) N
    format_string = '('//adjustl(trim(format_string))//'E24.16)'

    ! Open file for writing.
    call new_unit(iunit)
    open(unit=iunit, file=trim(filename), status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    
    do i = 1,N
       write(iunit,format_string) Tji(i,:)
    end do
    
    ! Close file
    close(iunit)
    write(*,*) 'Done.'
    return

    ! Report error if file could not be opened.
10  continue
    write(*,*) 'module transitionmatrix : writeTransitionMatrix: Error opening file ', trim(filename)
    stop

  end subroutine writeTransitionMatrix

  !=============================================================================================
  ! Compute MLE of the transition matrix at a single lag time from trajectory data.
  !
  ! Note that this method does not make optimal use of the log weights for optimal precision.
  ! TODO: 
  ! - Use quad-precision or logs to better estimate T_ji when there are very small and very large trajectory weights?
  !=============================================================================================
  function computeTransitionMatrix(ntrajectories, trajectories, nstates, stateassignments, tau, &
       use_overlapping, use_timereversed, ignore_empty_states)
    use utilities, only : generateUniformInteger

    ! Parameters.
    integer, intent(in) :: ntrajectories
      ! The number of trajectories.
    type(trajectory_t), dimension(ntrajectories), intent(in) :: trajectories
      ! The input trajectories.
    integer, intent(in) :: nstates             ! the number of states
    integer, dimension(:), intent(in) :: stateassignments ! the current state assignments
    integer, intent(in) :: tau                 ! the lag time to use for estimating the transition matrix 
    logical, intent(in), optional :: use_overlapping     ! if .TRUE., overlapping trajectory segments will be used in estimate (default .TRUE.)
    logical, intent(in), optional :: use_timereversed    ! if .TRUE., time-reversed trajectories will be used in estimate (default .TRUE.)
    logical, intent(in), optional :: ignore_empty_states ! if .TRUE., empty states cause the diagonal element of the corresponding state to be 1
    
    ! Return value.
    real(dp), dimension(nstates,nstates) :: computeTransitionMatrix
    
    ! Local variables.
    integer :: i, j, k, t, t0, n               ! state and loop variables
    real(dp) :: weight ! weight
    real(dp), dimension(nstates,nstates) :: N_ji  ! numerator for transition count
    real(dp) :: column_sum
    real(dp), dimension(ntrajectories) :: trajectory_weights ! (unnormalized) trajectory weights
    real(dp) :: max_log_weight
    integer :: tau_skip
    logical :: use_overlapping_
    logical :: use_timereversed_
    logical :: ignore_empty_states_
    real(dp), dimension(nstates,nstates) :: Tji   ! Tji(j,i) is the probability of transitioning occupying state j at time tau if in state i at time zero.
    integer :: trajectory_index
    integer :: trajectory_start

    ! Set default settings.
    use_overlapping_ = .TRUE.
    if(present(use_overlapping)) use_overlapping_ = use_overlapping
    use_timereversed_ = .TRUE.
    if(present(use_timereversed)) use_timereversed_ = use_timereversed   

    ! Set increment.
    if(use_overlapping_) then
       tau_skip = 1
    else
       tau_skip = tau
    end if

    ! Set whether empty states should be ignored.
    ignore_empty_states_ = .FALSE.
    if(present(ignore_empty_states)) ignore_empty_states_ = ignore_empty_states

    ! Perform sanity checks.
    if(any((stateassignments < 1) .or. (stateassignments > nstates))) then
       write(*,*) 'computeTransitionMatrix: stateassignments has elements out of the range 1 to ', nstates
       stop
    end if

    ! Compute trajectory weights.
    max_log_weight = maxval(trajectories(:)%log_weight)
    trajectory_weights = exp(trajectories(:)%log_weight - max_log_weight)
    
    ! Initialize matrix of (weighted) transition counts.
    N_ji(:,:) = 0

    ! Loop over trajectories to estimate transition probabilities.
    do trajectory_index = 1,ntrajectories
       ! Determine trajectory weight and length.
       weight = trajectory_weights(trajectory_index)
       T = trajectories(trajectory_index)%length
       trajectory_start = trajectories(trajectory_index)%start

       do t0 = 1,(T-tau),tau_skip
          ! Determine initial and final states.
          i = stateassignments(trajectory_start + t0 - 1)
          j = stateassignments(trajectory_start + t0 + tau - 1)

          ! Increment weight.
          N_ji(j,i) = N_ji(j,i) + weight
       end do
    end do

    if(use_timereversed_) then
       ! Enforce detailed balance.       
       N_ji = N_ji + transpose(N_ji)
    end if
    
    ! Compute transition matrix elements.
    do i = 1,nstates
       ! Compute column sum.
       column_sum = sum(N_ji(:,i))

       ! Error check.
       if(column_sum == 0.0) then
          if(ignore_empty_states_) then
             ! Ensure the diagonal element will be 1.
             N_ji(i,i) = 1
             column_sum = 1
          else
             ! Throw an error
             write(*,*) 'computeTransitionMatrix: column_sum for column ', i, ' is zero'
             stop
          end if
       end if

       Tji(:,i) = N_ji(:,i) / column_sum
    end do
    
    ! Return transition matrix.
    computeTransitionMatrix = Tji
    
  end function computeTransitionMatrix

  !=============================================================================================
  ! Compute symmetric weighted transition count matrix.
  ! Detailed balance is enforced by using time-reversed trajectories, and overlapping trajectory segments are used.
  ! N_ji is the (unnormalized) probability of observing a trajectory connecting states i and j.
  ! This matrix can be used for efficient lumping without the need to recompute transition probabilities.
  !=============================================================================================
  subroutine computeTransitionCountMatrix(ntrajectories, trajectories, nstates, stateassignments, tau, N_ji, &
       use_timereversed)

    ! Parameters.
    integer, intent(in) :: ntrajectories
      ! The number of trajectories.
    type(trajectory_t), dimension(ntrajectories), intent(in) :: trajectories
      ! The input trajectories.
    integer, intent(in) :: nstates             
      ! the number of states
    integer, dimension(:), intent(in) :: stateassignments 
      ! the current state assignments
    integer, intent(in) :: tau                 
      ! the lag time to use for estimating the transition counts
    real(dp), dimension(nstates,nstates), intent(out) :: N_ji  
      ! N_ji(j,i) is the symmetrized weighted transition count matrix
    logical :: use_timereversed
      ! if .true., will add count from time-reversed trajectories
    
    ! Local variables.
    integer :: i, j, k, t, t0, n               ! state and loop variables
    real(dp) :: weight ! weight
    real(dp), dimension(ntrajectories) :: trajectory_weights ! (unnormalized) trajectory weights
    real(dp) :: max_log_weight
    logical :: use_overlapping_
    integer :: trajectory_index
    integer :: trajectory_start

    ! Compute trajectory weights.
    max_log_weight = maxval(trajectories(:)%log_weight)
    trajectory_weights(:) = exp(trajectories(:)%log_weight - max_log_weight)

    ! Initialize matrix of weighted transition counts.
    N_ji(:,:) = 0

    ! Loop over trajectories to compute transition counts.
    do trajectory_index = 1,ntrajectories       
       ! Determine trajectory weight and length.
       weight = trajectory_weights(trajectory_index)
       T = trajectories(trajectory_index)%length
       trajectory_start = trajectories(trajectory_index)%start

       do t0 = 1,(T-tau)

          ! Determine initial and final states.
          i = stateassignments(trajectory_start + t0 - 1)
          j = stateassignments(trajectory_start + t0 + tau - 1)

          ! Accumulate symmetric count matrix.
          N_ji(j,i) = N_ji(j,i) + weight
          if (use_timereversed) then
             N_ji(i,j) = N_ji(i,j) + weight
          end if
       end do
    end do
    
  end subroutine computeTransitionCountMatrix

  !=============================================================================================
  ! Compute symmetrized transition matrix and equilibrium probabilities.
  ! Detailed balance is enforced by using time-reversed trajectories, and overlapping trajectory segments are used.
  !
  ! Note that this method does not make optimal use of the log weights for optimal precision.
  ! TODO: Use quad-precision or logs to better estimate T_ji?
  !=============================================================================================
  subroutine computeSymmetrizedTransitionMatrix(ntrajectories, trajectories, nstates, stateassignments, tau, Tsym_ji, ev)

    ! Parameters.
    integer, intent(in) :: ntrajectories
      ! The number of trajectories.
    type(trajectory_t), dimension(ntrajectories), intent(in) :: trajectories
      ! The input trajectories.
    integer, intent(in) :: nstates             ! the number of states
    integer, dimension(:), intent(in) :: stateassignments ! the current state assignments
    integer, intent(in) :: tau                 ! the lag time to use for estimating the transition matrix 
    real(dp), dimension(nstates,nstates), intent(out) :: Tsym_ji   ! Tsym_ji(j,i) is the symmetrized transition matrix
    real(dp), dimension(nstates), intent(out), optional :: ev  ! eigenvector of Tsym_ji with unit eigenvalue
    
    ! Local variables.
    integer :: i, j, k, t, t0, n               ! state and loop variables
    real(dp) :: weight ! weight
    real(dp), dimension(nstates,nstates) :: N_ji  
    real(dp), dimension(nstates) :: D_i           
    real(dp), dimension(ntrajectories) :: trajectory_weights ! (unnormalized) trajectory weights
    real(dp) :: max_log_weight
    integer :: tau_skip
    logical :: use_overlapping_
    integer :: trajectory_index
    integer :: trajectory_start
    logical :: use_timereversed

    ! For now, we must use time-reversed trajectories as well to ensure eigenvalues are real.
    use_timereversed = .true.
    
    ! Compute symmetric weighted transition count matrix.
    call computeTransitionCountMatrix(ntrajectories, trajectories, nstates, stateassignments, tau, N_ji, use_timereversed)
    
    ! Compute column sums (equivalent to row sums).
    D_i = sum(N_ji,2)

    ! Check to ensure no D_i are zero.
    if(any(D_i <= 0.0)) then
       write(*,*) 'transitionmatrix.f90: computeSymmetrizedTransitionMatrix: some D_i are zero!'       
       write(*,*) 'D_i = '
       write(*,*) D_i
       stop
    end if

    ! Take square root.
    D_i = sqrt(D_i)

    ! Compute symmetrized transition matrix.
    forall(i = 1:nstates, j = 1:nstates)
       Tsym_ji(j,i) = N_ji(j,i) / (D_i(i) * D_i(j))
    end forall

    ! If desired, compute equilibrium probabilities.
    if(present(ev)) ev = d_i
    
  end subroutine computeSymmetrizedTransitionMatrix

  !=============================================================================================
  ! Compute lumped transition matrix from transition count matrix
  !=============================================================================================
  subroutine computeLumpedTransitionMatrix(nmicrostates, N_ji, nmacrostates, T_ba, micro_to_macro)  
    
    ! Parameters.
    integer, intent(in) :: nmicrostates
      ! number of microstates
    real(dp), dimension(nmicrostates,nmicrostates), intent(in) :: N_ji
      ! microstates transition count matrix
    integer, intent(in) :: nmacrostates
      ! number of macrostates
    real(dp), dimension(nmacrostates,nmacrostates), intent(out) :: T_ba
      ! lumped macrostate transition matrix
    integer, dimension(nmicrostates), intent(in) :: micro_to_macro
      ! micro_to_macro(i) is the (compacted) macrostate associated with (compacted) microstate i

    ! Local variables.
    integer :: i, j
      ! (compacted) microstate indices
    integer :: a, b
      ! (compacted) macrostate indices    
    integer, dimension(nmacrostates) :: D_a 
      ! column sums

    ! First, compute lumped transition counts.
    T_ba(:,:) = 0.0
    do i = 1,nmicrostates
       do j = 1,nmicrostates
          a = micro_to_macro(i)
          b = micro_to_macro(j)
          T_ba(b,a) = T_ba(b,a) + N_ji(j,i)
       end do
    end do

    ! Normalize column sums.
    D_a = sum(T_ba,1)
    do a = 1,nmacrostates
       do b = 1,nmacrostates
          T_ba(b,a) = T_ba(b,a) / D_a(a)
       end do
    end do

  end subroutine computeLumpedTransitionMatrix

  !=============================================================================================
  ! Extract the lower triangle of the provided symmetrix matrix into a vector, row-wise.
  !=============================================================================================
  function lowerTriangle(N, Tsym_ji) result(packed_matrix)

    ! Arguments.
    integer, intent(in) :: N
    real(dp), dimension(N,N), intent(in) :: Tsym_ji

    ! Return values.
    real(dp), dimension(N*(N-1)/2) :: packed_matrix
    
    ! Local variables.
    integer :: index, i, j

    index = 1
    do j = 1,N
       do i = 1,j
          packed_matrix(index) = Tsym_ji(j,i)
          index = index + 1
       end do
    end do

  end function lowerTriangle

  !=============================================================================================
  ! Compute specified number of eigenvalues and eigenvectors (if desired) of the transition
  ! matrix, imposing detailed balance to ensure real ew and ev.
  !=============================================================================================
  subroutine transitionMatrixEigs(ntrajectories, trajectories, nstates, stateassignments, tau, neigs, ew, ev_left, ev_right)
    !use eispack, only : rs, rspp ! for symmetric eigenproblem
    !use eigs_trlan, only : eigs ! for symmetric eigenproblem
    use eigs_lapack, only : eigs ! for symmetric eigenproblem, using LAPACK
    !use eigs_arpack, only : eigs ! for symmetric eigenproblem
    use utilities, only : find

    ! Parameters.
    integer, intent(in) :: ntrajectories
      ! The number of trajectories.
    type(trajectory_t), dimension(ntrajectories), intent(in) :: trajectories
      ! The input trajectories.
    integer, intent(in) :: nstates
      ! number of states
    integer, dimension(:), intent(in) :: stateassignments
    integer, intent(in) :: tau
      ! lag time to use in computing transition matrix
    integer, intent(in) :: neigs
      ! number of dominant eigenvalues to compute
    real(dp), dimension(neigs), intent(out) :: ew
      ! ew(k) is the kth-largest eigenvalue of the transition matrix
    real(dp), dimension(nstates,neigs), intent(out), optional :: ev_left, ev_right
      ! ev(k,i) is component i of eigenvector k, corresponding to eigenvalue k

    ! Local variables.
    real(dp), dimension(nstates,nstates) :: N_ji
      ! Symmetrized transition count matrix.
    real(dp), dimension(nstates) :: D_i
      ! square root of state probabilities
    real(dp), dimension(nstates,nstates) :: Tsym_ji
      ! Transition matrix.
    real(dp), dimension(nstates,neigs) :: ev
      ! Storage for eigenvectors.
    integer :: i, j, k
      ! indices
    integer, dimension(:), pointer :: nonempty_indices
      ! indices of rows/columns that are not empty
    integer :: nnonempty
    logical :: use_timereversed

    ! For now, we must use time-reversed trajectories to ensure detailed balance and real eigenvalues.
    use_timereversed = .true.

    ! Compute symmetric weighted transition count matrix.
    call computeTransitionCountMatrix(ntrajectories, trajectories, nstates, stateassignments, tau, N_ji, use_timereversed)

    ! Compute square roots of column sums (equivalent to row sums).  This will be the stationary eigenvector.
    D_i = sqrt(sum(N_ji,2))

    ! If rows or columns are empty, compact the matrix.
    if(any(D_i <= 0) .or. any(isnan(D_i))) then

       if(present(ev_left) .or. present(ev_right)) then
          ! Throw an error, because we will screw things up if we try to compute eigenvectors for reduced matrix.
          write(*,*) 'transitionMatrixEigs: Some states have no samples, and eigenvectors were requested. Halting.'
          write(*,*) 'D_i = '
          write(*,*) D_i
          stop
       end if

       ! Find nonempty rows and columns.
       nonempty_indices => find( (D_i > 0.0) .and. .not. isnan(D_i) )
       nnonempty = size(nonempty_indices,1)

       if(debug) write(*,*) 'transitionMatrixEigs: eliminating empty rows/columns to get ', nnonempty, ' / ', nstates

       ! Throw an error if we've thrown out more states than we want eigenvalues for.
       if(nnonempty < neigs) then
          write(*,*) 'transitionMatrixEigs: Want ', neigs, ' ew but only ', nnonempty, ' nonempty rows/columns.'
          write(*,*) 'nonempty_indices ='
          write(*,*) nonempty_indices
          stop
       end if
       
       ! Compute symmetrized transition matrix.
       forall(i = 1:nnonempty, j = 1:nnonempty)
          Tsym_ji(j,i) = N_ji(nonempty_indices(j),nonempty_indices(i)) / (D_i(nonempty_indices(i)) * D_i(nonempty_indices(j)))
       end forall
       deallocate(nonempty_indices)       

       ! Sanity check.
       if(any(isnan(Tsym_ji(1:nnonempty,1:nnonempty)))) then
          write(*,*) '(a) Tsym_ji contains NaN!'
          write(19,*) N_ji
          write(20,*) Tsym_ji
          stop
       end if

       ! Compute eigenvalues.
       call eigs(nnonempty, Tsym_ji(1:nnonempty,1:nnonempty), neigs, ew)

    else

       ! Compute symmetrized transition matrix.
       forall(i = 1:nstates, j = 1:nstates)
          Tsym_ji(j,i) = N_ji(j,i) / (D_i(i) * D_i(j))
       end forall

       ! Sanity check.
       if(any(isnan(Tsym_ji))) then
          write(*,*) '(b) Tsym_ji contains NaN!'
          write(19,*) N_ji
          write(20,*) Tsym_ji
          stop
       end if

       ! Compute eigenvalues and eigenvectors.
       if(present(ev_left) .or. present(ev_right)) then
          call eigs(nstates, Tsym_ji, neigs, ew, evec = ev)
       else
          call eigs(nstates, Tsym_ji, neigs, ew)
       end if

       ! Compute left and right eigenvectors, as necessary.
       do k = 1,neigs
          if(present(ev_left)) ev_left(:,k) = ev(:,k) / sqrt(d_i)
          if(present(ev_right)) ev_right(:,k) = ev(:,k) * sqrt(d_i)
       end do

    end if

  end subroutine transitionMatrixEigs

  !=============================================================================================
  ! Estimate MLE estimate of timescales.
  !=============================================================================================
  subroutine computeTimescales(ntrajectories, trajectories, nstates, stateassignments, ntau, ntimescales, &
       tau_unit, taus, timescales)

    ! Parameters.
    integer, intent(in) :: ntrajectories
      ! The number of trajectories.
    type(trajectory_t), dimension(ntrajectories), intent(in) :: trajectories
      ! The input trajectories.
    integer, intent(in) :: nstates
      ! number of states
    integer, dimension(:), intent(in) :: stateassignments
    integer, intent(in) :: ntau
      ! number of tau values
    integer, intent(in) :: ntimescales
      ! number of timescales to compute (stationary one is ignored)
    real(sp), intent(in) :: tau_unit
      ! unit time associated with one sampling interval
    integer, dimension(ntau), intent(in) :: taus
      ! Tau values at which to compute the timescales and uncertainties.
    real(sp), dimension(ntimescales,ntau), intent(out) :: timescales
      ! timescales to compute
      ! timescales(k,i) is the timescale associated with k slowest process, lag time tau(i)

    ! Local variables.
    integer :: tau_index
      ! indices
    integer :: tau 
      ! current lag time
    integer :: neigs
      ! Number of eigenvalues to compute.
    real(dp), dimension(ntimescales+1) :: mu_k
      ! Storage for computed eigenvalues.

    ! Compute number of desired eigenvalues.
    neigs = ntimescales + 1

    ! Sanity check on arguments.
    if(ntimescales >= nstates) then
       write(*,*) 'computeTimescales: ntimescales must be from 0 to nstates-1.'
       stop
    end if

    ! Compute desired eigenvalues for all lag times.
    !$omp parallel do default(shared) private(tau, mu_k)
    do tau_index = 1,ntau
       ! Get lag time to compute transition matrix at.
       tau = taus(tau_index)

       ! Compute eigenvalues of transition matrix.
       call transitionMatrixEigs(ntrajectories, trajectories, nstates, stateassignments, tau, neigs, ew = mu_k)

       ! Compute timecales from eigenvalues.
       timescales(:,tau_index) = - (tau_unit * tau) / log(mu_k(2:neigs))
    end do
    !$omp end parallel do

  end subroutine computeTimescales

  !=============================================================================================
   ! Estimate timescales and uncertainties by bootstrap analysis.
  !=============================================================================================
  subroutine computeTimescalesBootstrap(ntrajectories, trajectories, nstates, stateassignments, ntau, ntimescales, &
       tau_unit, taus, timescales, confidence_interval, timescales_lower, timescales_upper, ntrials)

    use utilities, only : generateUniformIntegers, compute_confidence_bounds

    ! Parameters.
    integer, intent(in) :: ntrajectories
      ! The number of trajectories.
    type(trajectory_t), dimension(ntrajectories), intent(in) :: trajectories
      ! The input trajectories.
    integer, intent(in) :: nstates
      ! number of states
    integer, dimension(:), intent(in) :: stateassignments
    integer, intent(in) :: ntau
      ! number of tau values
    integer, intent(in) :: ntimescales
      ! number of timescales to compute (stationary one is ignored)
    real(sp), intent(in) :: tau_unit
      ! unit time associated with one sampling interval
    integer, dimension(ntau), intent(in) :: taus
      ! Tau values at which to compute the timescales and uncertainties.
    real(sp), dimension(ntimescales,ntau), intent(out) :: timescales
      ! timescales to compute
      ! timescales(k,i) is the timescale associated with k slowest process, lag time tau(i)
    real(sp), intent(in) :: confidence_interval
      ! The confidence interval (in [0,1]) to use for computing lower and upper confidence bounds.
    real(sp), dimension(ntimescales,ntau), intent(out) :: timescales_lower, timescales_upper
      ! Lower and upper confidence intervals
    integer, intent(in) :: ntrials
      ! number of bootstrap trials

    ! Local variables.
    integer :: trial
      ! index of bootstrap trial
    real(sp), dimension(ntimescales,ntau) :: timescales_mle
    real(sp), dimension(ntrials,ntimescales,ntau) :: timescales_samples
      ! samples from bootstrap sampling
    integer, dimension(ntrajectories) :: trajectory_indices
      ! trajectory indices of bootstrap sample
    integer :: k, t

    ! Compute MLE estimate of timescales.
    call computeTimescales(ntrajectories, trajectories, nstates, stateassignments, ntau, ntimescales, &
       tau_unit, taus, timescales_mle)

    ! Generate bootstrap samples.
    !!$omp parallel do default(shared) private(trajectory_indices)
    !!$omp parallel do default(private) shared(ntrajectories,trajectories,nstates,stateassignments,ntau,ntimescales) &
    !!$omp shared(tau_unit,taus,timescales_samples)
    do trial = 1, ntrials
       write(*,*) 'trial = ', trial, '/', ntrials

       ! Generate bootstrap sample of trajectory indices.
       call generateUniformIntegers(ntrajectories, trajectory_indices)
       
       ! Compute timescales of this sample.
       call computeTimescales(ntrajectories, trajectories(trajectory_indices), nstates, stateassignments, ntau, ntimescales, &
            tau_unit, taus, timescales_samples(trial,:,:))
    end do
    !!$omp end parallel do

    ! Compute confidence intervals.
    timescales = timescales_mle    
    do k = 1,ntimescales
       do t = 1,ntau
          ! timescales_uncertainties(k,t) = std(timescales_samples(:,k,t))
          call compute_confidence_bounds(timescales_samples(:,k,t), confidence_interval, &
               timescales_lower(k,t), timescales_upper(k,t))
       end do
    end do

  end subroutine computeTimescalesBootstrap

end module transitionmatrix
