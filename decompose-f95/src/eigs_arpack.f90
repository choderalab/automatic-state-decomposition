!=============================================================================================
! A simple interface to ARPACK for finding a few eigenvalues and eigenvectors of a real symmetric matrix.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
! - Add support for sparse specification of matrix A.
!=============================================================================================
module eigs_arpack
  use numeric_kinds ! for precision

  implicit none

  private
  public :: eigs
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  logical, parameter :: debug = .false.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

contains

  !=============================================================================================
  ! Compute the NEV larges eigenvalues of the real symmetric matrix A.
  !=============================================================================================
  subroutine eigs(nrow, A, nev, eval, evec, starting_guess)
    ! Parmeters.
    integer, intent(in) :: nrow
      ! Dimension of matrix A.
    real(dp), dimension(nrow,nrow), intent(in) :: A
      ! Real symmetric matrix to compute eigenvalue decomposition of.
    integer, intent(in) :: nev
      ! Number of eigenvalues desired.
    real(dp), dimension(nev), intent(out), optional :: eval
      ! Eigenvalues
    real(dp), dimension(nrow,nev), intent(out), optional :: evec
    real(dp), dimension(nrow), intent(in), optional :: starting_guess

    ! Local variables.
    real(dp), dimension(nrow) :: eval_
      ! Storage for eigenvalues.
    real(dp), dimension(nrow,nrow) :: evec_
      ! Storage for eigenvectors.

    real(dp), dimension(nrow) :: x, y
      ! Temporary vectors.

    ! Local arrays
    real(dp), dimension(3*nrow) :: workd
    real(dp), dimension(nrow) :: resid
    integer, dimension(11) :: iparam, ipntr
    
    real(dp), dimension(:), allocatable :: workl
    real(dp), dimension(:,:), allocatable :: d
    real(dp), dimension(:,:), allocatable :: v
    logical, dimension(:), allocatable :: select
    
    character :: bmat*1, which*2
    integer :: ido, n, ncv, lworkl, info, ierr, j, nx, ishfts, maxitr, mode1, nconv, ldv
    logical :: rvec
    real(dp) :: tol
    real(dp) :: sigma
    
    ! Constants.
    real(dp), parameter :: zero = 0.0D+0
    
    ! Set dimensions for this problem.
    n = nrow
    ldv = nrow

    ! Set length of ARNOLDI factorization.
    ncv = 2 * nev ! recommended choice
    bmat  = 'I'
    which = 'LA'

    ! Set the workspace size.
    lworkl = ncv*(ncv+8)
    allocate(workl(lworkl),d(ncv,2),v(nrow,ncv),select(ncv))    

    ! Set stopping rules.
    tol = zero

    ! Specify that a random starting vector should be chosen
    info = 0

    ! Set reverse communication parameter
    ido = 0

    ! Specify the algorithm mode
    ishfts = 1
    maxitr = 300 
    mode1 = 1

    iparam(1) = ishfts                
    iparam(3) = maxitr                  
    iparam(7) = mode1

    ! M A I N   L O O P (Reverse communication loop)
    do
       ! Repeatedly call the routine DSAUPD and take 
       ! actions indicated by parameter IDO until    
       ! either convergence is indicated or maxitr   
       ! has been exceeded.                          
       call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
            ncv, v, ldv, iparam, ipntr, workd, workl, &
            lworkl, info )

       if (ido .eq. -1 .or. ido .eq. 1) then
          ! Perform matrix-vector multiplication.
          ! TODO: Can replace this with sparse matrix-vector multiply.
          
          ! Extract vector from workd.
          x = workd(ipntr(1):(ipntr(1)+nrow-1))
          
          ! Apply matrix.
          y = matmul(A, x)
          
          ! Store result.
          workd(ipntr(2):(ipntr(2)+nrow-1)) = y
       else
          ! Exit loop.
          exit
       end if
    end do

    ! Either we have convergence or there is an error.                            
    if ( info .lt. 0 ) then
       ! Error message.
       print *, ' '
       print *, ' Error with _saupd, info = ', info
       print *, ' Check documentation in _saupd '
       print *, ' '
       stop
    end if
    
    ! No fatal errors occurred.                 
    ! Post-Process using DSEUPD.                
    
    ! Computed eigenvalues may be extracted.
    !
    ! Eigenvectors may be also computed now if desired (indicated by rvec = .true.)
    ! TODO: Turn off computation of eigenvectors if not needed.
    rvec = .true.
    
    call dseupd ( rvec, 'All', select, d, v, ldv, sigma, &
         bmat, n, which, nev, tol, resid, ncv, v, ldv, &
         iparam, ipntr, workd, workl, lworkl, ierr )
    
    !         %----------------------------------------------%
    !         | Eigenvalues are returned in the first column |
    !         | of the two dimensional array D and the       |
    !         | corresponding eigenvectors are returned in   |
    !         | the first NCONV (=IPARAM(5)) columns of the  |
    !         | two dimensional array V if requested.        |
    !         | Otherwise, an orthogonal basis for the       |
    !         | invariant subspace corresponding to the      |
    !         | eigenvalues in D is returned in V.           |
    !         %----------------------------------------------%
    
    if ( ierr .ne. 0) then
       ! Error.
       print *, ' '
       print *, ' Error with _seupd, info = ', ierr
       print *, ' Check the documentation of _seupd. '
       print *, ' '
       stop
    end if
    
    if(debug) then
       nconv =  iparam(5)
       do j = 1,nconv
          ! Compute the residual norm || A x - lambda x || for the NCONV accurately computed eigenvalues and eigenvectors.
          ! iparam(5) indices how many are accurate to the requested tolerance
          d(j,2) = sqrt(sum( (matmul(A, v(:,j)) - d(j,1) * v(:,j))**2 )) / abs(d(j,1))
       end do
       
       ! Display computed residuals
       call dmout(6, nconv, 2, d, ncv, -6, 'Ritz values and relative residuals')
       
       ! Print additional convergence information.
       if ( info .eq. 1) then
          print *, ' '
          print *, ' Maximum number of iterations reached.'
          print *, ' '
       else if ( info .eq. 3) then
          print *, ' ' 
          print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
          print *, ' '
       end if
       
       print *, ' '
       print *, ' _SSIMP '
       print *, ' ====== '
       print *, ' '
       print *, ' Size of the matrix is ', n
       print *, ' The number of Ritz values requested is ', nev
       print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
       print *, ' What portion of the spectrum: ', which
       print *, ' The number of converged Ritz values is ', nconv 
       print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
       print *, ' The number of OP*x is ', iparam(9)
       print *, ' The convergence criterion is ', tol
       print *, ' '       
    end if

    ! Provide eigenvalues and eigenvectors as desired.
    if(present(eval)) eval = d(:,1)
    if(present(evec)) evec = v(:,1:nev)

    ! Free memory.
    deallocate(workl,d,v,select)

  end subroutine eigs
end module eigs_arpack
