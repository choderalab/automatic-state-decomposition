!=============================================================================================
! A simple interface to LAPACK for finding a few eigenvalues and eigenvectors of a real symmetric matrix.
!
! This routine is horribly inefficient because *all* eigenvalues are found.
! It is provided just for testing purposes.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
! - Add support for sparse matrices.
!=============================================================================================
module eigs_lapack
  use numeric_kinds ! for precision

  implicit none

  private
  public :: eig, eigs
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  logical, parameter :: debug = .TRUE.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

contains

  !=============================================================================================
  ! Compute the NED larges eigenvalues of the real symmetric matrix A.
  !=============================================================================================
  subroutine eigs_dsyevx(nrow, A, ned, eval, evec, status)
    ! Parmeters.
    integer, intent(in) :: nrow
      ! Dimension of matrix A.
    real(dp), dimension(nrow,nrow), intent(in) :: A
      ! Real symmetric matrix to compute eigenvalue decomposition of.
    integer, intent(in) :: ned
      ! Number of eigenvalues desired.
    real(dp), dimension(ned), intent(out), optional :: eval
      ! Eigenvalues
    real(dp), dimension(nrow,ned), intent(out), optional :: evec
      ! Eigenvectors, stored as columns.
    integer, intent(out), optional :: status
      ! indicates failure if negative.  If failure and status is not specified, execution is stopped.

    ! Local variables.
    real(dp), dimension(nrow) :: eval_
    real(dp), dimension(nrow,nrow) :: evec_
    character :: jobz
      ! argument for dsyev: 'N' for only ew, 'V' for ev too
    character, parameter :: uplo = 'U'
      ! 'U' means upper triangle of matrix is stored by dsyev
    integer :: lwork
      ! size of work array
    real(dp), dimension(:), allocatable :: work
      ! workspace
    integer, dimension(5*nrow) :: iwork
      ! integer workspace
    integer :: info
      ! exit information
    integer :: il, iu
      ! Lower and upper indices of eigenvalues to be found.
    real(dp) :: vl, vu
      ! Upper and lower bounds on eigenvalues to be found (not used).
    real(dp), parameter :: abstol = 0.0
      ! Absolute tolerance to use in determining eigenvalues, or 0 if a suitable tolerance is to be computed by LAPACK.
    integer :: neval_found
      ! The number of eigenvalues found by LAPACK call.
    character(1) :: range
      ! Character identifying whether given index range il..iu of eigenvalues is to be found (if 'I') or given range (if 'V').
    real(dp), dimension(:,:), allocatable :: A_copy
      ! Copy of the matrix A, destroyed by dsyevr.
    integer, dimension(nrow) :: ifail

    status = 0

    ! Instruct dsyevr to find specified index range of eigenvalues.
    range = 'I'
    il = nrow-ned+1
    iu = nrow
    
    ! Compute only eigenvalues by default.
    jobz = 'N'
    ! If evec provided, compute eigenvectors too.
    if(present(evec)) jobz = 'V' 

    ! Make a copy of the matrix A which will be destroyed by dsyevr.
    allocate(A_copy(nrow,nrow))
    A_copy = A

    ! Inquire about optimal workspace size, lwork and liwork.
    lwork = -1
    allocate(work(1))
    call dsyevx(jobz, range, uplo, nrow, A_copy, nrow, vl, vu, il, iu, abstol, neval_found, eval_, evec_, nrow, &
         work, lwork, iwork, ifail, info)
    lwork = int(work(1))
    deallocate(work)

    ! Allocate workspace.
    allocate(work(lwork))

    ! Perform eigenvalue decomposition.
    call dsyevx(jobz, range, uplo, nrow, A_copy, nrow, vl, vu, il, iu, abstol, neval_found, eval_, evec_, nrow, &
         work, lwork, iwork, ifail, info)

    ! dsyev returns ew in ascending order, so flip.
    eval_(1:ned) = eval_(ned:1:-1)
    evec_(1:nrow,1:ned) = evec_(1:nrow,ned:1:-1)

    if(info /= 0) then
       ! DEBUG
       write(20,*) A
       write(*,*) 'eigs failure. info = ', info
       write(*,*) 'info = ', info
       write(*,*) 'eval_ = ', eval_
       stop

       if(present(status)) then
          ! Signal failure and return.
          status = -1
          return
       end if       

       ! Otherwise, halt execution.
       stop
    end if

    ! Deallocate workspace.
    deallocate(A_copy, work)

    ! Provide eigenvalues and eigenvectors as desired.
    if(present(eval)) eval = eval_(1:ned)
    if(present(evec)) evec = evec_(1:nrow,1:ned)

  end subroutine eigs_dsyevx

  !=============================================================================================
  ! Compute the NED larges eigenvalues of the real symmetric matrix A.
  !=============================================================================================
  subroutine eigs(nrow, A, ned, eval, evec)
    ! Parmeters.
    integer, intent(in) :: nrow
      ! Dimension of matrix A.
    real(dp), dimension(nrow,nrow), intent(in) :: A
      ! Real symmetric matrix to compute eigenvalue decomposition of.
    integer, intent(in) :: ned
      ! Number of eigenvalues desired.
    real(dp), dimension(ned), intent(out), optional :: eval
      ! Eigenvalues
    real(dp), dimension(nrow,ned), intent(out), optional :: evec
      ! Eigenvectors, stored as columns.

    ! Local variables.
    real(dp), dimension(nrow) :: eval_
    real(dp), dimension(nrow,nrow) :: evec_
    character :: jobz
      ! argument for dsyev: 'N' for only ew, 'V' for ev too
    character, parameter :: uplo = 'U'
      ! 'U' means upper triangle of matrix is stored by dsyev
    integer :: lwork
      ! size of work array
    real(dp), dimension(:), allocatable :: work
      ! workspace
    integer :: liwork
      ! size of integer work array
    integer, dimension(:), allocatable :: iwork
      ! integer workspace
    integer :: info
      ! exit information
    integer :: il, iu
      ! Lower and upper indices of eigenvalues to be found.
    real(dp) :: vl, vu
      ! Upper and lower bounds on eigenvalues to be found (not used).
    real(dp), parameter :: abstol = 0.0
      ! Absolute tolerance to use in determining eigenvalues, or 0 if a suitable tolerance is to be computed by LAPACK.
    integer :: neval_found
      ! The number of eigenvalues found by LAPACK call.
    character(1) :: range
      ! Character identifying whether given index range il..iu of eigenvalues is to be found (if 'I') or given range (if 'V').
    real(dp), dimension(:,:), allocatable :: A_copy
      ! Copy of the matrix A, destroyed by dsyevr.
    integer, dimension(2 * ned) :: isuppz
      ! The support of the eigenvectors.

    ! DEBUG
    if(any(isnan(A))) then
!    gfortran chokes on this -- the following might work with it...    
!    if(any(.not. (A <= 0.0 .or. A > 0.0))) then
       write(*,*) 'A contains NaN'
       write(18,*) A
       stop
    end if

    ! Instruct dsyevr to find specified index range of eigenvalues.
    range = 'I'
    il = nrow-ned+1
    iu = nrow
    
    ! Compute only eigenvalues by default.
    jobz = 'N'
    ! If evec provided, compute eigenvectors too.
    if(present(evec)) jobz = 'V' 

    ! Make a copy of the matrix A which will be destroyed by dsyevr.
    allocate(A_copy(nrow,nrow))
    A_copy = A

    ! Inquire about optimal workspace size, lwork and liwork.
    lwork = -1
    liwork = -1
    allocate(work(1), iwork(1))
    call dsyevr(jobz, range, uplo, nrow, A_copy, nrow, vl, vu, il, iu, abstol, neval_found, eval_, evec_, nrow, &
         isuppz, work, lwork, iwork, liwork, info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)

    ! Allocate workspace.
    allocate(work(lwork), iwork(liwork))

    ! Perform eigenvalue decomposition.
    call dsyevr(jobz, range, uplo, nrow, A_copy, nrow, vl, vu, il, iu, abstol, neval_found, eval_, evec_, nrow, &
         isuppz, work, lwork, iwork, liwork, info)

    ! dsyev returns ew in ascending order, so flip.
    eval_(1:ned) = eval_(ned:1:-1)
    evec_(1:nrow,1:ned) = evec_(1:nrow,ned:1:-1)

    if(info /= 0) then
       ! DEBUG
       write(20,*) A
       write(*,*) 'eigs failure. info = ', info
       write(*,*) 'info = ', info
       write(*,*) 'eval_ = ', eval_
       stop
    end if

    ! Deallocate workspace.
    deallocate(A_copy, work, iwork)

    ! Provide eigenvalues and eigenvectors as desired.
    if(present(eval)) eval = eval_(1:ned)
    if(present(evec)) evec = evec_(1:nrow,1:ned)

  end subroutine eigs

  !=============================================================================================
  ! Compute all eigenvalues and eigenvectors of the real symmetric matrix A.
  !=============================================================================================
  subroutine eig(nrow, A, eval, evec)
    ! Parmeters.
    integer, intent(in) :: nrow
      ! Dimension of matrix A.
    real(dp), dimension(nrow,nrow), intent(in) :: A
      ! Real symmetric matrix to compute eigenvalue decomposition of.
    real(dp), dimension(nrow), intent(out) :: eval
      ! Eigenvalues, in descending order.
    real(dp), dimension(nrow,nrow), intent(out) :: evec
      ! evec(:,i) is eigenvector i corresponding to eval(i)

    ! Local variables.
    character :: jobz
      ! argument for dsyev: 'N' for only ew, 'V' for ev too
    character, parameter :: uplo = 'U'
      ! 'U' means upper triangle of matrix is stored by dsyev
    integer :: lwork
      ! size of work array
    real(dp), dimension(:), allocatable :: work
      ! workspace
    integer :: info
      ! exit information
    
    ! Compute eigenvectors too.
    jobz = 'V' 

    ! Copy A into evec_
    evec = A

    ! Inquire about optimal workspace size.
    lwork = -1
    allocate(work(1))
    call dsyev(jobz, uplo, nrow, evec, nrow, eval, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)

    ! Allocate workspace.
    allocate(work(lwork))

    ! Perform eigenvalue decomposition.
    call dsyev(jobz, uplo, nrow, evec, nrow, eval, work, lwork, info)

    ! dsyev returns ew in ascending order, so flip.
    eval = eval(nrow:1:-1)
    evec = evec(1:nrow,nrow:1:-1)

    if(info /= 0) then
       ! DEBUG
       write(*,*) 'info = ', info
       write(*,*) 'eval = ', eval
       stop
    end if

    ! Deallocate workspace.
    deallocate(work)

  end subroutine eig

end module eigs_lapack
