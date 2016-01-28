!=============================================================================================
! A simple interface to for the TRLan package for finding a few eigenvalues and eigenvectors of
! a real symmetric matrix.
!
! TRLan can be obtained here:
!
! http://crd.lbl.gov/~kewu/trlan.html
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
! - Add sparse matrix-vector support.
!=============================================================================================
module eigs_trlan
  use numeric_kinds ! for precision

  implicit none

  private
  public :: eigs
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  logical, parameter :: debug = .FALSE.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

  real(dp), dimension(:,:), allocatable :: A_copy  ! symmetric matrix for multiplication

contains

  !=============================================================================================
  ! Matrix-vector operation for dense symmetric matrix.
  !=============================================================================================
  subroutine dense_sym_mv_op(nrow, ncol, xin, ldx, yout, ldy)
    implicit none

    ! Parameters.
    integer, intent(in) :: nrow 
    integer, intent(in) :: ncol 
    integer, intent(in) :: ldx
    integer, intent(in) :: ldy
    real(dp), dimension(nrow,ncol), intent(in) :: xin
    real(dp), dimension(nrow,ncol), intent(out) :: yout

    ! Local variables.
    integer :: i, j, ioff, joff
    
    ! Multiply matrix A times all column vectors desired.
    ! TODO: Use optimized BLAS?
    yout = matmul(A_copy, xin)

  end subroutine dense_sym_mv_op

  !=============================================================================================
  ! Compute the NED larges eigenvalues of the real symmetric matrix A.
  !=============================================================================================
  subroutine eigs(nrow, A, ned, eval, evec, starting_guess)
    use trl_info ! TRLan
    use trl_interface ! TRLan
    
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
    real(dp), dimension(nrow), intent(in), optional :: starting_guess

    ! Local variables.
    type(trl_info_t) :: info
      ! TRLan data structure.
    integer :: maxlan
      ! Maximum Lanczos basis size.
    integer, parameter :: lohi = 1
      ! Set to 1 to compute largest eigenvalues.
    real(dp), dimension(ned) :: eval_
    real(dp), dimension(nrow,ned) :: evec_
    integer :: i

    ! Make local copy of matrix A.
    allocate(A_copy(nrow,nrow))
    A_copy = A

    ! Choose the maximum Lanczos basis size.
    maxlan = ned + min(6, ned)
    maxlan = maxlan * 2
    
    ! Initialize info data structure.
    ! For info on choices, see: http://crd.lbl.gov/~kewu/ps/trlan_.html#SEC5
    call trl_init_info(info, nrow, maxlan, lohi, ned)

    ! Let TRLan know we will provide an initial starting vector guess.
    call trl_set_iguess(info, 1, -1)

    ! Initialize eigenvalues and eigenvectors.
    eval_ = 0
    evec_ = 0

    eval_(1) = 1.0
    evec_(1:nrow,1) = 1.0
    ! Use starting guess, if provided.
    if(present(starting_guess)) evec_(1:nrow,1) = starting_guess
    do i = 2,ned
       evec_(i,i) = 1
    end do

    write(*,*) evec_(1:10,1)

    ! Call TRLAN to compute the eigenvalues and eigenvectors.
    call trlan(dense_sym_mv_op, info, nrow, ned, eval_, evec_, nrow)

    ! DEBUG
    call trl_print_info(info, 2 * nrow**2)

    ! Provide eigenvalues and eigenvectors as desired.
    if(present(eval)) eval = eval_
    if(present(evec)) evec = evec_

    ! Free memory.
    deallocate(A_copy)

  end subroutine eigs

end module eigs_trlan
