!=============================================================================================
! Constants for describing the precision.
!=============================================================================================
module numeric_kinds
  implicit none
  public

  ! Named constants for 4, 2, and 1 byte integers.
  integer, parameter ::                              &
       i4b = selected_int_kind(9),                   & ! 4-byte integer
       i2b = selected_int_kind(4),                   & ! 2-byte integer
       i1b = selected_int_kind(2)                      ! 1-byte integer

  ! Named constants for single, double, and quadruple precision reals.
  integer, parameter ::                              & ! single precision (4-byte floating-point)
       sp = kind(1.0),                               & ! double precision (8-byte floating-point)
       dp = selected_real_kind(2*precision(1.0_sp)), & ! quad precision (16-byte floating-point)
       qp = selected_real_kind(2*precision(1.0_dp))  

  !=============================================================================================
  ! isnan replacement function for gfortran
  !=============================================================================================
  interface isnan
     module procedure isnan_sp, isnan_dp
  end interface

contains

  !=============================================================================================
  ! isnan replacement function for gfortran
  !=============================================================================================
  elemental function isnan_sp(x)
    real(sp), intent(in) :: x
    logical :: isnan_sp
    
    isnan_sp = .false.
    if (.not. (x <= 0.0 .or. x > 0.0)) isnan_sp = .true.

  end function isnan_sp

  elemental function isnan_dp(x)
    real(dp), intent(in) :: x
    logical :: isnan_dp
    
    isnan_dp = .false.
    if (.not. (x <= 0.0 .or. x > 0.0)) isnan_dp = .true.

  end function isnan_dp

end module numeric_kinds
