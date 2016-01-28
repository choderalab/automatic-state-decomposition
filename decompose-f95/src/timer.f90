!=============================================================================================
! Utility timing functions.
!
! Written by John Chodera.
!=============================================================================================
module timer
  use numeric_kinds ! for precision

  implicit none

  private
  public :: resetTimer, readTimer

  !=============================================================================================
  ! MODULE DATA
  !=============================================================================================

  integer :: count_start
  integer :: count_rate
  integer :: count_max

contains

  !=============================================================================================
  ! Reset the timer.
  !=============================================================================================
  subroutine resetTimer

    ! Get the start time.
    call system_clock(count_start, count_rate, count_max)       

  end subroutine resetTimer

  !=============================================================================================
  ! Return the seconds elapsed since the time was last reset.
  !=============================================================================================
  function readTimer() result(elapsed_seconds)
    
    ! Return value.
    real(dp) :: elapsed_seconds

    ! Parameters.
    integer :: count_stop

    ! Stop timer.
    call system_clock(count_stop, count_rate, count_max)
    elapsed_seconds = real(count_stop - count_start,dp)/real(count_rate,dp)
    
  end function readTimer

end module timer
