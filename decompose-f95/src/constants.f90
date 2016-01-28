!=============================================================================================
! Useful constants used by several modules.
!=============================================================================================
module constants
  implicit none
  public

  integer, parameter :: MAX_FILENAME_LENGTH = 4096  ! Maximum length of filenames
  integer, parameter :: MAX_LINE_LENGTH = 10*1024   ! Maximum length of a line of a file.

end module constants
