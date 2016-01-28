!=============================================================================================
! Test driver for file reading.
!=============================================================================================
program main
  use numeric_kinds
  use constants ! global constants
  use pdbio, only : new_unit 

  implicit none

  ! Constants.
  character(len=*), parameter :: filename = '../data/pentaalanine/ala5.data'
  
  ! Local variables.
  integer :: iunit
  character(len=MAX_LINE_LENGTH) :: line, token
  real(sp), dimension(108) :: array
  integer :: iostat
  integer :: line_index
  integer :: numbers_read, next, line_length

  ! Open file
  call new_unit(iunit)
  open(unit=iunit, file=filename, status='OLD', err=11, position='REWIND', &
       form='FORMATTED', action='READ')

  line_index = 0
  do
     ! Read a line
     line_index = line_index + 1
     read(unit=iunit, fmt='(a)', iostat=iostat, end=20) line     

     ! TODO: Terminate loop if EOF encountered?  Here or after parsing?
     ! if(iostat < 0) exit
     if(mod(line_index, 1000) == 0) write(*,*) 'read line ', line_index
!     write(*,*) trim(line)

     ! Process into tokens.
     numbers_read = 0
     next = 1
     line_length = len(trim(line))
     do
        ! Get delimited token.
        call gettext(line, token, next)
 
        ! Parse token
        !write(*,*) 'token: ', trim(token), ' next = ', next, ' len = ', line_length
        numbers_read = numbers_read + 1
        read(token,'(F16.8)') array(numbers_read)

        if(next .ge. line_length) exit        
     end do
     
     !write(*,'(108(F6.3,1X))') array(1:108)
  end do

  ! End of file has been reached.
20 continue
  
  ! Close file.
  close(iunit)

  stop

  ! Report error if file could not be opened.
11 continue
  write(*,*) 'Error opening file ', filename
  stop
  
end program main

