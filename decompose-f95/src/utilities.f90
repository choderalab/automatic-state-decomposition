!=============================================================================================
! Utility functions and subroutines.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
! TODO:
! - Have all calls for uniform random deviates go through one random number subroutine here.
! - Make that call $omp single protected.
!=============================================================================================
module utilities
  use numeric_kinds ! for precision

  implicit none

  private
  public :: setRandomSeed, histogram, generateUniform01, generateUniformInteger, generateUniformIntegers, &
       generateUniqueRandomIndices, sort, find, logSum, mean, std, unique, &
       compute_confidence_bounds, trace

  !=============================================================================================
  ! GENERIC INTERFACES
  !=============================================================================================
  interface mean
     module procedure mean_sp, mean_dp
  end interface

  interface std
     module procedure std_sp, std_dp
  end interface

contains

  !=============================================================================================
  ! Set the random number seed based on the current time and date.
  !=============================================================================================
  subroutine setRandomSeed

    ! Local variables.
    integer :: n
    ! seed size
    integer :: i
    ! loop variable
    integer, dimension(:), allocatable :: seed
    ! seed array
    integer, dimension(8) :: idate
    ! current date info

    ! Query seed array.
    call random_seed(size = n)

    ! Allocate seed array.
    allocate(seed(n))

    ! Retrieve seed array.
    call random_seed(get = seed)

    ! Initialize the random number seed with as many numbers as we have, starting with the
    ! most rapidly varying number -- idate(8) -- first.
    call date_and_time(values = idate)
    do i = 1, min(8, n)
       seed(i) = idate(8-i+1)
    end do

    ! Set seed
    call random_seed(put = seed)

    ! Deallocate
    deallocate(seed)

  end subroutine setRandomSeed

  !=============================================================================================
  ! Return the histogram count.
  !=============================================================================================
  pure function histogram(vector, nbins)
    ! Arguments.
    integer, dimension(:), intent(in) :: vector 
               ! the vector of entries to tally the histogram for, which must be in 1..nbins
    integer, intent(in) :: nbins           ! the number of bins in the histogram
    integer, dimension(nbins) :: histogram ! the resulting histogram

    ! Local variables.
    integer :: index
    integer :: bin_index

    ! Construct histogram count.
    histogram = 0
    do index = 1,size(vector)
       bin_index = vector(index)
       histogram(bin_index) = histogram(bin_index) + 1
    end do

  end function histogram

  !=============================================================================================
  ! Generate a uniform real on [0,1].
  !=============================================================================================
  function generateUniform01()
    
    ! Return values.
    real(dp) :: generateUniform01

    ! Generate random number.
    call random_number(generateUniform01)    

  end function generateUniform01

  !=============================================================================================
  ! Generate a uniform integer on the range 1...N.
  !=============================================================================================
  function generateUniformInteger(N)
    
    ! Parameters.
    integer, intent(in) :: N ! the maximum index to generate

    ! Return values.
    integer :: generateUniformInteger

    ! Local variables.
    real(dp) :: r

    ! Choose an index uniformly over range i..N
    r = generateUniform01()
    generateUniformInteger = floor(r * dble(N)) + 1

  end function generateUniformInteger

  !=============================================================================================
  ! Generate uniform integers on the range 1...N.
  !=============================================================================================
  subroutine generateUniformIntegers(N, x)
    
    ! Parameters.
    integer, intent(in) :: N ! the maximum index to generate
    integer, dimension(:), intent(out) :: x ! the vector to fill with integers

    ! Local variables.
    real(dp) :: r
    integer :: i

    ! Choose an index uniformly over range i..N
    ! Ensure only one thread at a time calls the random number generator.
    do i = 1,size(x,1)
       r = generateUniform01()
       x(i) = floor(r * dble(N)) + 1
    end do

  end subroutine generateUniformIntegers

  !=============================================================================================
  ! Exchange two numbers.
  !=============================================================================================  
  subroutine exchange(a, b)
    
    ! Parameters.
    integer, intent(inout) :: a, b
      ! numbers to be exchanged

    ! Local variables.
    integer :: tmp

    tmp = a
    a = b
    b = tmp

  end subroutine exchange
  
  !=============================================================================================
  ! Generate a list of K unique indices in the specified range 1...N.
  !
  ! The number K of indices to generate is determined by the dimension of 'indices'.
  !
  ! If provided, UNUSED_INDICES will be filled with the indices from 1...N not contained in INDICES.
  !
  ! WARNING: This method requires O(N) storage for scratch space.  This may be inefficient for K << N.
  ! Note that the random seed should be set if the desired sequence is to change.
  !=============================================================================================
  subroutine generateUniqueRandomIndices(N, K, indices, unused_indices)
    
    ! Parameters.
    integer, intent(in) :: N ! the maximum index to generate
    integer, intent(in) :: K ! number of unique indices to generate
    integer, dimension(K), intent(out) :: indices ! storage for the K unique generated indices
    integer, dimension(N-K), intent(out), optional :: unused_indices ! unused indices from 1...N not contained in INDICES

    ! Local variables.
    integer, dimension(N) :: allindices_permuted ! A list of unused indices.
    real(sp) :: r ! random number
    integer :: i, j ! loop indices
    integer :: tmpindex ! temporary index
    
    ! :TODO: Ensure that K <= N and that size(indices) >= K, or else report an error.
    
    ! Start by assigning the elements of permuation to the integers 1..N.
    allindices_permuted = (/ (i, i = 1,N) /)
    
    ! Generate K unique indices by choosing uniformly over the unselected indices, storing them in the first elements of allindices_permuted.
    do i = 1,K
       ! Choose an index uniformly over range i..N
       r = generateUniform01()
       j = floor(r * dble(N-i+1)) + i

       ! Swap index i with index j.
       call exchange(allindices_permuted(i), allindices_permuted(j))
    end do
    
    ! Recover the first K unique indices for the caller.
    indices(1:K) = allindices_permuted(1:K)

    ! If specified, provide unused indices.
    if(present(unused_indices)) then
       unused_indices = allindices_permuted((K+1):N)
    end if
    
  end subroutine generateUniqueRandomIndices

  !=============================================================================================
  ! Return a pointer to an array containing the indices of all the .true. elements in the given mask. 
  !=============================================================================================
  function find(mask) result(indices)
    
    ! Parameters.
    logical, dimension(:), intent(in) :: mask
      ! The logical array.
    
    ! Return value.
    integer, pointer, dimension(:) :: indices
      ! The returned integer array contains all the indices of mask that were .true.

    ! Local variables.
    integer :: nindices
      ! Number of indices to store.
    integer :: i, j

    ! Count number of indices.
    nindices = count(mask)

    ! Allocate storage.
    allocate(indices(nindices))
    
    ! Get indices.
    indices = pack((/ (i, i = 1,size(mask,1)) /), mask )

  end function find

  !=============================================================================================
  ! Sort the given array, arranging elements in ascending order.
  !=============================================================================================
  subroutine sort(a)
    ! Parameters.
    real(sp), dimension(:), intent(inout) :: a
      ! The array to be sorted.

    ! Call quicksort on the entire array span.
    call quicksort(a, 1, size(a,1))

    return

  contains

    ! Recursive quicksort subroutine.
    recursive subroutine quicksort(a, left, right)

      ! Parameters.
      real(sp), dimension(:), intent(inout) :: a 
      ! The array to be sorted.
      integer, intent(in) :: left, right
      ! Left and rightmost indices of span to be sorted.
      
      ! Local variables.
      integer :: pivotNewIndex
      
      if (left < right) then
         ! Select a new pivot point.
         pivotNewIndex = partition(a, left, right)
         
         ! Quicksort left and right of the pivot point.
         call quicksort(a, left, pivotNewIndex - 1)
         call quicksort(a, pivotNewIndex + 1, right)
      endif
      
    end subroutine quicksort

    ! Swap specified elements of an array.
    subroutine swap(a, i, j)
      
      ! Parameters.
      real(sp), dimension(:), intent(inout) :: a
      integer, intent(in) :: i, j

      ! Local variables.
      real(sp) :: tmp
      
      ! Swap
      tmp = a(i)
      a(i) = a(j)
      a(j) = tmp

    end subroutine swap

    ! In-place partitioning between indices left and right, inclusive.
    integer function partition(a, left, right)

      ! Parameters.
      real(sp), dimension(:), intent(inout) :: a
      integer, intent(in) :: left, right
      
      ! Local variables.
      integer :: pivotIndex
      real(sp):: pivotValue
      integer :: storeIndex
      integer :: i

      ! Choose a pivot index and get its value.
      pivotIndex = left
      pivotValue = a(pivotIndex)

      ! Move pivot to the end.
      call swap(a, pivotIndex, right)
      
      ! Store left index.
      storeIndex = left
      
      do i = left, right-1
         if(a(i) .le. pivotValue) then
            call swap(a, storeIndex, i)
            storeIndex = storeIndex + 1
         end if
      end do
      
      ! Move pivot to its final place.
      call swap(a, right, storeIndex)
      
      ! Return new pivot index.
      partition = storeIndex
    end function partition
  end subroutine sort

  !=============================================================================================
  ! Compute confidence bounds.
  !=============================================================================================
  subroutine compute_confidence_bounds(x_n, confidence_interval, lower, upper)

    ! Parameters.
    real(sp), dimension(:), intent(in) :: x_n
      ! x_n(n) is the nth sample from a random sample.
    real(sp), intent(in) :: confidence_interval
      ! The width of the desired symmetric confidence interval.
    real(sp), intent(out) :: lower, upper
      ! Lower and upper confidence bounds.

    ! Local variables.
    integer :: N
      ! Number of samples
    real(sp), dimension(size(x_n)) :: x_sorted
      ! Sorted copy of x.
    integer :: lower_index, upper_index
      ! Indices to take lower and upper confidence bounds from.

    ! Determine number of samples.
    N = size(x_n,1)

    ! Sanity check on requested confidence interval.
    if(confidence_interval <= 0.0 .or. confidence_interval > 1.0) then
       write(*,*) 'utilities.f90: confidence_interval: requested confidence interval ', confidence_interval, &
            ' exceeds allowed range'
       stop
    end if
    ! Sanity check on number of samples.
    ! TODO: Ensure N is large enough to give decent estimate of confidence interval.
    if(N < 1) then
       write(*,*) 'utilities.f90: confidence_interval: Too few samples (N = ', N, ') for confidence estimate.'
       stop
    end if
    
    ! Sort x in ascending order.
    x_sorted = x_n
    call sort(x_sorted)

    ! Determine indices to take lower and upper confidence bounds from.
    lower_index = ceiling( (N-1) * (0.5 - confidence_interval/2.0) ) + 1
    upper_index = floor( (N-1) * (0.5 + confidence_interval/2.0) ) + 1

    ! Report confidence interval.
    lower = x_sorted(lower_index)
    upper = x_sorted(upper_index)

  end subroutine compute_confidence_bounds

  !=============================================================================================
  ! Return the unique values in ARRAY (in order of appearance), and if desired, the compacted indices of ARRAY.
  !
  ! The result is returned in a pointer, which must be deallocated after use.
  !=============================================================================================
  function unique(array, compacted)
    
    ! Parameters.
    integer, dimension(:), intent(in) :: array
      ! The array whose unique elements are to be found.
    integer, dimension(size(array,1)), intent(out), optional :: compacted
      ! Storage for compacted translated indices that run from 1...nunique.
    
    ! Return value.
    integer, dimension(:), pointer :: unique
      ! An array containing the unique elements of ARRAY, in order of appearance.

    ! Local variables.
    integer, dimension(size(array,1)) :: unique_entries
    integer, dimension(size(array,1)) :: compacted_entries
    integer :: nunique
    integer :: i, j
    integer :: N
    integer :: element
    logical :: found

    ! Get size of array.
    N = size(array,1)

    ! Go through list of elements.
    nunique = 0
    do i = 1,N
       ! Look for the element of array in our list of unique states so far.
       found = .false.
       do j = 1,nunique
          if(array(i) == unique_entries(j)) then
             found = .true.
             compacted_entries(i) = j
             exit
          end if
       end do

       ! If not found, add it to the list.
       if(.not. found) then
          nunique = nunique + 1
          unique_entries(nunique) = array(i)
          compacted_entries(i) = nunique
       end if
    end do
    
    ! Allocate storage for result.
    allocate(unique(nunique))
    unique(1:nunique) = unique_entries(1:nunique)

    ! Store compacted entries if desired.
    if(present(compacted)) compacted = compacted_entries
    
  end function unique

  !=============================================================================================
  ! Compute the log of a sum of terms whose logarithms are provided.
  !=============================================================================================
  function logSum(arg_i)

    ! Arguments
    real(dp), dimension(:), intent(in) :: arg_i ! the logs of the terms to be summed
    real(dp) :: logSum                          ! the log of the sum of the exponentials of arg_i
    
    ! Local variables
    real(dp) :: max_arg

    ! Compute the maximum argument.
    max_arg = maxval(arg_i)

    ! Compute the log sum
    logSum = log( sum( exp(arg_i - max_arg) ) ) + max_arg
        
  end function logSum

  !=============================================================================================
  ! Compute the mean.
  !=============================================================================================
  function mean_sp(x) result(mean)
    
    ! Arguments
    real(sp), dimension(:), intent(in) :: x

    ! Return values.
    real(sp) :: mean

    mean = sum(x) / real(size(x,1),sp)

  end function mean_sp

  function mean_dp(x) result(mean)
    
    ! Arguments
    real(dp), dimension(:), intent(in) :: x

    ! Return values.
    real(dp) :: mean

    mean = sum(x) / real(size(x,1),dp)

  end function mean_dp
   
  !=============================================================================================
  ! Compute the sample standard deviation along the specified dimension.
  !=============================================================================================
  function std_sp(x) result(std)
    ! Arguments.
    real(sp), dimension(:), intent(in) :: x
    
    ! Local variables.
    real(sp) :: std

    ! Compute std.
    std = sqrt( sum( (x - mean(x))**2 ) / real(size(x,1)-1,sp) )
  end function std_sp

  function std_dp(x) result(std)
    ! Arguments.
    real(dp), dimension(:), intent(in) :: x
    
    ! Local variables.
    real(dp) :: std

    ! Compute std.
    std = sqrt( sum( (x - mean(x))**2 ) / real(size(x,1)-1,dp) )
  end function std_dp

  !=============================================================================================
  ! Compute the trace of a matrix.
  !=============================================================================================
  function trace(Tji)
    ! Arguments.
    real(dp), dimension(:,:), intent(in) :: Tji

    ! Return variables
    real(dp) :: trace

    ! Local variables.
    integer :: i

    trace = 0.0
    do i = 1,size(Tji,1)
       trace = trace + Tji(i,i)
    end do

  end function trace


end module utilities
