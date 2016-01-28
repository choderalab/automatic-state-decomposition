!=============================================================================================
! State-decomposition code.
! Written by John D. Chodera and Nina Singhal, 2005.
!=============================================================================================

!=============================================================================================
! Include files.
!=============================================================================================

! Constants defining the precision.
!include 'precision.f90'

! MPI inlcude files
!include 'mpi_constants.f90'
!include 'mpi1.f90'

! Least-squares RMSD and alignment.
!include 'ls_rmsd.f90'

! A distance metric must be defined.
!include 'distance_metric.f90'

!=============================================================================================
! A module to apply the weighted histogram analysis method (WHAM) to a number of independent
! canonical simulations, independent simulated tempering simulations, or parallel tempering
! simulations.  All simulations must be of the same length.
! NOTE: All 'real' variables should be compiled to double precision.
!=============================================================================================
module decompose
  implicit none

  private
  public :: initialize, finalize, testRmsd
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================

  integer, parameter, public :: snapshotReal = selected_real_kind(15,90) ! TODO: Make 4-byte float
  integer, parameter, public :: distanceReal = selected_real_kind(15,90) ! TODO: Make 4-byte float
  integer, parameter, public :: longReal = selected_real_kind(15,90)
  integer, parameter, public :: longInt = selected_int_kind(10)

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  ! Snapshot metadata.
  type SnapshotMetadataType
     integer id          ! unique identifier of this snapshot, compact starting from 1
     integer node        ! nodeId of node the snapshot is stored on
     integer localIndex  ! the local index of the snapshot on the node it is stored on
     integer microstate  ! the index of the microstate cluster the snapshot currently belongs to
     integer macrostate  ! the index of the macrostate the snapshot currently belongs to
  end type SnapshotMetadataType

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

  integer :: nsnapshots ! total number of snapshots
  integer :: numberOfAtoms ! total number of atoms per snapshot
  integer :: nodeId ! nodeId of this node (starting from 0)
  integer :: totalNodes ! number of processors
  type(SnapshotMetadataType), dimension(:), pointer :: snapshotMetadata
  real, dimension(:,:), pointer :: snapshot ! snapshot_ik(i,k) is the kth component of snapshot i

contains

  !=============================================================================================
  ! Generate a list of K unique indices in the specified range 1...N.
  !=============================================================================================
  subroutine generateUniqueRandomIndices(N, indices)
    
    ! Parameters.
    integer, intent(in) :: N ! the maximum index to generate
    integer, dimension(:), intent(out) :: indices ! storage for the K unique generated indices

    ! Local variables.
    integer :: K ! the number of indices to generate
    double precision :: r ! random number
    logical :: isUnique ! flag to indicate whether or not the chosen index is unique
    integer :: i, j ! loop indices
    
    ! Determine K from size of indices argument.
    K = size(indices)

    ! :TODO: Ensure that N >= K or else report an error.

    ! Generate K unique indices from 1..N using a simple algorithm.
    do i = 1,K
       isUnique = .false.
       do while (.not. isUnique)

          ! Choose an index uniformly over 1..N
          call random_number(r)
          indices(i) = floor(r * N) + 1
          
          ! Check to see if this index has not already been selected.
          isUnique = .true.
          do j = 1,(i-1)
             if(indices(i) == indices(j)) then
                isUnique = .false.
                exit
             end if
          end do

       end do       
    end do
    
  end subroutine generateUniqueRandomIndices
    
  !=============================================================================================
  ! Read the data.
  !=============================================================================================
  subroutine initialize(filename, nsnapshots_, numberOfAtoms_, nodeId_, totalNodes_)

    ! Parameters.
    character(len=*), intent(in) :: filename   ! name of binary-formatted snapshot data store
    integer, intent(in) :: nsnapshots_         ! the number of snapshots in the file
    integer, intent(in) :: numberOfAtoms_      ! the number of atoms
    integer, intent(in) :: nodeId_           ! nodeId of MPI process (counting from 0)
    integer, intent(in) :: totalNodes_ ! the number of processors

    ! Constants.
    integer, parameter :: unit = 2 ! unit to use for file access

    ! Local variables.
    integer :: i
    integer :: nsnapshotsThisNode
    integer :: snapshotRecordLength ! the size of a snapshot record in the binary snapshot store

    ! Store the node ID and number of processes.
    nsnapshots = nsnapshots_
    numberOfAtoms = numberOfAtoms_
    nodeId = nodeId_
    totalNodes = totalNodes_
   
    ! Create a master list of snapshots and their locations.
    allocate(snapshotMetadata(nsnapshots))
    do i = 1,nsnapshots
       snapshotMetadata(i)%id = i
       snapshotMetadata(i)%node = mod(i-1,totalNodes)
       snapshotMetadata(i)%localIndex = i/totalNodes + 1
       snapshotMetadata(i)%microstate = 1
    end do

    ! Count number of snapshots to be stored on this node.
    nsnapshotsThisNode = 0
    do i = 1,nsnapshots
       if(snapshotMetadata(i)%node == nodeId) then
          nsnapshotsThisNode = nsnapshotsThisNode + 1
       end if       
    end do

    ! Allocate storage for snapshots to be stored locally.
    write(*,*) 'Allocating storage for ', nsnapshotsThisNode, ' snapshots on node ', nodeId
    allocate(snapshot(nsnapshotsThisNode, numberOfAtoms*3))

    ! Load those snapshots assigned to our node.
    ! Determine record length of one snapshot entry.
    inquire(iolength=snapshotRecordLength) snapshot(1,:)
    ! Open the file for read-only, random access.
    open(unit, file=filename, status='old', access='direct', form='unformatted', action='read', recl=snapshotRecordLength)
    ! Read those snapshots that are assigned to this node into the local slot in the snapshot store.
    do i = 1,nsnapshots
       if( snapshotMetadata(i)%node == nodeId ) then
          read(unit, rec=snapshotMetadata(i)%id, err=1) snapshot(snapshotMetadata(i)%localIndex,:)
       end if
    end do
    ! Close the file.    
    close(unit)
    return

1   write(*,*) 'node ', nodeId, ' encountered a problem reading snapshot ', i, ' record ', snapshotMetadata(i)%id

  end subroutine initialize

  !=============================================================================================
  ! Test RMSD routines.
  !=============================================================================================
  subroutine testRmsd()
    use ls_rmsd ! least-squared RMSD using quaternions

    ! Local variables.
    integer :: i
    integer :: option
    logical :: calc_g
    double precision :: error
    double precision, dimension(3,numberOfAtoms) :: coord1, coord2, g
    double precision, dimension(3) :: x_center, y_center
    integer, dimension(2) :: reshape_array
    double precision, dimension(3,3) :: U

    ! Set options.
    option = 0
    calc_g = .false.

    ! Construct reshape array.
    reshape_array(1) = 3
    reshape_array(2) = numberOfAtoms
    
    do i = 1,100
       ! Copy coordinates of successive snapshots.
       coord1 = reshape(snapshot(i,:), reshape_array )
       coord2 = reshape(snapshot(i+1,:), reshape_array )
       
       ! Compute the rmsd between snapshots.
       ! call rmsd(numberOfAtoms, coord1, coord2, option, U, x_center, y_center, error, calc_g, g)    
       error = rmsd(numberOfAtoms, coord1, coord2)
       write(*,*) error
    end do
    
  end subroutine testRmsd

  !=============================================================================================
  ! Align one snapshot to another using least-squares rigid-body superposition.
  !=============================================================================================
  subroutine alignSnapshot(snapshot, referenceSnapshot)
    ! Parameters.
    real(snapshotReal), dimension(:) :: snapshot ! The snapshot to be aligned to the reference.
    real(snapshotReal), dimension(:) :: referenceSnapshot ! The reference snapshot to which snapshot is to be aligned.
    
    ! Local variables.
    integer :: numberOfAtoms
    
    ! Determine number of atoms.
    numberOfAtoms = size(snapshot) / 3

    ! TODO: Perform rigid-body superposition.
    

  end subroutine alignSnapshot
  
  !=============================================================================================
  ! Compute distances.
  !=============================================================================================
  function computeDistancesToGenerators(numberOfAtoms, snapshot, numberOfGenerators, generators)
    ! Parameters.
    integer, intent(in) :: numberOfAtoms
    real(snapshotReal), dimension(3,numberOfAtoms), intent(in) :: snapshot
    integer, intent(in) :: numberOfGenerators
    real(snapshotReal), dimension(numberOfGenerators, 3, numberOfAtoms), intent(in) :: generators
    real(distanceReal), dimension(numberOfGenerators) :: computeDistancesToGenerators
    
    ! Local variables.
    integer :: generatorIndex
    
    forall(generatorIndex = 1:numberOfGenerators)
       computeDistancesToGenerators(generatorIndex) = distance(snapshot, generators(generatorIndex,:,:))
    end forall
    
  end function computeDistances

  !=============================================================================================
  ! Assign snapshots to generators in parallel.
  !=============================================================================================
  subroutine assignSnapshotsToGenerators(snapshots, generators, shapshotAssignments)
    ! Parameters.
    type(SnapshotType), dimension(:,:), intent(in) :: snapshots  ! The snapshots to be assigned.
    type, dimension(:,:), intent(in) :: generators ! The generators to assign them to.
    integer, dimension(:), intent(out) :: snapshotAssignments ! the
    
    ! Local variables.
    integer :: numberOfSnapshots
    integer :: numberOfGenerators
    integer :: snapshotIndex
    integer :: generatorIndex
    real(distancePrecision) :: generatorDistance
    real(distancePrecision) :: closestGeneratorDistance
    integer :: closestGeneratorIndex

    ! Get dimensions.
    numberOfSnapshots = size(snapshots,1)
    numberOfGenerators = size(generators,1)
    
    ! Assign each snapshot to the closest generator.
    do snapshotIndex = 1,numberOfSnapshots
       ! Determine closest generator.
       closestGeneratorIndex = 1
       closestGeneratorDistance = distance(snapshots(snapshotIndex,:), generators(1,:))
       do generatorIndex = 2,numberOfGenerators
          generatorDistance = distance(snapshots(snapshotIndex,:), generators(generatorIndex,:))
          if(generatorDistance < closestGeneratorDistance) then
             closestGeneratorIndex = generatorIndex
             closestGeneratorDistance = generatorDistance
          end if
       end do
       
    end do

  end subroutine assignSnapshotsToGenerators
  

  !=============================================================================================
  ! Clean up allocated data.
  !=============================================================================================
  subroutine finalize

    ! Free snapshot table.
    deallocate(snapshotMetadata)
    ! Free snapshot storage.
    deallocate(snapshot)

  end subroutine finalize

end module decompose

!=============================================================================================
! Main driver.
!=============================================================================================
program main
  use mpi_constants
  use mpi1
  use decompose

  implicit none

  ! Constants.
  character(len=*), parameter :: snapshotFilename = 'trpzip2-100ps.f90trj' ! snapshot file
  integer, parameter :: nsnapshots = 32600
  integer, parameter :: numberOfAtoms = 145

  ! Local variables.
  integer count
  real data(0:99)
  integer dest
  integer i
  integer ierr
  integer totalNodes
  integer nodeId
  integer status(MPI_Status_size)
  integer tag
  real value(200)

  ! Initialize MPI.
  call MPI_Init(ierr)
    
  ! Determine this process's rank.
  call MPI_Comm_Rank(MPI_COMM_WORLD, nodeId, ierr)

  ! Get number of processors.
  call MPI_Comm_size (MPI_COMM_WORLD, totalNodes, ierr)

  ! DEBUG.
  write(*,*) 'initializing snapshots...'
  call initialize(snapshotFilename, nsnapshots, numberOfAtoms, nodeId, totalNodes)
  write(*,*) 'done.'
  
  call testRmsd()
  
  !  Process 0 expects to receive as much as 200 real values, from any source.
!  if ( nodeId == 0 ) then
!     tag = 55
!     call MPI_Recv(value, 200, MPI_REAL, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)
!     
!     write ( *, '(a,i1,a,i1)' ) 'P:', nodeId, ' Got data from processor ', status(MPI_SOURCE)
!     
!     call MPI_Get_count(status, MPI_REAL, count, ierr)
!!     
!     write ( *, '(a,i1,a,i3,a)' ) 'P:', nodeId, ' Got ', count, ' elements.'
!     
!     write ( *, '(a,i1,a,g14.6)' ) 'P:', nodeId, ' value(5) = ', value(5)
!     
!     !  Process 1 sends 100 real values to process 0.
!  else if ( nodeId == 1 ) then
!     
!     write ( *, '(a)' ) ' '
!!     write ( *, '(a,i1,a)' ) 'P:', nodeId, ' - setting up data to send to process 0.'
!     
!     do i = 0, 99
!        data(i) = real ( i )
!     end do
!     
!     dest = 0
!     tag = 55
!     call MPI_Send ( data, 100, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr )
!     
!  else
!     
!     write ( *, '(a)' ) ' '
!     write ( *, '(a,i1,a)' ) 'P:', nodeId, ' - MPI has no work for me!'
!     
!  end if
  
 ! Clean up.

  call MPI_Finalize ( ierr )

  if ( nodeId == 0 ) then
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end program main
