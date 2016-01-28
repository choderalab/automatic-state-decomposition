!=============================================================================================
! A module defining the distance metric to be used in conformational clustering of snapshots.
! This module is user-defined.  The user can choose whichever metric is most appropriate.
!
! This module currently uses the LS-RMSD algorithm of Vageli Coutsias and Chaok Seok, which
! uses quaternions, as implemented in the ls_rmsd module.
!=============================================================================================
module distance_metric
  use DefinePrecision
  use ls_rmsd
  implicit none

  private
  public :: initialize_distance, distance

  !=============================================================================================
  ! Private module data.
  !=============================================================================================
  integer :: numberOfAtoms             ! The number of atoms.
  integer, dimension(:), allocatable :: atomIndices ! The atom indices to include in evaluating the distance.
  real(snapshotReal), dimension(:,:), allocatable :: coordinateArray1, coordinateArray2 ! Scratch arrays.  TODO: pointers?

contains

  !=============================================================================================
  ! Initialize distance calculation with desired atomlist.
  !
  ! TODO: This interface should eventually be generalized to pass some generic user data,
  ! from which the atom indices to use can be computed.
  !=============================================================================================
  subroutine initialize_distance(numberOfAtoms_, atomIndices_)
    ! Parameters.
    integer, intent(in) :: numberOfAtoms_ ! The number of atoms.
    integer, intent(in), dimension(:) :: atomIndices_ ! The indices of atoms to include in evaluating the distance.

    ! Local variables.
    
    ! Store number of atoms.
    numberOfAtoms = numberOfAtoms_

    ! Store copy of atom indices.
    allocate(atomIndices(size(atomIndices_,1)))
    
    ! Allocate coordinate arrays.
    allocate(coordinateArray1(3,size(atomIndices)),coordinateArray2(3,size(atomIndices)))       

  end subroutine initialize_distance

  !=============================================================================================
  ! Return the distance between two snapshots.
  ! This implementation returns the least-squared RMSD for the initialized atom selection.
  !=============================================================================================
  function distance(snapshot1, snapshot2)

    ! Parameters.
    real(distanceReal) :: distance
    real(snapshotReal), dimension(:,:), intent(in) :: snapshot1, snapshot2 ! The snapshots to compute distance between, 3xN.

    ! Local variables.
    integer :: numberOfAtoms
    ! double precision, dimension(:,:) :: coord1, coord2
    integer, dimension(2) :: reshape_array

    ! Construct reshape array.
    reshape_array(1) = 3
    reshape_array(2) = numberOfAtoms

    ! Transform coordinates into proper format.
    coordinateArray1 = reshape(snapshot1, reshape_array )
    coordinateArray2 = reshape(snapshot2, reshape_array )

    ! Compute the LS-RMSD between two snapshots.
    distance = rmsd(numberOfAtoms, coordinateArray1, coordinateArray2)

  end function distance

end module distance_metric
