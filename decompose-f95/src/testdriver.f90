!=============================================================================================
! TEST CODE
!=============================================================================================
module tests

  use numeric_kinds ! for precision
  use decompose

  implicit none

  public

contains

  !---------------------------------------------------------------------------------------------
  ! Test RMSD code by computing RMSD between consecutive snapshots.
  !---------------------------------------------------------------------------------------------
  subroutine testRmsd()
    use kabsch_rmsd, only : kabsch_ls_rmsd => ls_rmsd ! least-squared RMSD
    use theobald_rmsd, only : theobald_ls_rmsd => ls_rmsd ! least-squared RMSD
    use kabsch_rmsd, only : ls_align
    use pdbio

    ! Local variables.
    integer :: i, n, k ! Loop variable
    real(sp) :: rmsd ! rmsd between two snapshots
    real(dp) :: avgrmsd ! cumulative average
    real(dp) :: count
    integer, parameter :: snapshot1 = 1, snapshot2 = 2
    real(sp) :: kabsch_rmsd_val, theobald_rmsd_val
    character(len=MAX_FILENAME_LENGTH) :: filename
      ! Filename to write to.
    character(len=59), dimension(1) :: remarks
      ! optional remarks to put in header

    rmsd_mode = allatom_rmsd

    ! Compute RMSD between two snapshots, writing to PDB for independent verification.
    kabsch_rmsd_val = kabsch_ls_rmsd(natoms, snapshots(snapshot1,:,:), snapshots(snapshot2,:,:))
    theobald_rmsd_val = theobald_ls_rmsd(natoms, snapshots(snapshot1,:,:), snapshots(snapshot2,:,:))
   
    write(*,'(A20,F16.8)') 'kabsch_rmsd = ', kabsch_rmsd_val
    write(*,'(A20,F16.8)') 'theobald_rmsd = ', theobald_rmsd_val

    call ls_align(natoms, snapshots(snapshot1,:,:), snapshots(snapshot2,:,:))
    
    ! Write PDB files.
    call pdb_write('structure1.pdb', pdbatoms, snapshots(snapshot1:snapshot1,:,:))    
    call pdb_write('structure2.pdb', pdbatoms, snapshots(snapshot2:snapshot2,:,:))    

!    ! Compute RMSD between all consecutive pairs of snapshots.
!    avgrmsd = 0.0
!    count = 0.0
!    do i = 1,nsnapshots
!       if(mod(i,10000) .eq. 0) write(*,*) i
!
!       ! Compute the rmsd between successive snapshots.
!       rmsd = ls_rmsd(natoms, snapshots(i,:,:), snapshots(i+1,:,:))
!
!       ! Accumulate average.
!       avgrmsd = avgrmsd + rmsd
!       count = count + 1.0
!    end do
!
!    avgrmsd = avgrmsd / count
!    write(*,*) 'avgrmsd = ', avgrmsd

  end subroutine testRmsd

end module tests

!=============================================================================================
! Main test driver program.
!=============================================================================================

program main
  use decompose ! for state space decomposition
  use utilities, only : setRandomSeed ! to set the random seed from the current time
  use tests

  implicit none

  ! Constants.
  character(len=*), parameter :: filename = 'input.xml' ! XML control filename
    ! TODO: Replace with a command-line argument.

  ! Local variables.
  
  ! Initialize random number generator.
  call setRandomSeed()

  ! Initialize.
  call initialize(filename)

  ! Test RMSD code.
  call testRMSD()
  
  ! Test split code.
!  call testSplit()

  ! Run algorithm.
!  call iterate()

  ! Clean up.
  call finalize()

  ! Terminate
  stop
end program main
