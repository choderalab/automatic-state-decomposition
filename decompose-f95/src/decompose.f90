!=============================================================================================
! Automatic state decomposition algorithm for the construction of Markov / master equation
! models.
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California.  All Rights Reserved.
!
! This program is free software; you can redistribute it and/or modify it under the terms of
! the GNU General Public License as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this program;
! if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! Boston, MA  02110-1301, USA.
!=============================================================================================
! TODO:
! - Add an option to XML control file to NOT perform least-squared alignment before computation of RMSD.
!   This could be useful for certain systems.
! - Compute state normalized fluctuation autocorrelation functions within decompose after each iteration, 
!   and report on estimate of independent transitions out of each microstate and macrostate.
!   This correlation time is also an upper bound on the internal equilibration time of each state.
! - Use estimate of independent transitions out of each state to lump microstates that have been split too small.
! - Compute lower bound estimate on internal equilibration time from transition matrix constructed from restriction to each macrostate.
! - Also write eigenvalues mu, for reference, and theoretical and achieved metastabilities as a function of iteration.
! - Write current and optimal metastability of splitting and lumping stages as a function of iteration, for convenience.
! - For MCSA lumping, can we speed things up by only computing diagonal elements of lumped transition matrix?
! - Replace eigensolver for eigs with efficient method like ARPACK using sparse matrices and specialized sparse matrix-vector multiplication routines?  Will need to change transition matrix computation routines to generate sparse matrices instead of full ones.
! - Split off K-medoid routine for easy reuse.
! - Change from using 'rmsd_mode' to allocating/defining 'rmsd_atomindices(:)' to use in compute_distance at the beginning of each iteration?
!   The current method appears to involve creating many temporary arrays, which may incur a performance hit.
! - Eliminate temporary copies in distance calculation for efficiency.
! - Can we pass snapshots directly, rather than indices of snapshots?
! - Improve efficiency of openmp parallelization of updateGenerators code.
! - Generalize splitting and lumping subroutines to allow the user to focus on particular states to split more finely?
!=============================================================================================
! PUBLIC SUBROUTINES:
! * initialize
!   Parses control file, initializes the various data structures, loads trajectories.
!
! * iterate
!   Performs cycles of splitting and lumping to refine the initial state decomposition.
! 
! * finalize
!   Frees allocated data.
!=============================================================================================

!=============================================================================================
! A module to store simulation data and perform decomposition operations on state space.
!=============================================================================================
module decompose
  use numeric_kinds ! for precision
  use pdbio, only : pdbatom_t ! PDB atom type, for storage of PDB file metadata
  use transitionmatrix, only : trajectory_t ! Trajectory type
  use constants ! global constants

  implicit none

  ! TODO: Sort out permissions.
  !private
  public :: initialize, iterate, finalize, readStateAssignments, writeTimescales
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================
  
  logical, parameter :: debug = .TRUE. ! Flag to enable printing of debugging information.
  integer, parameter :: MAX_PCDATA_LENGTH = 1024

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

  ! Molecular topology data.
  integer :: natoms ! total number of atoms per snapshot
  integer, dimension(:), pointer :: ca_atomindices ! alpha carbon atom indices
  integer, dimension(:), pointer :: heavy_atomindices ! heavy atom and polar hydrogen atom indices (excluding symmetry-related atoms)
  integer, dimension(:), pointer :: bb_atomindices ! backbone atom atomindices 
  integer, dimension(:), pointer :: user_atomindices ! user-specified atom indices
  integer, dimension(:), pointer :: ls_align_atomindices ! atom indices to use when aligning structures for writing PDBs

  ! Reference PDB data for reading and writing of PDB files.
  type(pdbatom_t), dimension(:), allocatable :: pdbatoms ! atom information from PDB file so that PDB files can be written
  real(sp), dimension(:,:), allocatable :: reference_snapshot ! reference snapshot for aligning all snapshots upon writing

  ! Snapshot (conformation) storage data.
  integer :: nsnapshots ! total number of snapshots
  real(sp), dimension(:,:,:), allocatable :: snapshots ! snapshots(i,n,k) is dimension k of atom n of snapshot i

  ! Trajectory (temporal connectivity of snapshots) data.
  ! NOTE: Trajectories are stored as contiguous sets of conformations in 'snapshots'.
  ! See transitionmatrix.f90 for a definition of type(trajectory_t).
  integer :: ntrajectories ! number of trajectories  
  type(trajectory_t), dimension(:), allocatable :: trajectories ! Trajectory data.

  ! Distance metric data.
  ! integer, dimension(:), allocatable :: rmsd_atomindices ! atom indices to use for computing RMSD
  integer, parameter :: allatom_rmsd = 1, ca_rmsd = 2, backbone_rmsd = 3, heavyatom_rmsd = 4, user_rmsd = 5
  integer :: rmsd_mode = allatom_rmsd ! choice of above mode for rmsd

  ! Run parameters.
  character(len = MAX_PCDATA_LENGTH) :: title ! the title of the system
  integer :: tau ! observation interval for computing transition matrix for lumping
  real(sp) :: tau_unit ! unit timestep (in ps or ns)
  integer :: niterations ! number of iterations of splitting and lumping to carry out   
  integer :: iteration ! current iteration of splitting and lumping

  integer :: resume_from_iteration ! specified iteration to resume from, or -1 if none specified

  integer :: target_microstate_size ! the target size (number of snapshots) that is used in determining how many states to split to
  integer :: max_microstates_per_macrostate_first_iteration ! the maximum number of microstates to split each macrostate into, or no limit if negative
  integer :: max_microstates_per_macrostate ! the maximum number of microstates to split each macrostate into, or no limit if negative
  integer :: rmsd_mode_first_iteration ! RMSD mode to use in first iteration
  integer :: rmsd_mode_subsequent_iterations ! RMSD mode to use in subsequent iterations
  logical :: full_kmedoid_update ! if .true., the exact data point that minimizes the intracluster variance is used -- otherwise, stochastic update with ntries = ngenerators is used.

  logical :: use_timereversed ! if .true., will use time-reversed trajectories
  integer :: target_nmacrostates ! target number of macrostates
  integer :: mcsa_ntrials ! number of MC/SA lumping optimization runs per lump step
  integer :: mcsa_nsteps ! number of steps per MC/SA lumping optimization run
  character(len=MAX_FILENAME_LENGTH) :: output_directory ! directory to write output files to

  integer :: max_tau ! maximum lag time to write timescales for
  integer :: tau_step ! lag time step to use
  integer :: max_ntimescales ! maximum number of timescales of interest (may be less)
  real(sp) :: timescales_confidence_interval ! confidence interval to use in computing confidence bounds on timescales
  integer :: timescales_bootstrap_ntrials ! number of bootstrap trials to use in computing timescales

  logical :: use_ev_initial_lumping  ! flag to set whether eigenvector scheme is to be used for initial lumping
  real(sp) :: generator_refinement_fraction ! fraction of snapshots from 0 to 1 to use for generator refinement
  integer :: generator_refinement_minsize ! minimum number of snapshots per microstate to use in generator refinement
  logical :: write_split_timescales ! .true. if we should write timescales each iteration for splitting
  logical :: write_merged_timescales ! .true. if we should write timescales each iteration for lumping
  logical :: write_merged_eigenvectors ! .true. if we should write merged eigenvectors after lumping

  logical :: write_split_pdbs ! .true. if we should write state PDB files each iteration after splitting
  integer :: split_pdb_nmodels ! number of models to write
  logical :: write_merged_pdbs ! .true. if we should write state PDB files each iteration after merging
  integer :: merged_pdb_nmodels ! number of models to write

  integer :: state_pdb_nmodels ! the number of models to be included in each state PDB file
  logical :: write_split_table ! .true. if a translation of lumped states -> split states is to be written
  logical :: write_lump_table ! .true. if a translation of split states -> lumped states is to be written
  real(sp) :: kmedoid_fraction ! fraction of data to use in k-medoid iterations
  integer :: kmedoid_iterations ! number of iterations of k-medoids to conduct (can be zero)

  ! State assignment data.
  integer :: nmicrostates ! number of microstates
  integer, dimension(:), allocatable :: microstate_assignments ! microstate assignments for each snapshot, each having the snapshot index of the generator
  integer, dimension(:), allocatable :: compacted_microstate_assignments ! microstate assignments for each snapshot, from 1...nmicrostates
  integer, dimension(:), pointer :: microstates ! unique ids of microstates

  integer :: nmacrostates ! number of macrostates
  integer, dimension(:), allocatable :: macrostate_assignments ! macrostate assignments for each snapshot, each having the snapshot index of a member of the state
  integer, dimension(:), allocatable :: compacted_macrostate_assignments ! macrostate assignments for each snapshot, from 1...nmacrostates
  integer, dimension(:), pointer :: macrostates ! unique ids of macrostates
  
contains

  !=============================================================================================
  ! Read reference PDB file, determining number of atoms, saving atom metadata, and storing reference snapshot.
  !
  ! After this call, the following module data is set: 
  !   natoms, pdbatoms, reference_snapshot
  !=============================================================================================
  subroutine readReferencePDB(reference_pdb_filename)

    use pdbio ! for PDB I/O

    ! Parameters.
    character(len=*), intent(in) :: reference_pdb_filename
      ! The filename of the reference PDB file.

    ! Local variables.
    integer :: nmodels 
      ! number of models in the PDB file
    real(sp), dimension(:,:,:), allocatable :: pdbmodels
      ! models read from the PDB file
    
    ! Query for number of atoms and models.
    call pdb_read(reference_pdb_filename, natoms=natoms, nmodels=nmodels)
    if(debug) print *, 'The reference PDB files contains ', natoms, ' atoms, and ', nmodels, ' models.'

    ! Allocate storage for atom records (metadata) and model coordinates.
    allocate(pdbatoms(natoms),pdbmodels(nmodels,natoms,3))

    ! Read atom records and model coordinates.
    call pdb_read(reference_pdb_filename, atoms=pdbatoms, models=pdbmodels)

    ! Store the first model as the reference snapshot.
    ! TODO: Average the models instead?
    allocate(reference_snapshot(natoms,3))
    reference_snapshot = pdbmodels(1,:,:)  

    ! Free model coordinates.
    deallocate(pdbmodels)

  end subroutine readReferencePDB

  !=============================================================================================
  ! Print out which atoms have been selected by the given mask.
  !=============================================================================================
  subroutine showMask(mask)
    
    ! Parameters.
    logical, dimension(natoms), intent(in) :: mask

    ! Local variables.
    integer :: n

    do n = 1,natoms
       if(mask(n)) then
          write(*,'(I5,1X,A3,1X,I5,1X,A4,1X,A1)') &
               pdbatoms(n)%serial, pdbatoms(n)%resName, pdbatoms(n)%resSeq, pdbatoms(n)%name, '*'
       else
          write(*,'(I5,1X,A3,1X,I5,1X,A4,1X,A1)') &
               pdbatoms(n)%serial, pdbatoms(n)%resName, pdbatoms(n)%resSeq, pdbatoms(n)%name, ' '
       end if
    end do

  end subroutine showMask

  !=============================================================================================
  ! Set up atom indices lists used for computing RMSD, LS-RMSD alignments, etc.
  !
  ! NOTE: This subroutine is fairly ugly because these sorts of selection operations are inelegant
  ! in Fortran 90.  Is there a better way to do this?
  !=============================================================================================
  subroutine setupAtomIndices()
    use utilities, only : find

    ! Local variables.
    logical, dimension(natoms) :: mask
    integer :: nindices
    integer :: n    

    ! alpha carbons (CA) or terminal methyl groups (CH3)
    mask = ((pdbatoms%name == " CA ") .or. (pdbatoms%name == " CH3"))
    ! Exclude ions.
    mask = mask .and. .not. (pdbatoms%name == "Cl- " .or. pdbatoms%name == "Na+ ")    
    ! Select.
    ca_atomindices => find(mask)
    
    if(debug) then
       print *, 'Selecting alpha carbons:'
       call showMask(mask)
       print *, 'atom indices: ', ca_atomindices
    end if

    ! backbone atoms (including H and O)
    mask = ((pdbatoms%name == " CA ") .or. (pdbatoms%name == " CH3") .or. (pdbatoms%name == " N  ") &
         .or. (pdbatoms%name == " C  ") .or. (pdbatoms%name == " O  ") .or. (pdbatoms %name == " H  "))
    ! Exclude ions.
    mask = mask .and. .not. (pdbatoms%name == "Cl- " .or. pdbatoms%name == "Na+ ")
    ! Select.
    bb_atomindices => find(mask)
    
    if(debug) then
       print *, 'Selecting backbone atoms:'
       call showMask(mask)
       print *, 'atom indices: ', bb_atomindices
    end if

    ! heavy atoms plus polar hydrogens
    ! select all nonhydrogens
    mask = ((pdbatoms%name(2:2) /= 'H'))
    ! Exclude ions.
    mask = mask .and. .not. (pdbatoms%name == "Cl- " .or. pdbatoms%name == "Na+ ")
    ! remove symmetry-related heavy atoms
    mask = mask .and. .not. &
         ( (pdbatoms%resName == "VAL" .and. pdbatoms%name == " CG ") &
         .or. (pdbatoms%resName == "LEU" .and. pdbatoms%name == " CD ") &
         .or. (pdbatoms%resName == "PHE" .and. (pdbatoms%name(1:3) == " CD" .or. pdbatoms%name(1:3) == " CE")) &
         .or. (pdbatoms%resName == "TYR" .and. (pdbatoms%name(1:3) == " CD " .or. pdbatoms%name(1:3) == " CE")) &
         .or. (pdbatoms%resName == "GLU" .and. (pdbatoms%name == " OD1" .or. pdbatoms%name == " OD2")) &
         .or. (pdbatoms%resName == "ASP" .and. (pdbatoms%name == " OG1" .or. pdbatoms%name == " OG2")) &
         .or. (pdbatoms%resName == "HIP" .and. (pdbatoms%name == " ND1" .or. pdbatoms%name == " NE2")) &
         .or. (pdbatoms%resName == "ARG" .and. (pdbatoms%name == " NH1" .or. pdbatoms%name == " NH2")) &
         .or. (pdbatoms%resSeq == pdbatoms(natoms)%resSeq .and. (pdbatoms%name == " O  " .or. pdbatoms%name == " OXT")) )
    ! add in polar hydrogens
    mask = mask .or. &
         ( (pdbatoms%resName == "HID" .and. pdbatoms%name == " HD1") &
         .or. (pdbatoms%resName == "HIE" .and. pdbatoms%name == " HE1") &
         .or. (pdbatoms%resName == "HIE" .and. pdbatoms%name == " HE1") &
         .or. (pdbatoms%resName == "SER" .and. pdbatoms%name == " HG ") &
         .or. (pdbatoms%resName == "THR" .and. pdbatoms%name == " HG ") &
         .or. (pdbatoms%resName == "TYR" .and. pdbatoms%name == " HH ") &
         .or. (pdbatoms%resName == "TRP" .and. pdbatoms%name == " HE1") )
    ! Select.
    heavy_atomindices => find(mask)

    if(debug) then
       print *, 'Selecting heavy atoms and polar hydrogens, excluding symmetry-related atoms:'
       call showMask(mask)
       print *, 'atom indices: ', heavy_atomindices
    end if

  end subroutine setupAtomIndices

  !=============================================================================================
  ! Read trajectories indexed in a trajectory-index file.
  !
  ! Trajectories can either be stored in AMBER netCDF trajectories or as multiple models in
  ! a PDB file.  The parameter trajectory_format should be set to either 'netCDF' or 'PDB'.
  !
  ! This subroutine sets the following module data:
  !   ntrajectories, trajectories, nsnapshots, snapshots
  !=============================================================================================
  subroutine readTrajectories(trajectory_index_filename, trajectory_format)

    use pdbio ! for automatic unit number determination and PDB file reading
    use netcdfio ! for netCDF trajectory reading
    use timer ! for timing info
    
    ! Parameters.
    character(len=*), intent(in) :: trajectory_index_filename
      ! name of trajectory index file
    character(len=*), intent(in) :: trajectory_format
      ! trajectory format -- one of 'netCDF' or 'PDB'

    ! Local variables.
    integer :: iunit
      ! unit number for reading file
    character(len=MAX_FILENAME_LENGTH) :: line
      ! an entire line from the PDB file
    character(len=MAX_FILENAME_LENGTH) :: trajectory_filename
      ! Filename of current trajectory
    integer :: trajectory_index
      ! Current trajectory index
    real(dp) :: log_weight
      ! Log weight of current trajectory
    integer :: trajectory_length
      ! Current length (in snapshots) of trajectory
    integer :: initial_snapshot, final_snapshot
      ! Indices of first and last snapshot from the trajectory in 'snapshots'
    type(netcdfio_t) :: nc
      ! netcdf file handle

    ! Open trajectory index file.
    if(debug) print *, 'opening trajectory index file ', trim(trajectory_index_filename)
    call resetTimer()
    call new_unit(iunit)
    open(unit=iunit, file=trajectory_index_filename, status='OLD', err=29, position='REWIND', &
         form='FORMATTED', action='READ')
    
    ! Read trajectory file once to determine number of trajectories.
    ntrajectories = 0
    do
       ! Read a line.
       read(unit=iunit, fmt='(a)', end=20, err=29) line

       ! Extract weight and trajectory filename.
       read(line, fmt = '(F16.8,1X,A)') log_weight, trajectory_filename

       ! Increment trajectory counter.
       ntrajectories = ntrajectories + 1
    end do
       
20  continue
    if(debug) print *, 'ntrajectories = ', ntrajectories

    ! Rewind file.
    rewind(unit=iunit)

    ! Allocate trajectory storage.
    allocate(trajectories(ntrajectories))

    ! Read and store trajectory file information, and count number of snapshots.
    trajectory_index = 1
    nsnapshots = 0
    do
       ! Read a line.
       read(unit=iunit, fmt='(a)', end=21, err=29) line

       ! Extract weight and trajectory filename.
       read(line, fmt = '(F16.8,1X,A)') log_weight, trajectory_filename

       ! Query trajectory file for length of trajectory.
       if(trim(trajectory_format) == 'netCDF') then
          nc = netcdf_open(trim(trajectory_filename), nframes=trajectory_length)
          call netcdf_close(nc)
       elseif(trim(trajectory_format) == 'PDB') then
          ! TODO: As a check, make sure number of atoms matches file.
          call pdb_read(trim(trajectory_filename), nmodels=trajectory_length)          
       else
          write(*,*) 'Trajectory format "', trim(trajectory_format), '" not supported.'
          stop
       end if

       ! Store trajectory start and stop indices.
       trajectories(trajectory_index)%filename = trajectory_filename
       trajectories(trajectory_index)%log_weight = log_weight
       trajectories(trajectory_index)%length = trajectory_length
       if (trajectory_index == 1) then
          trajectories(trajectory_index)%start = 1
       else
          trajectories(trajectory_index)%start = trajectories(trajectory_index-1)%start + trajectories(trajectory_index-1)%length
       end if

       ! Accumulate total number of snapshots.
       nsnapshots = nsnapshots + trajectory_length

       ! Increment trajectory counter.
       trajectory_index = trajectory_index + 1
    end do
       
21  continue
    ! Close file.
    close(iunit)
    if(debug) write(*,*) readTimer(), ' seconds elapsed.'

    ! Allocate storage for snapshots.
    if(debug) write(*,*) 'Allocating ', (nsnapshots*natoms*3*4/(1024.0*1024.0)), ' MB for storage for ', nsnapshots, ' snapshots...'
    call resetTimer()
    allocate(snapshots(nsnapshots,natoms,3))    

    ! Read and store snapshots.
    ! NOTE: This section could be parallelized if netCDF was thread-safe for opening and reading separate files.
    ! NOTE: We could potentially try opening all the trajectory files first, storing their 'nc' entries in an array, and then read them in parallel...
    if(debug) write(*,*) 'Reading ', nsnapshots,' snapshots from ', ntrajectories, ' trajectories...'
    do trajectory_index = 1,ntrajectories
       ! Extract relevant information from trajectory data structure.
       trajectory_filename = trajectories(trajectory_index)%filename
       initial_snapshot = trajectories(trajectory_index)%start
       final_snapshot = initial_snapshot + trajectories(trajectory_index)%length - 1

       ! Read snapshots from trajectory into their appropriate position in 'snapshots'.
       if(trim(trajectory_format) == 'netCDF') then
          nc = netcdf_open(trim(trajectory_filename))
          call netcdf_read(nc, snapshots(initial_snapshot:final_snapshot,:,:))
          call netcdf_close(nc)
       elseif(trim(trajectory_format) == 'PDB') then
          call pdb_read(trim(trajectory_filename), models = snapshots(initial_snapshot:final_snapshot,:,:))
       else
          write(*,*) 'Trajectory format "', trim(trajectory_format), '" not supported.'
          stop
       end if

    end do
    if(debug) write(*,*) readTimer(), ' seconds elapsed.'

    return

    ! ERROR in read.  Execution should only reach here upon error.
29  continue
    write(*,*) 'decompose: readTrajectories: Error opening ', trim(trajectory_index_filename)
    stop

  end subroutine readTrajectories

  !=============================================================================================
  ! Write trajectories in "old style" trajectory file format.
  !=============================================================================================
  subroutine writeTrajectoriesOldFormat(filename)

    use pdbio, only : new_unit

    ! Parameters.
    character(len=*), intent(in) :: filename

    ! Local variables.
    integer :: iunit
      ! unit number for reading file
    character(len=MAX_LINE_LENGTH) :: line, text
      ! an entire line from the PDB file
    integer :: trajectory_index
      ! Current trajectory index
    integer :: trajectory_length
      ! Current length (in snapshots) of trajectory
    integer :: initial_snapshot, final_snapshot
      ! Indices of first and last snapshot from the trajectory in 'snapshots'
    integer :: i
    integer :: snapshot_index

    ! Open trajectory index file.
    if(debug) print *, 'opening trajectories output ', trim(filename)
    call new_unit(iunit)
    open(unit=iunit, file=filename, status='REPLACE', err=39, position='REWIND', &
         form='FORMATTED', action='WRITE')

    ! Loop over trajectories.
    snapshot_index = 1
    do trajectory_index = 1,ntrajectories
       ! Construct line.
       write(text,'(I12)') trajectory_index
       write(iunit,'(A)',advance='no') trim(adjustl(text))


       write(iunit,'(A)',advance='no') ' 1'
              
       do i = 1,trajectories(trajectory_index)%length
          write(text,'(I12)') snapshot_index
          write(iunit,'(A)',advance='no') ' '//trim(adjustl(text))
          snapshot_index = snapshot_index + 1
       end do

       write(iunit,'(A)') ''
    end do

    ! Close file.
    close(iunit)

    return

    ! ERROR in read.  Execution should only reach here upon error.
39  continue
    write(*,*) 'decompose: writeTrajectoriesOldFormat: Error opening ', trim(filename)
    stop

  end subroutine writeTrajectoriesOldFormat

  !=============================================================================================
  ! If tag is present, set flag to .true., otherwise is .false.
  !=============================================================================================
  subroutine get_xml_flag(fxml, key, flag)
    use flib_xpath ! for XML files
    
    ! Parameters.
    type(xml_t), intent(inout) :: fxml
      ! XML file.    
    character(len=*), intent(in) :: key
      ! Key
    logical, intent(out) :: flag

    ! Local variables.
    integer  :: status
    character(len=MAX_PCDATA_LENGTH)  :: pcdata

    ! Assign default.
    flag = .false.

    ! Rewind.
    call rewind_xmlfile(fxml)
    
    ! Parse.
    call get_node(fxml, path="//"//trim(key), pcdata=pcdata, status=status)
    if(status >= 0) flag = .true.

    if(debug) print *, key, ' = ', flag
  
  end subroutine get_xml_flag

  !=============================================================================================
  ! Extract an integer value from XML file.
  !=============================================================================================
  subroutine get_xml_integer(fxml, key, value, default)
    use flib_xpath ! for XML files
    
    ! Parameters.
    type(xml_t), intent(inout) :: fxml
      ! XML file.    
    character(len=*), intent(in) :: key
      ! Key
    integer, intent(out) :: value
    integer, intent(in) :: default

    ! Local variables.
    integer  :: status
    character(len=MAX_PCDATA_LENGTH)  :: pcdata

    ! Assign default.
    value = default

    ! Rewind.
    call rewind_xmlfile(fxml)
    
    ! Parse.
    call get_node(fxml, path="//"//trim(key), pcdata=pcdata, status=status)
    if(status >= 0) read(pcdata,'(I16)') value

    if(debug) print *, key, ' = ', value
  
  end subroutine get_xml_integer

  !=============================================================================================
  ! Extract an integer value from XML file.
  !=============================================================================================
  function get_xml_integers(fxml, key, default, nmax) result(value)
    use flib_xpath ! for XML files
    
    ! Parameters.
    type(xml_t), intent(inout) :: fxml
      ! XML file.    
    character(len=*), intent(in) :: key
      ! Key
    integer, dimension(:), intent(in) :: default
    integer, intent(in) :: nmax
      ! maximum possible size for integer array    

    ! Return variable.
    integer, dimension(:), pointer :: value

    ! Local variables.
    integer  :: status
    character(len=MAX_PCDATA_LENGTH)  :: pcdata
    integer, dimension(nmax) :: temparray
    integer :: nelements

    ! Rewind.
    call rewind_xmlfile(fxml)
    
    ! Parse.
    call get_node(fxml, path="//"//trim(key), pcdata=pcdata, status=status)
    if(status >= 0) then
       nelements = 0
       call build_data_array(pcdata, temparray, nelements)
       allocate(value(nelements))
       value = temparray       
    else
       allocate(value(size(default,1)))
       value = default
    end if

    if(debug) print *, key, ' = ', value
  
  end function get_xml_integers

  !=============================================================================================
  ! Extract an integer value from XML file.
  !=============================================================================================
  subroutine get_xml_realsp(fxml, key, value, default)
    use flib_xpath ! for XML files
    
    ! Parameters.
    type(xml_t), intent(inout) :: fxml
      ! XML file.    
    character(len=*), intent(in) :: key
      ! Key
    real(sp), intent(out) :: value
    real(sp), intent(in) :: default

    ! Local variables.
    integer  :: status
    character(len=MAX_PCDATA_LENGTH)  :: pcdata

    ! Assign default.
    value = default

    ! Rewind.
    call rewind_xmlfile(fxml)
    
    ! Parse.
    call get_node(fxml, path="//"//trim(key), pcdata=pcdata, status=status)
    if(status >= 0) read(pcdata,'(F16.8)') value

    if(debug) print *, key, ' = ', value

  end subroutine get_xml_realsp

  !=============================================================================================
  ! Extract a string from an XML file.
  !=============================================================================================
  subroutine get_xml_string(fxml, key, value, default)
    use flib_xpath ! for XML files
    
    ! Parameters.
    type(xml_t), intent(inout) :: fxml
      ! XML file.    
    character(len=*), intent(in) :: key
      ! Key
    character(len=*), intent(out) :: value
    character(len=*), intent(in) :: default

    ! Local variables.
    integer :: status
    character(len=MAX_PCDATA_LENGTH) :: pcdata

    ! Assign default.
    value = default
    
    ! Rewind.
    call rewind_xmlfile(fxml)

    ! Parse.
    call get_node(fxml, path='//'//trim(key), pcdata=pcdata, status=status)
    if(status == 0) value = pcdata

    if(debug) print *, key, ' = ', trim(value)
  
  end subroutine get_xml_string

  function determine_rmsd_mode(string) result(rmsd_mode)

    ! Parameters.
    character(len=*), intent(in) :: string

    ! Return value.
    integer :: rmsd_mode
        
    if(string(1:2) .eq. 'ca') then
       rmsd_mode = ca_rmsd
    elseif(string(1:8) .eq. 'backbone') then
       rmsd_mode = backbone_rmsd
    elseif(string(1:) .eq. 'allatom') then
       rmsd_mode = allatom_rmsd
    elseif(string(1:9) .eq. 'heavyatom') then
       rmsd_mode = heavyatom_rmsd
    elseif(string(1:4) .eq. 'user') then
       rmsd_mode = user_rmsd
    else
       write(*,*) 'Unrecognized RMSD option: ', string
       stop
    end if

  end function determine_rmsd_mode

  !=============================================================================================
  ! Check to make sure no states are assigned to be zero.
  !=============================================================================================
  subroutine checkstates(stateassignments)
    ! Paramters.
    integer, dimension(:), intent(in) :: stateassignments

    ! Local variables.
    integer :: i
    
    if(any(stateassignments==0)) then
       write(*,*) 'some states assigned state 0!'
       do i = 1,size(stateassignments,1)
          write(19,'(I8,X,I8)') i, stateassignments(i)
       end do
       stop
    end if
    
  end subroutine checkstates

  !=============================================================================================
  ! Initialize the algorithm.
  !
  ! Read reference PDB data, setup atom indices lists, read snapshot data.
  !=============================================================================================
  subroutine initialize(control_filename, nocheck)

    use netcdfio ! for netCDF snapshot I/O
    use utilities, only : unique
    use flib_xpath ! for XML files
    
    ! Parameters.
    character(len=*), intent(in) :: control_filename   ! name of XML control file
    logical, intent(in), optional :: nocheck ! if .TRUE., no checking of input states is done

    ! Local variables.
    character(len=MAX_FILENAME_LENGTH) :: filename
      ! Temporary filename storage.
    type(xml_t) :: fxml
      ! XML file.
    type(dictionary_t) :: attributes
      ! storage for attributes associated with a node
    integer  :: status
      ! integer indicating success or failure of xmlf90 call
    character(len=MAX_PCDATA_LENGTH)  :: pcdata
      ! storage for node data
    character(len=MAX_PCDATA_LENGTH)  :: trajectory_format
      ! trajectory format type
    logical :: no_timereversed
      ! If .true., time-reversed trajectories will not be included.

    ! Open XML control file.
    if(debug) write(*,*) 'Opening XML control file ', trim(control_filename)
    call open_xmlfile(control_filename, fxml, status)

    ! DEBUG.
    ! call enable_debug(sax=.false.)

    ! Get title.
    call get_xml_string(fxml, 'title', title, '')

    ! Set target microstate sizes.
    call get_xml_integer(fxml, 'target_microstate_size', target_microstate_size, 100)
    call get_xml_integer(fxml, 'max_microstates_per_macrostate_first_iteration', max_microstates_per_macrostate_first_iteration, -1)
    call get_xml_integer(fxml, 'max_microstates_per_macrostate', max_microstates_per_macrostate, -1)

    ! Set splitting parameters.
    call get_xml_realsp(fxml, 'kmedoid_fraction', kmedoid_fraction, 1.0)
    if(kmedoid_fraction <= 0 .or. kmedoid_fraction > 1) then
       write(*,*) 'kmedoid_fraction must be in interval (0,1].'
       stop
    end if
    call get_xml_integer(fxml, 'kmedoid_iterations', kmedoid_iterations, 5)
    call get_xml_string(fxml, 'rmsd_mode_first_iteration', pcdata, 'ca')
    rmsd_mode_first_iteration = determine_rmsd_mode(pcdata)
    call get_xml_string(fxml, 'rmsd_mode_subsequent_iterations', pcdata, 'heavyatom')
    rmsd_mode_subsequent_iterations = determine_rmsd_mode(pcdata)
    call get_xml_flag(fxml, 'full_kmedoid_update', full_kmedoid_update)

    ! Set lumping parameters.
    call get_xml_integer(fxml, 'tau', tau, 1)
    call get_xml_realsp(fxml, 'tau_unit', tau_unit, 1.0)

    call get_xml_flag(fxml, 'use_ev_initial_lumping', use_ev_initial_lumping)

    call get_xml_integer(fxml, 'mcsa_ntrials', mcsa_ntrials, 20)
    call get_xml_integer(fxml, 'mcsa_nsteps', mcsa_nsteps, 20000)
    call get_xml_integer(fxml, 'target_nmacrostates', target_nmacrostates, 10)

    ! Set output directory.
    call get_xml_string(fxml, 'output_directory', output_directory, 'output')

    ! Set writing of splitting and lumping tables.
    write_split_table = .true.
    write_lump_table = .true.

    ! Read reference PDB file, storing atom metadata (for atom indices generation, PDB writing),
    ! and storing reference snapshot (for aligning states written to PDB files).
    call get_xml_string(fxml, 'reference_pdb_filename', filename, 'reference.pdb')
    if(debug) write(*,*) 'Reading reference PDB file ', trim(filename)
    call readReferencePDB(filename)

    ! Set up atom indices.
    call setupAtomIndices()

    ! Set up user-specified RMSD atom indices.
    user_atomindices => get_xml_integers(fxml, 'user_atomindices', heavy_atomindices, natoms)

    ! Set up atom indices for alignment.
    ls_align_atomindices => get_xml_integers(fxml, 'ls_align_atomindices', bb_atomindices, natoms)

    ! Read MD trajectories, setting up both trajectory and snapshot information.    
    call get_xml_string(fxml, 'trajectory_index_filename', filename, 'trajectory-index')
    call get_xml_string(fxml, 'trajectory_format', trajectory_format, 'netCDF')
    call readTrajectories(trim(filename), trim(trajectory_format))    

    ! Write which snapshots are in which trajectories if instructed to do so.
    call get_xml_string(fxml, 'write_trajectories', filename, '')
    if(len_trim(filename) > 0) call writeTrajectoriesOldFormat(filename)

    ! Set max_tau and tau_step for writing timescales.
    call get_xml_flag(fxml, 'write_split_timescales', write_split_timescales)
    call get_xml_flag(fxml, 'write_merged_timescales', write_merged_timescales)
    call get_xml_flag(fxml, 'write_merged_eigenvectors', write_merged_eigenvectors)
    call get_xml_integer(fxml, 'max_ntimescales', max_ntimescales, target_nmacrostates-1)
    call get_xml_integer(fxml, 'tau_step', tau_step, 1)
    call get_xml_integer(fxml, 'max_tau', max_tau, ceiling((maxval(trajectories%length)-1) / 2.0))
    call get_xml_realsp(fxml, 'timescales_confidence_interval', timescales_confidence_interval, 0.68)
    call get_xml_integer(fxml, 'timescales_bootstrap_ntrials', timescales_bootstrap_ntrials, 40)

    ! Determine whether time-reversibility is to be exploited.
    call get_xml_flag(fxml, 'no_timereversed', no_timereversed)
    use_timereversed = .not. no_timereversed    
    if(no_timereversed) then
       ! Disable all options that have to do with computing timescales, because current code assumes
       ! detailed balance is satisfied by using time-reversed trajectories.
       write(*,*) 'Trajectories are not time-reversible, so disabling computation of eigenvalues/vectors.'       
       write_split_timescales = .false.
       write_merged_timescales = .false.
       write_merged_eigenvectors = .false.
    end if

    ! Ensure that max_tau is shorter than the longest trajectory lengths.
    if(max_tau > maxval(trajectories%length)-1) then
       write(*,*) 'max_tau is longer than longest trajectory -- adjusting to a more sensible value.'
       max_tau = ceiling((maxval(trajectories%length)-1) / 2.0)
    end if

    ! Set state PDB ensemble options.
    write_split_pdbs = .false.
    split_pdb_nmodels = 30
    write_merged_pdbs = .false.
    merged_pdb_nmodels = 30
    call rewind_xmlfile(fxml)
    call get_node(fxml, path="//write_split_pdbs", pcdata=pcdata, attributes=attributes, status=status)
    if(status >= 0) then
       write_split_pdbs = .true.
       call get_value(attributes, 'nmodels', pcdata, status)
       if(status >= 0) read(pcdata,'(I8)') split_pdb_nmodels
       if(debug) write(*,*) 'Will write split PDBs with nmodels = ', split_pdb_nmodels
    end if
    call rewind_xmlfile(fxml)
    call get_node(fxml, path="//write_merged_pdbs", pcdata=pcdata, attributes=attributes, status=status)
    if(status >= 0) then
       write_merged_pdbs = .true.
       call get_value(attributes, 'nmodels', pcdata, status)
       if(status >= 0) read(pcdata,'(I8)') merged_pdb_nmodels
       if(debug) write(*,*) 'Will write merged PDBs with nmodels = ', merged_pdb_nmodels
    end if

    ! Allocate storage for microstate and macrostate assignments.
    allocate(microstate_assignments(nsnapshots),macrostate_assignments(nsnapshots))
    allocate(compacted_microstate_assignments(nsnapshots),compacted_macrostate_assignments(nsnapshots))
    
    ! Nullify microstates.
    nmicrostates = 0
    microstates => null()

    ! Initialize by assigning all snapshots to one macrostate.
    macrostate_assignments(:) = 1

    ! Start from iteration 1 (unless overridden).
    iteration = 1
    
    ! Get number of iterations to execute.
    call get_xml_integer(fxml, 'niterations', niterations, 10)

    ! Determine initialization of iteration and state assignments.
    call rewind_xmlfile(fxml)
    call get_node(fxml, path="//initialize", pcdata=pcdata, attributes=attributes, status=status)
    if(status >= 0) then
       ! Determine iteration to start from.
       call get_value(attributes, 'iteration', pcdata, status)
       if(status >= 0) read(pcdata,'(I8)') iteration

       ! Determine statassignments to start from.
       call get_value(attributes, 'initial_stateassignments', filename, status)
       if(status >= 0) then
          ! Attempt to load macrostate stateassignments.
          if(debug) write(*,*) 'Attempting to resume from end of iteration ', iteration
          call readStateassignments(filename, macrostate_assignments)          
       end if
    end if

    ! Determine macrostates and compacted macrostates.
    macrostates => unique(macrostate_assignments, compacted=compacted_macrostate_assignments)
    nmacrostates = size(macrostates,1)
    if(debug) write(*,*) 'There are ', nmacrostates, ' macrostates:'
    if(debug) write(*,*) macrostates
    
    ! DEBUG    
    if(.not. present(nocheck)) then
       call checkstates(macrostate_assignments)
       call checkstates(compacted_macrostate_assignments)
    end if

    ! Close XML control file.
    call close_xmlfile(fxml)

  end subroutine initialize

  !=============================================================================================
  ! Clean up allocated data.
  !=============================================================================================
  subroutine finalize

    ! Free snapshot storage.
    deallocate(snapshots)

    ! Free trajectory storage.
    deallocate(trajectories)

    ! Free atomindices.
    deallocate(ca_atomindices, bb_atomindices, heavy_atomindices)

    ! Free PDB atoms.
    deallocate(pdbatoms)

    ! Free state assignments.
    deallocate(microstate_assignments, macrostate_assignments)

  end subroutine finalize

  !=============================================================================================
  ! Compute the distance between two snapshots.
  !=============================================================================================
  function computeDistance(snapshot_index_1, snapshot_index_2)

    use theobald_rmsd, only : ls_rmsd ! for LS-RMSD using fast Newton iteration quaternion-based method
    !use kabsch_rmsd, only : ls_rmsd ! for LS-RMSD using Kabsch eigensolver method

    use pdbio
    
    ! Parameters.
    integer, intent(in) :: snapshot_index_1
      ! Index of first snapshot.
    integer, intent(in) :: snapshot_index_2
      ! Index of second snapshot.    

    ! Return value.
    real(sp) :: computeDistance
      ! The RMSD between the two snapshots.

    ! Compute snapshot RMSD.
    select case(rmsd_mode)
    case(heavyatom_rmsd)
       computeDistance = ls_rmsd(size(heavy_atomindices,1), snapshots(snapshot_index_1,heavy_atomindices,:), &
            snapshots(snapshot_index_2,heavy_atomindices,:))
    case(backbone_rmsd)
       computeDistance = ls_rmsd(size(bb_atomindices,1), snapshots(snapshot_index_1,bb_atomindices,:), &
            snapshots(snapshot_index_2,bb_atomindices,:))
    case(ca_rmsd)
       computeDistance = ls_rmsd(size(ca_atomindices,1), snapshots(snapshot_index_1,ca_atomindices,:), &
            snapshots(snapshot_index_2,ca_atomindices,:))
    case(user_rmsd)
       computeDistance = ls_rmsd(size(user_atomindices,1), snapshots(snapshot_index_1,user_atomindices,:), &
            snapshots(snapshot_index_2,user_atomindices,:))
    case default
       ! Compute all-atom RMSD by default.
       computeDistance = ls_rmsd(natoms, snapshots(snapshot_index_1,:,:), snapshots(snapshot_index_2,:,:))
    end select
    
  end function computeDistance

  !=============================================================================================
  ! Assign the specified snapshots to closest generators.
  !=============================================================================================
  subroutine assignSnapshotsToGenerators(nsnapshot_indices, snapshot_indices, ngenerators, generator_indices, stateassignments)

    ! Parameters.
    integer, intent(in) :: nsnapshot_indices
    integer, dimension(nsnapshot_indices), intent(in) :: snapshot_indices
    integer, intent(in) :: ngenerators
    integer, dimension(ngenerators), intent(in) :: generator_indices
    integer, dimension(nsnapshots), intent(out) :: stateassignments

    ! Local variables.
    integer :: i, j ! Loop variables.
    real(sp), dimension(ngenerators) :: generator_distances ! distances to all generators

    ! Assign each snapshot to the closest generator.
    if(debug) write(*,*) 'Assigning ', nsnapshot_indices, ' snapshots to generators...'
    !$omp parallel do default(shared) private(generator_distances)
    do i = 1, nsnapshot_indices
       if(debug .and. mod(i,1000).eq.0) write(*,*) 'snapshot ', i, ' / ', nsnapshots

       ! Compute distances to all generators.
       do j = 1, ngenerators
          generator_distances(j) = computeDistance(snapshot_indices(i), generator_indices(j))
       end do

       ! Assign state to the snapshot index of the closet generator.
       stateassignments(snapshot_indices(i)) = generator_indices(minloc(generator_distances,1))
    end do
    !$omp end parallel do
    if(debug) write(*,*) 'Done.'

    ! Sanity check.
    if(any(stateassignments(snapshot_indices) == 0)) then
       write(*,*) 'assignSnapshotsToGenerators: some stateassignments are zero!'
       stop
    end if

  end subroutine assignSnapshotsToGenerators

  !=============================================================================================
  ! Compute the intrastate variance of a set of snapshots from a given generator.
  !=============================================================================================
  function computeStateVariance(nindices, indices, generator_index) result(variance)
    
    ! Parameters.
    integer, intent(in) :: nindices
      ! number of snapshots in the state
    integer, dimension(:), intent(in) :: indices
      ! indices of snapshots belonging to the state
    integer, intent(in) :: generator_index 
      ! index of generator
      
    ! Return value.
    real(dp) :: variance

    ! Local variables.
    integer :: i
    integer :: snapshot_index
    
    ! Compute variance.
    variance = 0.0
    !$omp parallel do default(shared) reduction(+:variance)
    do i = 1,nindices
       variance = variance + computeDistance(generator_index, indices(i))**2
    end do
    !$omp end parallel do
    variance = variance / real(nindices-1,dp)

  end function computeStateVariance

  !=============================================================================================
  ! Use a stochastic procedure to attempt to find a better generator for a given state.
  ! A number of trials are attempted with randomly-chosen generators, attempting to reduce the intrastate variance.
  ! If ntries == nindices, then all data points in the state are considered to find the one that minimizes the intrastate variance.
  !=============================================================================================
  subroutine updateGenerator(nindices, indices, generator_index, ntries)

    use utilities, only : generateUniformInteger

    ! Parameters.
    integer, intent(in) :: nindices
      ! the number of indices in the state
    integer, dimension(nindices), intent(in) :: indices
      ! indices of the snapshots
    integer, intent(inout) :: generator_index
      ! the initial generator index, replaced by the new one if a snapshot with smaller variance is found
    integer, intent(in) :: ntries
      ! number of tries to improve generator

    ! Local variables.
    integer :: index
    integer :: i, j
    integer :: try
    integer :: trial_generator_index
    real(sp) :: variance, trial_variance

    if(ntries < nindices) then
       ! Use stochastic algorithm to try to update generator.

       ! Compute variance using initial generator_index.
       variance = computeStateVariance(nindices, indices, generator_index)

       do try = 1,ntries
          ! Select a snapshot at random from this state.
          trial_generator_index = indices(generateUniformInteger(nindices))
          
          ! Compute variance about this proposed generator.
          trial_variance = computeStateVariance(nindices, indices, trial_generator_index)
          
          ! Assign this to be the new generator if variance has decreased.
          if(trial_variance < variance) then
             variance = trial_variance
             generator_index = trial_generator_index
          end if
       end do

    else
       ! Try all members of state.
       do i = 1,nindices
          trial_generator_index = indices(i)

          ! Compute variance about this proposed generator.
          trial_variance = computeStateVariance(nindices, indices, trial_generator_index)          
          
          ! Assign this to be the new generator if variance has decreased.
          if(trial_variance < variance) then
             variance = trial_variance
             generator_index = trial_generator_index
          end if

       end do
    end if

  end subroutine updateGenerator

  !=============================================================================================
  ! Update generator positions using stochastic update scheme.
  !=============================================================================================
  subroutine updateGenerators(nsnapshot_indices, snapshot_indices, ngenerators, generator_indices, stateassignments)

    use utilities, only : find, generateUniformInteger

    ! Parameters.
    integer, intent(in) :: nsnapshot_indices
    integer, dimension(nsnapshot_indices), intent(in) :: snapshot_indices
    integer, intent(in) :: ngenerators
    integer, dimension(ngenerators), intent(out) :: generator_indices
    integer, dimension(nsnapshots), intent(in) :: stateassignments

    ! Local variables.
    integer :: i, j ! Loop variables.
    integer :: generator_index
    real(sp), dimension(ngenerators) :: generator_distances ! distances to all generators
    integer, dimension(:), pointer :: indices ! indices of snapshots belonging to current generator
    integer :: nindices ! number of indices

    do i = 1,ngenerators
       ! Find snapshots belonging to this microstate.
       indices => find(stateassignments(snapshot_indices) == generator_indices(i))
       indices = snapshot_indices(indices)
       nindices = size(indices,1)
       write(*,*) 'calling updateGenerator for generator ', i, ' (', nindices, ' members)'
       
       ! Try to find another member of this state that minimizes the intrastate variance about this generator.              
       if(full_kmedoid_update) then
          ! Use full update, considering all snapshots in state.
          call updateGenerator(nindices, indices, generator_indices(i), nindices)
       else
          ! Use stochastic update with ntries = ngenerators.
          call updateGenerator(nindices, indices, generator_indices(i), ngenerators)
       end if

       ! Clean up.
       deallocate(indices)
    end do

  end subroutine updateGenerators
  
  !=============================================================================================
  ! Split a set of snapshots into a specified number of subsets by spatial decomposition of the 
  ! region into a set of convex polytopes defined by proximity to particular 'generators'.
  ! In effect, a Vornoi decomposition of conformation space (in some space in which RMSD is the metric) is produced.
  !
  ! How the generators are arrived at is, in principle, arbitrary.
  ! This implementation uses a few iterations of an iterative k-medoid algorithm.
  !=============================================================================================
  subroutine splitState(nsnapshot_indices, snapshot_indices, ngenerators, stateassignments, generators)

    use utilities, only : generateUniqueRandomIndices, histogram, find

    ! Parameters.
    integer :: nsnapshot_indices
      ! Length of snapshot_indices
    integer, dimension(nsnapshot_indices), intent(in) :: snapshot_indices
      ! The list of indices of snapshots to split.
    integer, intent(in) :: ngenerators
      ! The number of microstates to split this state into.
    integer, dimension(nsnapshots), intent(out) :: stateassignments
      ! The state assignments (indices of generator snapshtots) for each snapshot in the list.
    integer, dimension(ngenerators), intent(out), optional :: generators
      ! Snapshot indices of the generators of the convex polytopes.
      ! NOTE that generators = unique(stateassignments), though not necessarily in that order.

    ! Local variables.
    integer :: i, j 
      ! Loop variables.
    integer, dimension(ngenerators) :: generator_indices 
      ! snapshot indices of generators
    integer :: iteration 
      ! iteration number of generator refinement
    integer :: nselected_snapshots 
      ! number of selected snapshots (for fast generator update)
    integer, dimension(:), allocatable :: selected_snapshot_indices 
      ! selected snapshots (for fast generator update)
    
    ! Sanity check on input parameters.  Ensure nsnapshots > ngenerators.
    if(ngenerators > nsnapshot_indices) then
       write(*,*) 'split: error: ngenerators > nsnapshot_indices'
       stop
    end if

    ! Choose how many snapshots to select for rapid determination of generators.
    nselected_snapshots = ceiling(nsnapshot_indices * kmedoid_fraction)
    if(nselected_snapshots < ngenerators) nselected_snapshots = nsnapshot_indices

    ! Determine snapshots to use in generator refinement.
    allocate(selected_snapshot_indices(nselected_snapshots))
    if(nselected_snapshots < nsnapshot_indices) then
       ! Select snapshots.
       if(debug) write(*,*) 'Choosing ', nselected_snapshots, ' of ', nsnapshot_indices, ' snapshots to refine generators...'
       call generateUniqueRandomIndices(nsnapshot_indices, nselected_snapshots, selected_snapshot_indices)    
       selected_snapshot_indices = snapshot_indices(selected_snapshot_indices)
    else
       ! Use all snapshots.  No need to select subset.
       selected_snapshot_indices = snapshot_indices
    end if

    ! Initiate algorithm by choosing a random set of ngenerators snapshots from the state, without replacement.
    if(debug) write(*,*) 'Choosing a random set of ', ngenerators, ' snapshots for generators...'
    call generateUniqueRandomIndices(nselected_snapshots, ngenerators, generator_indices)
    generator_indices = selected_snapshot_indices(generator_indices)  

    ! Iterate to refine generators.
    do iteration = 1,kmedoid_iterations
       ! Assign all snapshots to closest generator.
       write(*,*) 'Assigning snapshots to generators...'
       call assignSnapshotsToGenerators(nselected_snapshots, selected_snapshot_indices, &
            ngenerators, generator_indices, stateassignments)

       ! Update generators.
       write(*,*) 'Updating generators...'
       call updateGenerators(nselected_snapshots, selected_snapshot_indices, &
            ngenerators, generator_indices, stateassignments)
    end do

    ! Assign all snapshots to closest generator.
    call assignSnapshotsToGenerators(nsnapshot_indices, snapshot_indices, ngenerators, generator_indices, stateassignments)
    call checkstates(stateassignments(snapshot_indices))

    ! Provide generator indices if desired.
    if(present(generators)) generators = generator_indices

    ! Clean up.
    deallocate(selected_snapshot_indices)
    
  end subroutine splitState

  !=============================================================================================
  ! The transition matrix evaluation criteria to be maximized when trying to determine
  ! the optimal lumping of microstates into macrostates.
  !=============================================================================================
  function evaluationCriteria(Tji)
    
    use utilities, only : trace

    ! Parameters.
    real(dp), dimension(:,:) :: Tji
      ! Lumped transition matrix to be evaluated.
      ! Tji(j,i) is the transition probability of finding the system in state j a time tau later,
      ! given that it initially started in state i.

    ! Return value.
    real(sp) :: evaluationCriteria
      ! The criteria to be maximized over various possible lumpings of microstates into macrostates.

    ! Local variables.
    integer :: i

    ! In this implementation, we use the trace (which equals the sum of the eigenvalues) as the
    ! criteria to be maximized.  This is the sum of the self-transition probabilities, and tries
    ! to find the most metastable decomposition possible.    
    evaluationCriteria = trace(Tji)

    ! Sanity check.
    if(isnan(evaluationCriteria)) then
       write(*,*) 'evaluationCritera: trace is Nan'
       write(*,*) 'Tji = '
       do i = 1,nmacrostates
          write(*,'(6F8.5)') Tji(i,:)
       end do       
       stop
    end if

  end function evaluationCriteria

  !=============================================================================================
  ! Optimize a lumping of microstates into macrostates by using an MC/SA algorithm to maximize an objective.
  !=============================================================================================
  subroutine lump_mcsa(nmacrostates, best_criteria, best_micro_to_macro, N_ji)

    use utilities, only : generateUniformInteger, generateUniform01, histogram, generateUniqueRandomIndices
    use transitionmatrix, only : computeTransitionMatrix, computeLumpedTransitionMatrix

    ! Parameters.
    integer, intent(in) :: nmacrostates
      ! The number of macrostates to lump to.
    real(sp), intent(out) :: best_criteria
      ! The largest value of the objective function found.
    integer, dimension(nmicrostates), intent(inout) :: best_micro_to_macro
      ! The best microstate-to-macrostate lumping found.
    real(dp), dimension(nmicrostates,nmicrostates), intent(in) :: N_ji
      ! microstate transition count matrix

    ! Local variables.
    real(dp), dimension(nmacrostates,nmacrostates) :: T_ba
      ! Lumped transition matrix.
      ! T_ba(b,a) is the transition probability of finding the system in macrostate b a time tau later,
      ! given that it initially started in state a.
    integer, dimension(nmicrostates) :: micro_to_macro, micro_to_macro_proposed
      ! translation from microstate to macrostate index
    real(sp) :: criteria, criteria_proposed
      ! current and proposed objective function values of current and proposed lumpings
    integer, dimension(nmicrostates) :: indices
      ! indices storage
    integer :: i, j, n
      ! loop indices
    integer :: trial
      ! current MC/SA trial (max is ntrials)
    integer :: step
      ! current step in MC/SA trial (max is nsteps)
    integer, dimension(nmacrostates) :: macrostate_histogram
      ! histogram counts for macrostates, to make sure they don't end up empty
    integer :: microstate, macrostate
    logical :: accept
    real(sp) :: beta
    real(sp) :: r

    ! Take initial lumping from input best_micro_to_macro.
    micro_to_macro = best_micro_to_macro

    ! Determine how many microstates are lumped into each macrostate.
    macrostate_histogram = histogram(micro_to_macro, nmacrostates)
    
    ! Evaluate criteria function for initial lumping and store.
    ! T_ba = computeTransitionMatrix(ntrajectories, trajectories, nmacrostates, micro_to_macro(compacted_microstate_assignments), tau)
    call computeLumpedTransitionMatrix(nmicrostates, N_ji, nmacrostates, T_ba, micro_to_macro)  
    criteria = evaluationCriteria(T_ba)
    
    ! Store as best.
    best_criteria = criteria
    best_micro_to_macro = micro_to_macro
    
    if(debug) write(*,*) 'initial criteria: ', best_criteria
    
    ! Perform MC/SA run.
    do step = 1,mcsa_nsteps
       if(debug .and. mod(step,1000)==0) write(*,*) 'step ', step, ' of ', mcsa_nsteps, ': ', criteria, best_criteria
       
       ! Propose a move by randomly assigning a random microstate to a random macrostate.
       micro_to_macro_proposed = micro_to_macro
       microstate = generateUniformInteger(nmicrostates)
       macrostate = generateUniformInteger(nmacrostates)
       micro_to_macro_proposed(microstate) = macrostate

       ! Reject immediately if no change in lumping.
       if(micro_to_macro_proposed(microstate) == micro_to_macro(microstate)) cycle

       ! Reject immediately if one of the macrostates contains no microstates.
       if(any(histogram(micro_to_macro_proposed, nmacrostates) == 0)) cycle

       ! DEBUG
       !write(*,*) 'histogram = ', histogram(micro_to_macro_proposed(compacted_microstate_assignments), nmacrostates)
       
       ! Evaluate transition matrix and criteria.
       !T_ba = computeTransitionMatrix(ntrajectories, trajectories, nmacrostates, micro_to_macro_proposed(compacted_microstate_assignments), tau)
       call computeLumpedTransitionMatrix(nmicrostates, N_ji, nmacrostates, T_ba, micro_to_macro_proposed)  
       criteria_proposed = evaluationCriteria(T_ba)

       ! DEBUG
       !write(*,*) ' criteria: ', criteria_proposed

       ! Accept or reject based on Metropolis criteria.
       accept = .false.
       if(criteria_proposed >= criteria) then
          ! Always accept if criteria is increased.
          accept = .true.

          ! Check if criteria is optimal.
          if(criteria > best_criteria) then
             ! Store as best.
             best_criteria = criteria
             best_micro_to_macro = micro_to_macro_proposed
          end if
       else
          ! Set annealing temperature.
          beta = real(step,sp)
          ! Generate uniform variate.
          r = generateUniform01()
          ! Apply Metropolis test.
          if(r < exp(beta * (criteria_proposed - criteria))) accept = .true.
       end if

       ! Process acceptance.
       if(accept) then
          criteria = criteria_proposed             
          micro_to_macro = micro_to_macro_proposed
       end if

    end do

  end subroutine lump_mcsa

  !=============================================================================================
  ! Construct an initial random partitioning.
  !=============================================================================================
  subroutine random_lumping(nmicrostates, nmacrostates, micro_to_macro)

    use utilities, only : generateUniqueRandomIndices

    ! Arguments.
    integer, intent(in) :: nmicrostates
      ! Number of microstates.
    integer, intent(in) :: nmacrostates
      ! Number of macrostates.
    integer, dimension(nmicrostates), intent(out) :: micro_to_macro
      ! micro_to_macro(microstate) = macrostate which microstate belongs to.

    ! Local variables.
    integer :: i
      ! Loop variable.
    integer, dimension(nmicrostates) :: indices
      ! indices storage

    ! Pick an initial lumping of microstates to form nmacrostates macrostates, ensuring that each
    ! macrostate contains at least one microstate.    
    call generateUniqueRandomIndices(nmicrostates, nmicrostates, indices)
    do i = 1,nmicrostates
       micro_to_macro(indices(i)) = mod(i, nmacrostates) + 1
    end do

  end subroutine random_lumping

  !=============================================================================================
  ! Construct an initial guess at the partitioning using eigenvector-based PCCA-like algorithm.
  !=============================================================================================
  subroutine ev_lumping(nmicrostates, nmacrostates, mu_k, u_ik, micro_to_macro)
    use utilities, only : mean, find, unique, histogram

    ! Arguments.
    integer, intent(in) :: nmicrostates
      ! Number of microstates.
    integer, intent(in) :: nmacrostates
      ! Number of macrostates.
    integer, dimension(nmicrostates), intent(out) :: micro_to_macro
      ! micro_to_macro(microstate) = macrostate which microstate belongs to.
    real(dp), dimension(nmacrostates), intent(in) :: mu_k
      ! mu_k(k) is the kth largest eigenvalue of the transition matrix
    real(dp), dimension(nmicrostates,nmacrostates), intent(in) :: u_ik
      ! u_ik(i,k) is component i of eigenvector k.

    ! Local variables.
    integer :: i, j, k
      ! Loop variables.
    integer, dimension(:), pointer :: indices, indices2
      ! indices storage
    integer :: current_nmacrostates
      ! Current number of macrostates.
    real(dp), dimension(nmicrostates) :: ev
    real(sp), dimension(nmacrostates) :: deviation
      ! deviation of ev from mean over this macrostate
    real(dp) :: mean_ev
    integer :: macrostate_to_split    

    ! Start with all microstates in one macrostate.
    micro_to_macro(1:nmicrostates) = 1
    current_nmacrostates = 1

    ! Split using other eigenvectors.
    do k = 2,nmacrostates
       ! Get current eigenvector.
       ev = u_ik(:,k)

       ! Find the macrostate with the largest sum of differences between the mean eigenvector component and all of its members (L1 norm?).
       do i = 1,current_nmacrostates
          ! Get indices of microstates in macrostate i.
          indices => find(micro_to_macro .eq. i)
          
          ! Compute L1 norm of difference.
          mean_ev = mean(ev(indices))
          deviation(i) = sum(abs(ev(indices) - mean_ev))          
          
          ! Clean up
          deallocate(indices)
       end do
       macrostate_to_split = maxloc(deviation(1:current_nmacrostates),1)

       ! Split the state at the mean.
       indices => find(micro_to_macro .eq. macrostate_to_split)
       mean_ev = mean(ev(indices))
       indices2 => find(ev(indices) > mean_ev)
       micro_to_macro(indices(indices2)) = current_nmacrostates + 1
       current_nmacrostates = current_nmacrostates + 1
       deallocate(indices, indices2)       
    end do

    ! Check to make sure no macrostates are empty.
    if(any(histogram(micro_to_macro,nmacrostates) == 0)) then
       write(*,*) 'ev_lumping has generated an empty macrostate!'
       write(*,*) 'micro_to_macro = '
       write(*,*) micro_to_macro
       stop
    end if

  end subroutine ev_lumping

  !=============================================================================================
  ! Lump microstates into macrostates by kinetic clustering using a MC/SA algorithm.
  !=============================================================================================
  subroutine lump()
    use utilities, only : generateUniformInteger, generateUniform01, histogram, generateUniqueRandomIndices, unique
    use transitionmatrix, only : computeTransitionMatrix, transitionMatrixEigs, computeTransitionCountMatrix

    ! Parameters.
    
    ! Constants.

    ! Local variables.
    integer :: i, j, n
    integer :: trial
      ! current MC/SA trial (max is ntrials)
    integer :: best_trial
    real(sp) :: criteria
    integer, dimension(nmicrostates) :: micro_to_macro
    integer, dimension(nmicrostates) :: ev_micro_to_macro
      ! eigenvector-based lumping initial guess
    real(sp), dimension(:), allocatable :: trial_criteria
    integer, dimension(:,:), allocatable :: trial_micro_to_macro
    real(sp) :: criteria_upper_bound
    real(dp), dimension(nmacrostates) :: mu_k
      ! Dominant eigenvalues of the transition matrix.
    real(dp), dimension(nmicrostates,nmacrostates) :: u_ik
      ! Dominant eigenvectors of the transition matrix.
    integer :: status
    real(dp), dimension(nmicrostates,nmicrostates) :: N_ji
      ! symmetrized transition count matrix to use for efficiently computing lumped transition matrix

    write(*,*) 'Histogram of microstate populations:'
    write(*,'(100I6)') histogram(compacted_microstate_assignments, nmicrostates)
    write(*,*) 'total = ', sum(histogram(compacted_microstate_assignments, nmicrostates))
    if(debug) write(*,*) 'Lumping ', nmicrostates, ' microstates into ', nmacrostates, ' macrostates...'

    ! Compute eigenvalues and eigenvectors of transition matrix.
    if(debug) write(*,*) 'Computing eigenvalues and eigenvectors of transition matrix...'
    call transitionMatrixEigs(ntrajectories, trajectories, nmicrostates, compacted_microstate_assignments, &
         tau, nmacrostates, mu_k, ev_left = u_ik)
    if(debug) write(*,*) 'mu_k = ', mu_k

    ! Determine upper bound on metastability from eigenvalues.
    criteria_upper_bound = sum(mu_k)

    if(use_ev_initial_lumping) then
       ! Compute eigenvector-based initial lumping.
       write(*,*) 'Computing eigenvector-based lumping to use as initial guess...'
       call ev_lumping(nmicrostates, nmacrostates, mu_k, u_ik, ev_micro_to_macro)
    end if

    ! Compute symmetrized transition count matrix to use for efficiently computing transition matrix.
    call computeTransitionCountMatrix(ntrajectories, trajectories, nmicrostates, compacted_microstate_assignments, tau, N_ji, &
         use_timereversed)
    
    ! Perform trials.
    allocate(trial_criteria(mcsa_ntrials), trial_micro_to_macro(mcsa_ntrials, nmicrostates))
    !$omp parallel do default(shared) private(criteria,micro_to_macro)
    do trial = 1, mcsa_ntrials       
       if(debug) write(*,*) 'mcsa trial ', trial, ' of ', mcsa_ntrials

       ! Select initial lumping.
       if(use_ev_initial_lumping) then
          ! Use ev-based lumping.
          micro_to_macro = ev_micro_to_macro
       else
          call random_lumping(nmicrostates, nmacrostates, micro_to_macro)          
       end if

       ! Perform one pass of MC/SA to obtain a guess at the best criteria and lumping.
       call lump_mcsa(nmacrostates, criteria, micro_to_macro, N_ji)
       
       ! Store this lumping.
       ! TODO: Compact this into above call.
       trial_criteria(trial) = criteria
       trial_micro_to_macro(trial,:) = micro_to_macro
    end do
    !$omp end parallel do

    ! Determine best and convergence.
    write(*,*) 'trial_criteria'
    write(*,*) trial_criteria
    best_trial = maxloc(trial_criteria,1)
    criteria = trial_criteria(best_trial)
    micro_to_macro = trial_micro_to_macro(best_trial,:)

!    write(*,*) 'optimal lumping = '
!    write(*,'(I4)') micro_to_macro
    write(*,*) count(trial_criteria == criteria), ' / ', mcsa_ntrials, ' trials found the optimum.'
    write(*,*) 'criteria is ', criteria, ' / ', criteria_upper_bound

    ! Use this lumping to generate compacted macrostates assignments.
    compacted_macrostate_assignments = micro_to_macro(compacted_microstate_assignments)

    ! Determine microstate ids that go with macrostates.
    deallocate(macrostates, stat=status)
    allocate(macrostates(nmacrostates))
    do i = 1,nmicrostates
       macrostates(micro_to_macro(i)) = microstates(i)
    end do
    
    ! Map compacted to expanded macrostate ids.
    macrostate_assignments = macrostates(compacted_macrostate_assignments)

    deallocate(trial_criteria, trial_micro_to_macro)
    
  end subroutine lump

  !=============================================================================================  
  ! Write state assignments.
  !=============================================================================================
  subroutine writeStateassignments(filename, nstates, stateassignments)
    use pdbio, only : new_unit
    
    ! Parameters.
    character(len=*), intent(in) :: filename
      ! File to write stateassignments to.
    integer, intent(in) :: nstates
      ! Number of states.
    integer, dimension(:), intent(in) :: stateassignments
      ! Individual state assignments.

    ! Local variables.
    integer :: nsnapshots
    integer :: snapshot_index
    integer :: iunit
      ! unit to use for file management
    character(len=MAX_LINE_LENGTH) :: line, text

    ! Determine list size.
    nsnapshots = size(stateassignments,1)
    
    ! DEBUG
    if(debug) write(*,*) 'Writing state assignments to ', trim(filename), '...'
    
    ! Open file for writing.
    call new_unit(iunit)
    open(unit=iunit, file=filename, status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    
    ! Write.
    do snapshot_index = 1,nsnapshots
       ! Write line.
       write(iunit,'(I8,X,I8)') snapshot_index, stateassignments(snapshot_index)
    end do

    ! Close file.
    close(unit=iunit)
    if(debug) write(*,*) 'Done.'
    
    return

    ! Report error if file could not be opened.
10  continue
    write(*,*) 'decompose: writeStateassignments: Error opening file ', trim(filename)
    stop

  end subroutine writeStateassignments

  !=============================================================================================  
  ! Read state assignments.
  !=============================================================================================
  subroutine readStateassignments(filename, stateassignments)
    use pdbio, only : new_unit
    
    ! Parameters.
    character(len=*), intent(in) :: filename
      ! File to write stateassignments to.
    integer, dimension(:), intent(out) :: stateassignments
      ! Individual state assignments.

    ! Local variables.
    integer :: nsnapshots
    integer :: snapshot_index
    integer :: iunit
      ! unit to use for file management
    character(len=MAX_LINE_LENGTH) :: line, text

    ! Determine list size.
    nsnapshots = size(stateassignments,1)
    
    ! DEBUG
    if(debug) write(*,*) 'Reading state assignments from ', trim(filename), '...'
    
    ! Open file for reading.
    call new_unit(iunit)
    open(unit=iunit, file=trim(filename), status='OLD', err=10, position='REWIND', &
         form='FORMATTED', action='READ')    
    
    ! Read.
    read(iunit,'(9X,I8)') stateassignments

    ! Close file.
    close(unit=iunit)
    if(debug) write(*,*) 'Done.'
    
    return

    ! Report error if file could not be opened.
10  continue
    write(*,*) 'decompose: readStateassignments: Error opening file ', trim(filename)
    stop

  end subroutine readStateassignments
  
  !=============================================================================================  
  ! Write eigenvectors.
  !=============================================================================================  
  subroutine writeEigenvectors(prefix, nstates, compacted_stateassignments, tau, neigs)
    use pdbio, only : new_unit
    use transitionmatrix, only : transitionMatrixEigs
    use timer
    
    ! Parameters.
    character(len=*), intent(in) :: prefix
      ! File to write stateassignments to.
    integer, intent(in) :: nstates
      ! Number of states.
    integer, dimension(:), intent(in) :: compacted_stateassignments
      ! Individual state assignments
    integer, intent(in) :: tau
      ! Lag time to use
    integer, intent(in) :: neigs
      ! number of eigenvectors to write

    ! Local variables.
    integer :: nsnapshots
    integer :: snapshot_index
    integer :: iunit
      ! unit to use for file management
    character(len=MAX_LINE_LENGTH) :: line, text
    integer :: ntau 
      ! Number of lag times.
    integer, dimension(:), allocatable :: taus
    real(sp), dimension(:,:), allocatable :: timescales, timescales_lower, timescales_upper
      ! timescales(timescale_index, tau_index)
      ! lower and upper are confidence bounds
    integer :: i, k
      ! Loop variables.
    character(len=MAX_FILENAME_LENGTH) :: filename
      ! Filename to write to.
    integer :: tau_step_
      ! tau increment to use
    integer :: ntimescales
      ! number of timescales to use
    real(dp), dimension(nstates, neigs) :: ev_left, ev_right
    real(dp), dimension(neigs) :: mu_k

    if(debug) write(*,*) 'Computing eigenvalues and eigenvectors of transition matrix...'
    call transitionMatrixEigs(ntrajectories, trajectories, nstates, compacted_stateassignments, &
         tau, neigs, mu_k, ev_left = ev_left, ev_right = ev_right)
    if(debug) write(*,*) 'mu_k = ', mu_k

    ! Get a new file unit to use.
    call new_unit(iunit)

    ! Write eigenvectors.
    filename = trim(prefix)//'.ev_right'
    open(unit=iunit, file=filename, status='REPLACE', err=11, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    do k = 1,neigs
       line = ''
       do i = 1,nstates
          write(text,'(E16.8)') ev_right(i,k)
          line = trim(line)//' '//adjustl(text)
       end do

       ! Write to file.
       write(iunit,'(A)',advance='no') trim(line)
       write(iunit,'(A)') ''
    end do
    close(unit=iunit)

    ! Write eigenvectors.
    filename = trim(prefix)//'.ev_left'
    open(unit=iunit, file=filename, status='REPLACE', err=11, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    do k = 1,neigs
       line = ''
       do i = 1,nstates
          write(text,'(E16.8)') ev_left(i,k)
          line = trim(line)//' '//adjustl(text)
       end do

       ! Write to file.
       write(iunit,'(A)',advance='no') trim(line)
       write(iunit,'(A)') ''
    end do
    close(unit=iunit)

    if(debug) write(*,*) 'Done.'
    
    return

    ! Report error if file could not be opened.
11  continue
    write(*,*) 'decompose: writeEigenvectors: Error opening file ', trim(filename)
    stop

  end subroutine writeEigenvectors

  !=============================================================================================  
  ! Write timescales.
  !=============================================================================================
  subroutine writeTimescales(prefix, nstates, stateassignments, max_tau, tau_step)
    use pdbio, only : new_unit
    use transitionmatrix, only : computeTimescalesBootstrap
    use timer
    
    ! Parameters.
    character(len=*), intent(in) :: prefix
      ! File to write stateassignments to.
    integer, intent(in) :: nstates
      ! Number of states.
    integer, dimension(:), intent(in) :: stateassignments
      ! Individual state assignments
    integer, intent(in) :: max_tau
      ! Largest tau to use.  Must be shorter than longest trajectory.
    integer, intent(in), optional :: tau_step
      ! Increment to use in choosing taus.

    ! Local variables.
    integer :: nsnapshots
    integer :: snapshot_index
    integer :: iunit
      ! unit to use for file management
    character(len=MAX_LINE_LENGTH) :: line, text
    integer :: ntau 
      ! Number of lag times.
    integer, dimension(:), allocatable :: taus
    real(sp), dimension(:,:), allocatable :: timescales, timescales_lower, timescales_upper
      ! timescales(timescale_index, tau_index)
      ! lower and upper are confidence bounds
    integer :: i, tau_index
      ! Loop variables.
    character(len=MAX_FILENAME_LENGTH) :: filename
      ! Filename to write to.
    integer :: tau_step_
      ! tau increment to use
    integer :: ntimescales
      ! number of timescales to use

    if(debug) write(*,*) 'Computing timescales...'
    call resetTimer()

    ! Determine number of snapshots from stateassignments list.
    nsnapshots = size(stateassignments,1)

    ! Choose taus to evaluate at.
    tau_step_ = 1
    if(present(tau_step)) tau_step_ = tau_step
    ntau = int(max_tau / tau_step_)
    allocate(taus(ntau))
    taus = (/ (i, i = tau_step_,max_tau,tau_step_) /)

    ! Choose number of timescales based on number of states.
    ntimescales = min(nstates - 2, max_ntimescales)

    ! Allocate storage.
    allocate(timescales(ntimescales,ntau), timescales_lower(ntimescales,ntau), timescales_upper(ntimescales, ntau))

    ! Compute timescales for all taus.
    call computeTimescalesBootstrap(ntrajectories, trajectories, nstates, stateassignments, ntau, ntimescales, &
         tau_unit, taus, timescales, &
         timescales_confidence_interval, timescales_lower, timescales_upper, timescales_bootstrap_ntrials)

    ! Get a new file unit to use.
    call new_unit(iunit)

    ! Write lag times.    
    filename = trim(prefix)//'.lag_times'
    open(unit=iunit, file=filename, status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    do tau_index = 1,ntau
       write(iunit,*) taus(tau_index) * tau_unit
    end do
    close(unit=iunit)

    ! Write timescales.
    filename = trim(prefix)//'.timescales'
    open(unit=iunit, file=filename, status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    do tau_index = 1,ntau
       line = ''
       write(text,'(F12.3)') taus(tau_index) * tau_unit
       line = adjustl(text)//' '
       
       do i = 1,ntimescales
          write(text,'(E18.8)') timescales(i,tau_index)
          line = trim(line)//' '//adjustl(text)
       end do

       ! Write to file.
       write(iunit,'(A)',advance='no') trim(line)
       write(iunit,'(A)') ''
    end do
    close(unit=iunit)

    ! Write confidence bounds.
    filename = trim(prefix)//'.timescales_lower'
    open(unit=iunit, file=filename, status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    do tau_index = 1,ntau
       line = ''
       write(text,'(F12.3)') taus(tau_index) * tau_unit
       line = adjustl(text)//' '

       do i = 1,ntimescales
          write(text,'(E18.8)') timescales_lower(i,tau_index)
          line = trim(line)//' '//adjustl(text)
       end do

       ! Write to file.
       write(iunit,'(A)',advance='no') trim(line)
       write(iunit,'(A)') ''
    end do
    close(unit=iunit)

    filename = trim(prefix)//'.timescales_upper'
    open(unit=iunit, file=filename, status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')    
    do tau_index = 1,ntau
       line = ''
       write(text,'(F12.3)') taus(tau_index) * tau_unit
       line = adjustl(text)//' '

       do i = 1,ntimescales
          write(text,'(E18.8)') timescales_upper(i,tau_index)
          line = trim(line)//' '//adjustl(text)
       end do

       ! Write to file.
       write(iunit,'(A)',advance='no') trim(line)
       write(iunit,'(A)') ''
    end do
    close(unit=iunit)

    ! Clean up
    deallocate(taus,timescales,timescales_lower,timescales_upper)
    if(debug) write(*,*) readTimer(), ' seconds elapsed'

    return

    ! Report error if file could not be opened.
10  continue
    write(*,*) 'decompose: writeTimescales: Error opening file ', trim(filename)
    stop
    
  end subroutine writeTimescales

  !=============================================================================================
  ! Write transformation from initial states to final states to file.
  ! 
  ! For both splitting and lumping, the format is:
  !
  ! [lumped state] [splitstate1] [splitstate2] ... [splitstateN]
  !
  ! Note that, for splitting, the initial state before lumping is on the left and the states it was subsequently split into are on the right.
  ! For lumping, it is the opposite: The initial states that formed the lumped state are on the right, and the lumped state is on the left.
  !
  ! There is always a one-to-many relation. 
  !
  ! NOTE: We require ninitialstates > nfinalstates.
  !=============================================================================================
  subroutine writeStateTransformation(filename, lumped_stateassignments, split_stateassignments)
    use pdbio, only : new_unit
    use utilities, only : unique, find

    ! Parameters.
    character(len=*), intent(in) :: filename
      ! File to write stateassignments to.
    integer, dimension(:), intent(in) :: lumped_stateassignments, split_stateassignments
      ! Individual state assignments.
    
    ! Local variables.
    integer :: iunit
      ! unit to use for file management
    integer :: nsnapshots
    integer :: snapshot_index, state_index
    character(len=MAX_LINE_LENGTH) :: line, text   
    integer :: last_state
    integer, dimension(:), pointer :: snapshot_indices
      ! indices of snapshots matching the current state
    integer, dimension(:), pointer :: unique_lumped_states, unique_split_states
    integer :: state
      ! current state id
    integer :: nunique_split_states
    integer :: nlumpedstates
    integer :: i

    if(debug) write(*,*) 'Writing state transformation table...'

    ! Get number of snapshots.
    nsnapshots = size(lumped_stateassignments,1)
    write(*,*) 'nsnapshots = ', nsnapshots

    ! Get a new file unit to use.
    call new_unit(iunit)

    ! Open file for formatted writing.
    open(unit=iunit, file=filename, status='REPLACE', err=11, position='REWIND', &
         form='FORMATTED', action='WRITE')    

    ! Make a list of unique lumped states.
    unique_lumped_states => unique(lumped_stateassignments)
    nlumpedstates = size(unique_lumped_states,1)

    ! List all final states corresponding to each unique initial state.
    do state_index = 1,nlumpedstates
       ! Get current state.
       state = unique_lumped_states(state_index)

       ! Find all snapshots in this state.
       snapshot_indices => find(lumped_stateassignments == state)

       ! Make a list of all unique final states.
       unique_split_states => unique(split_stateassignments(snapshot_indices))
       nunique_split_states = size(unique_split_states,1)

       ! Construct a line of text to write to file.
       line = ''
       write(text,*) state
       line = trim(line)//adjustl(text)//' :'

       do i = 1,nunique_split_states
          write(text,*) unique_split_states(i)
          line = trim(line)//' '//adjustl(text)
       end do

       ! Write to file.
       write(iunit,'(A)',advance='no') trim(line)
       write(iunit,'(A)') ''

       ! Clean up.
       deallocate(snapshot_indices,unique_split_states)       
    end do
    
    ! Close file.
    close(unit=iunit)

    ! Clean up.
    deallocate(unique_lumped_states)

    if(debug) write(*,*) 'Done.'

    return

    ! Report error if file could not be opened.
11  continue
    write(*,*) 'decompose: writeStateTransformation: Error opening file ', trim(filename)
    stop


  end subroutine writeStateTransformation

  !=============================================================================================
  ! Write PDB files containing an ensemble of configurations chosen at random from the state.
  !=============================================================================================
  subroutine writeStatePDBs(prefix, nstates, state_representatives, stateassignments, nmodels_to_write)

    use utilities, only : find, generateUniqueRandomIndices
    use pdbio, only : pdb_write
    use timer
    use kabsch_rmsd, only : ls_align

    ! Arguments.
    character(len=*), intent(in) :: prefix
      ! Prefix of filename to write stateassignments to.
    integer, intent(in) :: nstates
      ! Number of states.
    integer, dimension(nstates), intent(in) :: state_representatives
      ! Indices of state representatives.
    integer, dimension(:), intent(in) :: stateassignments
      ! Individual state assignments
    integer, intent(in) :: nmodels_to_write
      ! Number of models to write to PDB file.

    ! Local variables.
    integer :: state_index
      ! State index, from 1...nstates.
    integer :: state_representative
      ! Index of representative snapshot.
    integer, dimension(:), pointer :: members
      ! Indices of snapshots that belong to the current state under consideration.
    integer :: nmembers
      ! Length of index list.
    integer :: nmodels
      ! Number of models to actually write for this PDB file (may be less than nmodels_to_write
      ! if there are fewer members than desired models).
    integer, dimension(:), allocatable :: indices
      ! Unique list of snapshot indices of models to write.
    character(len=MAX_FILENAME_LENGTH) :: filename
      ! Filename to write to.
    character(len=59), dimension(1) :: remarks
      ! optional remarks to put in header
    real(dp) :: elapsed_seconds
      ! Number of seconds elapsed.
    real(sp), dimension(:,:,:), allocatable :: models
      ! Models that have been aligned to the state representative.
    real(sp), dimension(natoms,3) :: representative_snapshot
      ! The representative snapshot, aligned to the reference snapshot.
    integer :: model_index
      ! Loop index.
    
    if(debug) write(*,*) 'Writing states to PDB files...'
    
    ! Loop over all states.
    call resetTimer()
    do state_index = 1,nstates
       ! Get the index of the state representative snapshot.
       state_representative = state_representatives(state_index)

       ! Skip this state if the representative is set to 0.
       if (state_representative .eq. 0) cycle

       ! Store the state representative as the first model.
       representative_snapshot = snapshots(state_representative,:,:)

       ! Align the representative to the 'reference' PDB file using backbone indices.
       call ls_align(natoms, representative_snapshot, reference_snapshot, ls_align_atomindices)

       ! Build a list of all snapshots in this state.
       members => find(stateassignments == state_representative)
       nmembers = size(members,1)
       
       ! Determine how many models to write.
       ! If requested to write more models than members, make sure that only nmembers models are written.
       nmodels = nmodels_to_write
       if(nmembers < nmodels) nmodels = nmembers

       ! Choose models at random from list of snapshots belonging to state.
       allocate(indices(nmodels))
       call generateUniqueRandomIndices(nmembers, nmodels, indices)       
       indices = members(indices)

       ! Replace the first model with the state representative.
       indices(1) = state_representative

       ! Store and align models.
       allocate(models(nmodels,natoms,3))
       models = snapshots(indices,:,:)
       do model_index = 1,nmodels
          call ls_align(natoms, models(model_index,:,:), representative_snapshot, ls_align_atomindices)
       end do
       
       ! Construct filename as '${prefix}-state_id.pdb'
       write(filename,*) state_representative
       filename = trim(prefix)//'-'//trim(adjustl(filename))//'.pdb'
       
       ! Construct remarks text, indicating state_id and number of members.
       write(remarks,*) 'State', state_representative, ', ', nmembers, ' members'

       ! Write models to PDB file.
       call pdb_write(filename, pdbatoms, models, remarks=remarks)

       ! Deallocate members list.
       deallocate(members, indices, models)
    end do
    if(debug) write(*,*) readTimer(), ' seconds elapsed'

  end subroutine writeStatePDBs

  !=============================================================================================
  ! Split all macrostates into microstates.
  !
  ! Each macrostate, if above a certain threshold size, is split into a number of microstates
  ! determined by its size.  These microstates are assigned as state identities the indices of
  ! their generators.
  !=============================================================================================
  subroutine split()
    use utilities, only : find, unique

    ! Local variables.
    integer :: macrostate, compacted_macrostate_index
      ! macrostate index
    integer, dimension(:), pointer :: indices
      ! snapshot indices of macrostate members
    integer :: nindices
      ! number of members of the macrostate
    integer :: nmembers
      ! Number of members of this state.
    integer :: nmaxmicrostates
      ! Maximum possible number of microstates.
    integer :: target_nmicrostates
      ! Target number of microstates to split to.
    integer :: status
    integer :: counter

    counter = 0    
    ! Attempt to split each macrostate.
    do compacted_macrostate_index = 1,nmacrostates
       ! Get the unique macrostate id for this state.
       macrostate = macrostates(compacted_macrostate_index)

       ! Get snapshot indices of all members of this macrostate.
       indices => find(compacted_macrostate_assignments == compacted_macrostate_index)

       ! Determine number of snapshots in this macrostate.
       nindices = size(indices,1)
       counter = counter + nindices

       ! Determine how many microstates to split to by using the given target microstate size.
       target_nmicrostates = (nindices / target_microstate_size)

       ! Limit the number of microstates we can split an individual macrostate into, if this option has been specified.
       if(iteration == 1) then
          if(max_microstates_per_macrostate_first_iteration >= 0) &
               target_nmicrostates = min(target_nmicrostates, max_microstates_per_macrostate_first_iteration)
       else
          if(max_microstates_per_macrostate >= 0) target_nmicrostates = min(target_nmicrostates, max_microstates_per_macrostate)
       end if

       ! Don't do anything if the state is too small to be split.
       if(target_nmicrostates <= 1) then
          ! Assign microstates from macrostate.
          microstate_assignments(indices) = macrostate          
          ! Cycle.
          cycle
       end if
       
       if(debug) write(*,*) 'Splitting macrostate ', compacted_macrostate_index, ' of ', nmacrostates, ', id ', macrostate, &
            ' (', nindices, ' members) into ', target_nmicrostates, ' microstates'

       ! Split the macrostate into microstates, assigning them unique ids by their generator.
       call splitState(nindices, indices, target_nmicrostates, microstate_assignments)
       
       ! Free list of macrostate members.
       deallocate(indices)
    end do

    write(*,*) 'counter = ', counter
    indices => find(microstate_assignments == 0)
    write(*,*) 'find found ', size(indices,1), ' zero indices:'
    write(*,*) indices

    ! Generate compacted microstate assignments, where indices go from 1...nmicrostates.
    deallocate(microstates, stat=status)
    microstates => unique(microstate_assignments, compacted=compacted_microstate_assignments)
    nmicrostates = size(microstates,1)
    call checkstates(microstate_assignments)
    
  end subroutine split

  !=============================================================================================
  ! Generate and refine a state decomposition by iterative splitting and lumping.
  !=============================================================================================
  subroutine iterate()
    use timer ! for timing info

    ! Local variables.
    real(dp) :: elapsed_seconds
      ! Time elapsed during each call.
    character(len=MAX_FILENAME_LENGTH) :: filename
      ! Storage for constructing filenames for writing out progress files.

    ! Iterate to refine state decomposition.
    do iteration = iteration, niterations
       if(debug) write(*,*) 'State decomposition iteration ', iteration, ' of ', niterations, '...'

       !=============================================================================================
       ! Split to microstates.
       !=============================================================================================

       ! Determine the distance metric to use for this iteration.
       if(iteration == 1) then 
          rmsd_mode = rmsd_mode_first_iteration
       else
          rmsd_mode = rmsd_mode_subsequent_iterations
       end if
       if(debug) write(*,*) 'Using rmsd_mode = ', rmsd_mode
       
       ! Split all states.
       call resetTimer()
       call split()
       write(*,*) 'checking microstate assignments...'
       call checkstates(microstate_assignments)
       call checkstates(compacted_microstate_assignments)
       if(debug) write(*,*) readTimer(), ' seconds elapsed'
       
       ! Write split stateassignments.
       write(filename,*) iteration
       filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.stateassignments-split'
       call writeStateassignments(filename, nmicrostates, microstate_assignments)
       write(filename,*) iteration
       filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.stateassignments-split-compacted'
       call writeStateassignments(filename, nmicrostates, compacted_microstate_assignments)

       ! Write microstates.
       write(filename,*) iteration
       filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.microstates'
       call writeStateassignments(filename, nmicrostates, microstates)

       if(write_split_table) then
          ! TODO: Write translation from lumped -> split state ids.
          call resetTimer()
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.split'
          call writeStateTransformation(filename, macrostate_assignments, microstate_assignments)
          if(debug) write(*,*) readTimer(), ' seconds elapsed'                 
       end if

       if(write_split_timescales) then
          ! Write timescales.
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.timescales-split'
          call writeTimescales(filename, nmicrostates, compacted_microstate_assignments, max_tau)
       end if

       if(write_split_pdbs) then
          ! Write state PDBs.
          call resetTimer()
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.ensemble-split'
          call writeStatePDBs(filename, nmicrostates, microstates, microstate_assignments, split_pdb_nmodels)
          if(debug) write(*,*) readTimer(), ' seconds elapsed'       
       end if

       !=============================================================================================
       ! Lump microstates to macrostates.
       !=============================================================================================

       ! Lump to desired number of macrostates.
       call resetTimer()
       nmacrostates = target_nmacrostates
       call lump()
       call checkstates(macrostate_assignments)
       call checkstates(compacted_microstate_assignments)
       if(debug) write(*,*) readTimer(), ' seconds elapsed'       

       ! Write lumped state assignments.
       write(filename,*) iteration
       filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.stateassignments-merged'
       call writeStateassignments(filename, nmacrostates, macrostate_assignments)
       write(filename,*) iteration
       filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.stateassignments-merged-compacted'
       call writeStateassignments(filename, nmacrostates, compacted_macrostate_assignments)

       ! Write macrostates.
       write(filename,*) iteration
       filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.macrostates'
       call writeStateassignments(filename, nmacrostates, macrostates)

       if(write_lump_table) then
          ! TODO: Write translation from split -> lumped state ids.
          call resetTimer()
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.lump'
          call writeStateTransformation(filename, macrostate_assignments, microstate_assignments)
          if(debug) write(*,*) readTimer(), ' seconds elapsed'                 
       end if

       if(write_merged_timescales) then
          ! Write timescales.
          call resetTimer()
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.timescales-merged'
          call writeTimescales(filename, nmacrostates, compacted_macrostate_assignments, max_tau)
          if(debug) write(*,*) readTimer(), ' seconds elapsed'       
       end if

       if(write_merged_eigenvectors) then
          ! Write eigenvectors.
          call resetTimer()
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.merged'
          call writeEigenvectors(filename, nmacrostates, compacted_macrostate_assignments, tau, min(max_ntimescales+1,nmacrostates))
          if(debug) write(*,*) readTimer(), ' seconds elapsed'       
       end if

       if(write_merged_pdbs) then
          ! Write state PDBs.
          call resetTimer()
          write(filename,*) iteration
          filename = trim(output_directory)//'/'//trim(adjustl(filename))//'.ensemble-merged'
          call writeStatePDBs(filename, nmacrostates, macrostates, macrostate_assignments, merged_pdb_nmodels)
          if(debug) write(*,*) readTimer(), ' seconds elapsed'       
       end if
    end do
    
  end subroutine iterate

end module decompose
