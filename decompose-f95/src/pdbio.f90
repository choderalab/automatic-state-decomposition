!=============================================================================================
! A module for reading and writing PDB files containing multiple models.
!
! The Brookhaven PDB format is documented here:
! http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/format_guide_noframes.html
!
! Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
!
! Copyright (c) 2006 The Regents of the University of California. 
! All Rights Reserved.
!=============================================================================================
module pdbio

  use numeric_kinds ! for precision

  implicit none

  private
  public :: pdb_write, pdb_read, new_unit
 
  !=============================================================================================
  ! Constants.
  !=============================================================================================

  logical, parameter :: debug = .FALSE.

  !=============================================================================================
  ! Data types.
  !=============================================================================================

  ! PDB atom information (excluding coordinates)
  type, public :: pdbatom_t
     integer :: serial                    ! atom serial number within model (from 1..natoms)
     character(len=4) :: name             ! atom name
     character :: altLoc                  ! alternate location indicator
     character(len=3) :: resName          ! residue name
     character :: chainID                 ! chain identifier
     integer :: resSeq                    ! residue sequence number
     character :: iCode                   ! code for insertion of residues
     real(sp) :: x, y, z                  ! Cartesian coordinates (may be unused)
     real(sp) :: occupancy                ! occupancy
     real(sp) :: tempFactor               ! temperature factor     
     character(len=4) :: segID            ! segment identifier, left-justified
     character(len=2) :: element          ! element symbol, right-justified
     character(len=2) :: charge           ! charge on the atom
  end type pdbatom_t

  !=============================================================================================
  ! Private module data.
  !=============================================================================================

contains
  
  !*******************************************************************************
  !
  !! GET_UNIT returns a free FORTRAN unit number.
  !
  !  Discussion:
  !
  !    A "free" FORTRAN unit number is an integer between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer IUNIT.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5 and 6).
  !
  !    Otherwise, IUNIT is an integer between 1 and 99, representing a
  !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  !    are special, and will never return those values.
  !
  !*******************************************************************************
  subroutine new_unit(iunit)
    
    ! Parameters.
    integer, intent(out) :: iunit
    
    ! Local variables.
    integer :: i
    integer :: ios
    logical :: lopen

    ! Set inuit to 0 (free unit not found) by default.
    iunit = 0

    ! Try to find an unused unit number.
    do i = 1, 99       
       ! Exclude stdout and stderr.
       if ( i == 5 .or. i == 6 ) cycle

       ! Inquire if the unit is available.
       inquire ( unit = i, opened = lopen, iostat = ios )
       if ( ios == 0 .and. .not. lopen ) then
          ! Unit found -- return.
          iunit = i
          return
       end if
       
    end do
    
  end subroutine new_unit

  !=============================================================================================
  ! Write a specified set of models to a PDB file.
  !=============================================================================================    
  subroutine pdb_write(filename, atoms, models, remarks)

    ! Parameters.
    character(len=*), intent(in) :: filename
      ! name of file for PDB to be written to
    type(pdbatom_t), dimension(:), intent(in) :: atoms 
      ! auxiliary (non-coordinate) on all atoms in the molecule
    real(sp), dimension(:,:,:), intent(in) :: models
      ! models to be written to PDB file
      ! models(i,n,k) is coordinate k of atom n of model i
    character(len=59), dimension(:), intent(in), optional :: remarks
      ! optional remarks to put in header

    ! Local variables.
    integer :: nmodels
      ! the number of models to be written
    integer :: natoms 
      ! the number of atoms
    integer :: iunit
      ! unit to use for file management
    integer :: i, n
      ! loop indices
    type(pdbatom_t) :: atom
      ! temporary PDB atom container
    integer :: serial
      ! serial index of record, incremented after each ATOM and TER record written
    integer :: remarkNum
      ! index of remark

    ! Determine the number of models to be written.
    nmodels = size(models,1)
    
    ! Determine the number of atoms.
    natoms = size(models,2)        

    ! TODO: Check to make sure size(models,3) == 3.

    ! Open PDB file for writing.
    call new_unit(iunit)
    open(unit=iunit, file=filename, status='REPLACE', err=10, position='REWIND', &
         form='FORMATTED', action='WRITE')

    ! ---------------------------------------------------------------------------------
    ! Write remarks.
    !
    ! COLUMNS      DATA TYPE       FIELD          DEFINITION
    ! ---------------------------------------------------------------------------------
    ! 1 -  6      Record name     "REMARK"
    !
    ! 8 - 10      Integer         remarkNum      Remark number. It is not an error
    !                                            for remark n to exist in an entry
    !                                            when remark n-1 does not.
    !
    !12 - 70      LString         empty          Left as white space in first line of
    !                                            each new remark.
    ! ---------------------------------------------------------------------------------
    if(present(remarks)) then
       remarkNum = 1
       write(iunit, '(A6,1X,I3,1X,A59)') 'REMARK', remarkNum, ''
       do i = 1,size(remarks,1)
          write(iunit, '(A6,1X,I3,1X,A59)') 'REMARK', remarkNum, remarks(i)
       end do
       write(iunit, '(A6,1X,I3,1X,A59)') 'REMARK', remarkNum, ''
    end if

    ! Write models.
    serial = 1
    do i = 1,nmodels
       ! ---------------------------------------------------------------------------------
       ! Write model header.
       !
       ! COLUMNS       DATA TYPE      FIELD         DEFINITION
       ! ----------------------------------------------------------------------
       !  1 -  6       Record name    "MODEL "
       ! 11 - 14       Integer        serial        Model serial number.
       ! ---------------------------------------------------------------------------------
       write(iunit, '(A6,4X,I4)') 'MODEL ', i

       ! ---------------------------------------------------------------------------------
       ! Write atom records.
       !
       ! COLUMNS        DATA TYPE       FIELD         DEFINITION
       ! ---------------------------------------------------------------------------------
       !  1 -  6        Record name     "ATOM  "
       !  7 - 11        Integer         serial        Atom serial number.
       ! 13 - 16        Atom            name          Atom name.
       ! 17             Character       altLoc        Alternate location indicator.
       ! 18 - 20        Residue name    resName       Residue name.
       ! 22             Character       chainID       Chain identifier.
       ! 23 - 26        Integer         resSeq        Residue sequence number.
       ! 27             AChar           iCode         Code for insertion of residues.
       ! 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
       !                                              Angstroms.
       ! 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
       !                                              Angstroms.
       ! 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
       !                                              Angstroms.
       ! 55 - 60        Real(6.2)       occupancy     Occupancy.
       ! 61 - 66        Real(6.2)       tempFactor    Temperature factor.
       ! 73 - 76        LString(4)      segID         Segment identifier, left-justified.
       ! 77 - 78        LString(2)      element       Element symbol, right-justified.
       ! 79 - 80        LString(2)      charge        Charge on the atom.
       ! ---------------------------------------------------------------------------------
       do n = 1,natoms
          ! Get all data for this atom.
          atom = atoms(n)

          ! Fill in coordinates.
          atom%x = models(i,n,1)
          atom%y = models(i,n,2)
          atom%z = models(i,n,3)
                             
          ! Write atom record.
          write(iunit, '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)') 'ATOM  ', serial, atom%name, atom%altLoc, &
               atom%resName, atom%chainID, atom%resSeq, atom%icode, atom%x, atom%y, atom%z, atom%occupancy, atom%tempFactor, &
               atom%segID, atom%element, atom%charge
          
          ! Increment serial record counter.
          serial = serial + 1
       end do

       ! ---------------------------------------------------------------------------------
       ! Write terminator.
       !
       ! COLUMNS         DATA TYPE         FIELD        DEFINITION
       ! -------------------------------------------------------------------------
       !  1 -  6         Record name       "TER   "
       !  7 - 11         Integer           serial       Serial number.
       ! 18 - 20         Residue name      resName      Residue name.
       ! 22              Character         chainID      Chain identifier.
       ! 23 - 26         Integer           resSeq       Residue sequence number.
       ! 27              AChar             iCode        Insertion code.
       ! ---------------------------------------------------------------------------------
       write(iunit, '(A6,I5,6X,A3,1X,A1,I4,A1)') 'TER   ', serial, atom%resName, atom%chainID, atom%resSeq, atom%icode
       serial = serial + 1

       ! ---------------------------------------------------------------------------------
       ! Write end-of-model record.
       !
       ! COLUMNS         DATA TYPE        FIELD           DEFINITION
       ! ------------------------------------------------------------------
       ! 1 -  6         Record name      "ENDMDL"
       ! ---------------------------------------------------------------------------------
       write(iunit, '(A6)') 'ENDMDL'
    end do

    ! Close the PDB file.
    close(iunit)

    ! Return
    return
    
    ! Report error if file could not be opened.
10  continue
    write(*,*) 'pdbio: pdb_write: Error opening file ', filename
    stop

  end subroutine pdb_write

  !=============================================================================================
  ! Read models from a PDB file, or query for number of atoms or models.
  !
  ! Example usage to read models:
  !
  ! ! Query for number of atoms and models.
  ! call pdb_read(reference_pdb_filename, natoms=natoms, nmodels=nmodels)
  ! ! Allocate storage for atom records (metadata) and model coordinates.
  ! allocate(pdbatoms(natoms),pdbmodels(nmodels,natoms,3))
  ! ! Read atom records and model coordinates.
  ! call pdb_read(reference_pdb_filename, atoms=pdbatoms, models=pdbmodels)
  !=============================================================================================
  subroutine pdb_read(filename, atoms, models, natoms, nmodels)
    
    ! Parameters.
    character(*), intent(in) :: filename
      ! name of PDB file
    type(pdbatom_t), dimension(:), intent(out), optional :: atoms 
      ! auxiliary (non-coordinate) on all atoms in the molecule
    real(sp), dimension(:,:,:), intent(out), optional :: models
      ! models to be written to PDB file
      ! models(i,n,k) is coordinate k of atom n of model i
    integer, intent(out), optional :: natoms
      ! number of atoms per model in the PDB file
    integer, intent(out), optional :: nmodels
      ! number of models in the PDB file

    ! Local variables.
    integer :: iunit
      ! unit to use for file management
    type(pdbatom_t) :: atom
      ! temporary PDB atom container
    integer :: serial
      ! serial index of record, incremented after each ATOM and TER record written
    character(len=128) :: line
      ! an entire line from the PDB file
    integer :: lineNumber
      ! line number of PDB file
    integer :: iostat
      ! iostat after read
    integer :: currentModel
      ! index of current model
    integer :: firstAtomSerial
      ! first atom serial number in a model
    integer :: natomsThisModel
      ! number of atoms read for this model

    ! Open file for reading, ensuring it exists.
    call new_unit(iunit)
    open(unit=iunit, file=filename, status='OLD', err=11, position='REWIND', &
         form='FORMATTED', action='READ')
   
    ! Read and parse the file.
    lineNumber = 0
    currentModel = 0
    firstAtomSerial = -1
    natomsThisModel = 0
    do 
       ! Read a line.
       lineNumber = lineNumber + 1
       read(unit=iunit, fmt='(a)', iostat=iostat, end=20) line
       ! TODO: Terminate loop if EOF encountered?  Here or after parsing?
       ! if(iostat < 0) exit

       ! DEBUG
       if(debug) write(*,*) 'line ', lineNumber
              
       ! DEBUG
       if(debug .or. lineNumber>515000) write(*,'(I5,1X,A)') lineNumber, line

       ! Parse line, using first six characters to determine record type.
       select case(line(1:6))
       case ("MODEL ") ! New model.
          ! Increment model number.
          currentModel = currentModel + 1
          ! Reset the counter for the first atom in a model.
          firstAtomSerial = -1

       case ("ENDMDL") ! End of model.
          
       case ("ATOM  ") ! Parse atom record.
          ! Parse atom record.
          read(line, '(6x,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)') serial, atom%name, atom%altLoc, &
               atom%resName, atom%chainID, atom%resSeq, atom%icode, atom%x, atom%y, atom%z, atom%occupancy, atom%tempFactor, &
               atom%segID, atom%element, atom%charge          

          ! If no "MODEL...ENDMODEL" specified, assume this will be the only model.
          if(currentModel == 0) currentModel = 1

          ! Reset first atom serial number when first ATOM record in the model is read.
          if(firstAtomSerial == -1) firstAtomSerial = serial

          ! DEBUG: Write atom record.
          if(debug) then
             write(*, '(6x,a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)') 'ATOM  ', serial, atom%name, atom%altLoc, &
                  atom%resName, atom%chainID, atom%resSeq, atom%icode, atom%x, atom%y, atom%z, atom%occupancy, atom%tempFactor, &
                  atom%segID, atom%element, atom%charge          
          end if

          ! Determine atom number within this model
          natomsThisModel = serial - firstAtomSerial + 1
          ! Store that atom number in the atom's serial entry
          atom%serial = natomsThisModel          

          ! Store coordinates, if desired.
          if(present(models)) models(currentModel,natomsThisModel,:) = (/ atom%x, atom%y, atom%z /)

          ! Store atoms, if desired.
          if(present(atoms)) atoms(natomsThisModel) = atom

       case ("TER   ") ! Atom block termination record.
          
       case default
          ! Don't parse the line.
       end select


    end do

    ! End of file has been reached.
20  continue
    
    ! Close file.
    close(iunit)

    ! Determine number of atoms, if desired.
    if(present(natoms)) natoms = natomsThisModel

    ! Determine number of models, if desired.
    if(present(nmodels)) nmodels = currentModel

    if(debug) write(*,*) 'Done reading PDB file.'
    
    ! Return control.
    return

    ! Report error if file could not be opened.
11  continue
    write(*,*) 'pdbio: pdb_read: Error opening file ', filename
    stop
    
  end subroutine pdb_read

end module pdbio
