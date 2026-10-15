      MODULE CWSTWVCP
!
!svn $Id: waves_coupler.F 756 2008-09-14 20:18:28Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between SWAN   !
!  other coupled model(s) using the Model Coupling Toolkit (MCT).      !
!                                                                      !
!=======================================================================
!
!  Componenet model registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain decompositin descriptor datatype and assocoiated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field storage data types and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communitcations scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY : MCT_isend => isend
      USE m_Transfer, ONLY : MCT_irecv => irecv
      USE m_Transfer, ONLY : MCT_waitr => waitrecv
      USE m_Transfer, ONLY : MCT_waits => waitsend
!
      USE W3ODATMD, ONLY: SCREEN
      implicit none
!
!     PRIVATE
      PUBLIC :: INIT_WVCP
      PUBLIC :: initialize_wav_routers
      PUBLIC :: COAWST_CPL
      PUBLIC :: wav2ocn_coupling
      PUBLIC :: wavfocn_coupling
      PUBLIC :: wav2atm_coupling
      PUBLIC :: wavfatm_coupling
      PUBLIC :: finalize_wav_coupling
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapSWAN         ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)
      TYPE T_AttrVect_G
        TYPE(AttrVect) :: wav2ocn_AV            ! AttrVect variables
        TYPE(AttrVect) :: ocn2wav_AV
        TYPE(AttrVect) :: atm2wav_AV            ! AttrVec variables
        TYPE(AttrVect) :: wav2atm_AV            ! AttrVec variables
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)
      TYPE T_Router_O
        type(Router)   :: SWANtoROMS            ! Router variables
      END TYPE T_Router_O
      TYPE (T_Router_O), ALLOCATABLE :: Router_O(:,:)
      TYPE T_GSMapInterp_A
        TYPE(GlobalSegMap) :: GSMapWRF        ! GloabalSegMap variables
      END TYPE T_GSMapInterp_A
      TYPE (T_GSMapInterp_A), ALLOCATABLE :: GSMapInterp_A(:,:)
      TYPE T_Router_A
        type(Router)   :: SWANtoWRF           ! Router variables
      END TYPE T_Router_A
      TYPE (T_Router_A), ALLOCATABLE :: Router_A(:,:)
      CONTAINS
!***********************************************************************
!                                                                      *
      SUBROUTINE COAWST_CPL (first)
!                                                                      *
!***********************************************************************
!
      USE MCT_COUPLER_PARAMS
      USE W3GDATMD, ONLY: W3SETG
      USE WMMDATMD, ONLY: MDST, MDSE
      USE W3IDATMD, ONLY: W3SETI
      USE W3ADATMD, ONLY: W3SETA
      USE W3WDATMD, ONLY: W3SETW, VA
      USE W3UPDTMD, ONLY: W3ULEV, W3UICE, W3UIC1 !! cyy added
      USE W3IOGOMD, ONLY: W3OUTG_COAWST          !! cyy added
!
      IMPLICIT NONE
!
      INTEGER,  intent(in) :: first
      INTEGER   io, iw, ia, offset, run_couple
!
      LOGICAL   :: FLFRST
!
! cyy note: Re-organizing the sending/receiving order
! The updated order:
! CALL WAVFATM_COUPLING
! CALL WAVFOCN_COUPLING
! CALL W3ULEV
! CALL W3UICE
! CALL W3UIC1
! CALL W3OUTG_COAWST
! CALL WAV2ATM_COUPLING
! CALL WAV2OCN_COUPLING           
!
!     Call to get data from atm model.
      DO iw=1,Nwav_grids
        DO ia=1,Natm_grids
          run_couple=1
          IF (MOD(first, nWAVFATM(Natm_grids,1)).ne.0) run_couple=0
          IF (run_couple.eq.1) THEN
            CALL W3SETG ( iw, MDSE, MDST )
            CALL W3SETI ( iw, MDSE, MDST )
            CALL W3SETA ( iw, MDSE, MDST )
            CALL W3SETW ( iw, MDSE, MDST )
            CALL WAVFATM_COUPLING (iw, ia)
          ELSE
            GOTO 55
          END IF
        END DO
      END DO
  55  CONTINUE
!
!     Call to get data from ocn model.
      FLFRST = .TRUE.
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          run_couple=1
          IF (MOD(first, nWAVFOCN(Nwav_grids,1)).ne.0) run_couple=0
          IF (run_couple.eq.1) THEN
            CALL W3SETG ( iw, MDSE, MDST )
            CALL W3SETI ( iw, MDSE, MDST )
            CALL W3SETA ( iw, MDSE, MDST )
            CALL W3SETW ( iw, MDSE, MDST )
            CALL WAVFOCN_COUPLING (iw, io)
            CALL W3ULEV ( VA, VA )
            CALL W3UICE ( VA )
            CALL W3UIC1 ( FLFRST ) 
          ELSE
            GOTO 50
          END IF
        END DO
      END DO
  50  CONTINUE
!
!     Call the output routine to compute the mean properties.
      CALL W3OUTG_COAWST ( VA, .FALSE., .FALSE., .FALSE. )
!
!
!     Send data to atm model.
      DO iw=1,Nwav_grids
        DO ia=1,Natm_grids
          run_couple=1
          IF (MOD(first, nWAV2ATM(Natm_grids,1)).ne.0) run_couple=0
          IF (run_couple.eq.1) THEN
            CALL W3SETG ( iw, MDSE, MDST )
            CALL W3SETI ( iw, MDSE, MDST )
            CALL W3SETA ( iw, MDSE, MDST )
            CALL W3SETW ( iw, MDSE, MDST )
            CALL WAV2ATM_COUPLING (iw, ia)
          ELSE
            GOTO 45
          END IF
        END DO
      END DO
  45  CONTINUE
!
!     Send data to ocn model.
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          run_couple=1
          IF (MOD(first, nWAV2OCN(Nwav_grids,1)).ne.0) run_couple=0
          IF (run_couple.eq.1) THEN
            CALL W3SETG ( iw, MDSE, MDST )
            CALL W3SETI ( iw, MDSE, MDST )
            CALL W3SETA ( iw, MDSE, MDST )
            CALL W3SETW ( iw, MDSE, MDST )
            CALL WAV2OCN_COUPLING (iw, io)
          ELSE
            GOTO 40
          END IF
        END DO
      END DO
  40  CONTINUE
      RETURN
      END SUBROUTINE COAWST_CPL
      SUBROUTINE INIT_WVCP (ng)
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE MCT_COUPLER_PARAMS
      USE W3GDATMD, ONLY: NX, NY, NSEA, NK
!
      include 'mpif.h'
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: gsmsize, Nprocs
      integer :: i, j, io, ia, Isize, Jsize, Asize
      integer :: nRows, nCols, num_sparse_elems
      integer :: cid, cad
      integer :: IZ
      character (len=70)  :: nc_name
      character (len=20)  :: to_add
      character (len=1200) :: wostring
      character (len=1200) :: owstring
      character (len=5)   :: LSTCODE
      real :: cff
!     integer, dimension(2) :: src_grid_dims, dst_grid_dims
      integer, allocatable :: start(:), length(:)
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (WAV_COMM_WORLD, Nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nwav_grids))
        ALLOCATE(AttrVect_G(Nwav_grids))
      END IF
!
      WAVid=wavids(ng)
      IF (Nwav_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      WAV_COMM_WORLD,myids=wavids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      WAV_COMM_WORLD,WAVid)
      END IF
!
!  Initialize a Global Segment Map for non-haloed transfer of data for
!  SWAN. Determine non-haloed start and length arrays for this
!  processor. For now, this will set up a tiled exchange that is
!  not identical to thw ww3 tiling.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      IF (.not.allocated(start)) THEN
        allocate ( start(1) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(1) )
      END IF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
      gsmsize=length(1)
!
      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapSWAN, start,      &
     &                        length, 0, WAV_COMM_WORLD, WAVid)
      deallocate (start)
      deallocate (length)
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model.
!
      cad=LEN(wostring)
      DO i=1,cad
        wostring(i:i)=''
      END DO
      cid=1
!
      to_add='DISBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DISSURF'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DISWCAP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':HSIGN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':RTP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':ABOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':UBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DIRE'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
      to_add=':DIRN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLEN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLENP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':QB'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WDSPR'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WQP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':TAUOCE'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':TAUOCN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
! cyy added for wave radiation stress in CICE
      to_add=':TAUICEE'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':TAUICEN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      DO IZ=1,NK
        IF (IZ.LT.10) THEN
         WRITE(LSTCODE,"(A4,I1)") ":KS0",IZ
        ELSE
         WRITE(LSTCODE,"(A3,I2)") ":KS",IZ
        END IF
        to_add=LSTCODE
        cad=LEN_TRIM(to_add)
        write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
        cid=cid+cad
      END DO
      DO IZ=1,NK
        IF (IZ.LT.10) THEN
         WRITE(LSTCODE,"(A4,I1)") ":US0",IZ
        ELSE
         WRITE(LSTCODE,"(A3,I2)") ":US",IZ
        END IF
        to_add=LSTCODE
        cad=LEN_TRIM(to_add)
        write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
        cid=cid+cad
      END DO
      DO IZ=1,NK
        IF (IZ.LT.10) THEN
         WRITE(LSTCODE,"(A4,I1)") ":VS0",IZ
        ELSE
         WRITE(LSTCODE,"(A3,I2)") ":VS",IZ
        END IF
        to_add=LSTCODE
        cad=LEN_TRIM(to_add)
        write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
        cid=cid+cad
      END DO
!
!  Finalize and remove trailing spaces from the wostring
!  for the rlist.
!
      cad=LEN_TRIM(wostring)
      wostring=wostring(1:cad)
!
      CALL AttrVect_init(AttrVect_G(ng)%wav2ocn_AV,                     &
     &                   rList=TRIM(wostring),lsize=gsmsize)
      CALL AttrVect_zero(AttrVect_G(ng)%wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      cad=LEN(owstring)
      DO i=1,cad
        owstring(i:i)=''
      END DO
      cid=1
!
      to_add='DEPTH'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLEV'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELX'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELY'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':ZO'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':SIC'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':SIT'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad      
!
!  Finalize and remove trailing spaces from the owstring
!  for the rlist.
!
      cad=LEN_TRIM(owstring)
      owstring=owstring(1:cad)
!
      CALL AttrVect_init (AttrVect_G(ng)%ocn2wav_AV,                    &
     &                    rList=TRIM(owstring),lsize=gsmsize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2wav_AV)
!
!
!  Initialize attribute vector holding the export data code string of
!  the atmosphere model.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,            &
     &                         WAV_COMM_WORLD)
      CALL AttrVect_init (AttrVect_G(ng)%atm2wav_AV, rlist="U10:V10",   &
     &                    lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%atm2wav_AV)
!
! Initialize atribute vector holding wave data to atm.
!
      CALL AttrVect_init(AttrVect_G(ng)%wav2atm_AV,                     &
     &                   rList="HSIGN:WLENP:RTP",                       &
     &                   lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%wav2atm_AV)
!
      RETURN
      END SUBROUTINE INIT_WVCP
      SUBROUTINE INITIALIZE_WAV_ROUTERS
!
!=======================================================================
!                                                                      !
!  Initialize waves routers for wave model.                            !
!                                                                      !
!=======================================================================
!
      USE MCT_COUPLER_PARAMS
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: ng, iw, ia
!
!  Initialize MCT Routers.
!
      ALLOCATE(Router_O(Nwav_grids,Nocn_grids))
!
!  Initialize a router to the ocean model component.
!
      DO ng=1,Nocn_grids
        DO iw=1,Nwav_grids
          OCNid=ocnids(ng)
          CALL Router_init (OCNid, GlobalSegMap_G(iw)%GSMapSWAN,        &
     &                      WAV_COMM_WORLD, Router_O(iw,ng)%SWANtoROMS)
        END DO
      END DO
      ALLOCATE(Router_A(Nwav_grids,Natm_grids))
!
!  Initialize a router to the atmosphere model component.
!
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          ATMid=atmids(ia)
          CALL Router_init (ATMid, GlobalSegMap_G(iw)%GSMapSWAN,        &
     &                      WAV_COMM_WORLD, Router_A(iw,ia)%SWANtoWRF)
        END DO
      END DO
      RETURN
      END SUBROUTINE INITIALIZE_WAV_ROUTERS
      SUBROUTINE WAV2OCN_COUPLING (iw, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the OCEAN model:                                 !
!                                                                      !
!     * Wave direction (degrees)                                       !
!     * Significant wave height (m)                                    !
!     * Average wave length (m)                                        !
!     * Surface wave relative peak period (s)                          !
!     * Bottom wave period (s)                                         !
!     * Percent of breakig waves (nondimensional)                      !
!     * Wave energy dissipation (W/m2)                                 !
!     * Wave bottom orbital velocity (m/s)                             !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Bathymetry, bottom elevation (m)                               !
!     * Free-surface, water surface elevation (m)                      !
!     * Depth integrated u-momentum (m/s)                              !
!     * Depth integrated v-momentum (m/s)                              !
!                                                                      !
!=======================================================================
!
      USE MCT_COUPLER_PARAMS
      USE CONSTANTS, ONLY: PI
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF, NK, NTH
      USE W3ADATMD, ONLY: HS, PHIBBL, PHIOC, FP0, T0M1, UBA, ABA,       &
     &                    STK_STOKES, STU_STOKES, STV_STOKES,           &
     &                    THM, WLM, WLP, WBT, THS, QP,                  &
     &                    PHIBRKX, PHIBRKY, TAUOSX, TAUOSY,             &
     &                    PHICAPX, PHICAPY
      USE W3ADATMD, ONLY: TAUICE
      USE W3ADATMD, ONLY: QB
      USE W3WDATMD, ONLY: VA, UST, USTDIR, RHOAIR
      USE W3IOGOMD
!
      implicit none
      include 'mpif.h'
!
!  Imported variable declarations.
!
      integer :: Numcouple, iw, io
      integer :: IP, ISEA, IX, IY
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: i, id, j, gsmsize, ierr, indx, Tag, JSEA
      integer :: Istr, Iend, Jstr, Jend
      integer :: start, length, grdsize
      integer :: Isize, Jsize, Nprocs, ipp
      integer, pointer :: points(:)
      character (len=4) :: SCODE
      integer :: IZ
      real, pointer :: SND_BUF(:), RCV_BUF(:)
      real(m8) :: fac, cff
      real(m8), pointer :: avdata(:)
      real(m8), pointer :: DIRE(:)
      real(m8), pointer :: DIRN(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, Nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      allocate ( DIRE(gsmsize),stat=ierr )
      allocate ( DIRN(gsmsize),stat=ierr )
      avdata=0.0_m8
      DIRE=0.0_m8
      DIRN=0.0_m8
!
      grdsize=NX*NY
      allocate ( SND_BUF(grdsize),stat=ierr )
      SND_BUF=0.0
      allocate ( RCV_BUF(grdsize),stat=ierr )
      RCV_BUF=0.0
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           MyRank, points)
!
!  Determine grid tiling for exchanges.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
!
!  Call the output routine to compute the mean properties.
!  cyy note: masked the calling W3OUTG_COAWST in WAV2OCN_COUPLING
!      CALL W3OUTG_COAWST ( VA, .FALSE., .FALSE., .FALSE. )
!
!
!  Load WW3 data into MCT storage buffers.
!  The data is exported using ROMS definition for real kind m8=r8.
!
!-------------------------------------------------------------------
!  DISBOT: Dissipation Bottom
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=PHIBBL(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "DISBOT",avdata)
!-------------------------------------------------------------------
!  DISSURF: Dissipation surface
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        cff=PHIBRKX(i)**2+PHIBRKY(i)**2
        IF (cff.ge.10000.0) cff=0.
        SND_BUF(IP)=SQRT(cff)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "DISSURF",avdata)
!
!  DISWCAP: Dissipation white capping
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        cff=PHICAPX(i)**2+PHICAPY(i)**2
        IF (cff.ge.10000.0) cff=0.
        SND_BUF(IP)=SQRT(cff)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "DISWCAP",avdata)
!-------------------------------------------------------------------
!  HS: Signfiicant wave height
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=HS(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "HSIGN",avdata)
!-------------------------------------------------------------------
!  RTP: Peak surface period
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=FP0(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
!     fac=2.0_m8*PI
      fac=1.0_m8
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=fac/MAX(REAL(RCV_BUF(i),m8),0.001_m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "RTP",avdata)
!-------------------------------------------------------------------
!  ABOT: Mean bottom period
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=ABA(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "ABOT",avdata)
!-------------------------------------------------------------------
!  UBOT: bottom orbitral velocity
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=UBA(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "UBOT",avdata)
!-------------------------------------------------------------------
!  DIRE and DIRN: mean direction, break into 2 components
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0_m8
      RCV_BUF=0.0_m8
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=THM(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
!     fac=PI/180.0
      fac=1.0_m8
      DO IP=1,gsmsize
        cff=450.0_m8*PI/180.0_m8-(avdata(IP)+PI)
        DIRE(IP)=1.0_m8*SIN(cff*fac)
        DIRN(IP)=1.0_m8*COS(cff*fac)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                           "DIRE",DIRE)
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                           "DIRN",DIRN)
!-------------------------------------------------------------------
!  NEED DIRE and DIRN: Peak directions, break into 2 components
!-------------------------------------------------------------------
!  WLEN: mean wave length
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=WLM(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "WLEN",avdata)
!-------------------------------------------------------------------
!  WLENP: peak wave length
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=WLP(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "WLENP",avdata)
!-------------------------------------------------------------------
!  QB: percent wave break
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
!       SND_BUF(IP)=WBT(i)
        SND_BUF(IP)=QB(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "QB",avdata)
!-------------------------------------------------------------------
!  WDSPR: wave dir spread
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=THS(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "WDSPR",avdata)
!-------------------------------------------------------------------
!  WQP: wave spectrum peakedness
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=QP(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                             "WQP",avdata)
!-------------------------------------------------------------------
!  Tauocx:  Stress in east direction
!  cyy note: this is COAWST added variable, which is already in N/m2
!            and this is only for guard BULK_FLUXES && WAV2OCN_FLUXES
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!  Notice Tauocx and Tauocy are indixed ISEA.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        ISEA=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IP=(IY-1)*NX+IX
        cff=TAUOSX(i)
        IF (cff.le.-900.0) cff=0.
        SND_BUF(IP)=cff
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                           "TAUOCE",avdata)
!-------------------------------------------------------------------
!  Tauocy: Ustar in north direction
!  cyy note: this is COAWST added variable, which is already in N/m2
!            and this is only for guard BULK_FLUXES && WAV2OCN_FLUXES
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!  Notice Tauocx and Tauocy are indixed ISEA.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        ISEA=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IP=(IY-1)*NX+IX
        cff=TAUOSY(i)
        IF (cff.le.-900.0) cff=0.
        SND_BUF(IP)=cff
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                           "TAUOCN",avdata)
!-------------------------------------------------------------------
!  Tauicex:  Wave stress on sea ice in east direction (m2/s2)
!  cyy note : added for CICE
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!  Notice Tauicex and Tauicey are indixed ISEA.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        ISEA=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IP=(IY-1)*NX+IX
        cff=TAUICE(i,1)
        !IF (cff.le.-900.0) cff=0. ! cyy note: may not need this
        SND_BUF(IP)=cff
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                           "TAUICEE",avdata)
!-------------------------------------------------------------------
!  Tauicey:  Wave stress on sea ice in north direction (m2/s2)
!  cyy note: added for CICE
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!  Notice Tauocx and Tauocy are indixed ISEA.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        ISEA=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IP=(IY-1)*NX+IX
        cff=TAUICE(i,2)
        !IF (cff.le.-900.0) cff=0. ! cyy note: may not need this
        SND_BUF(IP)=cff
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,             &
     &                           "TAUICEN",avdata)
!-------------------------------------------------------------------
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
!  STK
      DO IZ=1,NK
        IF (IZ.LT.10) THEN
          WRITE(SCODE,"(A3,I1)") "KS0",IZ
        ELSE
          WRITE(SCODE,"(A2,I2)") "KS", IZ
        END IF
!
        SND_BUF=0.0
        DO i=1,NSEAL
          IP=(MyRank+1)+(i-1)*Nprocs
          IX     = MAPSF(IP,1)
          IY     = MAPSF(IP,2)
          IP=(IY-1)*NX+IX
          SND_BUF(IP)=STK_STOKES(i,IZ)
        END DO
!
!  Gather up all the data.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          avdata(IP)=REAL(RCV_BUF(i),m8)
        END DO
!
        CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,           &
     &                             SCODE,avdata)
      END DO
!
!  STU
      DO IZ=1,NK
        IF (IZ.LT.10) THEN
          WRITE(SCODE,"(A3,I1)") "US0",IZ
        ELSE
          WRITE(SCODE,"(A2,I2)") "US", IZ
        END IF
!
        SND_BUF=0.0
        DO i=1,NSEAL
          IP=(MyRank+1)+(i-1)*Nprocs
          IX     = MAPSF(IP,1)
          IY     = MAPSF(IP,2)
          IP=(IY-1)*NX+IX
          SND_BUF(IP)=STU_STOKES(i,IZ)
        END DO
!
!  Gather up all the data.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          avdata(IP)=REAL(RCV_BUF(i),m8)
        END DO
!
        CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,           &
     &                             SCODE,avdata)
      END DO
!
!  STV
      DO IZ=1,NK
        IF (IZ.LT.10) THEN
          WRITE(SCODE,"(A3,I1)") "VS0",IZ
        ELSE
          WRITE(SCODE,"(A2,I2)") "VS", IZ
        END IF
!
        SND_BUF=0.0
        DO i=1,NSEAL
          IP=(MyRank+1)+(i-1)*Nprocs
          IX     = MAPSF(IP,1)
          IY     = MAPSF(IP,2)
          IP=(IY-1)*NX+IX
          SND_BUF(IP)=STV_STOKES(i,IZ)
        END DO
!
!  Gather up all the data.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          avdata(IP)=REAL(RCV_BUF(i),m8)
        END DO
!
        CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2ocn_AV,           &
     &                             SCODE,avdata)
      END DO
!
!
!-----------------------------------------------------------------------
!  Send wave fields bundle to ocean model, ROMS.
!-----------------------------------------------------------------------
!
          Tag=io*100+0*10+iw
          CALL MCT_isend (AttrVect_G(iw)%wav2ocn_AV,                    &
     &                    Router_O(iw,io)%SWANtoROMS, Tag)
          CALL MCT_waits (Router_O(iw,io)%SWANtoROMS)
          IF (MyRank.EQ.0) THEN
            WRITE (SCREEN,36)' == WW3 grid ',iw,                        &
     &                       ' sent wave data to ROMS grid ', io
 36         FORMAT (a14,i2,a29,i2)
          END IF
          IF (MyError.ne.0) THEN
            WRITE (*,*)'coupling send fail swancplr, Error= ', MyError
            CALL FINALIZE_WAV_COUPLING(iw)
          END IF
      deallocate (avdata, points, DIRE, DIRN)
      deallocate (SND_BUF, RCV_BUF)
!
      RETURN
      END SUBROUTINE WAV2OCN_COUPLING
      SUBROUTINE WAV2ATM_COUPLING (iw, ia)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  atm and wave models. Currently, the following data streams are      !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the ATM model:                                   !
!                                                                      !
!     * Significant wave height (m)                                    !
!     * Surface wave relative peak period (s)                          !
!     * Surface wave length (m)                                        !
!                                                                      !
!  Fields imported from the ATM Model:                                 !
!                                                                      !
!     * Wind Speed (m/s)                                               !
!                                                                      !
!=======================================================================
!
      USE CONSTANTS, ONLY: PI
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF
      USE MCT_COUPLER_PARAMS
      USE W3ADATMD, ONLY: HS, PHIBBL, PHIOC, FP0, T0M1, UBA
      USE W3ADATMD, ONLY: WLM
!
      implicit none
      include 'mpif.h'
!
!  Imported variable declarations.
!
      integer :: Numcouple, iw, ia
      integer :: IP, IX, IY
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend, Asize
      integer :: Isize, Jsize, INDXG, Nprocs, OFFSET
      integer :: start, length, grdsize
      integer, pointer :: points(:)
      real, pointer :: SND_BUF(:), RCV_BUF(:)
      real(m8) :: fac
      real(m8), pointer :: avdata(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to WRF.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, Nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
!
      grdsize=NX*NY
      allocate ( SND_BUF(grdsize),stat=ierr )
      SND_BUF=0.0
      allocate ( RCV_BUF(grdsize),stat=ierr )
      RCV_BUF=0.0
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           MyRank, points)
!
!  Determine grid tiling for exchanges.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
!
!  Load WW3 data into MCT storage buffers.
!  The data is exported using WRF definition for real kind m8=r8.
!-------------------------------------------------------------------
!  HS: Signfiicant wave height
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=HS(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2atm_AV,             &
     &                             "HSIGN",avdata)
!-------------------------------------------------------------------
!  RTP: Peak surface period
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=FP0(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
!     fac=2.0_m8*PI cyy note: FP0 is in Hz, fac should be changed.
      fac=1.0_m8
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=fac/MAX(REAL(RCV_BUF(i),m8),0.001_m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2atm_AV,             &
     &                             "RTP",avdata)
!-------------------------------------------------------------------
!  WLEN: mean wave length
!
!  Fill wet parts of array SND_BUF that is NXxNY length.
!  The local variable is only 1:NSEAL(M) long.
!
      SND_BUF=0.0
      DO i=1,NSEAL
        IP=(MyRank+1)+(i-1)*Nprocs
        IX     = MAPSF(IP,1)
        IY     = MAPSF(IP,2)
        IP=(IY-1)*NX+IX
        SND_BUF(IP)=WLM(i)
      END DO
!
!  Gather up all the data.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Now extract the section of data from this tile
!  and fill the mct array.
!
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        avdata(IP)=REAL(RCV_BUF(i),m8)
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(iw)%wav2atm_AV,             &
     &                             "WLENP",avdata)
!-------------------------------------------------------------------
!  Send wave fields bundle to atm model, WRF.
!-----------------------------------------------------------------------
!
 35     FORMAT (a14,i2,a23,i5)
!
!  Send fields to atmosphere model.
!
        Tag=ia*10+iw
        CALL MCT_isend (AttrVect_G(iw)%wav2atm_AV,                      &
     &                  Router_A(iw,ia)%SWANtoWRF, Tag)
        CALL MCT_waits (Router_A(iw,ia)%SWANtoWRF)
        IF (MyRank.EQ.0) THEN
          WRITE (SCREEN,35) '== WW3 grid ',iw,' sent data to WRF grid ' &
     &                      ,ia
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupl fail swancplr, MyStatus= ', MyError
          CALL FINALIZE_WAV_COUPLING(iw)
        END IF
      deallocate (avdata)
      deallocate (SND_BUF, RCV_BUF)
!
      RETURN
      END SUBROUTINE WAV2ATM_COUPLING
      SUBROUTINE WAVFOCN_COUPLING (ng, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the OCEAN model:                                 !
!                                                                      !
!     * Wave direction (degrees)                                       !
!     * Significant wave height (m)                                    !
!     * Average wave length (m)                                        !
!     * Surface wave relative peak period (s)                          !
!     * Bottom wave period (s)                                         !
!     * Percent of breakig waves (nondimensional)                      !
!     * Wave energy dissipation (W/m2)                                 !
!     * Wave bottom orbital velocity (m/s)                             !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Bathymetry, bottom elevation (m)                               !
!     * Free-surface, water surface elevation (m)                      !
!     * Depth integrated u-momentum (m/s)                              !
!     * Depth integrated v-momentum (m/s)                              !
!                                                                      !
!=======================================================================
!
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF, NK, NTH
      USE W3SERVMD
      USE W3IDATMD
      USE MCT_COUPLER_PARAMS
!     USE W3WDATMD, ONLY: WLV
!     USE W3IDATMD, ONLY: WLEV
      USE W3ADATMD, ONLY: CX, CY
!
      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
      integer :: ng, io
!
!  Local variable declarations.
!
      integer :: IP, IX, IY, grdsize
      integer :: MyStatus, MyError, MySize, MyRank, Nprocs
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend, start, length
      integer :: Isize, Jsize, INDXG, OFFSET, IAPROC
!  cyy added: exchange timer
      integer :: timer_error, timer_nprocs
      double precision :: timer_start, timer_end, wait_local
      double precision :: wait_min, wait_max, wait_sum, wait_mean
      !double precision, external :: MPI_WTIME
!
      real :: cff, cffmin, cffmax
      real :: cff2
      real, parameter              :: Large = 1.0E+20
      real, dimension(2)           :: range
      real, pointer                :: SND_BUF(:), RCV_BUF(:)
      real(m8), pointer            :: avdata(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, Nprocs, MyError)
      IAPROC=MyRank+1
      grdsize=NX*NY
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
      allocate ( SND_BUF(grdsize),stat=ierr )
      SND_BUF=0.
      allocate ( RCV_BUF(grdsize),stat=ierr )
      RCV_BUF=0.
!
!-----------------------------------------------------------------------
!  Receive from ROMS: Depth, Water Level, VELX, and VELY.
!-----------------------------------------------------------------------
!
!  Schedule receiving field from ocean model.
!
      Tag=io*100+0*10+ng
      ! cyy added: exchange timer
      timer_start=MPI_WTIME()
      CALL MCT_irecv (AttrVect_G(ng)%ocn2wav_AV,                        &
     &                Router_O(ng,io)%SWANtoROMS, Tag)
!
!     Wait to make sure the OCN data has arrived.
!
      CALL MCT_waitr (AttrVect_G(ng)%ocn2wav_AV,                        &
     &                 Router_O(ng,io)%SWANtoROMS)
      ! cyy added: exchange timer
      timer_end=MPI_WTIME()
      wait_local=timer_end-timer_start
!      CALL MPI_REDUCE (wait_local, wait_min, 1,                         &
!     &                 MPI_DOUBLE_PRECISION, MPI_MIN, 0,                &
!     &                 WAV_COMM_WORLD, timer_error)
!      CALL MPI_REDUCE (wait_local, wait_max, 1,                         &
!     &                 MPI_DOUBLE_PRECISION, MPI_MAX, 0,                &
!     &                 WAV_COMM_WORLD, timer_error)
      CALL MPI_REDUCE (wait_local, wait_sum, 1,                         &
     &                 MPI_DOUBLE_PRECISION, MPI_SUM, 0,                &
     &                 WAV_COMM_WORLD, timer_error)
!
      IF (MyRank.EQ.0) THEN
        WRITE (SCREEN,35) ' == WW3 grid ',ng,                           &
     &                    ' recv data from ROMS grid ', io
        ! cyy added: exchange timer
        wait_mean=wait_sum/DBLE(Nprocs)
        WRITE (SCREEN,'(A,1X,F12.6,A)')                                 &
     &  ' Exchange timer [WW3 <- ROMS] wait mean:', wait_mean, ' s'
      END IF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling fail ww3cplr, MyStatus= ', MyError
        CALL FINALIZE_WAV_COUPLING(ng)
      END IF
 35   FORMAT (a14,i2,a26,i2)
!
! Compute local non-halo data size.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
!
 40   FORMAT (a36,1x,2(1pe14.6))
!
!  Bottom elevation....................................................
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "DEPTH",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoWW3 Min/Max DEPTH   (m):     ',       &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array Depth.
!
        IF (MyRank.eq.0) THEN
          write(SCREEN,*) 'Depth from ROMS to WW3 not coded yet.'
        END IF
!
!  Water surface elevation............................................
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "WLEV",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoWW3 Min/Max WLEV    (m):     ',       &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array WLEV.
!
        DO i=1,NSEA
          IX     = MAPSF(i,1)
          IY     = MAPSF(i,2)
          IP=(IY-1)*NX+IX
          IF (io.eq.1) THEN
            WLEV(IX,IY)=RCV_BUF(IP)
          ELSE
            WLEV(IX,IY)=WLEV(IX,IY)+RCV_BUF(IP)
          END IF
        END DO
!
!  Depth-integrated u-velocity........................................
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "VELX",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) ' ROMStoWW3 Min/Max VELX    (ms-1):    ',    &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array VELX.
!
       DO i=1,NSEA
         IX     = MAPSF(i,1)
         IY     = MAPSF(i,2)
         IP=(IY-1)*NX+IX
         IF (io.eq.1) THEN
           CX(i)=RCV_BUF(IP)
         ELSE
           CX(i)=CX(i)+RCV_BUF(IP)
         END IF
       END DO
!
!  Depth-integrated v-velocity........................................
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "VELY",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) ' ROMStoWW3 Min/Max VELY    (ms-1):    ',    &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array VELY.
!
       DO i=1,NSEA
         IX     = MAPSF(i,1)
         IY     = MAPSF(i,2)
         IP=(IY-1)*NX+IX
         IF (io.eq.1) THEN
           CY(i)=RCV_BUF(IP)
         ELSE
           CY(i)=CY(i)+RCV_BUF(IP)
         END IF
       END DO
!
!  Bottom roughness...................................................
!
        CALL AttrVect_exportRAttr(AttrVect_G(ng)%ocn2wav_AV,            &
     &                            "ZO",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoWW3 Min/Max BottZ0  (m):     ',       &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array Z0.
!
        IF (MyRank.eq.0) THEN
          write(SCREEN,*) 'Z0 from ROMS to WW3 not coded yet.'
        END IF
!
!  Sea ice concentration..............................................
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "SIC",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoWW3 Min/Max SIC     (-):     ',       &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array ICEI.
!
        DO i=1,NSEA
          IX     = MAPSF(i,1)
          IY     = MAPSF(i,2)
          IP=(IY-1)*NX+IX
          IF (io.eq.1) THEN
            ICEI(IX,IY)=RCV_BUF(IP)
          ELSE
            ICEI(IX,IY)=ICEI(IX,IY)+RCV_BUF(IP)
          END IF
        END DO
!
!  Sea ice thickness.................................................
!
        CALL AttrVect_exportRAttr (AttrVect_G(ng)%ocn2wav_AV,           &
     &                             "SIT",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        SND_BUF=0.0
        IP=0
        DO i=start,start+length-1
          IP=IP+1
          range(1)=MIN(range(1),REAL(avdata(IP)))
          range(2)=MAX(range(2),REAL(avdata(IP)))
          SND_BUF(i)=REAL(avdata(IP))
        END DO
        CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,               &
                           MPI_MIN, WAV_COMM_WORLD, MyError)
        CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,               &
                           MPI_MAX, WAV_COMM_WORLD, MyError)
        IF (MyRank.eq.0) THEN
          write(SCREEN,40) 'ROMStoWW3 Min/Max SIT     (m):     ',       &
     &                      cffmin, cffmax
        END IF
!
!  now scatter data to all nodes.
!
        CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                   &
     &                     MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array ICEP1.
!
        DO i=1,NSEA
          IX     = MAPSF(i,1)
          IY     = MAPSF(i,2)
          IP=(IY-1)*NX+IX
          IF (io.eq.1) THEN
            ICEP1(IX,IY)=RCV_BUF(IP)
          ELSE
            ICEP1(IX,IY)=ICEP1(IX,IY)+RCV_BUF(IP)
          END IF
        END DO
!
      deallocate (avdata, SND_BUF, RCV_BUF)
!
      RETURN
      END SUBROUTINE WAVFOCN_COUPLING
      SUBROUTINE WAVFATM_COUPLING (iw, ia)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!                                                                      !
!  Fields imported from the ATM Model:                                 !
!                                                                      !
!     * Wind E and N (m/s)                                             !
!                                                                      !
!=======================================================================
!
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF
      USE W3SERVMD
      USE W3IDATMD
      USE MCT_COUPLER_PARAMS
      USE W3ADATMD, ONLY: UA, UD
!
      implicit none
      include 'mpif.h'
!
!  Imported variable declarations.
!
      integer :: iw, ia
!
!  Local variable declarations.
!
      integer :: IP, IX, IY, grdsize
      integer :: MyStatus, MyError, MySize, MyRank, Nprocs
      integer :: i, id, j, gsmsize, ierr, indx, Tag
      integer :: Istr, Iend, Jstr, Jend, start, length
      integer :: Isize, Jsize, INDXG, OFFSET, IAPROC
      real :: cff, cffmin, cffmax
      real, parameter    ::  Large = 1.0E+20
      real, dimension(2) :: range
      real, pointer      :: SND_BUF(:), RCV_BUF(:)
      real, pointer      :: WND_U10(:), WND_V10(:)
      real(m8), pointer  :: avdata(:)
!  cyy added: exchange timer
      integer :: timer_error, timer_nprocs
      double precision :: timer_start, timer_end, wait_local
      double precision :: wait_min, wait_max, wait_sum, wait_mean
      !double precision, external :: MPI_WTIME
!
!-----------------------------------------------------------------------
!  Get wind data from atm.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, Nprocs, MyError)
      IAPROC=MyRank+1
      grdsize=NX*NY
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(iw)%GSMapSWAN,          &
     &                           WAV_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
      MyError=0
      allocate ( SND_BUF(grdsize),stat=ierr )
      SND_BUF=0.
      allocate ( RCV_BUF(grdsize),stat=ierr )
      RCV_BUF=0.
      allocate ( WND_U10(grdsize),stat=ierr )
      WND_U10=0.
      allocate ( WND_V10(grdsize),stat=ierr )
      WND_V10=0.
!
!-----------------------------------------------------------------------
!  RCV the data from WRF.
!-----------------------------------------------------------------------
!
 35   FORMAT (a14,i2,a24,i2)
!
!  Receive fields from atmosphere model.
!
      Tag=0*100+ia*10+iw
      ! cyy added: exchange timer
      timer_start=MPI_WTIME()
      CALL MCT_irecv (AttrVect_G(iw)%atm2wav_AV,                        &
     &                Router_A(iw,ia)%SWANtoWRF, Tag)
!     Wait to make sure the WRF data has arrived.
      CALL MCT_waitr (AttrVect_G(iw)%atm2wav_AV,                        &
     &                Router_A(iw,ia)%SWANtoWRF)
      ! cyy added: exchange timer
      timer_end=MPI_WTIME()
      wait_local=timer_end-timer_start
!      CALL MPI_REDUCE (wait_local, wait_min, 1,                         &
!     &                 MPI_DOUBLE_PRECISION, MPI_MIN, 0,                &
!     &                 WAV_COMM_WORLD, timer_error)
!      CALL MPI_REDUCE (wait_local, wait_max, 1,                         &
!     &                 MPI_DOUBLE_PRECISION, MPI_MAX, 0,                &
!     &                 WAV_COMM_WORLD, timer_error)
      CALL MPI_REDUCE (wait_local, wait_sum, 1,                         &
     &                 MPI_DOUBLE_PRECISION, MPI_SUM, 0,                &
     &                 WAV_COMM_WORLD, timer_error)
      IF (MyRank.EQ.0) THEN
        WRITE (SCREEN,35) ' == WW3 grid ',iw,                           &
     &                    ' recv data from WRF grid', ia
	! cyy added: exchange timer
        wait_mean=wait_sum/DBLE(Nprocs)
        WRITE (SCREEN,'(A,1X,F12.6,A)')                                 &
     &  ' Exchange timer [WW3 <- WRF] wait mean:', wait_mean, ' s'
      END IF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling fail swancplr, MyStatus= ', MyError
        CALL FINALIZE_WAV_COUPLING(iw)
      END IF
!
! Compute local non-halo data size.
!
      Jsize=INT(NY/Nprocs)
      IF (MyRank.eq.Nprocs-1) THEN
        Jsize=NY-Jsize*(Nprocs-1)
      ENDIF
      start=(MyRank*INT(NY/Nprocs))*NX+1
      length=Jsize*NX
!
!
!  U-wind.
!
      CALL AttrVect_exportRAttr(AttrVect_G(iw)%atm2wav_AV,"U10",        &

     &                          avdata,gsmsize)
      range(1)= Large
      range(2)=-Large
      SND_BUF=0.0
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        range(1)=MIN(range(1),REAL(avdata(IP)))
        range(2)=MAX(range(2),REAL(avdata(IP)))
        SND_BUF(i)=REAL(avdata(IP))
      END DO
      CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,                 &
                         MPI_MIN, WAV_COMM_WORLD, MyError)
      CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,                 &
                         MPI_MAX, WAV_COMM_WORLD, MyError)
      IF (MyRank.eq.0) THEN
        write(SCREEN,40) ' WRFtoWW3 Min/Max U10     (ms-1):    ',       &
     &                      cffmin, cffmax
      END IF
!
!  now scatter data to all nodes.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array WND_U10 as temporary for now.
!
      DO i=1,grdsize
        IF (ia.eq.1) THEN
          WND_U10(i)=RCV_BUF(i)
        ELSE
          WND_U10(i)=WND_U10(i)+RCV_BUF(i)
        END IF
      END DO
!     DO i=1,NSEA
!       IX     = MAPSF(i,1)
!       IY     = MAPSF(i,2)
!       IP=(IY-1)*NX+IX
!       WND_SPD(i)=RCV_BUF(IP)
!     END DO
!
!  V-wind.
!
      CALL AttrVect_exportRAttr(AttrVect_G(iw)%atm2wav_AV,"V10",        &
     &                          avdata,gsmsize)
      range(1)= Large
      range(2)=-Large
      SND_BUF=0.0
      IP=0
      DO i=start,start+length-1
        IP=IP+1
        range(1)=MIN(range(1),REAL(avdata(IP)))
        range(2)=MAX(range(2),REAL(avdata(IP)))
       SND_BUF(i)=REAL(avdata(IP))
      END DO
      CALL MPI_ALLREDUCE(range(1), cffmin, 1, MPI_REAL,                 &
                         MPI_MIN, WAV_COMM_WORLD, MyError)
      CALL MPI_ALLREDUCE(range(2), cffmax, 1, MPI_REAL,                 &
                         MPI_MAX, WAV_COMM_WORLD, MyError)
      IF (MyRank.eq.0) THEN
        write(SCREEN,40) ' WRFtoWW3 Min/Max V10     (ms-1):    ',       &
     &                    cffmin, cffmax
      END IF
!
!  now scatter data to all nodes.
!
      CALL MPI_ALLREDUCE(SND_BUF, RCV_BUF, grdsize,                     &
     &                   MPI_REAL, MPI_SUM, WAV_COMM_WORLD, MyError)
!
!  Scatter to array WND_V10 as temporary for now.
!
      DO i=1,grdsize
        IF (ia.eq.1) THEN
          WND_V10(i)=RCV_BUF(i)
        ELSE
          WND_V10(i)=WND_V10(i)+RCV_BUF(i)
        END IF
      END DO
!
!  Now we need to combine wnd speed and dir and scatter.
!
      DO i=1,NSEA
        IX     = MAPSF(i,1)
        IY     = MAPSF(i,2)
        IP=(IY-1)*NX+IX
        IF (ia.eq.1) THEN
          UA(i)=SQRT(WND_U10(IP)**2+WND_V10(IP)**2+0.000001)
          UD(i)=ATAN2(WND_V10(IP),WND_U10(IP))
        ELSE
          UA(i)=UA(i)+SQRT(WND_U10(IP)**2+WND_V10(IP)**2+0.000001)
          UD(i)=UD(i)+ATAN2(WND_V10(IP),WND_U10(IP))
        END IF
      END DO
 40     FORMAT (a36,1x,2(1pe14.6))
!
      deallocate (avdata)
      deallocate (SND_BUF, RCV_BUF)
      deallocate (WND_U10, WND_V10)
!
      RETURN
      END SUBROUTINE WAVFATM_COUPLING
      SUBROUTINE FINALIZE_WAV_COUPLING(iw)
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
      USE MCT_COUPLER_PARAMS
!
!  Local variable declarations.
!
      integer :: iw, io, ia, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      DO io=1,Nocn_grids
        CALL Router_clean (Router_O(iw,io)%SWANtoROMS, MyError)
      END DO
      CALL AttrVect_clean (AttrVect_G(iw)%wav2ocn_AV, MyError)
      DO ia=1,Natm_grids
        CALL Router_clean (Router_A(iw,ia)%SWANtoWRF, MyError)
      END DO
      CALL AttrVect_clean (AttrVect_G(iw)%atm2wav_AV, MyError)
      CALL GlobalSegMap_clean (GlobalSegMap_G(iw)%GSMapSWAN, MyError)
      END SUBROUTINE FINALIZE_WAV_COUPLING
      END MODULE CWSTWVCP
