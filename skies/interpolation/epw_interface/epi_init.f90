!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                  !!!
!!!           Module for initializing EPI EPW facilities.            !!!
!!!           Major part of the code in this file is                 !!!
!!!           copied directly from EPW software                      !!!
!!!           (e.g. ephwann_shuffle.f90 and other src files)         !!!
!!!           to set up all the required global variables.           !!!
!!!                                                                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE epi_initialize
    !
    USE kinds,              ONLY : DP
    USE constants_epw,      ONLY : zero, czero, eps8
    !
    USE mp,                 ONLY : mp_bcast
    USE mp_global,          ONLY : mp_startup
    USE mp_pools,           ONLY : npool
    USE mp_images,          ONLY : nimage
    USE parallel_include
    USE mp_world,           ONLY : mpime, mp_world_start, world_comm, mp_world_end
    USE io_epw,             ONLY : epw_read
    USE control_epw,        ONLY : wannierize
    USE epwcom,             ONLY : vme, eig_read, lifc, etf_mem, nbndsub, use_ws, &
                                   epwread, epbread, lpolar, lphase, system_2d, &
                                   longrange, nkc1, nkc2, nkc3, nqc1, nqc2, nqc3
    USE io_global,          ONLY : ionode_id, meta_ionode
    USE io_files,           ONLY : tmp_dir, prefix
    USE io_var,             ONLY : crystal, iunepmatwp, iunepmatwp2
    USE control_flags,      ONLY : gamma_only
    !
    USE modes,              ONLY : nmodes
    USE pwcom,              ONLY : nelec
    USE noncollin_module,   ONLY : noncolin
    USE cell_base,          ONLY : at, bg, alat, omega
    USE ions_base,          ONLY : nat, tau, ityp, amass
    USE elph2,              ONLY : chw, epmatwp, eps_rpa, &
                                   wf, qrpl
    USE wigner,             ONLY : wigner_seitz_wrap
    !
    USE wan2bloch,          ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch, &
                                   ephwan2blochp, ephwan2bloch, vmewan2bloch, &
                                   dynifc2blochf, vmewan2blochp
    USE iso_c_binding
    !
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! Contains the state of the opened file
    INTEGER :: ierr
    !! Error index when reading/writing a file
    !
    INTEGER :: nrr_k
    !! Number of WS vectors for the electrons
    INTEGER :: nrr_q
    !! Number of WS vectors for the phonons
    INTEGER :: nrr_g
    !! Number of WS vectors for the electron-phonons
    !
    REAL(KIND = DP), ALLOCATABLE :: irvec_r(:, :)
    !! Wigner-Size supercell vectors, store in real instead of integer
    INTEGER, ALLOCATABLE :: irvec_k(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, ALLOCATABLE :: irvec_q(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE :: irvec_g(:, :)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for
    !electron-phonon
    INTEGER, ALLOCATABLE :: ndegen_k(:, :, :)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons
    !grid
    INTEGER, ALLOCATABLE :: ndegen_q(:, :, :)
    !! Wigner-Seitz weights for the phonon grid that depend on atomic
    !positions $R
    !+ \tau(nb) - \tau(na)$
    INTEGER, ALLOCATABLE :: ndegen_g(:, :, :)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R - \tau(na)$
    REAL(KIND = DP), ALLOCATABLE :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE :: wslen_g(:)
    !! real-space length for electron-phonons, in units of alat
    !
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :, :)
    !! e-p matrix  in el wannier - fine Bloch phonon grid
    COMPLEX(KIND = DP), ALLOCATABLE :: uf(:, :)
    !! Rotation matrix for phonons
    REAL(KIND = DP), ALLOCATABLE :: w(:)
    !
    INTEGER :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER :: dims2
    !! Dims is either nat if use_ws or 1 if not
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    INTEGER :: direct_io_factor
    !! Type of IOLength
    INTEGER*8 :: unf_recl
    !! Record length
    INTEGER :: lrepmatw
    !! record length while reading file
    REAL(KIND = DP) :: dummy(3)
    !! Dummy variable
    !
    !
    INTERFACE
    SUBROUTINE fill_nbands(nbnd) bind(c, name='fillNbands')
      !
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER(c_int), INTENT(IN) :: nbnd
      !
    END SUBROUTINE fill_nbands
    END INTERFACE
    !
    !
    INTERFACE
    SUBROUTINE fill_nelec(nelec) bind(c, name='fillNelec')
      !
      USE iso_c_binding
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: nelec
      !
    END SUBROUTINE fill_nelec
    END INTERFACE
    !
    !
    INTERFACE
    SUBROUTINE fill_nmodes(nmds) bind(c, name='fillNmodes')
      !
      USE iso_c_binding
      IMPLICIT NONE
      INTEGER(c_int), INTENT(IN) :: nmds
      !
    END SUBROUTINE fill_nmodes
    END INTERFACE
    !
    !
    INTERFACE
    SUBROUTINE fill_lattice(lat_const, unit_cell_vol, coords) bind(c, name='fillLattice')
      !
      USE iso_c_binding
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: lat_const
      REAL(c_double), INTENT(IN) :: unit_cell_vol
      REAL(c_double), INTENT(IN) :: coords(3, 3)
      !
    END SUBROUTINE fill_lattice
    END INTERFACE
    !
    !
    CONTAINS
    !
    !
    SUBROUTINE epi_init(rank) bind(C, name='epiInit')
      !
        IMPLICIT NONE
        INTEGER(c_int), INTENT(IN) :: rank
        INTEGER(c_int) :: use_ws_int
        INTEGER(c_int) :: nkc1, nkc2, nkc3, nqc1, nqc2, nqc3
        !
        REAL(KIND = DP), ALLOCATABLE :: w_centers(:, :)
        !REAL(KIND = DP) :: eF
        !eF = 7.778!12.3606
        !! Wannier centers
        !
        gamma_only = .FALSE.
        !
        meta_ionode = (rank == 0)
        !CALL epw_readin()
        vme = 'wannier'
        lifc = .FALSE.
        etf_mem = 1
        eig_read = .FALSE.
        epwread = .TRUE.
        epbread = .FALSE.
        wannierize = .FALSE.
        longrange = .FALSE.

        lpolar = .FALSE.
        lphase = .FALSE.
        system_2d  = .FALSE.
        qrpl = .FALSE.
        use_ws = .FALSE.
        nkc1 = 6
        nkc2 = 6
        nkc3 = 6
        nqc1 = 6
        nqc2 = 6
        nqc3 = 6
        use_ws_int = .FALSE.
        ! !
        ALLOCATE(w_centers(3, nbndsub), STAT = ierr)
        IF (ierr /= 0) CALL errore('epi_init', 'Error allocating w_centers', 1)
        !
        CALL read_crystal_info(w_centers)
        CALL fill_lattice(alat, omega, at)
        CALL epw_read(nrr_k, nrr_q, nrr_g) !!! need to be fixed from inside in EPW !!!
        !
        CALL fill_nelec(nelec)
        CALL fill_nbands(nbndsub)
        CALL fill_nmodes(nmodes)
        CALL allocate_epi_stuff_arrays()
        !
        IF (use_ws) THEN
                ! Use Wannier-centers to contstruct the WS for electonic part and el-ph part
                ! Use atomic position to contstruct the WS for the phonon part
                dims  = nbndsub
                dims2 = nat
                CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q, irvec_g, &
                                       ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                                       w_centers, dims, tau, dims2)
        ELSE
                ! Center the WS at Gamma for electonic part, the phonon part and el-ph part
                dims  = 1
                dims2 = 1
                dummy(:) = [0.0, 0.0, 0.0]
                CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q, irvec_g, &
                                       ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                                       dummy, dims, dummy, dims2)
        ENDIF
        !! This is simply because dgemv take only real number (not integer)
        ALLOCATE(irvec_r(3, nrr_k), STAT = ierr)
        IF (ierr /= 0) CALL errore('epi_init', 'Error allocating irvec_r', 1)
        irvec_r = REAL(irvec_k, KIND = DP)
        ! Open the ephmatwp file here
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        filint = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
        INQUIRE(IOLENGTH = direct_io_factor) dummy(1)
        unf_recl = direct_io_factor * INT(lrepmatw, KIND = KIND(unf_recl))
        IF (unf_recl <= 0) CALL errore('epw_init', 'wrong record length', 3)
        OPEN(iunepmatwp, FILE = TRIM(ADJUSTL(filint)), IOSTAT = ierr, FORM='unformatted', &
            STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
        IF (ierr /= 0) CALL errore('epw_init', 'error opening ' // TRIM(filint), 1)
    !
    END SUBROUTINE epi_init
    !
    !
    SUBROUTINE allocate_epi_stuff_arrays()
      !
      ! The same as in ephwann_shuffle.f90
      ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('allocate_epi_stuff_arrays', 'Error allocating epmatwef', 1)
      ALLOCATE(uf(nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('allocate_epi_stuff_arrays', 'Error allocating uf', 1)
      ALLOCATE(w(nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('allocate_epi_stuff_arrays', 'Error allocating w2', 1)
      epmatwef(:, :, :, :) = czero
      uf(:, : ) = czero
      !
    END SUBROUTINE allocate_epi_stuff_arrays
    !
    !
    SUBROUTINE deallocate_epi_stuff_arrays() bind(C, name='epiFinalize')
      !
      DEALLOCATE(epmatwef, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_epi_stuff_arrays', 'Error deallocating epmatwef', 1)
      DEALLOCATE(uf, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_epi_stuff_arrays', 'Error deallocating uf', 1)
      DEALLOCATE(w, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_epi_stuff_arrays', 'Error deallocating w2', 1)
      DEALLOCATE(tau, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_epi_stuff_arrays', 'Error deallocating tau', 1)
      DEALLOCATE(ityp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_epi_stuff_arrays', 'Error deallocating ityp', 1)
      DEALLOCATE(irvec_r, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_epi_stuff_arrays', 'Error deallocating irvec_r', 1)
      !
      RETURN
      !
    END SUBROUTINE deallocate_epi_stuff_arrays
    !
    !
    SUBROUTINE read_crystal_info(w_centers)
      IMPLICIT NONE
      REAL(KIND = DP), INTENT(INOUT) :: w_centers(3, nbndsub)
      !! Wannier centers
      !
      IF (mpime == ionode_id) THEN
        !
        OPEN(UNIT = crystal, FILE = 'crystal.fmt', STATUS = 'old', IOSTAT = ios)
        IF (ios /= 0) CALL errore('read_crystal_info', 'error opening crystal.fmt', crystal)
        READ(crystal, *) nat
        READ(crystal, *) nmodes
        READ(crystal, *) nelec
        READ(crystal, *) at
        READ(crystal, *) bg
        READ(crystal, *) omega
        READ(crystal, *) alat
        ALLOCATE(tau(3, nat), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_crystal_info', 'Error allocating tau', 1)
        READ(crystal, *) tau
        READ(crystal, *) amass
        ALLOCATE(ityp(nat), STAT = ierr)
        IF (ierr /= 0) CALL errore('read_crystal_info', 'Error allocating ityp', 1)
        READ(crystal, *) ityp
        READ(crystal, *) noncolin
        READ(crystal, *) w_centers
        !
      ENDIF ! mpime == ionode_id
      CALL mp_bcast(nat      , ionode_id, world_comm)
      IF (mpime /= ionode_id) ALLOCATE(ityp(nat))
      CALL mp_bcast(nmodes   , ionode_id, world_comm)
      CALL mp_bcast(nelec    , ionode_id, world_comm)
      CALL mp_bcast(at       , ionode_id, world_comm)
      CALL mp_bcast(bg       , ionode_id, world_comm)
      CALL mp_bcast(omega    , ionode_id, world_comm)
      CALL mp_bcast(alat     , ionode_id, world_comm)
      IF (mpime /= ionode_id) ALLOCATE(tau(3, nat))
      CALL mp_bcast(tau      , ionode_id, world_comm)
      CALL mp_bcast(amass    , ionode_id, world_comm)
      CALL mp_bcast(ityp     , ionode_id, world_comm)
      CALL mp_bcast(noncolin , ionode_id, world_comm)
      CALL mp_bcast(w_centers, ionode_id, world_comm)
      IF (mpime == ionode_id) THEN
        CLOSE(crystal)
      ENDIF 
      !
    END SUBROUTINE read_crystal_info
    !
END MODULE epi_initialize
