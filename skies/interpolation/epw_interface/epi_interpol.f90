!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                  !!!
!!!           Module for EPI EPW interpolation facilities            !!!
!!!           Major part of the code in this file is                 !!!
!!!           copied directly from EPW software                      !!!
!!!           (e.g. ephwann_shuffle.f90 and other src files)         !!!
!!!           to reuse the existing subroutines in the               !!!
!!!           transport spectral function calculation.               !!!
!!!                                                                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE epi_interpol
    !
    USE kinds,              ONLY : DP
    USE constants_epw,      ONLY : zero, czero, eps8, twopi, ci, two
    !
    USE mp,                 ONLY : mp_bcast
    USE mp_global,          ONLY : mp_startup, world_comm
    USE io_epw,             ONLY : epw_read
    USE control_epw,        ONLY : wannierize
    USE epwcom,             ONLY : vme, eig_read, lifc, etf_mem, nbndsub, &
                                   use_ws, nqc1, nqc2, nqc3, nkc1, nkc2, nkc3, &
                                   epwread, epbread
    USE io_global,          ONLY : ionode_id
    USE io_files,           ONLY : prefix, tmp_dir
    USE io_var,             ONLY : crystal
    USE control_flags,      ONLY : gamma_only

    USE modes,              ONLY : nmodes
    USE pwcom,              ONLY : nelec
    USE noncollin_module,   ONLY : noncolin
    USE cell_base,          ONLY : at, bg, alat, omega
    USE ions_base,          ONLY : nat, tau, ityp, amass
    USE elph2,              ONLY : chw, epmatwp, eps_rpa, wf
    USE wigner,             ONLY : wigner_seitz_wrap
    !
    USE wan2bloch,          ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch, &
                                   ephwan2blochp, ephwan2bloch, vmewan2bloch, &
                                   dynifc2blochf, vmewan2blochp
    !
    USE epi_initialize
    USE iso_c_binding
    !
    IMPLICIT NONE
    !
    CONTAINS 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!        These subroutines are parts of the elphwan_shuffle.f90 file        !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !
    ! adaptor functions to use in iso_c_binding
    !
    ! PHONON FREQUENCIES
    ! IT'S WRONG TO CALL THIS SUBROUTINE IN MANY THREADS (uf and epmatwef are shared vars)!
    SUBROUTINE interpolate_eigen_freq_at(kx, ky, kz, frs) bind(c, name='interpEigenFreqAt')
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz
      REAL(c_double), INTENT(INOUT) :: frs(nmodes)
      !! output value of frequencies at given q
      !
      REAL(KIND = DP), DIMENSION(3) :: q
      !! q-point in crystal coordinates
      REAL(KIND = DP) :: w2(nmodes)
      !! Interpolated phonon squared frequency
      INTEGER :: nu, mu
      INTEGER :: na
      !
      !! Counter on atom
      !
      uf(:, :) = czero
      epmatwef(:, :, :, :) = czero
      !
      q = [kx, ky, kz]
      !
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, q, uf, w2)
      !ELSE
        !CALL errore('lifc is .TRUE.')
        !CALL dynifc2blochf(nmodes, rws, nrws, q, uf, w2)
      ENDIF
      !
      DO nu = 1, nmodes
        IF (w2(nu) > -eps8) THEN
          frs(nu) =  DSQRT(ABS(w2(nu)))
        ELSE
          frs(nu) = -DSQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      !
      w(:) = frs(:)
      !
      ! ...then take into account the mass factors and square-root the frequencies...
      !
      DO nu = 1, nmodes
        DO mu = 1, nmodes
          na = (mu - 1) / 3 + 1
          uf(mu, nu) = uf(mu, nu) / DSQRT(amass(ityp(na)))
        ENDDO
      ENDDO
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
      ! --------------------------------------------------------------
      !
      CALL ephwan2blochp_(nmodes, q, irvec_g, ndegen_g, nrr_g, uf, epmatwef, nbndsub, nrr_k, dims, dims2)
      !
    END SUBROUTINE interpolate_eigen_freq_at
    !
    !
    SUBROUTINE interpolate_eigen_freq1D_at(kx, ky, kz, frs) bind(c, name='interpEigenFreq1DAt')
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz
      REAL(c_double), INTENT(INOUT) :: frs(nmodes)
      !! output value of frequencies at given q
      !
      REAL(KIND = DP), DIMENSION(3) :: q
      !! q-point in crystal coordinates
      REAL(KIND = DP) :: w2(nmodes)
      !! Interpolated phonon squared frequency
      INTEGER :: nu
      !
      uf(:, :) = czero
      epmatwef(:, :, :, :) = czero
      !
      q = [kx, ky, kz]
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, q, uf, w2)
      !ELSE
        !CALL errore('lifc is .TRUE.')
        !CALL dynifc2blochf(nmodes, rws, nrws, q, uf, w2)
      ENDIF
      !
      DO nu = 1, nmodes
        IF (w2(nu) > -eps8) THEN
          frs(nu) =  DSQRT(ABS(w2(nu)))
        ELSE
          frs(nu) = -DSQRT(ABS(w2(nu)))
        ENDIF
      ENDDO
      !
    END SUBROUTINE interpolate_eigen_freq1D_at
    !
    ! ELECTRON ENERGIES
    SUBROUTINE interpolate_eigen_value_at(kx, ky, kz, e) bind(c, name='interpEigenValueAt')
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz
      !! k-point for interpolation at
      REAL(c_double), INTENT(INOUT) :: e(nbndsub)
      !! output value of eigvals at given k
      INTEGER(c_int) :: iw
      !! Counter on bands when use_ws == .TRUE.
      INTEGER(c_int) :: iw2
      !! Counter on bands when use_ws == .TRUE.
      INTEGER(c_int) :: ir
      !! Counter for WS loop
      COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
      !! Rotation matrix, fine mesh, points k
      COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
      !! Used to store $e^{2\pi r \cdot k}$ exponential
      REAL(KIND = DP) :: rdotk(nrr_k)
      !! $r\cdot k$
      !
      REAL(c_double), DIMENSION(3) :: k
      !
      k = [kx, ky, kz]
      ! k is assumed to be in crys coord
      !
      cufkk(:, :)          = czero
      cfac(:, :, :)        = czero
      rdotk(:)             = zero
      !
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, k, 1, 0.0_DP, rdotk, 1)
      !
      IF (use_ws) THEN
        DO iw = 1, dims
          DO iw2 = 1, dims
            DO ir = 1, nrr_k
              IF (ndegen_k(ir, iw2, iw) > 0) THEN
                cfac(ir, iw2, iw)  = EXP(ci * rdotk(ir))  / ndegen_k(ir, iw2, iw)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        cfac(:, 1, 1)  = EXP(ci * rdotk(:))  / ndegen_k(:, 1, 1)
      ENDIF
      !
      ! ------------------------------------------------------
      ! hamiltonian : Wannier -> Bloch
      ! ------------------------------------------------------
      !
      CALL hamwan2bloch(nbndsub, nrr_k, cufkk, e(:), chw, cfac, dims)
      !
    END SUBROUTINE interpolate_eigen_value_at
    !
    !
    ! ELECTRON VELOCITIES
    SUBROUTINE interpolate_velocity_at(kx, ky, kz, v, dir) bind(c, name='interpVelocityAt')
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz
      !! k-point for interpolation at
      REAL(c_double), INTENT(INOUT) :: v(nbndsub)
      !! output value of velocities at given k
      INTEGER(c_int), INTENT(IN) :: dir
      !! cartesian direction of velocity (must be one of 1, 2 or 3)
      !
      REAL(KIND = DP) :: e(nbndsub)
      !! output value of eigvals at given k
      REAL(KIND = DP) :: eks(nbndsub)
      !! dummy variable for calling vmewan2bloch (eig_read = .FALSE.)
      INTEGER :: I
      !! index over bands
      COMPLEX(KIND = DP) :: vmef(3, nbndsub, nbndsub)
      !
      COMPLEX(KIND = DP) :: cufkk(nbndsub, nbndsub)
      !! Rotation matrix, fine mesh, points k
      COMPLEX(KIND = DP) :: cfac(nrr_k, dims, dims)
      !! Used to store $e^{2\pi r \cdot k}$ exponential
      !
      REAL(c_double), DIMENSION(3) :: k
      !
      k = [kx, ky, kz]
      !
      CALL eigenvalues_preapre_at_k(k, cufkk, cfac, e)
      !
      ! ------------------------------------------------------
      !  velocity: Wannier -> Bloch
      ! ------------------------------------------------------
      !
      vmef(:, :, :) = czero
      eks(:) = zero
      !
      CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :), e, eks, chw, cfac, dims)
      DO I = 1, nbndsub
        v(I) = REAL(vmef(dir, I, I))
      END DO
      !
    END SUBROUTINE interpolate_velocity_at
    !
    !
    ! ELECTRON-PHONON INTERACTION
    SUBROUTINE interpolate_gmatrix_at(kx, ky, kz, qx, qy, qz, nu, nini, nfin, g) bind(c, name="interpGmatAtq")
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz, qx, qy, qz
      INTEGER(c_int), INTENT(IN) :: nu, nini, nfin
      REAL(c_double), INTENT(INOUT) :: g
      !
      REAL(KIND = DP),  DIMENSION(3) :: k
      !! k-point in crystal coordinates
      REAL(KIND = DP),  DIMENSION(3) :: q
      !! q-point in crystal coordinates
      REAL(c_double), DIMENSION(nbndsub) :: e
      !! output value of eigvals at given k
      COMPLEX(KIND = DP), DIMENSION(nbndsub, nbndsub) :: cufkk
      !! Rotation matrix, fine mesh, for points k
      COMPLEX(KIND = DP), DIMENSION(nbndsub, nbndsub) :: cufkq
      !! Rotation matrix, fine mesh, for points k+q
      COMPLEX(KIND = DP), DIMENSION(nrr_k, dims, dims) :: cfac
      !! Used to store $e^{2\pi r \cdot k}$ exponential
      COMPLEX(KIND = DP), DIMENSION(nrr_k, dims, dims) :: cfacq
      !! Used to store $e^{2\pi r \cdot k+q}$ exponential
      COMPLEX(KIND = DP) epmatf(nbndsub, nbndsub, nmodes)
      !! e-p matrix  in smooth Bloch basis, fine mesh
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
      ! --------------------------------------------------------------
      !
      k = [kx, ky, kz]
      q = [qx, qy, qz]
      !
      epmatf(:, :, :)      = czero
      !
      CALL eigenvalues_preapre_at_k(k, cufkk, cfac, e)
      CALL eigenvalues_preapre_at_k(k + q, cufkq, cfacq, e)
      CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, dims)
      g = ABS(epmatf(nfin + 1, nini + 1, nu + 1)) / SQRT(two * w(nu + 1))
      !  
    END SUBROUTINE interpolate_gmatrix_at
    !
    !
    ! EPh MATRIX ELEMENTS
    SUBROUTINE interpolate_gmatrix2_at(kx, ky, kz, qx, qy, qz, nu, nini, nfin, g2) bind(c, name="interpGmat2Atq")
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz, qx, qy, qz
      INTEGER(c_int), INTENT(IN) :: nu, nini, nfin
      REAL(c_double), INTENT(INOUT) :: g2
      !
      COMPLEX(KIND = DP) :: g
      ! matrix element
      REAL(KIND = DP),  DIMENSION(3) :: k
      !! k-point in crystal coordinates
      REAL(KIND = DP),  DIMENSION(3) :: q
      !! q-point in crystal coordinates
      REAL(c_double), DIMENSION(nbndsub) :: e
      !! output value of eigvals at given k
      COMPLEX(KIND = DP), DIMENSION(nbndsub, nbndsub) :: cufkk
      !! Rotation matrix, fine mesh, for points k
      COMPLEX(KIND = DP), DIMENSION(nbndsub, nbndsub) :: cufkq
      !! Rotation matrix, fine mesh, for points k+q
      COMPLEX(KIND = DP), DIMENSION(nrr_k, dims, dims) :: cfac
      !! Used to store $e^{2\pi r \cdot k}$ exponential
      COMPLEX(KIND = DP), DIMENSION(nrr_k, dims, dims) :: cfacq
      !! Used to store $e^{2\pi r \cdot k+q}$ exponential
      COMPLEX(KIND = DP) epmatf(nbndsub, nbndsub, nmodes)
      !! e-p matrix  in smooth Bloch basis, fine mesh
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
      ! --------------------------------------------------------------
      !
      k = [kx, ky, kz]
      q = [qx, qy, qz]
      !
      epmatf(:, :, :)      = czero
      !
      CALL eigenvalues_preapre_at_k(k, cufkk, cfac, e)
      CALL eigenvalues_preapre_at_k(k + q, cufkq, cfacq, e)
      CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, dims)
      g = epmatf(nfin + 1, nini + 1, nu + 1)
      g2 = 0.5E0_DP * (REAL(g) * REAL(g) + AIMAG(g) * AIMAG(g)) / w(nu + 1)
      ! 
    END SUBROUTINE interpolate_gmatrix2_at
    !
    !
    ! helper subroutine to use in other subroutines
    SUBROUTINE eigenvalues_preapre_at_k(k, cufkk, cfac, e)
      !
      IMPLICIT NONE
      REAL(KIND = DP), INTENT(IN), DIMENSION(3) :: k
      !! k-point in crystal coordinates
      REAL(c_double), INTENT(INOUT), DIMENSION(nbndsub) :: e
      !! output value of eigvals at given k
      COMPLEX(KIND = DP), INTENT(INOUT), DIMENSION(nbndsub, nbndsub) :: cufkk
      !! Rotation matrix, fine mesh, for points k
      COMPLEX(KIND = DP), INTENT(INOUT), DIMENSION(nrr_k, dims, dims) :: cfac
      !! Used to store $e^{2\pi r \cdot k}$ exponential
      !
      !
      REAL(KIND = DP), DIMENSION(nrr_k) :: rdotk
      !! $r\cdot k$
      INTEGER(c_int) :: iw
      !! Counter on bands when use_ws == .TRUE.
      INTEGER(c_int) :: iw2
      !! Counter on bands when use_ws == .TRUE.
      INTEGER(c_int) :: ir
      !! Counter for WS loop
      !
      !
      e(:)                 = zero
      cufkk(:, :)          = czero
      rdotk(:)             = zero
      cfac(:, :, :)        = czero
      !
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, k, 1, 0.0_DP, rdotk, 1)
      !
      IF (use_ws) THEN
        DO iw = 1, dims
          DO iw2 = 1, dims
            DO ir = 1, nrr_k
              IF (ndegen_k(ir, iw2, iw) > 0) THEN
                cfac(ir, iw2, iw) = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        cfac(:, 1, 1) = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
      ENDIF
      !
      ! ------------------------------------------------------
      ! hamiltonian : Wannier -> Bloch
      ! ------------------------------------------------------
      !
      CALL hamwan2bloch(nbndsub, nrr_k, cufkk, e, chw, cfac, dims)
      !
    END SUBROUTINE eigenvalues_preapre_at_k
    !
    !
    ! This special version does much recalculations just because its for drawing a figure
    SUBROUTINE interpolate_gmatrix1D_at(kx, ky, kz, qx, qy, qz, nu, nini, nfin, g) bind(c, name="interpGmat1DAtq")
      !
      IMPLICIT NONE
      REAL(c_double), INTENT(IN) :: kx, ky, kz, qx, qy, qz
      INTEGER(c_int), INTENT(IN) :: nu, nini, nfin
      REAL(c_double), INTENT(INOUT) :: g
      !
      REAL(KIND = DP),  DIMENSION(3) :: k
      !! k-point in crystal coordinates
      REAL(KIND = DP),  DIMENSION(3) :: q
      !! q-point in crystal coordinates
      REAL(KIND = DP), DIMENSION(nbndsub) :: e
      !! output value of eigvals at given k
      COMPLEX(KIND = DP), DIMENSION(nbndsub, nbndsub) :: cufkk
      !! Rotation matrix, fine mesh, for points k
      COMPLEX(KIND = DP), DIMENSION(nbndsub, nbndsub) :: cufkq
      !! Rotation matrix, fine mesh, for points k+q
      COMPLEX(KIND = DP), DIMENSION(nrr_k, dims, dims) :: cfac
      !! Used to store $e^{2\pi r \cdot k}$ exponential
      COMPLEX(KIND = DP), DIMENSION(nrr_k, dims, dims) :: cfacq
      !! Used to store $e^{2\pi r \cdot k+q}$ exponential
      REAL(KIND = DP) :: om(nmodes)
      !! output value of frequencies at given q
      COMPLEX(KIND = DP) epmatf(nbndsub, nbndsub, nmodes)
      !! e-p matrix  in smooth Bloch basis, fine mesh
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
      ! --------------------------------------------------------------
      !
      k = [kx, ky, kz]
      q = [qx, qy, qz]
      !
      epmatf(:, :, :)      = czero
      !
      CALL interpolate_eigen_freq_at(qx, qy, qz, om)        ! --> uf, epmatwef
      CALL eigenvalues_preapre_at_k(k, cufkk, cfac, e)
      CALL eigenvalues_preapre_at_k(k + q, cufkq, cfacq, e)          ! --> cufkq
      !
      CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, dims)
      g = ABS(epmatf(nfin + 1, nini + 1, nu + 1)) / SQRT(two * om(nu + 1))
      !  
    END SUBROUTINE interpolate_gmatrix1D_at
    !
    !
    SUBROUTINE ephwan2blochp_(nmodes, xxq, irvec_g, ndegen_g, nrr_g, cuf, epmatf, nbnd, nrr_k, dims, nat)
      !!
      !!   Adopted from similar subroutine in wan2bloch.f90
      !!
      USE kinds,            ONLY : DP
      USE epwcom,           ONLY : etf_mem, use_ws
      USE elph2,            ONLY : epmatwp
      USE constants_epw,    ONLY : twopi, ci, czero, cone
      USE io_var,           ONLY : iunepmatwp
      USE io_epw,           ONLY : rwepmatw
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: nmodes
      !! Total number of modes
      INTEGER, INTENT(in) :: nrr_g
      !! Number of phononic WS points
      INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
      !! Coordinates of WS points
      INTEGER, INTENT(in) :: dims
      !! Is equal to the number of Wannier function if use_ws == .TRUE. Is equal to 1 otherwise.
      INTEGER, INTENT(in) :: nat
      !! Is equal to the number of atoms if use_ws == .TRUE. or 1 otherwise
      INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
      !! Number of degeneracy of WS points
      INTEGER, INTENT(in) :: nbnd
      !! Number of bands
      INTEGER, INTENT(in) ::  nrr_k
      !! Number of electronic WS points
      REAL(KIND = DP), INTENT(in) :: xxq(3)
      !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
      COMPLEX(KIND = DP), INTENT(in) :: cuf(nmodes, nmodes)
      !! e-p matrix in Wanner representation
      COMPLEX(KIND = DP), INTENT(out) :: epmatf(nbnd, nbnd, nrr_k, nmodes)
      !! e-p matrix in Bloch representation, fine grid
      !
      ! Local variables
      INTEGER :: ir
      !! Real space WS index
      INTEGER :: iw
      !! Wannier function index
      INTEGER :: irn
      !! Combined WS and atom index
      INTEGER :: ir_start
      !! Starting ir for this cores
      INTEGER :: ir_stop
      !! Ending ir for this pool
      INTEGER :: na
      !! Atom index
      INTEGER :: imode
      !! Number of modes
      INTEGER :: diff
      !! Difference between starting and ending on master core
      INTEGER :: ierr
      !! Error status
      REAL(KIND = DP) :: rdotk
      !! Exponential for the FT
      COMPLEX(KIND = DP) :: eptmp(nbnd, nbnd, nrr_k, nmodes)
      !! Temporary matrix to store el-ph
      COMPLEX(KIND = DP) :: cfac(dims, nat, nrr_g)
      !! Factor for the FT
      COMPLEX(KIND = DP), ALLOCATABLE :: epmatw(:, :, :, :)
      !! El-ph matrix elements
      !
      ir_start = 1
      ir_stop = nrr_g * nmodes
      IF (use_ws) ir_stop = nrr_g * nat
      diff = ir_stop - ir_start
      !
      eptmp(:, :, :, :) = czero
      cfac(:, :, :) = czero
      !
      IF (use_ws) THEN
        DO irn = ir_start, ir_stop
          ir = (irn - 1) / nat + 1
          na = MOD(irn - 1, nat) +1
          rdotk = twopi * DOT_PRODUCT(xxq, DBLE(irvec_g(:, ir)))
          DO iw = 1, dims
            IF (ndegen_g(iw, ir, na) > 0) &
              cfac(iw, na, ir) = EXP(ci * rdotk) / DBLE(ndegen_g(iw, ir, na))
          ENDDO
        ENDDO
      ELSE
        DO irn = ir_start, ir_stop
          ir = (irn - 1) / nmodes + 1
          rdotk = twopi * DOT_PRODUCT(xxq, DBLE(irvec_g(:, ir)))
          cfac(1, 1, ir) = EXP(ci * rdotk) / DBLE(ndegen_g(1, ir, 1))
        ENDDO
      ENDIF
      !
      IF (etf_mem == 0) then
        !
        IF (use_ws) THEN
          DO irn = ir_start, ir_stop
            ir = (irn - 1) / nat + 1
            na = MOD(irn - 1, nat) + 1
            DO iw = 1, dims
              CALL ZAXPY(nbnd * nrr_k * 3, cfac(iw, na, ir), epmatwp(iw, :, :, 3 * (na - 1) + 1:3 * na, ir), 1, &
                   eptmp(iw, :, :, 3 * (na - 1) + 1:3 * na), 1)
            ENDDO
          ENDDO
        ELSE ! use_ws
          DO irn = ir_start, ir_stop
            ir = (irn - 1) / nmodes + 1
            imode = MOD(irn-1, nmodes) + 1
            CALL ZAXPY(nbnd * nbnd * nrr_k, cfac(1, 1, ir), epmatwp(:, :, :, imode, ir), 1, eptmp(:, :, :, imode), 1)
          ENDDO
        ENDIF
        !
      ELSE ! etf_mem != 0
        !
        IF (use_ws) THEN
          ALLOCATE(epmatw(nbnd, nbnd, nrr_k, nmodes), STAT = ierr)
          IF (ierr /= 0) CALL errore('ephwan2blochp_', 'Error allocating epmatw', 1)
          DO irn = ir_start, ir_stop
            ir = (irn - 1) / nat + 1
            na = MOD(irn - 1, nat) + 1
            CALL rwepmatw(epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
            DO iw = 1, dims
              CALL ZAXPY(nrr_k * 3 * nbnd, cfac(iw, na, ir), epmatw(iw, :, :, 3 * (na - 1) + 1:3 * na), 1, &
                   eptmp(iw, :, :, 3 * (na - 1) + 1:3 * na), 1)
            ENDDO
          ENDDO ! irn
        ELSE ! use_ws
          ALLOCATE(epmatw(nbnd, nbnd, nrr_k, nmodes), STAT = ierr)
          IF (ierr /= 0) CALL errore('ephwan2blochp_', 'Error allocating epmatw', 1)
          !
          DO irn = ir_start, ir_stop
            ir = (irn - 1) / nmodes + 1
            imode = MOD(irn - 1, nmodes) + 1
            CALL rwepmatw(epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
            CALL ZAXPY(nbnd * nbnd * nrr_k, cfac(1, 1, ir), &
                epmatw(:, :, :, imode), 1, eptmp(:, :, :, imode), 1)
          ENDDO
        ENDIF ! use_ws
        DEALLOCATE(epmatw, STAT = ierr)
        IF (ierr /= 0) CALL errore('ephwan2blochp_', 'Error deallocating epmatw', 1)
      ENDIF ! etf_mem
      !
      Call ZGEMM('n', 'n', nbnd * nbnd * nrr_k, nmodes, nmodes, cone, eptmp, &
                  nbnd * nbnd * nrr_k, cuf, nmodes, czero, epmatf, nbnd * nbnd * nrr_k)
      !
    END SUBROUTINE ephwan2blochp_
    !
END MODULE epi_interpol
