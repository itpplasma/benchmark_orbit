program test_field_representation_error
    !> Compare Bmod from:
    !>   Pass 0: direct Fourier summation from wout.nc (reference)
    !>   Pass 1: VMEC splines        (magfie with VMEC field type)
    !>   Pass 2: Boozer splines      (magfie with BOOZER field type)
    !>   Pass 3: canonical splines   (magfie with CANFLUX field type)
    !>
    !> Passes 1-3 are compared against the Fourier reference from pass 0.
    !>
    !> Usage: test_field_representation_error.x [wout.nc] [ns_s] [ns_tp] [multharm]
    !>                                          [n_r] [n_th] [n_ph]
    !> Match ns_s, ns_tp, multharm to the values used in simple.in.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi
    use simple, only: init_vmec
    use new_vmec_stuff_mod, only: nper, multharm, ns_s, ns_tp
    use vector_potentail_mod, only: ns
    use magfie_sub, only: magfie, init_magfie, VMEC, BOOZER, CANFLUX
    use boozer_sub, only: get_boozer_coordinates, vmec_to_boozer, boozer_to_vmec
    use get_can_sub, only: get_canonical_coordinates, vmec_to_can, can_to_vmec
    use spline_vmec_sub, only: splint_vmec_data
    use netcdf

    implicit none

    character(len=1024) :: wout_file, arg
    integer :: n_r, n_th, n_ph, n_points
    integer :: i_r, i_th, i_ph, k
    real(dp) :: fper, r, theta, phi, s, phi_period
    real(dp) :: x(3), bmod, sqrtg
    real(dp) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: theta_c, phi_c, theta_B, phi_B, err

    real(dp), allocatable :: bmod_ref(:), sqrtg_vmec(:)
    real(dp), allocatable :: bmod_vmec_arr(:), bmod_boozer_arr(:), bmod_can_arr(:)
    real(dp), allocatable :: R_arr_phys(:), Z_arr_phys(:)
    real(dp), allocatable :: r_arr(:), theta_arr(:), phi_arr(:)

    real(dp) :: max_err_vmec_spl, max_err_boozer, max_err_can
    real(dp) :: sum_err_vmec_spl, sum_err_boozer, sum_err_can, sum_weight
    real(dp) :: vol_avg_err_vmec_spl, vol_avg_err_boozer, vol_avg_err_can

    ! R, Z from each spline representation (via inverse coordinate transform)
    real(dp), allocatable :: R_vmec_spl(:), Z_vmec_spl(:)
    real(dp), allocatable :: R_boozer_spl(:), Z_boozer_spl(:)
    real(dp), allocatable :: R_can_spl(:), Z_can_spl(:)
    real(dp) :: max_Rerr_vmec, max_Zerr_vmec
    real(dp) :: max_Rerr_boozer, max_Zerr_boozer
    real(dp) :: max_Rerr_can, max_Zerr_can
    ! Workspace for splint_vmec_data
    real(dp) :: A_phi_w, A_theta_w, dA_phi_w, dA_theta_w, aiota_w
    real(dp) :: Rw, Zw, aw, dRds, dRdt, dRdp, dZds, dZdt, dZdp, dlds, dldt, dldp
    real(dp) :: theta_v_back, phi_v_back

    ! ── command line ─────────────────────────────────────────────────────────
    wout_file = 'wout.nc'
    ns_s     = 5
    ns_tp    = 5
    multharm = 5
    n_r  = 10
    n_th = 20
    n_ph = 10

    if (command_argument_count() >= 1) call get_command_argument(1, wout_file)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg); read(arg, *) ns_s
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read(arg, *) ns_tp
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg); read(arg, *) multharm
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg); read(arg, *) n_r
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg); read(arg, *) n_th
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg); read(arg, *) n_ph
    end if

    n_points = n_r * n_th * n_ph
    allocate(bmod_ref(n_points), sqrtg_vmec(n_points))
    allocate(bmod_vmec_arr(n_points), bmod_boozer_arr(n_points), bmod_can_arr(n_points))
    allocate(R_arr_phys(n_points), Z_arr_phys(n_points))
    allocate(R_vmec_spl(n_points), Z_vmec_spl(n_points))
    allocate(R_boozer_spl(n_points), Z_boozer_spl(n_points))
    allocate(R_can_spl(n_points), Z_can_spl(n_points))
    allocate(r_arr(n_points), theta_arr(n_points), phi_arr(n_points))

    print '(A)',           '================================================='
    print '(A)',           ' Field representation error vs VMEC Fourier series'
    print '(A)',           '================================================='
    print '(A,A)',         ' VMEC file : ', trim(wout_file)
    print '(A,I4,A,I4,A,I4)', ' Spline    : ns_s=', ns_s, '  ns_tp=', ns_tp, '  multharm=', multharm
    print '(A,I5,A,I5,A,I5)', ' Grid      : n_r=', n_r, '  n_th=', n_th, '  n_ph=', n_ph

    ! ── initialise all SIMPLE field representations ───────────────────────────
    call init_vmec(trim(wout_file), ns_s, ns_tp, multharm, fper)

    print '(A,I4,A,I4)', ' VMEC surfaces: ns=', ns, '  nper=', nper

    call get_boozer_coordinates
    call get_canonical_coordinates

    phi_period = twopi / real(nper, dp)

    ! Build grid coordinate arrays (same for all passes)
    k = 0
    do i_r = 1, n_r
        r = 0.1_dp + 0.8_dp * real(i_r - 1, dp) / real(max(n_r - 1, 1), dp)
        do i_th = 1, n_th
            theta = twopi * real(i_th - 1, dp) / real(n_th, dp)
            do i_ph = 1, n_ph
                phi = phi_period * real(i_ph - 1, dp) / real(n_ph, dp)
                k = k + 1
                r_arr(k) = r; theta_arr(k) = theta; phi_arr(k) = phi
            end do
        end do
    end do

    ! ── pass 0: Fourier reference — evaluate B from R,Z,lambda Fourier series ─
    print '(A)', ' Evaluating B from R,Z,lambda Fourier series ...'
    call evaluate_rz_bmod(trim(wout_file), n_points, r_arr, theta_arr, phi_arr, &
                          bmod_ref, R_arr_phys, Z_arr_phys)

    ! ── pass 1: VMEC splines ─────────────────────────────────────────────────
    call init_magfie(VMEC)
    max_err_vmec_spl = 0.0_dp; sum_err_vmec_spl = 0.0_dp; sum_weight = 0.0_dp
    max_Rerr_vmec = 0.0_dp; max_Zerr_vmec = 0.0_dp
    do k = 1, n_points
        s = r_arr(k)**2
        x = [s, theta_arr(k), phi_arr(k)]
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bmod_vmec_arr(k) = bmod
        sqrtg_vmec(k) = abs(sqrtg)
        call splint_vmec_data(s, theta_arr(k), phi_arr(k), &
            A_phi_w, A_theta_w, dA_phi_w, dA_theta_w, aiota_w, &
            Rw, Zw, aw, dRds, dRdt, dRdp, dZds, dZdt, dZdp, dlds, dldt, dldp)
        R_vmec_spl(k) = Rw;  Z_vmec_spl(k) = Zw
        max_Rerr_vmec = max(max_Rerr_vmec, abs(Rw - R_arr_phys(k)) / max(R_arr_phys(k), 1.0e-30_dp))
        max_Zerr_vmec = max(max_Zerr_vmec, abs(Zw - Z_arr_phys(k)) / max(Rw, 1.0e-30_dp))
        err = abs(bmod - bmod_ref(k)) / max(bmod_ref(k), 1.0e-30_dp)
        max_err_vmec_spl = max(max_err_vmec_spl, err)
        sum_err_vmec_spl = sum_err_vmec_spl + err * sqrtg_vmec(k)
        sum_weight       = sum_weight + sqrtg_vmec(k)
    end do
    vol_avg_err_vmec_spl = sum_err_vmec_spl / max(sum_weight, 1.0e-30_dp)

    ! ── pass 2: Boozer splines ───────────────────────────────────────────────
    call init_magfie(BOOZER)
    max_err_boozer = 0.0_dp; sum_err_boozer = 0.0_dp; sum_weight = 0.0_dp
    max_Rerr_boozer = 0.0_dp; max_Zerr_boozer = 0.0_dp
    do k = 1, n_points
        s = r_arr(k)**2
        call vmec_to_boozer(s, theta_arr(k), phi_arr(k), theta_B, phi_B)
        x = [s, theta_B, phi_B]
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bmod_boozer_arr(k) = bmod
        ! Roundtrip: Boozer → VMEC → R,Z
        call boozer_to_vmec(s, theta_B, phi_B, theta_v_back, phi_v_back)
        call splint_vmec_data(s, theta_v_back, phi_v_back, &
            A_phi_w, A_theta_w, dA_phi_w, dA_theta_w, aiota_w, &
            Rw, Zw, aw, dRds, dRdt, dRdp, dZds, dZdt, dZdp, dlds, dldt, dldp)
        R_boozer_spl(k) = Rw;  Z_boozer_spl(k) = Zw
        max_Rerr_boozer = max(max_Rerr_boozer, abs(Rw - R_arr_phys(k)) / max(R_arr_phys(k), 1.0e-30_dp))
        max_Zerr_boozer = max(max_Zerr_boozer, abs(Zw - Z_arr_phys(k)) / max(Rw, 1.0e-30_dp))
        err = abs(bmod - bmod_ref(k)) / max(bmod_ref(k), 1.0e-30_dp)
        max_err_boozer = max(max_err_boozer, err)
        sum_err_boozer = sum_err_boozer + err * sqrtg_vmec(k)
        sum_weight     = sum_weight + sqrtg_vmec(k)
    end do
    vol_avg_err_boozer = sum_err_boozer / max(sum_weight, 1.0e-30_dp)

    ! ── pass 3: canonical (CANFLUX) splines ──────────────────────────────────
    call init_magfie(CANFLUX)
    max_err_can = 0.0_dp; sum_err_can = 0.0_dp; sum_weight = 0.0_dp
    max_Rerr_can = 0.0_dp; max_Zerr_can = 0.0_dp
    do k = 1, n_points
        s = r_arr(k)**2
        call vmec_to_can(s, theta_arr(k), phi_arr(k), theta_c, phi_c)
        x = [s, theta_c, phi_c]
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bmod_can_arr(k) = bmod
        ! Roundtrip: canonical → VMEC → R,Z
        call can_to_vmec(s, theta_c, phi_c, theta_v_back, phi_v_back)
        call splint_vmec_data(s, theta_v_back, phi_v_back, &
            A_phi_w, A_theta_w, dA_phi_w, dA_theta_w, aiota_w, &
            Rw, Zw, aw, dRds, dRdt, dRdp, dZds, dZdt, dZdp, dlds, dldt, dldp)
        R_can_spl(k) = Rw;  Z_can_spl(k) = Zw
        max_Rerr_can = max(max_Rerr_can, abs(Rw - R_arr_phys(k)) / max(R_arr_phys(k), 1.0e-30_dp))
        max_Zerr_can = max(max_Zerr_can, abs(Zw - Z_arr_phys(k)) / max(Rw, 1.0e-30_dp))
        err = abs(bmod - bmod_ref(k)) / max(bmod_ref(k), 1.0e-30_dp)
        max_err_can = max(max_err_can, err)
        sum_err_can = sum_err_can + err * sqrtg_vmec(k)
        sum_weight  = sum_weight + sqrtg_vmec(k)
    end do
    vol_avg_err_can = sum_err_can / max(sum_weight, 1.0e-30_dp)

    ! ── report ───────────────────────────────────────────────────────────────
    print '(A)', ''
    print '(A,I8)', ' Total grid points : ', n_points
    print '(A)', ''
    print '(A)', ' Relative |Bmod_rep - Bmod_Fourier| / Bmod_Fourier'
    print '(A)', ' ---------------------------------------------------'
    print '(A)', '                         Maximum     Volume-averaged'
    print '(A,2ES18.6)', ' VMEC splines  :  ', max_err_vmec_spl, vol_avg_err_vmec_spl
    print '(A,2ES18.6)', ' Boozer        :  ', max_err_boozer,   vol_avg_err_boozer
    print '(A,2ES18.6)', ' Canonical     :  ', max_err_can,      vol_avg_err_can
    print '(A)', ''

    print '(A)', ' R/Z consistency: |R_spl - R_Fourier|/R_Fourier  and  |Z_spl - Z_Fourier| [cm]'
    print '(A)', ' ------------------------------------------------------------------------------'
    print '(A)', '                    max(dR/R)       max(dZ/R)'
    print '(A,2ES18.6)', ' VMEC splines  :  ', max_Rerr_vmec,   max_Zerr_vmec
    print '(A,2ES18.6)', ' Boozer (rt)   :  ', max_Rerr_boozer, max_Zerr_boozer
    print '(A,2ES18.6)', ' Canonical (rt):  ', max_Rerr_can,    max_Zerr_can
    print '(A)', ' (rt) = roundtrip via inverse coordinate transform'
    print '(A)', ''

    if (max_err_vmec_spl > 1.0e-3_dp) then
        print '(A)', ' WARNING: VMEC spline error > 0.1% -- consider more VMEC surfaces (ns)'
    end if
    if (max_err_boozer > max_err_vmec_spl*10.0_dp .or. max_err_can > max_err_vmec_spl*10.0_dp) then
        print '(A)', ' WARNING: Boozer/canonical error >> VMEC spline error --'
        print '(A)', '          consider increasing multharm in simple.in'
    end if
    print '(A)', '================================================='

    ! ── dump all data to file for post-processing ─────────────────────────────
    call dump_field_data('field_representation_data.dat', n_points, &
                         n_r, n_th, n_ph, phi_period, &
                         r_arr, theta_arr, phi_arr, R_arr_phys, Z_arr_phys, &
                         bmod_ref, bmod_vmec_arr, bmod_boozer_arr, bmod_can_arr, &
                         R_vmec_spl, Z_vmec_spl, R_boozer_spl, Z_boozer_spl, &
                         R_can_spl, Z_can_spl)

    deallocate(bmod_ref, sqrtg_vmec, r_arr, theta_arr, phi_arr)
    deallocate(bmod_vmec_arr, bmod_boozer_arr, bmod_can_arr)
    deallocate(R_arr_phys, Z_arr_phys)
    deallocate(R_vmec_spl, Z_vmec_spl, R_boozer_spl, Z_boozer_spl, R_can_spl, Z_can_spl)

contains

    !> Write all field data to a columnar ASCII file for post-processing.
    !> Columns: r theta phi R[cm] Z[cm] B_Fourier B_VMEC B_Boozer B_Canonical  (B in Gauss)
    subroutine dump_field_data(filename, npts, nr, nth, nph, phiper, &
                               r_in, theta_in, phi_in, Rcyl_four, Zcyl_four, &
                               bmod_four, bmod_vmec, bmod_boozer, bmod_can, &
                               Rcyl_vmec, Zcyl_vmec, Rcyl_boozer, Zcyl_boozer, &
                               Rcyl_can, Zcyl_can)
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: npts, nr, nth, nph
        real(dp),         intent(in) :: phiper
        real(dp),         intent(in) :: r_in(npts), theta_in(npts), phi_in(npts)
        real(dp),         intent(in) :: Rcyl_four(npts), Zcyl_four(npts)
        real(dp),         intent(in) :: bmod_four(npts), bmod_vmec(npts)
        real(dp),         intent(in) :: bmod_boozer(npts), bmod_can(npts)
        real(dp),         intent(in) :: Rcyl_vmec(npts),   Zcyl_vmec(npts)
        real(dp),         intent(in) :: Rcyl_boozer(npts), Zcyl_boozer(npts)
        real(dp),         intent(in) :: Rcyl_can(npts),    Zcyl_can(npts)
        integer :: ipt

        open(unit=42, file=trim(filename), status='replace')
        ! Grid metadata on first line for Python post-processing
        write(42, '(A,I0,A,I0,A,I0,A,ES18.8)') &
            '# n_r=', nr, ' n_th=', nth, ' n_ph=', nph, ' phi_period=', phiper
        write(42, '(A)') &
            '# r  theta  phi  R_four  Z_four' // &
            '  B_Fourier  B_VMEC  B_Boozer  B_Can' // &
            '  R_VMEC  Z_VMEC  R_Boozer  Z_Boozer  R_Can  Z_Can'
        do ipt = 1, npts
            write(42, '(15ES18.8)') r_in(ipt), theta_in(ipt), phi_in(ipt), &
                                    Rcyl_four(ipt), Zcyl_four(ipt), &
                                    bmod_four(ipt), bmod_vmec(ipt), &
                                    bmod_boozer(ipt), bmod_can(ipt), &
                                    Rcyl_vmec(ipt),   Zcyl_vmec(ipt), &
                                    Rcyl_boozer(ipt), Zcyl_boozer(ipt), &
                                    Rcyl_can(ipt),    Zcyl_can(ipt)
        end do
        close(42)
        print '(A,A)', ' Written: ', trim(filename)

    end subroutine dump_field_data


    !> Evaluate |B| at each point (r=sqrt(s), theta, phi) directly from the
    !> VMEC R(s,theta,phi) and Z(s,theta,phi) Fourier representation using the
    !> magnetic field ansatz:
    !>
    !>   B^theta = -dA_phi/ds / sqg,   B^phi = dA_theta/ds / sqg
    !>   |B|^2   = B^theta * B_theta + B^phi * B_phi
    !>
    !> R, Z, lambda and all partial derivatives are evaluated by direct Fourier
    !> summation (rmnc, zmns, lmns from wout.nc). Radial derivatives use finite
    !> differences between adjacent full-mesh VMEC surfaces. The vector potential
    !> derivatives dA_phi/ds and dA_theta/ds are obtained from the existing
    !> 1-D spline (highly accurate) via splint_vmec_data.
    subroutine evaluate_rz_bmod(filename, npts, r_in, theta_in, phi_in, bmod_out, R_out, Z_out)
        use spline_vmec_sub, only: splint_vmec_data, compute_field_components
        use new_vmec_stuff_mod, only: vmec_RZ_scale

        character(len=*), intent(in)  :: filename
        integer,          intent(in)  :: npts
        real(dp),         intent(in)  :: r_in(npts), theta_in(npts), phi_in(npts)
        real(dp),         intent(out) :: bmod_out(npts), R_out(npts), Z_out(npts)

        ! Fourier data read from wout.nc
        integer :: ncid, varid, dimid
        integer :: ns_nc, mnmax_nc, lasym_int
        logical :: lasym
        real(dp), allocatable :: xm(:), xn(:)
        real(dp), allocatable :: rmnc_nc(:,:), zmns_nc(:,:), lmns_nc(:,:)
        real(dp), allocatable :: rmns_nc(:,:), zmnc_nc(:,:), lmnc_nc(:,:)
        real(dp), allocatable :: phi_nc(:)        ! toroidal flux at each surface
        real(dp), allocatable :: s_nc(:)          ! normalised toroidal flux
        real(dp), allocatable :: almns(:,:), almnc(:,:)  ! full-mesh lambda

        ! Per-point workspace
        real(dp) :: s_pt, rho_pt, rho_lo, rho_hi, w1, w2, drho_inv
        real(dp) :: phase, cph, sph
        real(dp) :: R_lo, R_hi, Z_lo, Z_hi, lam_lo, lam_hi
        real(dp) :: dR_dt_lo, dR_dt_hi, dZ_dt_lo, dZ_dt_hi, dl_dt_lo, dl_dt_hi
        real(dp) :: dR_dp_lo, dR_dp_hi, dZ_dp_lo, dZ_dp_hi, dl_dp_lo, dl_dp_hi
        real(dp) :: R, Z, alam, dR_ds, dZ_ds, dl_ds
        real(dp) :: dR_dt, dZ_dt, dl_dt, dR_dp, dZ_dp, dl_dp
        ! Vector potential from existing 1-D spline (used as-is)
        real(dp) :: A_phi_spl, A_theta_spl, dA_phi_ds, dA_theta_ds, aiota_spl
        real(dp) :: R_spl, Z_spl, alam_spl
        real(dp) :: dR_ds_spl, dR_dt_spl, dR_dp_spl
        real(dp) :: dZ_ds_spl, dZ_dt_spl, dZ_dp_spl
        real(dp) :: dl_ds_spl, dl_dt_spl, dl_dp_spl
        ! Field components
        real(dp) :: sqg, Bctr_th, Bctr_ph, Bcov_s, Bcov_th, Bcov_ph
        integer :: i, ipt, is_lo, mn

        ! ── read wout.nc ─────────────────────────────────────────────────────
        call check_nc(nf90_open(trim(filename), nf90_nowrite, ncid))

        call check_nc(nf90_inq_dimid(ncid, 'radius', dimid))
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=ns_nc))
        call check_nc(nf90_inq_dimid(ncid, 'mn_mode', dimid))
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=mnmax_nc))

        call check_nc(nf90_inq_varid(ncid, 'lasym__logical__', varid))
        call check_nc(nf90_get_var(ncid, varid, lasym_int))
        lasym = (lasym_int == 1)

        allocate(xm(mnmax_nc), xn(mnmax_nc))
        allocate(rmnc_nc(mnmax_nc, ns_nc), zmns_nc(mnmax_nc, ns_nc))
        allocate(lmns_nc(mnmax_nc, ns_nc))
        allocate(phi_nc(ns_nc))

        call check_nc(nf90_inq_varid(ncid, 'xm',   varid)); call check_nc(nf90_get_var(ncid, varid, xm))
        call check_nc(nf90_inq_varid(ncid, 'xn',   varid)); call check_nc(nf90_get_var(ncid, varid, xn))
        call check_nc(nf90_inq_varid(ncid, 'rmnc', varid)); call check_nc(nf90_get_var(ncid, varid, rmnc_nc))
        call check_nc(nf90_inq_varid(ncid, 'zmns', varid)); call check_nc(nf90_get_var(ncid, varid, zmns_nc))
        call check_nc(nf90_inq_varid(ncid, 'lmns', varid)); call check_nc(nf90_get_var(ncid, varid, lmns_nc))
        call check_nc(nf90_inq_varid(ncid, 'phi',  varid)); call check_nc(nf90_get_var(ncid, varid, phi_nc))

        if (lasym) then
            allocate(rmns_nc(mnmax_nc, ns_nc), zmnc_nc(mnmax_nc, ns_nc), lmnc_nc(mnmax_nc, ns_nc))
            call check_nc(nf90_inq_varid(ncid, 'rmns', varid)); call check_nc(nf90_get_var(ncid, varid, rmns_nc))
            call check_nc(nf90_inq_varid(ncid, 'zmnc', varid)); call check_nc(nf90_get_var(ncid, varid, zmnc_nc))
            call check_nc(nf90_inq_varid(ncid, 'lmnc', varid)); call check_nc(nf90_get_var(ncid, varid, lmnc_nc))
        end if

        call check_nc(nf90_close(ncid))

        ! Apply same unit scaling as vmecin (vmecinm_m.f90): m -> cm
        ! fac_r = fac_r0 * vmec_RZ_scale = 100 * vmec_RZ_scale
        ! Lambda (lmns) is dimensionless — no scaling needed.
        rmnc_nc = rmnc_nc * 1.0e2_dp * vmec_RZ_scale
        zmns_nc = zmns_nc * 1.0e2_dp * vmec_RZ_scale
        if (lasym) then
            rmns_nc = rmns_nc * 1.0e2_dp * vmec_RZ_scale
            zmnc_nc = zmnc_nc * 1.0e2_dp * vmec_RZ_scale
        end if

        ! Normalised toroidal flux at each full-mesh surface
        allocate(s_nc(ns_nc))
        s_nc = abs(phi_nc) / max(abs(phi_nc(ns_nc)), 1.0e-30_dp)

        ! Convert lambda from half-mesh (lmns) to full-mesh (almns),
        ! matching the convention in vmecin (vmecinm_m.f90):
        !   almns(:,0) = 0
        !   almns(:,i) = 0.5*(lmns(:,i+1) + lmns(:,i))  for i=1..ns-2
        !   almns(:,ns-1) = lmns(:,ns-1) + 0.5*(lmns(:,ns-1) - lmns(:,ns-2))
        allocate(almns(mnmax_nc, 0:ns_nc-1))
        ! vmecin uses 0-based lmns(:,j), so lmns(:,j) = lmns_nc(:,j+1).
        ! vmecin formula: almns(:,i) = 0.5*(lmns(:,i+1)+lmns(:,i))
        !                            = 0.5*(lmns_nc(:,i+2)+lmns_nc(:,i+1))
        almns(:, 0) = 0.0_dp
        do i = 1, ns_nc - 2
            almns(:, i) = 0.5_dp * (lmns_nc(:, i+2) + lmns_nc(:, i+1))
        end do
        almns(:, ns_nc-1) = lmns_nc(:, ns_nc) + 0.5_dp * &
                            (lmns_nc(:, ns_nc) - lmns_nc(:, ns_nc-1))
        if (lasym) then
            allocate(almnc(mnmax_nc, 0:ns_nc-1))
            almnc(:, 0) = 0.0_dp
            do i = 1, ns_nc - 2
                almnc(:, i) = 0.5_dp * (lmnc_nc(:, i+2) + lmnc_nc(:, i+1))
            end do
            almnc(:, ns_nc-1) = lmnc_nc(:, ns_nc) + 0.5_dp * &
                                 (lmnc_nc(:, ns_nc) - lmnc_nc(:, ns_nc-1))
        end if

        ! ── evaluate at each query point ──────────────────────────────────────
        do ipt = 1, npts
            s_pt  = r_in(ipt) * r_in(ipt)
            rho_pt = r_in(ipt)

            ! Find bounding full-mesh surfaces in s
            is_lo = 1
            do i = 1, ns_nc - 1
                if (s_nc(i) <= s_pt) is_lo = i
            end do
            is_lo = min(is_lo, ns_nc - 1)

            ! Work in rho = sqrt(s) space, matching splint_vmec_data which
            ! builds splines in rho and applies dX_ds = 0.5*dX_drho/rho (line 593-595).
            ! s-space finite differences are inaccurate near the axis (Ds -> 0).
            rho_lo   = sqrt(s_nc(is_lo))
            rho_hi   = sqrt(s_nc(is_lo+1))
            drho_inv = 1.0_dp / max(rho_hi - rho_lo, 1.0e-30_dp)
            w2 = (rho_pt - rho_lo) * drho_inv
            w2 = max(0.0_dp, min(1.0_dp, w2))
            w1 = 1.0_dp - w2

            ! Initialise accumulators for both bounding surfaces
            R_lo=0; R_hi=0; Z_lo=0; Z_hi=0; lam_lo=0; lam_hi=0
            dR_dt_lo=0; dR_dt_hi=0; dZ_dt_lo=0; dZ_dt_hi=0
            dl_dt_lo=0; dl_dt_hi=0
            dR_dp_lo=0; dR_dp_hi=0; dZ_dp_lo=0; dZ_dp_hi=0
            dl_dp_lo=0; dl_dp_hi=0

            ! Direct Fourier summation at the two bounding surfaces
            do mn = 1, mnmax_nc
                phase = xm(mn)*theta_in(ipt) - xn(mn)*phi_in(ipt)
                cph = cos(phase)
                sph = sin(phase)

                R_lo  = R_lo  + rmnc_nc(mn, is_lo)   * cph
                R_hi  = R_hi  + rmnc_nc(mn, is_lo+1) * cph
                Z_lo  = Z_lo  + zmns_nc(mn, is_lo)   * sph
                Z_hi  = Z_hi  + zmns_nc(mn, is_lo+1) * sph
                lam_lo = lam_lo + almns(mn, is_lo-1)  * sph
                lam_hi = lam_hi + almns(mn, is_lo)    * sph

                ! Angular derivatives: d/dtheta, d/dphi
                dR_dt_lo  = dR_dt_lo  - xm(mn)*rmnc_nc(mn, is_lo)   * sph
                dR_dt_hi  = dR_dt_hi  - xm(mn)*rmnc_nc(mn, is_lo+1) * sph
                dZ_dt_lo  = dZ_dt_lo  + xm(mn)*zmns_nc(mn, is_lo)   * cph
                dZ_dt_hi  = dZ_dt_hi  + xm(mn)*zmns_nc(mn, is_lo+1) * cph
                dl_dt_lo  = dl_dt_lo  + xm(mn)*almns(mn, is_lo-1)   * cph
                dl_dt_hi  = dl_dt_hi  + xm(mn)*almns(mn, is_lo)     * cph

                dR_dp_lo  = dR_dp_lo  + xn(mn)*rmnc_nc(mn, is_lo)   * sph
                dR_dp_hi  = dR_dp_hi  + xn(mn)*rmnc_nc(mn, is_lo+1) * sph
                dZ_dp_lo  = dZ_dp_lo  - xn(mn)*zmns_nc(mn, is_lo)   * cph
                dZ_dp_hi  = dZ_dp_hi  - xn(mn)*zmns_nc(mn, is_lo+1) * cph
                dl_dp_lo  = dl_dp_lo  - xn(mn)*almns(mn, is_lo-1)   * cph
                dl_dp_hi  = dl_dp_hi  - xn(mn)*almns(mn, is_lo)     * cph

                if (lasym) then
                    R_lo  = R_lo  + rmns_nc(mn, is_lo)   * sph
                    R_hi  = R_hi  + rmns_nc(mn, is_lo+1) * sph
                    Z_lo  = Z_lo  + zmnc_nc(mn, is_lo)   * cph
                    Z_hi  = Z_hi  + zmnc_nc(mn, is_lo+1) * cph
                    lam_lo = lam_lo + almnc(mn, is_lo-1) * cph
                    lam_hi = lam_hi + almnc(mn, is_lo)   * cph
                    dR_dt_lo  = dR_dt_lo  + xm(mn)*rmns_nc(mn, is_lo)   * cph
                    dR_dt_hi  = dR_dt_hi  + xm(mn)*rmns_nc(mn, is_lo+1) * cph
                    dZ_dt_lo  = dZ_dt_lo  - xm(mn)*zmnc_nc(mn, is_lo)   * sph
                    dZ_dt_hi  = dZ_dt_hi  - xm(mn)*zmnc_nc(mn, is_lo+1) * sph
                    dl_dt_lo  = dl_dt_lo  - xm(mn)*almnc(mn, is_lo-1)   * sph
                    dl_dt_hi  = dl_dt_hi  - xm(mn)*almnc(mn, is_lo)     * sph
                    dR_dp_lo  = dR_dp_lo  - xn(mn)*rmns_nc(mn, is_lo)   * cph
                    dR_dp_hi  = dR_dp_hi  - xn(mn)*rmns_nc(mn, is_lo+1) * cph
                    dZ_dp_lo  = dZ_dp_lo  + xn(mn)*zmnc_nc(mn, is_lo)   * sph
                    dZ_dp_hi  = dZ_dp_hi  + xn(mn)*zmnc_nc(mn, is_lo+1) * sph
                    dl_dp_lo  = dl_dp_lo  + xn(mn)*almnc(mn, is_lo-1)   * sph
                    dl_dp_hi  = dl_dp_hi  + xn(mn)*almnc(mn, is_lo)     * sph
                end if
            end do

            ! Interpolate to query point and form radial derivatives (d/ds)
            R    = w1*R_lo   + w2*R_hi
            Z    = w1*Z_lo   + w2*Z_hi
            alam = w1*lam_lo + w2*lam_hi
            dR_dt = w1*dR_dt_lo + w2*dR_dt_hi
            dZ_dt = w1*dZ_dt_lo + w2*dZ_dt_hi
            dl_dt = w1*dl_dt_lo + w2*dl_dt_hi
            dR_dp = w1*dR_dp_lo + w2*dR_dp_hi
            dZ_dp = w1*dZ_dp_lo + w2*dZ_dp_hi
            dl_dp = w1*dl_dp_lo + w2*dl_dp_hi

            ! Get dA_phi_ds and dA_theta_ds from the existing 1-D spline — these
            ! are smooth radial functions and well-represented by the spline.
            call splint_vmec_data(s_pt, theta_in(ipt), phi_in(ipt), &
                                  A_phi_spl, A_theta_spl, dA_phi_ds, dA_theta_ds, &
                                  aiota_spl, &
                                  R_spl, Z_spl, alam_spl, &
                                  dR_ds_spl, dR_dt_spl, dR_dp_spl, &
                                  dZ_ds_spl, dZ_dt_spl, dZ_dp_spl, &
                                  dl_ds_spl, dl_dt_spl, dl_dp_spl)

            ! Compute field components using the Fourier R,Z,lambda
            call compute_field_components(R, dR_ds_spl, dR_dt, dR_dp, &
                                          dZ_ds_spl, dZ_dt, dZ_dp, &
                                          dA_theta_ds, dA_phi_ds, &
                                          dl_ds_spl, dl_dt, dl_dp, &
                                          sqg, Bctr_th, Bctr_ph, &
                                          Bcov_s, Bcov_th, Bcov_ph)

            bmod_out(ipt) = sqrt(abs(Bctr_th*Bcov_th + Bctr_ph*Bcov_ph))
            R_out(ipt)    = R
            Z_out(ipt)    = Z
        end do

        deallocate(xm, xn, rmnc_nc, zmns_nc, lmns_nc, phi_nc, s_nc, almns)
        if (lasym) deallocate(rmns_nc, zmnc_nc, lmnc_nc, almnc)

    end subroutine evaluate_rz_bmod


    subroutine check_nc(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            print *, 'NetCDF error: ', trim(nf90_strerror(status))
            error stop 'NetCDF read failed in test_field_representation_error'
        end if
    end subroutine check_nc

end program test_field_representation_error
