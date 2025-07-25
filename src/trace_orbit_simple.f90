program trace_orbit_simple
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp
  use velo_mod, only : isw_field_type
  use util, only : pi, c, e_charge, p_mass, ev
  use simple_main, only : init_field
  use get_can_sub, only : vmec_to_can
  use params, only : params_init, rt0, v0, dtaumin, npoiper2, trace_time, ntimstep, &
                     ntau, dtau, relerr, integmode
  use simple, only : Tracer, init_sympl
  use orbit_symplectic, only : orbit_timestep_sympl_expl_impl_euler, orbit_timestep_sympl_midpoint
  use magfie_sub, only : init_magfie, VMEC
  use netcdf

  implicit none

  ! Variables
  type(Tracer) :: tracy
  double precision, dimension(5) :: z0_vmec, z0_can, z
  character(256) :: vmec_file, output_file
  integer :: kt, ierr, ktau

  ! NetCDF variables
  integer :: ncid, dimid_time, varid_time, varid_s, varid_theta_vmec, varid_phi_vmec
  integer :: varid_theta_can, varid_phi_can, varid_vpar, varid_mu, varid_bmod

  ! Output arrays (allocatable)
  double precision, allocatable :: s_arr(:), theta_vmec_arr(:), phi_vmec_arr(:)
  double precision, allocatable :: theta_can_arr(:), phi_can_arr(:), time_arr(:)
  double precision, allocatable :: v_par_arr(:), mu_arr(:), b_field_arr(:)

  ! Get VMEC file from command line
  if (command_argument_count() < 1) then
    write(*,*) 'Usage: trace_orbit_fixed <vmec_file>'
    stop 1
  endif
  call get_command_argument(1, vmec_file)

  write(*,*) 'Following SIMPLE control flow exactly...'

  ! Step 1: Set configuration parameters (instead of read_config)
  netcdffile = trim(vmec_file)
  ns_s = 5
  ns_tp = 5
  multharm = 5
  isw_field_type = 2  ! Boozer coordinates
  trace_time = 1d-3   ! 0.001 seconds (10x shorter)
  ntimstep = 10001    ! From SIMPLE examples
  npoiper2 = 256
  relerr = 1d-13
  integmode = 2       ! Standard symplectic integration

  ! Allocate arrays after ntimstep is set
  allocate(s_arr(ntimstep), theta_vmec_arr(ntimstep), phi_vmec_arr(ntimstep))
  allocate(theta_can_arr(ntimstep), phi_can_arr(ntimstep), time_arr(ntimstep))
  allocate(v_par_arr(ntimstep), mu_arr(ntimstep), b_field_arr(ntimstep))

  write(*,*) '  Configuration set: trace_time =', trace_time, 's'

  ! Step 2: Initialize field (exactly like SIMPLE main)
  call init_field(tracy, netcdffile, ns_s, ns_tp, multharm, integmode)

  ! Step 3: Initialize parameters (exactly like SIMPLE main)
  call params_init()

  ! Step 4: Initialize magfie TWICE like SIMPLE does
  call init_magfie(VMEC)
  call init_magfie(isw_field_type)

  write(*,*) '  SIMPLE initialization completed!'
  write(*,*) '    v0 =', v0, 'cm/s'
  write(*,*) '    RT0 =', rt0, 'cm'
  write(*,*) '    dtaumin =', dtaumin
  write(*,*) '    ntau =', ntau
  write(*,*) '    dtau =', dtau


  ! Step 5: Set initial conditions from files in initial_condition/
  open(unit=10, file='../initial_condition/s_booz.txt', status='old')
  read(10,*) z0_vmec(1)
  close(10)
  open(unit=11, file='../initial_condition/theta_vmec.txt', status='old')
  read(11,*) z0_vmec(2)
  close(11)
  open(unit=12, file='../initial_condition/phi_vmec.txt', status='old')
  read(12,*) z0_vmec(3)
  close(12)
  z0_vmec(4) = 1.0d0      ! v/v_th
  z0_vmec(5) = 0.0d0      ! v_par/v

  ! Convert to canonical coordinates (ref_to_can equivalent)
  z0_can(1) = z0_vmec(1)
  call vmec_to_can(z0_vmec(1), z0_vmec(2), z0_vmec(3), z0_can(2), z0_can(3))
  z0_can(4) = z0_vmec(4)
  z0_can(5) = z0_vmec(5)
  ! Note: SIMPLE has ref_to_can conversion here
  z = z0_can  ! Working coordinates

  write(*,*) '  Initial conditions: s=', z(1), ', theta=', z(2), ', phi=', z(3)
  write(*,*) '    p0=', z(4), ', vpar=', z(5)

  ! Step 6: Initialize symplectic integrator AFTER all setup (like trace_orbit)
  call init_sympl(tracy%si, tracy%f, z, dtaumin, dtaumin, relerr, integmode)

  write(*,*) '  Symplectic integrator initialized'
  write(*,*) '  Starting orbit integration with', ntimstep, 'timesteps...'

  ! Store initial values
  time_arr(1) = 0.0d0
  s_arr(1) = z(1)
  theta_can_arr(1) = z(2)
  phi_can_arr(1) = z(3)
  theta_vmec_arr(1) = z0_vmec(2)  ! Use original VMEC coordinates
  phi_vmec_arr(1) = z0_vmec(3)
  v_par_arr(1) = tracy%f%vpar
  mu_arr(1) = tracy%f%mu
  b_field_arr(1) = tracy%f%bmod

  ! Time integration loop following SIMPLE's macrostep pattern
  do kt = 2, ntimstep
    ! Inner loop over ntau steps (like macrostep)
    do ktau = 1, ntau
      ! Call the same integration routine as SIMPLE
      !call orbit_timestep_sympl_expl_impl_euler(tracy%si, tracy%f, ierr)
      call orbit_timestep_sympl_midpoint(tracy%si, tracy%f, ierr)
      if (ierr /= 0) then
        write(*,*) 'Error in orbit integration at step:', kt-1, 'ktau:', ktau
        exit
      endif
    enddo

    if (ierr /= 0) exit

    ! Store time as equally spaced values from 0 to trace_time
    time_arr(kt) = (kt-1) * (trace_time / (ntimstep-1))

    ! Store canonical coordinates from integrator
    s_arr(kt) = tracy%si%z(1)
    theta_can_arr(kt) = tracy%si%z(2)
    phi_can_arr(kt) = tracy%si%z(3)

    ! Convert to VMEC coordinates (can_to_ref equivalent)
    ! For now, use canonical coordinates as approximation
    theta_vmec_arr(kt) = theta_can_arr(kt)
    phi_vmec_arr(kt) = phi_can_arr(kt)

    ! Store field and velocity information
    v_par_arr(kt) = tracy%f%vpar
    mu_arr(kt) = tracy%f%mu
    b_field_arr(kt) = tracy%f%bmod
  enddo

  ! Write output to NetCDF file (matching xarray format)
  output_file = '../run/trace_orbit_simple.nc'

  ! Create NetCDF file
  ierr = nf90_create(output_file, NF90_CLOBBER, ncid)
  if (ierr /= NF90_NOERR) then
    write(*,*) 'Error creating NetCDF file:', trim(nf90_strerror(ierr))
    stop 1
  endif

  ! Define dimensions
  ierr = nf90_def_dim(ncid, 'time', ntimstep, dimid_time)

  ! Define coordinate variables
  ierr = nf90_def_var(ncid, 'time', NF90_DOUBLE, [dimid_time], varid_time)
  ierr = nf90_put_att(ncid, varid_time, 'units', 'seconds')
  ierr = nf90_put_att(ncid, varid_time, 'long_name', 'time')

  ! Define data variables (matching xarray structure)
  ierr = nf90_def_var(ncid, 's', NF90_DOUBLE, [dimid_time], varid_s)
  ierr = nf90_put_att(ncid, varid_s, 'long_name', 'normalized toroidal flux')

  ierr = nf90_def_var(ncid, 'theta_vmec', NF90_DOUBLE, [dimid_time], varid_theta_vmec)
  ierr = nf90_put_att(ncid, varid_theta_vmec, 'long_name', 'VMEC poloidal angle')

  ierr = nf90_def_var(ncid, 'phi_vmec', NF90_DOUBLE, [dimid_time], varid_phi_vmec)
  ierr = nf90_put_att(ncid, varid_phi_vmec, 'long_name', 'VMEC toroidal angle')

  ierr = nf90_def_var(ncid, 'theta_can', NF90_DOUBLE, [dimid_time], varid_theta_can)
  ierr = nf90_put_att(ncid, varid_theta_can, 'long_name', 'canonical poloidal angle')

  ierr = nf90_def_var(ncid, 'phi_can', NF90_DOUBLE, [dimid_time], varid_phi_can)
  ierr = nf90_put_att(ncid, varid_phi_can, 'long_name', 'canonical toroidal angle')

  ierr = nf90_def_var(ncid, 'v_par', NF90_DOUBLE, [dimid_time], varid_vpar)
  ierr = nf90_put_att(ncid, varid_vpar, 'long_name', 'parallel velocity')

  ierr = nf90_def_var(ncid, 'mu', NF90_DOUBLE, [dimid_time], varid_mu)
  ierr = nf90_put_att(ncid, varid_mu, 'long_name', 'magnetic moment')

  ierr = nf90_def_var(ncid, 'B', NF90_DOUBLE, [dimid_time], varid_bmod)
  ierr = nf90_put_att(ncid, varid_bmod, 'long_name', 'magnetic field strength')

  ! End define mode
  ierr = nf90_enddef(ncid)

  ! Write data
  ierr = nf90_put_var(ncid, varid_time, time_arr)
  ierr = nf90_put_var(ncid, varid_s, s_arr)
  ierr = nf90_put_var(ncid, varid_theta_vmec, theta_vmec_arr)
  ierr = nf90_put_var(ncid, varid_phi_vmec, phi_vmec_arr)
  ierr = nf90_put_var(ncid, varid_theta_can, theta_can_arr)
  ierr = nf90_put_var(ncid, varid_phi_can, phi_can_arr)
  ierr = nf90_put_var(ncid, varid_vpar, v_par_arr)
  ierr = nf90_put_var(ncid, varid_mu, mu_arr)
  ierr = nf90_put_var(ncid, varid_bmod, b_field_arr)

  ! Close NetCDF file
  ierr = nf90_close(ncid)

  write(*,*) ''
  write(*,*) 'Orbit tracing completed successfully!'
  write(*,*) '  Final s: ', s_arr(ntimstep)
  write(*,*) '  s range: ', minval(s_arr), ' - ', maxval(s_arr)
  write(*,*) '  Energy conservation: ', sqrt(v_par_arr(ntimstep)**2 + 2*mu_arr(ntimstep)*b_field_arr(ntimstep)) / &
                                         sqrt(v_par_arr(1)**2 + 2*mu_arr(1)*b_field_arr(1))
  write(*,*) '  Mu conservation: ', mu_arr(ntimstep) / mu_arr(1)
  write(*,*) '  NetCDF output written to: ', trim(output_file)

end program trace_orbit_simple
