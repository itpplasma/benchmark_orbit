program trace_orbit
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp
  use velo_mod, only : isw_field_type
  use parmot_mod, only : rmu, ro0
  use util, only : pi, c, e_charge, p_mass, ev
  use spline_vmec_sub
  use get_can_sub, only: vmec_to_can, can_to_vmec
  use orbit_symplectic, only : orbit_timestep_sympl_expl_impl_euler, SymplecticIntegrator
  use simple, only : Tracer, init_sympl
  use simple_main, only : init_field
  use params, only : params_init, rt0, v0, dtaumin, npoiper2
  use field_can_mod, only : FieldCan
  use netcdf

  implicit none

  ! Parameters
  double precision, parameter :: E_alpha = 3.5d6  ! eV
  integer, parameter :: n_d = 4   ! mass number (alpha)
  integer, parameter :: n_e = 2   ! charge number (alpha)
  integer, parameter :: nt = 10000  ! number of timesteps
  
  ! Variables
  type(Tracer) :: tracy
  type(SymplecticIntegrator) :: si
  type(FieldCan) :: f
  double precision, dimension(5) :: z0_vmec, z0_can
  double precision :: rlarm, bmod_ref, facE_al
  integer :: kt, ierr
  character(256) :: vmec_file, output_file
  
  ! NetCDF variables
  integer :: ncid, dimid_time, varid_time, varid_s, varid_theta_vmec, varid_phi_vmec
  integer :: varid_theta_can, varid_phi_can, varid_vpar, varid_mu, varid_bmod
  integer :: varid_r_cyl, varid_z_cyl
  
  ! Output arrays
  double precision, dimension(nt) :: s_arr, theta_vmec_arr, phi_vmec_arr
  double precision, dimension(nt) :: theta_can_arr, phi_can_arr, time_arr
  double precision, dimension(nt) :: v_par_arr, mu_arr, b_field_arr
  double precision, dimension(nt) :: r_cyl_arr, z_cyl_arr

  ! Get VMEC file from command line
  if (command_argument_count() < 1) then
    write(*,*) 'Usage: trace_orbit <vmec_file>'
    stop 1
  endif
  call get_command_argument(1, vmec_file)

  ! Set parameters exactly like SIMPLE examples
  netcdffile = trim(vmec_file)
  ns_s = 3
  ns_tp = 3
  multharm = 3
  isw_field_type = 2  ! Boozer coordinates
  
  ! Initialize field and parameters
  call init_field(tracy, netcdffile, ns_s, ns_tp, multharm, 2)  ! mode 2 = Boozer
  call params_init()
  
  ! Reference field strength (from SIMPLE examples)
  bmod_ref = 10000.0d0  ! Gauss
  facE_al = 1.0d0
  
  ! Calculate particle parameters
  rlarm = v0 * n_d * p_mass * c / (n_e * e_charge * bmod_ref)
  ro0 = rlarm
  rmu = 1d8  ! Large inverse relativistic temperature
  
  ! Initial conditions: s=0.5, theta=0.5, phi=0.5, v/v_th=1.0, v_par/v=0.5
  z0_vmec(1) = 0.5d0      ! s
  z0_vmec(2) = 0.5d0      ! theta_vmec
  z0_vmec(3) = 0.5d0      ! phi_vmec  
  z0_vmec(4) = 1.0d0      ! v/v_th
  z0_vmec(5) = 0.5d0      ! v_par/v
  
  ! Convert VMEC to canonical coordinates
  z0_can = z0_vmec
  call vmec_to_can(z0_vmec(1), z0_vmec(2), z0_vmec(3), z0_can(2), z0_can(3))
  
  ! Initialize symplectic integrator
  call init_sympl(tracy%si, tracy%f, z0_can, dtaumin, dtaumin, 1d-13, tracy%integmode)
  
  write(*,*) 'Orbit tracing with SIMPLE Fortran'
  write(*,*) '  VMEC file: ', trim(vmec_file)
  write(*,*) '  Field mode: ', isw_field_type, ' (Boozer coordinates)'
  write(*,*) '  Initial s, theta, phi: ', z0_vmec(1:3)
  write(*,*) '  v_par/v: ', z0_vmec(5)
  write(*,*) '  Timesteps: ', nt
  write(*,*) '  dtaumin: ', dtaumin
  write(*,*) '  RT0: ', rt0, ' cm'
  write(*,*) '  v0: ', v0, ' cm/s'
  
  ! Store initial values
  time_arr(1) = 0.0d0
  s_arr(1) = z0_vmec(1)
  theta_vmec_arr(1) = z0_vmec(2)
  phi_vmec_arr(1) = z0_vmec(3)
  theta_can_arr(1) = z0_can(2)
  phi_can_arr(1) = z0_can(3)
  v_par_arr(1) = tracy%f%vpar
  mu_arr(1) = tracy%f%mu
  b_field_arr(1) = tracy%f%bmod
  
  ! Time integration loop
  do kt = 2, nt
    ! Advance one timestep
    call orbit_timestep_sympl_expl_impl_euler(tracy%si, tracy%f, ierr)
    if (ierr /= 0) then
      write(*,*) 'Error in orbit integration at step:', kt-1
      exit
    endif
    
    ! Store time
    time_arr(kt) = (kt-1) * dtaumin
    
    ! Store canonical coordinates
    s_arr(kt) = tracy%si%z(1)
    theta_can_arr(kt) = tracy%si%z(2)
    phi_can_arr(kt) = tracy%si%z(3)
    
    ! Convert canonical to VMEC coordinates
    call can_to_vmec(s_arr(kt), theta_can_arr(kt), phi_can_arr(kt), &
                     theta_vmec_arr(kt), phi_vmec_arr(kt))
    
    ! Store field and velocity information
    v_par_arr(kt) = tracy%f%vpar
    mu_arr(kt) = tracy%f%mu
    b_field_arr(kt) = tracy%f%bmod
  enddo
  
  ! Write output to NetCDF file (matching xarray format)
  output_file = 'build/orbit_trace_fortran.nc'
  
  ! Create NetCDF file
  ierr = nf90_create(output_file, NF90_CLOBBER, ncid)
  if (ierr /= NF90_NOERR) then
    write(*,*) 'Error creating NetCDF file:', trim(nf90_strerror(ierr))
    stop 1
  endif
  
  ! Define dimensions
  ierr = nf90_def_dim(ncid, 'time', nt, dimid_time)
  
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
  write(*,*) 'Orbit tracing completed!'
  write(*,*) '  Final s: ', s_arr(nt)
  write(*,*) '  s range: ', minval(s_arr), ' - ', maxval(s_arr)
  write(*,*) '  Energy conservation: ', sqrt(v_par_arr(nt)**2 + 2*mu_arr(nt)*b_field_arr(nt)) / &
                                         sqrt(v_par_arr(1)**2 + 2*mu_arr(1)*b_field_arr(1))
  write(*,*) '  Mu conservation: ', mu_arr(nt) / mu_arr(1)
  write(*,*) '  NetCDF output written to: ', trim(output_file)

end program trace_orbit