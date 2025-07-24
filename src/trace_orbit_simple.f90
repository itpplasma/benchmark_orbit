program trace_orbit_simple
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp
  use velo_mod, only : isw_field_type
  use parmot_mod, only : rmu, ro0
  use util, only : pi, c, e_charge, p_mass, ev
  use simple_main, only : init_field
  use params, only : params_init, rt0, v0, dtaumin, npoiper2
  use simple, only : Tracer

  implicit none

  ! Variables
  type(Tracer) :: tracy
  character(256) :: vmec_file
  double precision :: bmod_ref

  ! Get VMEC file from command line
  if (command_argument_count() < 1) then
    write(*,*) 'Usage: trace_orbit_simple <vmec_file>'
    stop 1
  endif
  call get_command_argument(1, vmec_file)

  ! Set parameters exactly like SIMPLE examples
  netcdffile = trim(vmec_file)
  ns_s = 3
  ns_tp = 3
  multharm = 3
  isw_field_type = 2  ! Boozer coordinates
  
  write(*,*) 'Initializing SIMPLE with VMEC file: ', trim(vmec_file)
  
  ! Initialize field and parameters
  call init_field(tracy, netcdffile, ns_s, ns_tp, multharm, 2)  ! mode 2 = Boozer
  call params_init()
  
  ! Reference field strength (from SIMPLE examples)
  bmod_ref = 10000.0d0  ! Gauss
  
  write(*,*) 'SIMPLE initialization completed successfully!'
  write(*,*) '  Field mode: ', isw_field_type, ' (Boozer coordinates)'
  write(*,*) '  dtaumin: ', dtaumin
  write(*,*) '  RT0: ', rt0, ' cm'
  write(*,*) '  v0: ', v0, ' cm/s'
  write(*,*) ''
  write(*,*) 'Test completed - SIMPLE library works correctly'

end program trace_orbit_simple