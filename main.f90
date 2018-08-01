 
program main

  use global_numbers

  implicit none


  integer i,Nxx,Ntt
  integer every_0Dt, every_1Dt


  Namelist /RHD_Input/ xmin, xmax, n_cores,&
       res_num, g_pts, nvarsHyd, nvarsRad, ngroups, &
       Nxx, courant, Ntt, &
       amplitude, x0, & 
       every_0Dt, every_1Dt, &
       gamma, Floor, Floorrad, &
       rho_iniL, rho_iniR, &
       p_iniL, p_iniR, &
       vx_iniL, vx_iniR, & 
       Er_iniL, Er_iniR, &
       cte_rad, opacity_a, opacity_s,close_type,rad_type,eos_type, &  
       Limiter, Riemann, num_det, initial_det, space_between_det,jobname    

  open (3, file='input.par', status = 'old' )
  read (3, nml = RHD_Input)
  close(3)

  zero  = 0.0d0
  third = 1.0d0/3.0d0
  half  = 0.5d0
  one   = 1.0d0
  two   = 2.0d0
  three = 3.0d0
  four  = 4.0d0
  five  = 5.0d0
  pii   = 4.0d0*atan(1.0d0)

  Nx = 2**(res_num-1)*Nxx
  Nt = 2**(res_num-1)*Ntt
  every_0D = 2**(res_num-1)*every_0Dt
  every_1D = 2**(res_num-1)*every_1Dt
  
  chunk_size = Nx/n_cores
  print *, "#################################################"
  print *, "#################################################"
  print *, "################## RRHD1D  v2 ###################"
  print *, "############# Now with multi-group ##############"
  print *, "#################################################"
  print *, "#################################################"

  
  call check_parameters

  call allocate
  call evolve


  print *, '=================================='
  print *, '=CAFE+RAD Code in 1D has finished='
  print *, '=================================='

end program main

