 
module global_numbers

  real(kind=8) xmin, xmax, courant
  real(kind=8) dx, t, dt,pressradguess
  real(kind=8) amplitude, x0,rtsafe

  real(kind=8) Nf,fmin,fmax,deltaf
  integer Nx, Nt, g_pts, iGA
  integer every_0D, every_1D
  integer res_num, nvarsHyd,nvarsRad,ngroups

  integer method
  
  integer n_cores, chunk_size

  real(kind=8) zero,third,half,one,two,three,four,five,pii

  integer       num_det
  real(kind=8)  initial_det,space_between_det

  character(len=20) :: Limiter
  character(len=20) :: Riemann
  character(len=20) :: close_type,eos_type,rad_type,recv_type
  character*128 jobname

  !  ==================================================
  !  Relativistic 1D Riemann Problem

  real(kind=8) gamma, Floor,floorrad
  real(kind=8) rho_iniL, rho_iniR,Lum1,tao 
  real(kind=8) maxT, maxTr
  real(kind=8) p_iniL, p_iniR 
  real(kind=8) vx_iniL, vx_iniR, W_iniL, W_iniR 
  real(kind=8) Er_iniL,Er_iniR, Pr_iniL, Pr_iniR 
  real(kind=8) Fr_iniL,Fr_iniR  
  real(kind=8) Er_iniObsL,Er_iniObsR 
  real(kind=8) Fr_iniObsL,Fr_iniObsR 
  real(kind=8) Pr_iniObsL,Pr_iniObsR 
  real(kind=8) cte_rad,opacity_a,opacity_s

  ! courtant stuff 

  real(kind=8) lamb_max_x

end module global_numbers

