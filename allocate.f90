
subroutine allocate

  use arrays
  use global_numbers 

  implicit none
  
  !hydrodynamics variables
  allocate(x(-g_pts:Nx+g_pts))
  allocate(SS(-g_pts:Nx+g_pts))

  allocate(uHyd(1:nvarsHyd,-g_pts:Nx+g_pts))
  allocate(uHyd_p(1:nvarsHyd,-g_pts:Nx+g_pts))
  allocate(rhs_uHyd(1:nvarsHyd,-g_pts:Nx+g_pts))
  allocate(FluxHyd(1:nvarsHyd,-g_pts:Nx+g_pts))
  
  allocate(rho(-g_pts:Nx+g_pts))
  allocate(vx(-g_pts:Nx+g_pts))
  allocate(press(-g_pts:Nx+g_pts))
  allocate(gradPress(-g_pts:Nx+g_pts))
  allocate(W(-g_pts:Nx+g_pts))
  allocate(Re1(-g_pts:Nx+g_pts))
  
  allocate(Temp(-g_pts:Nx+g_pts,1:ngroups))

  allocate(ee(-g_pts:Nx+g_pts))
  allocate(hh(-g_pts:Nx+g_pts))
  allocate(CS(-g_pts:Nx+g_pts))
  allocate(rhou(-g_pts:Nx+g_pts),gam(-g_pts:Nx+g_pts))
  
  allocate(P_GUESS(-g_pts:Nx+g_pts),v_GUESS(-g_pts:Nx+g_pts),T_GUESS(-g_pts:Nx+g_pts,1:ngroups))
  allocate(P_STARK(-g_pts:Nx+g_pts),v_STARK(-g_pts:Nx+g_pts),T_STARK(-g_pts:Nx+g_pts,1:ngroups))
  allocate(P_STARK_p(-g_pts:Nx+g_pts),v_STARK_p(-g_pts:Nx+g_pts),T_STARK_p(-g_pts:Nx+g_pts,1:ngroups))
  
!radiation variables  
  allocate(uRad(1:nvarsRad,-g_pts:Nx+g_pts,1:ngroups))
  allocate(uRad_p(1:nvarsRad,-g_pts:Nx+g_pts,1:ngroups))
  allocate(rhs_uRad(1:nvarsRad,-g_pts:Nx+g_pts,1:ngroups))
  allocate(FluxRad(1:nvarsRad,-g_pts:Nx+g_pts,1:ngroups))
  allocate(sRad(1:nvarsRad,-g_pts:Nx+g_pts,1:ngroups))
  allocate(sHyd(1:nvarsHyd,-g_pts:Nx+g_pts,1:ngroups))

  allocate(Er(-g_pts:Nx+g_pts,1:ngroups),ErComovil(-g_pts:Nx+g_pts,1:ngroups))
  allocate(Fr(-g_pts:Nx+g_pts,1:ngroups),FrComovil(-g_pts:Nx+g_pts,1:ngroups))
  allocate(PrComovil(-g_pts:Nx+g_pts,1:ngroups))
  allocate(Pr(-g_pts:Nx+g_pts,1:ngroups),PrThick(-g_pts:Nx+g_pts,1:ngroups))
  allocate(Force_r(-g_pts:Nx+g_pts,1:ngroups))
  allocate(Temp_r(-g_pts:Nx+g_pts,1:ngroups))
  
  allocate(Lum(-g_pts:Nx+g_pts,1:ngroups),Lum_aux(-g_pts:Nx+g_pts,1:ngroups), Er_p(-g_pts:Nx+g_pts,1:ngroups))
  allocate(Det(1:num_det))
  allocate(f(-g_pts:Nx+g_pts,1:ngroups))
  allocate(zeta1(-g_pts:Nx+g_pts,1:ngroups))
  allocate(Re(-g_pts:Nx+g_pts,1:ngroups))
  allocate(chi_t(-g_pts:Nx+g_pts,1:ngroups),chi_s(-g_pts:Nx+g_pts,1:ngroups),chi(-g_pts:Nx+g_pts,1:ngroups))

  allocate(limx(-g_pts:Nx+g_pts,1:ngroups))
  allocate(optical_depth(-g_pts:Nx+g_pts,1:ngroups))
  
  allocate(freq(0:ngroups))
  allocate(phi1(-g_pts:Nx+g_pts,1:ngroups),phi2(-g_pts:Nx+g_pts,1:ngroups))
  allocate(psi(-g_pts:Nx+g_pts,1:ngroups),ar(-g_pts:Nx+g_pts,1:ngroups))
  
  allocate(Nxxx(-g_pts:Nx+g_pts,1:ngroups))
  allocate(PlanckFunction(-g_pts:Nx+g_pts,1:ngroups))
  
end subroutine allocate
