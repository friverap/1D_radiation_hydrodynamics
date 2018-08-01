
  module arrays

  implicit none
  !hydrodynamical variables
  real(kind=8), allocatable, dimension (:) :: x
  real(kind=8), allocatable, dimension (:) :: SS

  real(kind=8), allocatable, dimension (:,:) :: uHyd
  real(kind=8), allocatable, dimension (:,:) :: uHyd_p
  real(kind=8), allocatable, dimension (:,:) :: rhs_uHyd
  real(kind=8), allocatable, dimension (:,:) :: FluxHyd

  real(kind=8), allocatable, dimension (:) :: rho
  real(kind=8), allocatable, dimension (:) :: vx
  real(kind=8), allocatable, dimension (:) :: press
  real(kind=8), allocatable, dimension (:) :: gradPress
  real(kind=8), allocatable, dimension (:) :: W
  real(kind=8), allocatable, dimension (:) :: Re1
  
  real(kind=8), allocatable, dimension (:,:) :: Temp
  real(kind=8), allocatable, dimension (:,:) :: T_guess, T_stark, T_stark_p
  
  real(kind=8), allocatable, dimension (:) :: ee
  real(kind=8), allocatable, dimension (:) :: hh  
  real(kind=8), allocatable, dimension (:) :: CS
  real(kind=8), allocatable, dimension (:) :: rhou,gam,tauprime
  
  real(kind=8), allocatable, dimension (:) :: P_GUESS, v_guess   
  real(kind=8), allocatable, dimension (:) :: P_STARK, v_stark
  real(kind=8), allocatable, dimension (:) :: P_STARK_p, v_stark_p
  
! radiation variables  
  real(kind=8), allocatable, dimension (:,:,:) :: uRad
  real(kind=8), allocatable, dimension (:,:,:) :: uRad_p
  real(kind=8), allocatable, dimension (:,:,:) :: rhs_uRad
  real(kind=8), allocatable, dimension (:,:,:) :: FluxRad
  real(kind=8), allocatable, dimension (:,:,:) :: sRad
  real(kind=8), allocatable, dimension (:,:,:) :: sHyd

  real(kind=8), allocatable, dimension (:,:) :: Er,ErComovil
  real(kind=8), allocatable, dimension (:,:) :: Fr,FrComovil
  real(kind=8), allocatable, dimension (:,:) :: Pr,PrThick,PrComovil
  real(kind=8), allocatable, dimension (:,:) :: Force_r,chi_s,chi_t,chi
  real(kind=8), allocatable, dimension (:,:) :: Temp_r
  real(kind=8), allocatable, dimension (:,:) :: f,zeta1,limx,Re
  real(kind=8), allocatable, dimension (:,:) :: optical_depth
  
  real(kind=8), allocatable, dimension (:,:) :: Lum,Lum_aux,Er_p
  
  real(kind=8), allocatable, dimension (:) :: Det,freq
  
  real(kind=8), allocatable, dimension (:,:) :: phi1,phi2,ar,psi
  
  real(kind=8), allocatable, dimension (:,:) :: Nxxx,PlanckFunction
  end module arrays
