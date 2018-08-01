
  subroutine lpm(i,u_im1,u_i,u_ip1,u_ip2,uL,uR)

!  use arrays
  use global_numbers

  implicit none

  real(kind=8), intent(in) :: u_im1,u_i,u_ip1,u_ip2
  real(kind=8), intent(out) :: uL,uR
  integer, intent(in)  :: i
 
  real(kind=8) reconstructor
  real(kind=8) sigmaL,sigmaR
  real(kind=8) smaxL,sminL,smaxR,sminR
  real(kind=8) idx

  idx = one/dx

! ***********************
!      Notación
! ***********************
!   u_im1  =  u(i-1)
!   u_i    =  u(i)
!   u_ip1  =  u(i+1)
!   u_ip2  =  u(i+2)
! ***********************

  if (i.le.g_pts.and.i.ge.Nx+g_pts-1) then

  uL = u_i
  uR = u_ip1

  else

  smaxL = (u_ip1 - u_i)*idx
  sminL = (u_i   - u_im1)*idx
  sigmaL = reconstructor(sminL,smaxL)

  smaxR = (u_ip2 - u_ip1)*idx
  sminR = (u_ip1 - u_i )*idx
  sigmaR = reconstructor(sminR,smaxR)

  uL = u_i   + sigmaL*half*dx
  uR = u_ip1 - sigmaR*half*dx

  end if

  end subroutine lpm

