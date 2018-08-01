
  subroutine rhs
 
  use arrays
  use global_numbers

  implicit none

  integer i,j

  real(kind=8) idx

  idx  = 1.0D0/dx

  call recv_primitives
  call sources

  if (Riemann.eq."hlle") then
     call flux_hlle
  end if
  
!   call sources
!$OMP PARALLEL DO COLLAPSE(1) PRIVATE(j) SCHEDULE(STATIC,chunk_size)

  do i=-g_pts+1,Nx+g_pts-1
   do j=1,ngroups
     rhs_uHyd(1,i) = - (FluxHyd(1,i) - FluxHyd(1,i-1))*idx + sHyd(1,i,j)

     rhs_uHyd(2,i) = - (FluxHyd(2,i) - FluxHyd(2,i-1))*idx + sHyd(2,i,j)

     rhs_uHyd(3,i) = - (FluxHyd(3,i) - FluxHyd(3,i-1))*idx + sHyd(3,i,j)

     rhs_uRad(1,i,j) = - (FluxRad(1,i,j) - FluxRad(1,i-1,j))*idx + sRad(1,i,j)

     rhs_uRad(2,i,j) = - (FluxRad(2,i,j) - FluxRad(2,i-1,j))*idx + sRad(2,i,j)
   end do
  end do

!$OMP END PARALLEL DO
end subroutine rhs


