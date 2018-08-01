
subroutine check_parameters

  use arrays
  use global_numbers

  implicit none

  if ((Limiter.ne.'minmod_TVD').and. &
       (Limiter.ne.'MC_TVD')) then
     print *, '================================================================================'
     print *, 'The Limiter you selected --->',Limiter,'<--- is not available'
     print *, '================================================================================'
     stop
  end if

  if ((Riemann.ne.'hlle')) then
     print *, '================================================================================'
     print *, 'The Approximate Riemann solver you selected --->',Riemann,'<--- is not available'
     print *, '================================================================================'
     stop
  end if

end subroutine check_parameters
