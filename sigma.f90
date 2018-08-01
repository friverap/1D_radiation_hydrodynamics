
   real(kind=8) FUNCTION sigma(l,N,dh,u_im,u_i,u_ip)  

   implicit none 

   real(kind=8):: smax,smin,idh, dh
   real(kind=8):: u_im, u_i, u_ip
   integer l, N

   idh = 1.0d0/dh  
   
   if (l.eq.0.or.l.eq.N-1) then 

   sigma = 0.0d0 

   else 

   smax = (u_ip - u_i)*idh
   smin = (u_i - u_im)*idh

   sigma = 0.5D0*(sign(1.0D0,smin) + sign(1.0D0,smax))*min(abs(smin),abs(smax)) 

   end if

   end function sigma
