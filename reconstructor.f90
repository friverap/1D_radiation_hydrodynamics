
  real(kind=8) FUNCTION reconstructor(a,b)  

   use global_numbers   

   implicit none 

   real(kind=8)::a,b,c


   c = 0.5D0*(a + b)

   if (Limiter.eq.'minmod_TVD') then 

     reconstructor = 0.5D0*(sign(1.0D0,a) + sign(1.0D0,b))*min(abs(a),abs(b)) 

   else if (Limiter.eq.'MC_TVD') then 

     if (abs(a).lt.abs(b).and.2.0D0*abs(a).lt.abs(c).and.a*b.gt.0.0D0)  then
     reconstructor = 2.0D0*a
     else if (abs(b).lt.abs(a).and.2.0D0*abs(b).lt.abs(c).and.a*b.gt.0.0D0) then
     reconstructor = 2.0D0*b
     else if (abs(c).lt.2.0D0*abs(a).and.abs(c).lt.2.0D0*abs(b).and.a*b.gt.0.0D0) then
     reconstructor = c
     else
     reconstructor = 0.0D0
     end if

   end if

  end function reconstructor
