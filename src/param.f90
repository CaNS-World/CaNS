module mod_param
implicit none
public
   real(8), parameter :: pi = acos(-1.d0)
   real(8), parameter :: small = 1.d-9
   include 'setup.h90'
   logical, parameter, dimension(2,3) :: no_outflow = & 
       reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
                 .false.,.false.,   & ! no outflow in y lower,upper bound
                 .false.,.false./), & ! no outflow in z lower,upper bound
                 shape(no_outflow))
   character(len=100), parameter :: datadir = 'data/'
   real(8), parameter, dimension(2,3) :: rkcoeff = reshape((/ 32.d0/60.d0,  0.d0       , &
                                                              25.d0/60.d0, -17.d0/60.d0, &
                                                              45.d0/60.d0, -25.d0/60.d0/), shape(rkcoeff))
   real(8), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
end module mod_param
