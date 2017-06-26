module mod_param
implicit none
public
   include 'setup.h90'
   character(len=100), parameter :: datadir = 'data/'
   real(8), parameter :: pi = acos(-1.d0)
   real(8), parameter, dimension(2,3) :: rkcoeff = reshape((/ 32.d0/60.d0,  0.d0       , &
                                                              25.d0/60.d0, -17.d0/60.d0, &
                                                              45.d0/60.d0, -25.d0/60.d0/), shape(rkcoeff))
   real(8), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
end module mod_param
