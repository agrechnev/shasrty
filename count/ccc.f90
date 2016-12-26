! Routines for count
module ccc
 implicit none
 
 real(8), parameter :: pi=3.14159265359_8

contains

 !--------------------
 function tp2v(t,p) ! Convert theta, phi -> unit vector
  real(8), dimension(1:3) :: tp2v
 
  real(8), intent(in) :: t,p

  ! Code
  tp2v(1) = sin(t)*cos(p) 
  tp2v(2) = sin(t)*sin(p)  
  tp2v(3) = cos(t) 

 end function tp2v

 !-----------------------
 ! Check if the 2 unitary vectors are "equal"
 function match(s1,s2)

  logical :: match
  real(8), dimension(1:3), intent(in) :: s1,s2

  ! local data
  real(8) :: prod ! Scalar product
  real(8), parameter :: cutoff=1.0e-4

  ! code
  prod = s1(1)*s2(1) + s1(2)*s2(2) + s1(3)*s2(3)  ! Scalar product
  if (prod > 1.0_8) prod=1.0_8
  match = (acos(prod) < cutoff)

 end function match

end module ccc
