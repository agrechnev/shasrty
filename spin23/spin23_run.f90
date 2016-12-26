! Routines run_2spins and run_3spins
module spin23_run
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
 subroutine v2rtp(v,r,t,p) ! Convert vector -> r,theta,phi

  real(8), dimension(1:3), intent(in) :: v
  real(8), intent(out) :: r,t,p

  ! Code
  r=sqrt(v(1)*v(1)+ v(2)*v(2)+ v(3)*v(3)) ! Calculate r
  
  if (r < 1.0e-20_8) then  ! Too short vector
   t=0
   p=0
   return
  end if

  t=acos(v(3)/r) ! Calculate theta

  if (v(1)*v(1) < 1.0e-20) then ! Check for v(1)=0

   if (v(2) < 0) then
    p=2*pi/3
   else
    p=pi/2
   end if

   return
  end if
  
  ! Calculate phi
  p=atan(v(2)/v(1))
  
  if (v(1) < 0) p=p+pi
  
  if (p<0) p=p+2*pi
   
 end subroutine v2rtp

 !---------------
 ! Analyze a pair of spins
 subroutine  run_2spins(theta1,phi1,theta2,phi2)

  real(8), intent(in) :: theta1,phi1,theta2,phi2

  real(8) :: s1(1:3), s2(1:3), s(1:3)  ! Spins as vectors

  real(8) :: s1s2 ! Scalar product

  real(8) :: r,t,p ! s in polar coordinates

  ! The code

  s1=tp2v(theta1,phi1) ! Get 2 spins and their sum as vectors
  s2=tp2v(theta2,phi2)
  s=s1+s2
  
  s1s2=s1(1)*s2(1) + s1(2)*s2(2) + s1(3)*s2(3) ! Scalar product
  print *,"s1*s2=",s1s2
  print *,"angle = ",acos(s1s2)

  call v2rtp(s,r,t,p)
  print *,"Vector sum: S = S1 + S2 : S,theta,phi = "
  print *,r,t,p
  
 end subroutine  run_2spins

 !---------------
 ! Analyze 3 spins
 subroutine  run_3spins(theta1,phi1,theta2,phi2,theta3,phi3)

  real(8), intent(in) :: theta1,phi1,theta2,phi2,theta3,phi3

  real(8) :: s1(1:3), s2(1:3), s3(1:3), s(1:3)  ! Spins as vectors

  real(8) :: r,t,p ! s in polar coordinates

  ! The code

  s1=tp2v(theta1,phi1) ! Get 3 spins and their sum as vectors
  s2=tp2v(theta2,phi2)
  s3=tp2v(theta3,phi3)
  s=s1+s2+s3
  
  call v2rtp(s,r,t,p)
  print *,"Vector sum: S = S1 + S2 + S3 : S,theta,phi = "
  print *,r,t,p
  
 end subroutine  run_3spins

end module spin23_run
