! Calculate the energy difference E_calc-E_theo
! for umbrella phase, \rho=2.5
!

program proce
 implicit none
 
 ! Imput data
 integer :: nump ! Number of data points

 ! More data
 real(8) :: h ! The field

 real(8) :: ecalc, etheo, de ! Calc, theoretic energy and their difference 

 real(8) :: dummy ! M on nput -- not used

 integer :: i ! Loop counter
 
 ! --- The code

 write (6,*) "Number of data points?"
 read (5,*) nump

 open(21,file='mag.dat',status='old')
 open(22,file='del.dat',status='replace')

 do i=1,nump
  read(21,*) h,dummy,ecalc ! Read h and ecalc
  
  ! Calculate etheo
  if (h> 9.8_8) then ! Spin flip phase
   etheo=3.25_8 - h
  else  ! Umbrella phase
   etheo=-1.65_8 - 2.5_8*h*h/49 
  end if

  ! Calculate the difference
  de=ecalc-etheo
  
  ! Write the result
  write (22,"(2f25.15)") h,de

 end do !i

 close(21)
 close(22)

end program proce
