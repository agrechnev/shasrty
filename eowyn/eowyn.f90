! Relax spins in a SS lattice in an external field


program eowyn
 use blonde
 implicit none

 integer :: shuf ! shuffle number
 

 integer :: hiter !Field iteration


 real(8) :: e_theo, m_theo ! Theoretical energy and magnetization for special cases

 real(8) :: beste62,bestm62 ! Best energy and mag_mom for the 6x2 calculation

 ! The code

 ! Read data
 open(21,file='eowyn.in',status='old')
 read (21,*) nx,ny,nz  ! Lattice size
 read (21,*) J1,J2,J3,Bu,Be  ! Hamiltonian parameters
 read (21,*) numh,H1,H2 ! Magnetic field
 read (21,*) maxit,lambda,tol_g,tol_e  ! Convergence parameters
 read (21,*) nshuffle                  ! # of shuffles
 read (21,*) mode_alpha, alpha1, alpha2 ! NEW! for trigen
 close(21)
 
 write (6,*) "Eowyn: I am no man !!!!"
 write (6,*) "  "

 if ((mod(nx,2) .ne. 0) .or. (mod(ny,2) .ne. 0)) then
  print *,"Error: nx,ny must be even"
  stop
 end if

 pstep = max(maxit/50,1) ! Set up pstep

 call blonde_init  ! Initialize everything

 
 lambda0=lambda ! Save lambda

 open(23,file='mag.dat',status='replace')
 do hiter=1,numh ! field loop
    
  if (numh>1) then
   H = H1 + (H2-H1)*(hiter-1)/(numh-1) ! Set field 
  else
   H=H1
  end if

  write (6,"(a)") "----------------------------------------------------------"
  write (6,*) "H=",H
  write (6,"(a)") "----------------------------------------------------------"

  do shuf=1,nshuffle ! Shuffle loop
    
    lambda=lambda0

    write (6,"(a)") "----------------------------------------------------------"
    write (6,*) "Shuffle : ",shuf
    write (6,"(a)") "----------------------------------------------------------"

  
    ! Data template
    if (numh==1) then
     write (6,"(a)") "iter : energy, mag_mom, energy_diff, maxgrad_t, maxgrad_p "
     write (6,"(a)") "----------------------------------------------------------"
    end if !numh==1

    call blonde_run(numh==1)

    ! Find the best energy value algorithm
    if ((shuf==1) .or. (energy < best_e-tol_e/2)) then
    !write (6,*) "Found !",energy,best_e,best_e-tol_e/2

    best_phi=phi
    best_theta=theta
    best_e=energy
    best_m=mag_mom
    best_shuf=shuf
    end if

  end do ! shuf

    ! Write the final data

  write (6,*) "Found on shuffle ",best_shuf
  write (6,*) "E=",best_e
  write (6,*) "M=",best_m

  write (23,"(3f20.14)") H,best_m,best_e 
  
  !------------------

  if (numh==1) then
    call blonde_write

    !call blonde_trigen

    ! Print the exact result, if available
    
    if (Bu>0) then  ! Easy axis uniaxial anisotropy parallel to H

      write (6,*) "Uniaxial anisotropy (parallel to H):"

      if (H< sqrt(Bu*(8*J1+Bu))) then ! Try Neel
       write (6,*) "Try Neel:"
       e_theo=-J1*2+J2/2
       write (6,"(a,f25.15)") "E=", e_theo
       write (6,"(a,f25.15)") "M=", 0.0_8
       if (e_theo<best_e-1.0e-10) write (6,*) "Neel is better !"
      end if ! Neel

      ! Try UUD
      write (6,*) "Try UUD:"
      e_theo=-J1*2/3-J2/6-H/3
      write (6,"(a,f25.15)") "E=", e_theo
      write (6,"(a,f25.15)") "M=", 1.0_8/3
      if (e_theo<best_e-1.0e-10) write (6,*) "UUD is better !"

      if (J2 .ge. 2*J1) then ! Try Dimer
        write (6,*) "Try Dimer:"
        e_theo=-J2/2
        write (6,"(a,f25.15)") "E=", e_theo
        write (6,"(a,f25.15)") "M=", 0.0_8
        if (e_theo<best_e-1.0e-10) write (6,*) "Dimer is better !"

      end if ! Dimer
      
      !if (abs((H-4*J1)/(2*J2-Bu)) < 1.0_8) then ! Try Y
      !  write (6,*) "Try Y:"
      !  write (6,"(a,f25.15)") "E=", 0.25_8*(- (H-4*J1)**2/(2*J2-Bu) + Bu - 2*H)
      !  write (6,"(a,f25.15)") "M=", 0.5_8*(H-4*J1)/(2*J2-Bu) + 0.5_8
      ! 
      !end if ! Y1


      if (J2 .le. J1) then
       if (H < 8*J1-Bu) then 
          write (6,*) "Try Spin Flop:"
          e_theo=     -J1*2+J2/2+Bu/2 - H*H*0.5_8/(8*J1-Bu)    
	  write (6,"(a,f25.15)") "E=", e_theo
	  write (6,"(a,f25.15)") "M=", H/(8*J1-Bu)
          if (e_theo<best_e-1.0e-10) write (6,*) "Spin floip is better !"
       else
	  write (6,*) "Try Spin Flip (FM):"
          e_theo=J1*2+J2/2-H
	  write (6,"(a,f25.15)") "E=", e_theo
	  write (6,"(a,f25.15)") "M=", 1.0_8
      if (e_theo<best_e-1.0e-10) write (6,*) "Spin flip is better !"
       end if ! H < ..

      else

      if (H < 2*(J1+J2)*(J1+J2)/J2-Bu) then 
	write (6,*) "Try Umbrella:"
        e_theo=-J1*J1/J2-J2/2+Bu/2 - H*H*0.5_8/(2*(J1+J2)*(J1+J2)/J2-Bu)
	write (6,"(a,f25.15)") "E=", e_theo
	write (6,"(a,f25.15)") "M=", H/(2*(J1+J2)*(J1+J2)/J2-Bu)
        if (e_theo<best_e-1.0e-10) write (6,*) "Umbrella is better !"
      else
	write (6,*) "Try Spin Flip (FM):"
        e_theo=J1*2+J2/2-H
	write (6,"(a,f25.15)") "E=", e_theo
	write (6,"(a,f25.15)") "M=", 1.0_8
      if (e_theo<best_e-1.0e-10) write (6,*) "Spin flip is better !"
      end if

      end if ! J2 .le. J1

      
    end if ! (Bu>0)

    
    if ((H>0) .and. (J3==0) .and. (Bu==0) .and. (Be==0)) then ! Try spin flip/flop/umbrella

      if ((H < 8*J1) .and. (J2 .le. J1)) then ! Spin flop

      e_theo = - J1*2 +J2/2 - H*H/(J1*16)
      write (6,*) "Try canted Neel (spin flop):"
      write (6,"(a,f25.15)") "E=", e_theo 
      write (6,"(a,f25.15)") "M=", H/(8*J1)
      if (e_theo<best_e-1.0e-10) write (6,*) "Spin flop is better !"

      else if ((J2>J1) .and. (H < 2*(J1+J2)*(J1+J2)/J2)) then ! Umbrella
      
      e_theo = -J1*J1/J2-J2/2 - H*H*J2/((J1+J2)*(J1+J2)*4)
    
      write (6,*) "Try Umbrella :"
      write (6,"(a,f25.15)") "E=", e_theo 
      write (6,"(a,f25.15)") "M=", H*J2/((J1+J2)*(J1+J2)*2)
      write (6,"(a,f20.10)") "theta=",acos(H*J2/((J1+J2)*(J1+J2)*2))
      write (6,"(a,f20.10,a,f20.10)") "phi=",acos(J1/J2),"=pi*",acos(J1/J2)/pi
      if (e_theo<best_e-1.0e-10) write (6,*) "Umbrella is better !"

      else ! Spin flip (FM phase)

      write (6,*) "Try FM phase (spin flip):"
      e_theo=J1*2+J2/2-H
      write (6,"(a,f25.15)") "E=", e_theo
      write (6,"(a,f25.15)") "M=", 1.0_8
      if (e_theo<best_e-1.0e-10) write (6,*) "Spin flip is better !"

      end if ! (H < 8*J1) ...

    end if ! (H>0) ...

    if ((H==0) .and. (J3==0) .and. (Bu==0) .and. (Be==0)) then  ! SS, H=0 
      
      if (J2<J1) then ! Neel
      write (6,*) "Neel phase :"
      write (6,"(a,2f25.15)") "E=", -J1*2+J2/2, 0.0
      else  ! "Spiral"
      write (6,*) "Spiral phase :"
      write (6,"(a,2f25.15)") "E=", -J1*J1/J2-J2/2, 0.0
      end if  ! J2 < J1

    else if ((J2==2*J1) .and. (J3==0) .and. (Bu==0) .and. (Be==0)) then
      
      ! Moliner's situation
      if (H<9*J1) then
      write (6,*) "Moliner unsaturated :"
      write (6,"(a,2f25.15)") "E:", -H*H/(18*J1) - 1.5_8*J1, H/(9*J1)
      else
      write (6,*) "Moliner spin flip :"
      write (6,"(a,2f25.15)") "E:", 3*J1-H, 1.0
      end if
    end if  ! Cases


   !-------- Print the 6x2 result if available
   if (nx .ge. 6) then ! Check because of array sizes
    nx=6 ! Good data for 6x2 .. hopefully
    ny=2
    nz=1
    nsite=12
    lambda=0.2
    maxit=10000
    beste62=0
    do shuf=1,20
     call blonde_run(.false.)
     if ((energy<beste62) .or. (shuf==1)) then 
      beste62=energy
      bestm62=mag_mom
     end if 
    end do !shuf
    
    write (6,*) "6x2 data:"
    write (6,*) "E=",beste62
    write (6,*) "M=",bestm62
    
    if (beste62<best_e-1.0e-9) write (6,*) "6x2 is better !"
    if (best_e<beste62-1.0e-9) write (6,*) "6x2 loses ...."
   end if !  nx .ge. 6

  end if ! numh==1

 end do ! hiter
 
 call blonde_finish ! Finish the calculation
end program eowyn
