! Create a nice xfig file for the SS lattice
! Kind of a truncated version of count
! Optional: plot the circles

program plotss
 implicit none

 integer, parameter :: nx=4
 integer, parameter :: ny=6

 integer,parameter :: max_xy=8000 ! Max coordinates in xfig plot - not used

 integer,parameter :: rr=500 ! radius

 integer :: stepx,stepy,offsx,offsy ! steps, offsets

 integer ::xxmax,yymax ! Max coordinates (max_xy - roundabout errows)

 integer :: xx,yy, xx2,yy2 ! Coords of line ends

 integer :: ix,iy ! Loop counters

 ! Start the code
 
 ! Data for xfig
 open(24,file='tail',status='replace')

 stepx=2000
 offsx=stepx
 stepy=2000
 offsy=stepy

 xxmax=(nx-1)*stepx+offsx
 yymax=(ny-1)*stepx+offsx

 ! Vertical lines
 do ix=1,nx
  xx=(ix-1)*stepx+offsx
  write (24,"(a)") "2 1 0 5 0 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
  write (24,"(a,4i7)") "     ",xx,offsy,xx,yymax
 end do

 ! Horizontal lines
 do iy=1,ny
  yy=(ny-iy)*stepy+offsy
  write (24,"(a)") "2 1 0 5 0 7 60 -1 -1 0.0 0 0 -1 0 0 2"  ! Line header
  write (24,"(a,4i7)") "     ",offsx,yy,xxmax,yy
 end do

! Diagonals
do iy=1,ny-1 ! Loop over "rows", minus the last one
  yy=(ny-iy)*stepy+offsy
  yy2=yy-stepy
 
  do ix=1,nx
   xx=(ix-1)*stepx+offsx
   
   if (mod(ix,2)==1) then ! Odd ix only
    xx2=xx+stepx*(2*mod(iy,2)-1)

    if ((xx2 .ge. offsx) .and. (xx2 .le. xxmax)) then ! Add diag line
      write (24,"(a)") "2 1 1 5 0 7 60 -1 -1 12.0 0 0 -1 0 0 2"  ! Line header
      write (24,"(a,4i7)") "     ",xx,yy,xx2,yy2
    end if

   end if  ! mod(ix,2)=1

  end do ! ix
 end do ! iy

 ! Circles

 do ix=1,nx   ! Loop over all spins
 do iy=1,ny
  ! Data for xfig
  xx=(ix-1)*stepx+offsx
  yy=(ny-iy)*stepy+offsy

  write (24,"(a,2i2,a,8i7)") "1 3 0 5 ",0,0," 50 -1 20 0.0 1 0.0 ",xx,yy,rr,rr,xx,yy,xx+rr,yy

 end do ! iy
 end do ! ix
 
 close(24)

 call system("cat head tail > temp.fig") 

end program plotss