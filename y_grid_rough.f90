subroutine y_grid_rough

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!   YGRID for roughness  !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Creates the geometry in the y-direction
!  ! y(j) = a(j-q) + b(j-q)**p
  
!  y is the physical coordinate
!  theta is the mapping coordinate (constant increments: dtheta)
!  Since (1) dtheta is constant and (2) the mapping between y and
!  theta is analytical,
!  the second order of the centered finite difference is preserved.
  
!  n   is the total number of point from wall to wall
!  q   is the index of the centerline (n/2) (y(q) = 0)
!  dyq is dymax/dymin (max stretching)
!  p   polynomial exponent
!  a and b are obtained such that dy(nb)/dy(q)=dyq and y(0)=-h_canopy

   use declaration
   implicit none
   integer j,iband,sumj
   real(8) aaa, bbb, aat, aab, dy1, dye, ccc, qqq
  
!   dny = dsty;
   nn  = N(3,nband)+1
 print*, "nn = ", nn
   qqq = nn/2d0
   aaa = 15d0*(dyq - 1d0)/((8d0*dyq + 7d0)*(qqq**ppp))
   bbb =-30d0*(dyq - 1d0)/((8d0*dyq + 7d0)*(qqq**(ppp-2)))
   ccc = 15d0*dyq/((8d0*dyq + 7d0)*qqq)
   dy1 = aaa*(-qqq)**4 + bbb*(-qqq)**2 + ccc

   dny = ceiling(posth/dy1)
   sumj = 0 
   do j = -dny, -1 
     sumj = sumj + j
   end do 
   ! Lineal (almost constant) stretching in the immersed boundary region (below y=-1)
   aab = (posth-dy1*dny)/sumj
  
   dtheta     = 1d0
   dthetai    = 1d0/dtheta
   ddthetavi  = 1.0d0/dtheta


  
   allocate(yu    (  -dny:nn+dny+1))
   allocate(dyu2i (3,-dny:nn+dny+1))
   allocate(dthdyu(  -dny:nn+dny+1))
  
   allocate(yv    (  -dny:nn+dny  ))
   allocate(dyv2i (3,-dny:nn+dny  ))
   allocate(dthdyv(  -dny:nn+dny  ))


 !!! Grid between tips

   do j=0,nn-1
     yv(j)= aaa*(j-qqq)**ppp/5d0 + bbb*(j-qqq)**(ppp-2)/3d0 + ccc*(j-qqq)
     yu(j+1)= aaa*(j+.5d0-qqq)**ppp/5d0 + bbb*(j+.5d0-qqq)**(ppp-2)/3d0 + ccc*(j+.5d0-qqq)
   
     dthdyv(j) = 1d0/(aaa*(j-qqq)**(ppp-1)+bbb*(j-qqq)**(ppp-3)+ccc)
     dthdyu(j+1) = 1d0/(aaa*(j+.5d0-qqq)**(ppp-1)+bbb*(j+.5d0-qqq)**(ppp-3)+ccc)

   end do

   if (dny == 0) then

      yu(0) = aaa*(-0.5d0-qqq)**ppp/5d0+bbb*(-0.5d0-qqq)**(ppp-2)/3d0+ccc*(-0.5d0-qqq)
      dthdyu(0) = 1d0/(aaa*(-.5d0-qqq)**(ppp-1)+bbb*(-.5d0-qqq)**(ppp-3)+ccc)
   
   else

!!!! Grid inside canopies

      yv(-dny)    = -1d0-posth
      yu(-dny)    = yv(-dny)-dy1*0.5d0-aab*(-dny-0.5d0)*0.5d0
      dthdyu(-dny)= 1d0/(dy1+aab*(-dny-1d0)) 
      dthdyv(-dny)= 1d0/(dy1+aab*(-dny-0.5d0)) 
      do j=-dny+1,-1
        yv(j)  = yv(j-1)+dy1+aab*(j-1)
        yu(j)= yu(j-1)+dy1+aab*(j-1.5d0)

        dthdyv(j)  = 1d0/(dy1+aab*(j-0.5d0)) 
        dthdyu(j)= 1d0/(dy1+aab*(j-1d0))
      end do
      yu(0) = yu(-1) +dy1+aab*(-1.5d0) 
      dthdyu(0) = 1d0/(dy1+aab*(-1d0))  

      do j=0,dny
        yv(j+nn)  = yv(j+nn-1)+dy1-aab*j
        yu(j+nn+1)= -yu(-j) 

        dthdyv(j+nn)  = 1d0/(dy1-aab*(j+0.5d0)) 
        dthdyu(j+nn+1)= dthdyu(-j)
      end do
   end if

print*, "dthdyu = ", dthdyu(0), " ", dthdyu(nn+1)
print*, "Shape u, v", shape(dthdyu), ", ", shape(dthdyv)


! Second-order finite difference coefficient for 1. u(j-1), 2. u(j), 3. u(j+1)

  dyub2 = ((yu(nyu21+1)-yu(nyu21))*(yu(nyu21+1)-yu(nyu21-1)))
  dyut2 = ((yu(nyu12)-yu(nyu12-1))*(yu(nyu12+1)-yu(nyu12-1)))

  dyvb2 = ((yv(nyv21+1)-yv(nyv21))*(yv(nyv21+1)-yv(nyv21-1)))
  dyvt2 = ((yv(nyv12)-yv(nyv12-1))*(yv(nyv12+1)-yv(nyv12-1)))

  j=-dny
  dyu2i(1,j)= 0d0
  dyu2i(2,j)=-2d0/(yu(j+1)-yu(j))**2
  dyu2i(3,j)= 1d0/(yu(j+1)-yu(j))**2

  dyv2i(1,j)= 0d0
  dyv2i(2,j)=-2d0/(yv(j+1)-yv(j))**2
  dyv2i(3,j)= 1d0/(yv(j+1)-yv(j))**2
  do j=-dny+1,nn+dny-1
    dyu2i(1,j)=2d0/((yu(j)-yu(j-1))*(yu(j+1)-yu(j-1)))
    dyu2i(2,j)=-2d0/((yu(j+1)-yu(j))*(yu(j)-yu(j-1)))
    dyu2i(3,j)=2d0/((yu(j+1)-yu(j))*(yu(j+1)-yu(j-1)))

    dyv2i(1,j)=2d0/((yv(j)-yv(j-1))*(yv(j+1)-yv(j-1)))
    dyv2i(2,j)=-2d0/((yv(j+1)-yv(j))*(yv(j)-yv(j-1)))
    dyv2i(3,j)=2d0/((yv(j+1)-yv(j))*(yv(j+1)-yv(j-1)))
  end do
  j=nn+dny
  dyu2i(1,j)=2d0/((yu(j)-yu(j-1))*(yu(j+1)-yu(j-1)))
  dyu2i(2,j)=-2d0/((yu(j+1)-yu(j))*(yu(j)-yu(j-1)))
  dyu2i(3,j)=2d0/((yu(j+1)-yu(j))*(yu(j+1)-yu(j-1)))

  dyv2i(1,j)= 1d0/(yv(j)-yv(j-1))**2
  dyv2i(2,j)=-2d0/(yv(j)-yv(j-1))**2
  dyv2i(3,j)= 0d0
  
  j=nn+dny+1
  dyu2i(1,j)= 1d0/(yu(j)-yu(j-1))**2
  dyu2i(2,j)=-2d0/(yu(j)-yu(j-1))**2
  dyu2i(3,j)= 0d0

  nn=N(3,nband)

  N(3,0)=-dny
  N(4,0)=-dny

  do iband = 1,botband
    j = N(3,iband-1)
    h_ny(iband) = yu(0) + h_ny(iband)
    do while (yu(j)<h_ny(iband))
      j = j+1
    end do
    N(3,iband) = j
    N(3,nband-iband) = N(3,nband)-j
  end do
  N(3,nband)=nn+dny

  N(4,1) = N(3,1)+1 !Band 1-2
  N(4,2) = N(3,2) !Band 2-3 
  N(4,3) = N(3,3)+1 !Band 3-end 

  do iband=0,nband+1
    Ngal(3,iband)=N(3,iband)
    Ngal(4,iband)=N(4,iband)
  end do

  Ngal(4,1) = Ngal(3,1)+1 !Band 1-2
  Ngal(4,2) = Ngal(3,2) !Band 2-3 
  Ngal(4,3) = Ngal(3,3)+1 !Band 3-end 

!C! Ny - temp variable, needs removing...
  !Ny(1,:) - v grid points
  !Ny(2,:) - u/w grid points
  !Ny(3,:) - p grid points
  Ny(ugrid,:  ) = N(4,:)
  Ny(vgrid,:  ) = N(3,:)
  Ny(pgrid,0  ) = N(4,0)+1
  Ny(pgrid,1  ) = N(4,1)-1
  Ny(pgrid,2  ) = N(4,2)+1
  Ny(pgrid,3  ) = N(4,3)-1


open(10,file=trim(dirout)//'y_grid_rough.dat',form='unformatted',access='stream')
write(10) N, dny
write(10) yu, yv, dthdyu, dthdyv
write(10) dyu2i, dyv2i
close(10)


end subroutine
