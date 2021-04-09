subroutine y_grid_canopy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!   YGRID for canopies  !!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !Fairhall, C. T. (2019). The influence of superhydrophobic surfaces on near-wall turbulence

   !  Creates the geometry in the y-direction
   !  y(j) = a5(j-q) + b(j-q)**p
  
   !  y is the physical coordinate
   !  theta is the mapping coordinate (constant increments: dtheta)
   !  Since (1) dtheta is constant and (2) the mapping between y and
   !  theta is analytical,
   !  the second order of the centered finite difference is preserved.
  
   !  n   is the total number of point from wall to wall
   !  q   is the index of the centerline (n/2) (y(q) = 0)
   !  dyq is dymax/dymin (max stretching)
   !  p   polynomial exponent
   !  a5 and b are obtained such that dy(nb)/dy(q)=dyq and y(0)=-h_canopy

   use declaration
   implicit none
   integer j,iband, dny
   real(8) aaa, bbb, aat, aab, dy1, dye, ccc, qqq
   real(8) grid_tip, grid_mid
   real(8) a5, a3, a1, b7, b5, b3, b1, b0
   integer dny1

   dny = dsty !!grids in canopy
   dny1 = dsty1
   nn  = N(3,nband)+1 !!nn = ny-1; pts between canopy tips
   print*, "nn = ", nn
	
   qqq = nn/2d0
   aaa = 15d0*(dyq - 1d0)/((8d0*dyq + 7d0)*(qqq**ppp))
   bbb =-30d0*(dyq - 1d0)/((8d0*dyq + 7d0)*(qqq**(ppp-2)))
   ccc = 15d0*dyq/((8d0*dyq + 7d0)*qqq)
   dy1 = aaa*(-qqq)**4 + bbb*(-qqq)**2 + ccc
   aab = (posth-dy1*dny)/(dny**3)
  
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

      j = nn
      yv(j)= aaa*(j-qqq)**ppp/5d0 + bbb*(j-qqq)**(ppp-2)/3d0 + ccc*(j-qqq)
      yu(j+1)= aaa*(j+.5d0-qqq)**ppp/5d0 + bbb*(j+.5d0-qqq)**(ppp-2)/3d0 + ccc*(j+.5d0-qqq)
   
      dthdyv(j) = 1d0/(aaa*(j-qqq)**(ppp-1)+bbb*(j-qqq)**(ppp-3)+ccc)
      dthdyu(j+1) = 1d0/(aaa*(j+.5d0-qqq)**(ppp-1)+bbb*(j+.5d0-qqq)**(ppp-3)+ccc)
   else


      !! Grid inside canopies

      !!!!! Matlab Code to find [a5 a3 a1] !!!!!
      ! syms a5 a4 a3 a2 a1 a0 y x dny1 posth1 grid_tip
      ! y = a5*x^4 + a3*x^2 + a1;%y
      ! dy = diff(y,x);
      ! int_y = int(y,x)+a0;
      ! eqn(1) = subs(dy,x,-dny1)==0;
      ! eqn(2) = subs(dy,x,0)    ==0;%wall
      ! eqn(3) = subs(y,x,0) == grid_tip;%tip slope
      ! eqn(4) = subs(int_y,x,-dny1) == -posth1;
      ! S = solve(eqn,[a5 a3 a1 a0]);

      grid_tip = yv(1)-yv(0)

      if(dny1==0 .OR. dny1==dny .OR. posth1==0d0 .OR. posth1==posth) then

         a5 = -(15d0*(posth - grid_tip*dny)) / (7d0*dny**5d0)
         a3 = (30d0*(posth - grid_tip*dny)) / (7d0*dny**3d0)
         a1 = grid_tip
   
         yu(-dny)    = a5*(-dny-.5d0)**5/5d0 + a3*(-dny-.5d0)**3/3d0 + a1*(-dny-.5d0) - 1d0
         dthdyu(-dny)= 1d0/(a5*(-dny-.5d0)**4 + a3*(-dny-.5d0)**2 + a1)
      
         do j=-dny,-1
            yv(j)  = a5*j**5d0/5d0 + a3*j**3d0/3d0 + a1*j - 1d0
            yu(j+1)= a5*(j+.5d0)**5/5d0 + a3*(j+.5d0)**3/3d0 + a1*(j+.5d0) - 1d0
   
            dthdyv(j) = 1d0/(a5*j**4d0 + a3*j**2d0 +a1)
            dthdyu(j+1)= 1d0/(a5*(j+.5d0)**4 + a3*(j+.5d0)**2 + a1)
         end do
   
         do j=0,dny
            yv(j+nn)  = a5*j**5d0/5d0 + a3*j**3d0/3d0 + a1*j + 1d0
            yu(j+nn+1)= a5*(j+.5d0)**5/5d0 + a3*(j+.5d0)**3/3d0 + a1*(j+.5d0) + 1d0
   
            dthdyv(j+nn)  = 1d0/(a5*j**4d0 + a3*j**2d0 +a1)
            dthdyu(j+nn+1)= 1d0/(a5*(j+.5d0)**4 + a3*(j+.5d0)**2 + a1)
         end do

      else
         a5 = -(15d0*(posth1 - grid_tip*dny1)) / (7d0*dny1**5d0)
         a3 = (30d0*(posth1 - grid_tip*dny1)) / (7d0*dny1**3d0)
         a1 = grid_tip

         !upper half of bottom canopy
         do j=-dny,-1
            yv(j)       = a5*j**5d0/5d0 + a3*j**3d0/3d0 + a1*j - 1d0
            yu(j+1)     = a5*(j+.5)**5d0/5d0 + a3*(j+.5)**3d0/3d0 + a1*(j+.5) - 1d0

            dthdyv(j)   = 1d0/(a5*j**4d0 + a3*j**2d0 +a1)
            dthdyu(j+1) = 1d0/(a5*(j+.5)**4d0 + a3*(j+.5)**2d0 + a1)
         end do

         !upper half of top canopy
         do j=0,dny1
            yv(j+nn)      = a5*j**5d0/5d0 + a3*j**3d0/3d0 + a1*j + 1d0
            yu(j+nn+1)    = a5*(j+.5)**5d0/5d0 + a3*(j+.5)**3d0/3d0 + a1*(j+.5) + 1d0
            
            dthdyv(j+nn)   = 1d0/(a5*j**4d0 + a3*j**2d0 +a1)
            dthdyu(j+nn+1) = 1d0/(a5*(j+.5)**4d0 + a3*(j+.5)**2d0 + a1)
         end do

         grid_mid = yv(-dny1+1) - yv(-dny1);

         !!!!! Matlab Code to find [b7-b0] !!!!!
         ! syms a7 a5 a3 a1 a0 y x dny dny1 posth posth1 grid_mid
         ! y = a7*x^6 + a5*x^4 + a3*x^2 + a1;%y
         ! dy = diff(y,x);
         ! int_y = int(y,x)+a0;
         ! eqn(1) = subs(dy,x,-dny) ==0;%wall
         ! eqn(2) = subs(dy,x,-dny1)== 0;%tip slope
         ! eqn(3) = subs(int_y,x,-dny) == -posth;
         ! eqn(4) = subs(int_y,x,-dny1)== -posth1;
         ! eqn(5) = subs(y,x,-dny1) == grid_mid;%tip slope
         ! S = solve(eqn,[a7 a5 a3 a1 a0]);


         b7 = -(70d0*(posth - posth1 - dny*grid_mid + dny1*grid_mid))/&
         & ((dny - dny1)*(11d0*dny**6d0 + 11d0*dny**5d0*dny1 - 38d0*dny**4d0*dny1**2d0 - 38d0*dny**3d0*dny1**3d0 +&
         & 67d0*dny**2d0*dny1**4d0 + 11d0*dny*dny1**5d0 - 24*dny1**6d0))

         b5 = (105d0*(dny**2d0 + dny1**2d0)*(posth - posth1 - dny*grid_mid + dny1*grid_mid))/&
         & ((dny - dny1)*(11d0*dny**6d0 + 11d0*dny**5d0*dny1 - 38d0*dny**4d0*dny1**2d0 - 38d0*dny**3d0*dny1**3d0 +&
         & 67d0*dny**2d0*dny1**4d0 + 11d0*dny*dny1**5d0 - 24d0*dny1**6d0))

         b3 = (210d0*dny1**2d0*(dny**3d0*grid_mid - dny**2d0*posth + dny**2d0*posth1 - dny**2d0*dny1*grid_mid))/&
         & ((dny - dny1)*(11d0*dny**6d0 + 11d0*dny**5d0*dny1 - 38d0*dny**4d0*dny1**2d0 - 38d0*dny**3d0*dny1**3d0 +&
         & 67d0*dny**2d0*dny1**4d0 + 11d0*dny*dny1**5d0 - 24d0*dny1**6d0))

         b1 = (11d0*dny**7d0*grid_mid - 11d0*dny1**7d0*grid_mid - 35d0*dny1**6d0*posth + 35d0*dny1**6d0*posth1 +&
         & 49d0*dny**2d0*dny1**5d0*grid_mid - 49d0*dny**5d0*dny1**2d0*grid_mid + 105d0*dny**2d0*dny1**4d0*posth -&
         & 105d0*dny**2d0*dny1**4d0*posth1)/((dny - dny1)*(11d0*dny**6d0 + 11d0*dny**5d0*dny1 - 38d0*dny**4d0*dny1**2d0 -&
         & 38d0*dny**3d0*dny1**3d0 + 67d0*dny**2d0*dny1**4d0 + 11d0*dny*dny1**5d0 - 24d0*dny1**6d0))
         
         b0 = -(- 11d0*grid_mid*dny**7d0*dny1 + 11d0*posth1*dny**7d0 + 49d0*grid_mid*dny**5d0*dny1**3d0 -&
         & 49d0*posth1*dny**5d0*dny1**2d0 - 49d0*grid_mid*dny**3d0*dny1**5d0 + 105d0*posth1*dny**3d0*dny1**4d0 -&
         & 56d0*posth*dny**2d0*dny1**5d0 + 11d0*grid_mid*dny*dny1**7d0 - 35d0*posth1*dny*dny1**6d0 + 24d0*posth*dny1**7d0)/&
         & ((dny - dny1)*(11d0*dny**6d0 + 11d0*dny**5d0*dny1 - 38d0*dny**4d0*dny1**2d0 - 38d0*dny**3d0*dny1**3d0 +&
         & 67d0*dny**2d0*dny1**4d0 + 11d0*dny*dny1**5d0 - 24d0*dny1**6d0))
         
         !lower half of bottom canopy
         do j=-dny,-dny1-1
            yv(j)  = b7*j**7d0/7d0      + b5*j**5d0/5d0      + b3*j**3d0/3d0      + b1*j      + b0 - 1d0
            yu(j+1)= b7*(j+.5)**7d0/7d0 + b5*(j+.5)**5d0/5d0 + b3*(j+.5)**3d0/3d0 + b1*(j+.5) + b0 - 1d0
            
            dthdyv(j)   = 1d0/(b7*j**6d0      + b5*j**4d0      + b3*j**2d0      + b1)
            dthdyu(j+1) = 1d0/(b7*(j+.5)**6d0 + b5*(j+.5)**4d0 + b3*(j+.5)**2d0 + b1)
         end do

         j = -dny
         yu(-dny)    = b7*(j-.5)**7d0/7d0 + b5*(j-.5)**5d0/5d0 + b3*(j-.5)**3d0/3d0 + b1*(j-.5) + b0 - 1d0
         dthdyu(-dny)= 1d0/(b7*(j-.5)**6d0 + b5*(j-.5)**4d0 + b3*(j-.5)**2d0 + b1)
         
         !lower half of top canopy
         do j=dny1+1,dny
            yv(j+nn)      = b7*j**7d0/7d0      + b5*j**5d0/5d0      + b3*j**3d0/3d0      + b1*j      - b0 + 1d0
            yu(j+nn+1)    = b7*(j+.5)**7d0/7d0 + b5*(j+.5)**5d0/5d0 + b3*(j+.5)**3d0/3d0 + b1*(j+.5) - b0 + 1d0

            dthdyv(j+nn)  = 1d0/(b7*j**6d0 + b5*j**4d0 + b3*j**2d0 + b1)
            dthdyu(j+nn+1)= 1d0/(b7*(j+.5)**6d0 + b5*(j+.5)**4d0 + b3*(j+.5)**2d0 + b1)
         end do
      end if

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


   open(10,file=trim(dirout)//'y_grid_canopy.dat',form='unformatted',access='stream')
   write(10) N, dny
   write(10) yu, yv, dthdyu, dthdyv
   write(10) dyu2i, dyv2i
   close(10)


end subroutine
