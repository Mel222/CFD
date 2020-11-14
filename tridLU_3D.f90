!  ! ****** package of utilities for tridiagonal LU factorisation *****
!  !                                   RGM ETSI Aeronauticos 2008
!  ! *************************************************************


! *************************************************************************************************************** !
!
! Last modified by C 20/05/2016
!   - Set up for SHS
!   - Tridiag set up for FFT in middle of bands 1 and 3
!   - Split SolU and SolV as V BC can be solved in Fourier space
!
! TODO
!   - More efficient FFT
!   - DG needs to contatin option for non-zero v velocity at wall
!
! *************************************************************************************************************** !
!
! Contains: LUsolU    - Called from littleharsh for solveU : calls LU_build and LU_dec
!                       - Solves the (1-L) matrix for the U and W velocities
!                       - Decomposed into (1-Lxz)(1-y)
!                       - Contains FFT to solve wall in physical space and interface in Fourier space
!           LUsolV    - Called from littleharsh for solveV : calls LU_build and LU_dec
!                       - Solves the (1-L) matrix for the V velocities
!                       - Decomposed into (1-Lxz)(1-y)
!                       - No FFT all Fourier
!           LUsolP    - Called from littleharsh for solveP : calls LU_buildP and LU_decP
!                       - Solves the lplacian of the pressure (divergance of the gradient)
!           LUsol0    - Called in flowrate_corr (part of const. flow rate) : calls LU_build0 and LU_dec0
!                       - Solves mean mode to calculate pressure gradient for constant mass flow rate
!           LU_build  - Called from Lusol
!                       - Builds the (1-L) matrices. Different for u/w and v
!           LU_buildP - Called from start (with allocation)
!                       - Builds the DG matrix for the pressure
!                       - Only called once as no dt dependence
!           LU_build0 - Called from Lusol0
!                       - Builds LHS for LUsol0
!           LU_dec    - Called from LUsol
!                       - LU Decomposition
!           LU_decP   - Called from LUsolP
!                       - LU Decomposition
!                       - Only called once as no dt dependence
!           LU_dec0   - Called from LUsol0
!                       - LU Decomposition
!
! *************************************************************************************************************** !

subroutine LUsolU(x,grid,myid)
!----------------------------------------------------------------------*
!      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
!      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! INPUT m_j number of unknowns
!       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
!       a(2,:,i,k) inverse of the diagonal of L_j
!       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
!       (note that the diagonal of U_j is 1)
!       x(nystart:nyend,i,k)   rhs of the problem
! OUTPUT
!       all untouched but the solution 'x'
!----------------------------------------------------------------------*

  use declaration
  implicit none
  
  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr

  integer    :: i,k,j,column,iband,grid
  integer    :: iopt,myid
  type(cfield) x(sband:eband)
    
  real(8) , allocatable:: a(:,:,:,:)

  
  !a(diag,column,j,iband)
  allocate(a(3,maxval(columns_num),jlim(1,grid,2):jlim(2,grid,2),nband))

  do iband = sband,eband
    call LU_build(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,a)
    call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,a(:,:,:,:))
  enddo
  
! -Original   
  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      x(iband)%f(jlim(1,grid,iband),column)=x(iband)%f(jlim(1,grid,iband),column)*a(2,column,jlim(1,grid,iband),iband)
      do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
        x(iband)%f(j,column) = (x(iband)%f(j,column)-a(1,column,j,iband)*x(iband)%f(j-1,column))*a(2,column,j,iband)
      end do
      do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
        x(iband)%f(j,column) =  x(iband)%f(j,column)-a(3,column,j,iband)*x(iband)%f(j+1,column)
      end do
    end do
  enddo

  deallocate(a)
  
end subroutine

subroutine LUsolV(x,grid,myid)
!----------------------------------------------------------------------*
!      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
!      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! INPUT m_j number of unknowns
!       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
!       a(2,:,i,k) inverse of the diagonal of L_j
!       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
!       (note that the diagonal of U_j is 1)
!       x(nystart:nyend,i,k)   rhs of the problem
! OUTPUT
!       all untouched but the solution 'x'
!----------------------------------------------------------------------*

  use declaration
  implicit none
  
  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr

  integer    :: i,k,j,column,iband,grid
  integer    :: iopt,myid
  type(cfield) x(sband:eband)
    
  real(8) , allocatable:: a(:,:,:,:)

  
  !a(diag,column,j,iband)
  allocate(a(3,maxval(columns_num),jlim(1,grid,2):jlim(2,grid,2),nband))

  do iband = sband,eband
    call LU_build(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,a)
    call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,a(:,:,:,:))
  enddo
  
! -Original   
  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      x(iband)%f(jlim(1,grid,iband),column)=x(iband)%f(jlim(1,grid,iband),column)*a(2,column,jlim(1,grid,iband),iband)
      do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
        x(iband)%f(j,column) = (x(iband)%f(j,column)-a(1,column,j,iband)*x(iband)%f(j-1,column))*a(2,column,j,iband)
      end do
      do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
        x(iband)%f(j,column) =  x(iband)%f(j,column)-a(3,column,j,iband)*x(iband)%f(j+1,column)
      end do
    end do
  enddo

  deallocate(a)
  
end subroutine

subroutine LUsolP(x,myid,iband,nystart,nyend)
!----------------------------------------------------------------------*
!      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
!      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! INPUT m_j number of unknowns
!       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
!       a(2,:,i,k) inverse of the diagonal of L_j
!       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
!       (note that the diagonal of U_j is 1)
!       x(nystart:nyend,i,k)   rhs of the problem
! OUTPUT
!       all untouched but the solution 'x'
!----------------------------------------------------------------------*

  use declaration
  implicit none
  integer i,k,j,column,myid,iband
  integer nystart,nyend
  complex(8) x(nystart:nyend,columns_num(iband,myid))
   
  do column = 1,columns_num(iband,myid) 
    i = columns_i(column,iband,myid)
    k = columns_k(column,iband,myid)
    x(nystart,column) = x(nystart,column)*DG(iband)%f_dg(2,nystart,column)
    do j = nystart+1,nyend
      x(j,column) = (x(j,column)-DG(iband)%f_dg(1,j,column)*x(j-1,column))*DG(iband)%f_dg(2,j,column)
    end do
    do j = nyend-1,nystart,-1
      x(j,column) =  x(j,column)-DG(iband)%f_dg(3,j,column)*x(j+1,column)
    end do

  end do


end subroutine

subroutine LUsol0(x,nystart,nyend)
!----------------------------------------------------------------------*
!                  LUsol for real (0,0) mode
!----------------------------------------------------------------------*

  use declaration
  implicit none
  integer j
  integer nystart,nyend
  real(8) x(nystart:nyend)
  real(8), allocatable:: a(:,:)
  allocate(a(3,nystart:nyend))

  call LU_build0(nystart,nyend,a)
  call LU_dec0  (nystart,nyend,a)

  x(nystart) = x(nystart)*a(2,nystart)
  do j = nystart+1,nyend
    x(j) = (x(j)-a(1,j)*x(j-1))*a(2,j)
  end do
  do j = nyend-1,nystart,-1
    x(j) = x(j)-a(3,j)*x(j+1)
  end do

  deallocate(a)

end subroutine

subroutine LU_build(nystart,nyend,grid,myid,iband,a)
!-------------------------------------------------------!
!       specifies original values of a(1:3,j,i,k)       !
!-------------------------------------------------------!

  use declaration
  implicit none
  integer i,k,j,grid,myid,column,iband
  integer nystart,nyend
  real(8) k2x,k2z,beta
 
  !real(8) axz(maxval(columns_num),nband)
  
  !ay(BC,diag,j,iband)
  !real(8) ay(2,3,jlim(1,grid,2):jlim(2,grid,2),nband)
  
  !a(diag,column,j,iband)
  real(8) a(3,maxval(columns_num),jlim(1,grid,2):jlim(2,grid,2),nband)

  beta = bRK(kRK)*dt/Re
  

  !!!!!!!!!!!!!!!!!!!     u velocity:       !!!!!!!!!!!!!!!!!!!
  if (grid==ugrid) then
  
    do column = 1,columns_num(iband,myid) 
      i = columns_i(column,iband,myid)
      k = columns_k(column,iband,myid)
      k2x = k2F_x(i)
      k2z = k2F_z(k)

        a(1,column,nystart,iband) = 0d0
        a(2,column,nystart,iband) = 1d0  
        a(3,column,nystart,iband) = gridweighting(iband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
      do j = nystart+1,nyend-1
        a(1,column,j      ,iband) =    -beta* dyu2i(1,j)
        a(2,column,j      ,iband) = 1d0-beta*(dyu2i(2,j)+k2x+k2z)
        a(3,column,j      ,iband) =    -beta* dyu2i(3,j)
      end do
        a(1,column,nyend  ,iband) = gridweighting(iband,2)!0d0 !Free-shear -1 no-slip 0
        a(2,column,nyend  ,iband) = 1d0
        a(3,column,nyend  ,iband) = 0d0
        
    enddo
    
  !!!!!!!!!!!!!!!!!!!     v velocity:       !!!!!!!!!!!!!!!!!!!  
  else if (grid==vgrid) then
  
    do column = 1,columns_num(iband,myid) 
      i = columns_i(column,iband,myid)
      k = columns_k(column,iband,myid)
      k2x = k2F_x(i)
      k2z = k2F_z(k)
        a(1,column,nystart,iband) = 0d0
        a(2,column,nystart,iband) = 1d0  
        a(3,column,nystart,iband) = 0d0
      do j = nystart+1,nyend-1
        a(1,column,j      ,iband) =    -beta* dyv2i(1,j)
        a(2,column,j      ,iband) = 1d0-beta*(dyv2i(2,j)+k2x+k2z)
        a(3,column,j      ,iband) =    -beta* dyv2i(3,j)
      end do
        a(1,column,nyend  ,iband) = 0d0
        a(2,column,nyend  ,iband) = 1d0
        a(3,column,nyend  ,iband) = 0d0
        
    enddo
    
  end if
      

end subroutine
!  
subroutine LU_buildP(nystart,nyend,myid,iband,a)
!-------------------------------------------------------!
!       specifies original values of a(1:3,j,i,k)       !
!-------------------------------------------------------!

  use declaration
  implicit none
  type(rfield_dg)  a(sband:eband)
  integer column,i,k,j,myid,iband
  integer nystart,nyend
  real(8) k2x,k2z

  real(8) D_vec_b(nystart-1:nyend+1)
  real(8) G_vec_b(nystart-1:nyend+1)
  real(8) D_vec_t(nystart-1:nyend+1)
  real(8) G_vec_t(nystart-1:nyend+1)

  do j = nystart,nyend
      D_vec_b(j-1) = -ddthetavi*dthdyu(j)
      G_vec_b(j  ) = -ddthetavi*dthdyv(j)
      D_vec_t(j  ) =  ddthetavi*dthdyu(j)
      G_vec_t(j+1) =  ddthetavi*dthdyv(j)
  end do
  
  
  
  
  do column = 1,columns_num(iband,myid)
    i = columns_i(column,iband,myid)
    k = columns_k(column,iband,myid)
    
    !For exact wavenumbers
     k2x = k2F_x(i) 
     k2z = k2F_z(k)
    
    !For 2nd order centrered difference wavenumbers
!    k2x = k1F_x(i)*k1F_x(i) 
!    k2z = k1F_z(k)*k1F_z(k)
   
    a(iband)%f_dg(2,nystart,column) = (D_vec_t(nystart)*G_vec_b(nystart  )) + k2x + k2z
    a(iband)%f_dg(3,nystart,column) = (D_vec_t(nystart)*G_vec_t(nystart+1))
    a(iband)%f_dg(2,nyend,  column) = (D_vec_b(nyend-1)*G_vec_t(nyend    )) + k2x + k2z
    a(iband)%f_dg(1,nyend,  column) = (D_vec_b(nyend-1)*G_vec_b(nyend-1  ))
    do j = nystart+1,nyend-1
      a(iband)%f_dg(2,j,column) = (D_vec_b(j-1)*G_vec_t(j  )) + (D_vec_t(j)*G_vec_b(j)) + k2x + k2z
      a(iband)%f_dg(1,j,column) = (D_vec_b(j-1)*G_vec_b(j-1))
      a(iband)%f_dg(3,j,column) = (D_vec_t(j  )*G_vec_t(j+1))
    end do

  end do
  
  if(myid==0)then
    if(iband==midband)then 
      a(iband)%f_dg(1,nystart,1) = 0d0
      a(iband)%f_dg(2,nystart,1) = 1d0/((yu(1+1)-yu(1)))**2 !C!
      a(iband)%f_dg(3,nystart,1) = 0d0
    end if
  end if

! For modified wavenumbers, need BC on last pressure mode (as k2x=k2z=0)
!  do column = 1,columns_num(iband,myid)
!  i = columns_i(column,iband,myid)
!  k = columns_k(column,iband,myid)
!  if(abs(k1F_x(i)*k1F_x(i))<10e-10.and.abs(k1F_z(k)*k1F_z(k))<10e-10)then
!  a(iband)%f_dg(1,nystart,column) = 0d0
!  a(iband)%f_dg(2,nystart,column) = 1d0/((yu(1+1)-yu(1)))**2 !C!
!  a(iband)%f_dg(3,nystart,column) = 0d0
!  endif
!  enddo
  

  call LU_decP(nystart,nyend,columns_num(iband,myid),a(iband)%f_dg(:,nystart:nyend,1:columns_num(iband,myid)),iband,myid)

  
end subroutine

subroutine LU_build0(nystart,nyend,a)
!----------------------------------------------------------------------*
!                  LU_build for real (0,0) mode
!----------------------------------------------------------------------*

  use declaration
  implicit none
  integer j
  integer nystart,nyend
  real(8) a(3,nystart:nyend)
  real(8) beta

  !!!!!!!!!!!!!!!!!!     velocity     !!!!!!!!!!!!!!!!!!
  beta = bRK(kRK)*dt/Re
  a(1,nystart) = 0d0
  a(2,nystart) = 1d0
  a(3,nystart) = gridweighting(midband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
  do j = nystart+1,nyend-1
    a(1,j) =    -beta*dyu2i(1,j)
    a(2,j) = 1d0-beta*dyu2i(2,j)
    a(3,j) =    -beta*dyu2i(3,j)
  end do
  a(1,nyend) = gridweighting(midband,2)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
  a(2,nyend) = 1d0
  a(3,nyend) = 0d0

end subroutine

subroutine LU_dec(nystart,nyend,grid,myid,iband,a)
!----------------------------------------------------------------------*
!         GIVEN 'a_j' FOR j=nystart:nyend
!         PERFORMS a_j=L_j*U_j
!         WITH 'a_j' allmost TRIDIAGONAL
! INPUT  nystart:nyend matrix leading dimension
!        iopt selects original value of 'a'
!        a_j is given as a tridiagonal with
!        a(1,nystart:nyend-1) lower diagonal,
!        a(2,:) main diagonal,
!        a(3,nystart+1:nyend) upper diagonal
!        the matrix 'a_j' is assembeld as:
!  a_j(2,1) a_j(3,1)    0        0         0      .............       0
!  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
!     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
!     0     ...................................................       0
!     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
!     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
!
! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
!         a(2,nystart:nyend) inverse of the diagonal of L_j
!         a(3,nystart:nyend-1) upper diagonal of U_j
!         (notice that the diagonal of U_j is 1)
!----------------------------------------------------------------------*

  use declaration
  implicit none
  integer j,column,iband,myid,grid
  integer nystart,nyend
  real(8) a(3,maxval(columns_num),jlim(1,grid,2):jlim(2,grid,2),nband)

    
    do column = 1,columns_num(iband,myid)
      a(2,column,nystart,iband) = 1d0/a(2,column,nystart,iband)
      do j = nystart+1,nyend
	a(3,column,j-1,iband) = a(3,column,j-1,iband)*a(2,column,j-1,iband)
	a(2,column,j  ,iband) = 1d0/(a(2,column,j,iband)-a(1,column,j,iband)*a(3,column,j-1,iband))
      end do
    enddo
  

end subroutine
!  
subroutine LU_decP(nystart,nyend,columns,a,iband,myid)
!----------------------------------------------------------------------*
!         GIVEN 'a_j' FOR j=nystart:nyend
!         PERFORMS a_j=L_j*U_j
!         WITH 'a_j' allmost TRIDIAGONAL
! INPUT  nystart:nyend matrix leading dimension
!        iopt selects original value of 'a'
!        a_j is given as a tridiagonal with
!        a(1,nystart:nyend-1) lower diagonal,
!        a(2,:) main diagonal,
!        a(3,nystart+1:nyend) upper diagonal
!        the matrix 'a_j' is assembeld as:
!  a_j(2,1) a_j(3,1)    0        0         0      .............       0
!  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
!     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
!     0     ...................................................       0
!     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
!     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
!
! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
!         a(2,nystart:nyend) inverse of the diagonal of L_j
!         a(3,nystart:nyend-1) upper diagonal of U_j
!         (notice that the diagonal of U_j is 1)
!----------------------------------------------------------------------*

  use declaration
  implicit none

  integer j,column
  integer columns,iband,myid
  integer nystart,nyend
  real(8) a(3,nystart:nyend,columns)
  
  do column = 1,columns
    a(2,nystart,column) = 1d0/a(2,nystart,column)
    do j = nystart+1,nyend
      a(3,j-1,column) = a(3,j-1,column)*a(2,j-1,column)     
      a(2,j  ,column) = 1d0/(a(2,j,column) - a(1,j,column)*a(3,j-1,column))
    end do
  end do 

end subroutine

subroutine LU_dec0(nystart,nyend,a)
!----------------------------------------------------------------------*
!                  LU_dec for real (0,0) mode
!----------------------------------------------------------------------*

  use declaration
  implicit none
  integer j
  integer nystart,nyend
  real(8) a(3,nystart:nyend)

  a(2,nystart) = 1d0/a(2,nystart)
  do j = nystart+1,nyend
    a(3,j-1) = a(3,j-1)*a(2,j-1)
    a(2,j  ) = 1d0/(a(2,j)-a(1,j)*a(3,j-1))
  end do

end subroutine
