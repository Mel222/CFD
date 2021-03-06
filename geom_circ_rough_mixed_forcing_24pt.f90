subroutine boundary_circ_rough
! This function creates the geometry- tall, sparse filaments.
! It defines lists with the points at the forcing points (list_ib1),
!   

! ATTENTION, it raises an error if the number of points in x (Ngal(1,1)) and z (Ngal(2,1))
!   is not a multiple of the number of tiles (ntilex, ntilez)

! It is called by 'getbounds' in 'start.f90'

!                          |- Periodic unit -|
!    dloly    ***          ***               ***                         ***           ***           ***                    
!     .       ***          ***               ***                         ***           ***           ***                    
!    dsty      *            *                 *                           *             *             *                     
! Y   2        *            *                 *                           *             *             *                      
!     1        *            *                 *                           *             *             *                     
!   wall  *    *     *      *         *       *      *      *      *      *      *      *      *      *      *     *  
!                           1      2      3                           7      8     9(1)  10(2)  11(3)
!                                                      X,Z


!   ------------------------ ny22
!   \  /\  /\  /\  /\  /\  /
!    \/  \/  \/  \/  \/  \/
!   ------------------------ ny12
!
!
!
!   ------------------------ ny21
!    /\  /\  /\  /\  /\  /\
!   /  \/  \/  \/  \/  \/  \
!   ------------------------ ny11



  use declaration
  implicit none
  integer k_ele, i_ele, j,ilist,ilist_xz,ix,iz,shift
  integer nlist_xz, nlist_ib_bot_v, nlist_ib_bot_u  
  integer points_stem_v, points_stem_u, ind 
  integer, allocatable:: list_ib_bot_v(:,:), list_ib_top_v(:,:)  
  integer, allocatable:: list_ib_bot_u(:,:), list_ib_top_u(:,:)  
  real(8), allocatable:: list_ib_w_bot_u(:,:), list_ib_w_top_u(:,:)  
  real(8), allocatable:: list_ib_w_bot_v(:,:), list_ib_w_top_v(:,:)  
  integer, allocatable:: sortInd(:)  
  real(8) :: r, c_i, c_k, shift_stag
  integer dstx1, dstx2, dstz1, dstz2
  
  real(8), allocatable :: test_circ(:,:) 
  integer, allocatable :: list_ib_xz(:,:)  
  real(8), allocatable :: list_ib_w_xz(:,:)  
  logical :: rig, lef, top, bot
  integer :: interpp(4), i_ele_o, k_ele_o
  real(8) :: weicoef(2), zbound, xbound, circ


! Define the limits of the IB in the y direction   
  nyu11 = -dsty+1
  nyu21 = 0
  nyv11 = -dsty+1        ! Dont want to include wall point in IB list
  nyv21 = 0-1            !-1 To make the boundary align with u_grid

  nyu12 = Ngal(4,3)-dsty+1
  nyu22 = Ngal(4,3)

  nyv12 = Ngal(3,3)+1-dsty+1  !+1 To make the boundary align with u_grid
  nyv22 = Ngal(3,3)           ! Dont want to include wall point in IB list

print*, "nyv11, 21, 12, 22", nyv11, nyv21, nyv12, nyv22
print*, "nyu11, 21, 12, 22", nyu11, nyu21, nyu12, nyu22

! Check that the number of grid points in x and z is a multiple of the number of tiles 
  if (ntilex.eq.0) then
    write(*,*) 'ERROR: ntilex equals 0'
    stop
  end if
  if (ntilez.eq.0) then
    write(*,*) 'ERROR: ntilez equals 0'
    stop
  end if
  if (mod(ntilex,2).ne.0) then
    write(*,*) 'ERROR: ntilex is not even'
    stop
  end if
  if (mod(Ngal(1,1),ntilex)/=0) then
    write(*,*) 'ERROR: nx not multiple of ntilex'
    stop
  end if
  if (mod(Ngal(2,1),ntilez)/=0) then
    write(*,*) 'ERROR: nz not multiple of ntilez'
    stop
  end if

  dnx = Ngal(1,1)/ntilex  ! Width:  Number of points per tile the in streamwise direction
  dnz = Ngal(2,1)/ntilez  ! Width:  Number of points per tile the in spanwise   direction


! CREATE CIRCULAR CROSS SECTION IN I-K PLANE FOR ONE TILE
!
!      k
!      ^
!      |      ___          
!      |   _-     -_
!      |  -         -
!  c_k | *     *     *     
!      |  _         _
!      |   -       -
!      |      ---          
! 0.5 -|-----------------> i
!     0.5     c_i  

! READ IN POSITION OF CIRCLES 
  open(69,file='geom_circ_rough.txt',form='formatted')
  10 FORMAT(7X,D10.1)
!  20 FORMAT(7X,I10)
  read(69,10) r
  read(69,10) c_i
  read(69,10) c_k
  read(69,10) shift_stag
  close(69)
!  r = 4.5d0 ! Values for the hexagon 
!  c_i = 9d0
!  c_k = 9d0
!  shift_stag = 0.5d0

! LIMITS OF THE SQUARE WHICH CONTAINS THE CIRCLE
  dstx1 = floor(  c_i-r-1);
  dstx2 = ceiling(c_i+r+1);
  dstz1 = floor(  c_k-r-1);
  dstz2 = ceiling(c_k+r+1);

  if (dstx1 < 1 .or. dstz1 < 1) then
    write(*,*)'Circ too close to edge'
    stop
  end if

  if (dstx2 > dnx .or. dstz2 > dnz) then  
    write(*,*) 'Circ too close to edge'
    stop
  end if 

  allocate(test_circ(dstx1:dstx2, dstz1:dstz2))
  allocate(list_ib_xz(  6,(dstx2-dstx1)*(dstz2-dstz1)))  
  allocate(list_ib_w_xz(5,(dstx2-dstx1)*(dstz2-dstz1)))  

  do k_ele = dstz1, dstz2
    do i_ele = dstx1, dstx2
      test_circ(i_ele, k_ele) = circ(i_ele*1d0, k_ele*1d0, r, c_i, c_k)
    end do
  end do

  nlist_xz = 0 !num of points for one stem in one xz plane
  list_ib_xz = 0
  list_ib_w_xz = 0d0
 
  do k_ele = dstz1+1, dstz2-1
    do i_ele = dstx1+1, dstx2-1

      if (test_circ(i_ele, k_ele) .lt. 0d0) then
        ! Inside the circle  

        rig = test_circ(i_ele-1,k_ele  ) > 0d0
        lef = test_circ(i_ele+1,k_ele  ) > 0d0
        bot = test_circ(i_ele  ,k_ele+1) > 0d0
        top = test_circ(i_ele  ,k_ele-1) > 0d0 
      
        if (rig .or. lef .or. bot .or. top) then
          ! Forcing is needed
          interpp = 0
          weicoef = 0d0
          zbound = 0d0
          xbound = 0d0
          call lin_interp(i_ele*1d0,k_ele*1d0,interpp,weicoef,zbound,xbound,r,c_i,c_k)

          if (weicoef(1) < -1d0 .or. weicoef(2) < -1d0) then 
            ! Make this a solid point and force outside
            ! This is a solid point 
            nlist_xz = nlist_xz + 1
            call add_solid_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz),i_ele,k_ele)
       
            ! Force outside 
            i_ele_o = interpp(2)
            k_ele_o = interpp(1)
            interpp = 0
            weicoef = 0d0
            zbound = 0d0
            xbound = 0d0
            call lin_interp(i_ele_o*1d0,k_ele_o*1d0,interpp,weicoef,zbound,xbound,r,c_i,c_k)

            nlist_xz = nlist_xz + 1
            call add_forcing_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz), &
&           i_ele_o,k_ele_o,interpp,weicoef,xbound,zbound)

          else 
            ! Force inside 
            nlist_xz = nlist_xz + 1
            call add_forcing_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz), &
&           i_ele,k_ele,interpp,weicoef,xbound,zbound)

          end if 

        else
          ! This is a solid point
          nlist_xz = nlist_xz + 1
          call add_solid_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz),i_ele,k_ele)
        end if
      end if

      if (test_circ(i_ele, k_ele) == 0d0) then
        ! On the circle 
        ! This is a solid point
        nlist_xz = nlist_xz + 1
        call add_solid_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz),i_ele,k_ele)
      end if

    end do
  end do

  ! Add some additional force outside points
  i_ele = 8
  k_ele = 16
  interpp = 0
  weicoef = 0d0
  zbound = 0d0
  xbound = 0d0
  call lin_interp(i_ele*1d0,k_ele*1d0,interpp,weicoef,zbound,xbound,r,c_i,c_k)
  nlist_xz = nlist_xz + 1
  call add_forcing_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz), &
&           i_ele,k_ele,interpp,weicoef,xbound,zbound)
  i_ele = 16
  k_ele = 16
  interpp = 0
  weicoef = 0d0
  zbound = 0d0
  xbound = 0d0
  call lin_interp(i_ele*1d0,k_ele*1d0,interpp,weicoef,zbound,xbound,r,c_i,c_k)
  nlist_xz = nlist_xz + 1
  call add_forcing_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz), &
&           i_ele,k_ele,interpp,weicoef,xbound,zbound)
  i_ele = 16
  k_ele = 8
  interpp = 0
  weicoef = 0d0
  zbound = 0d0
  xbound = 0d0
  call lin_interp(i_ele*1d0,k_ele*1d0,interpp,weicoef,zbound,xbound,r,c_i,c_k)
  nlist_xz = nlist_xz + 1
  call add_forcing_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz), &
&           i_ele,k_ele,interpp,weicoef,xbound,zbound)
  i_ele = 8
  k_ele = 8
  interpp = 0
  weicoef = 0d0
  zbound = 0d0
  xbound = 0d0
  call lin_interp(i_ele*1d0,k_ele*1d0,interpp,weicoef,zbound,xbound,r,c_i,c_k)
  nlist_xz = nlist_xz + 1
  call add_forcing_point(list_ib_xz(1,nlist_xz),list_ib_w_xz(1,nlist_xz), &
&           i_ele,k_ele,interpp,weicoef,xbound,zbound)
  
  
  deallocate(test_circ)


! CREATE FIRST STEM MEL
!            _______
!           |       |
!           |__   __|
!              | |
!              | |
!              | |
!              |_|

  points_stem_u   = dsty*nlist_xz
  points_stem_v   = (dsty-1)*nlist_xz
  nlist_ib_bot_u   = ntilex*ntilez*points_stem_u      ! Number of points in all stems.
  nlist_ib_bot_v   = ntilex*ntilez*points_stem_v  

  allocate (list_ib_bot_u(  9, nlist_ib_bot_u))
  allocate (list_ib_w_bot_u(3, nlist_ib_bot_u))
  allocate (list_ib_top_u(  9, nlist_ib_bot_u))
  allocate (list_ib_w_top_u(3, nlist_ib_bot_u))
  allocate (list_ib_bot_v(  9, nlist_ib_bot_v))
  allocate (list_ib_w_bot_v(3, nlist_ib_bot_v))
  allocate (list_ib_top_v(  9, nlist_ib_bot_v))
  allocate (list_ib_w_top_v(3, nlist_ib_bot_v))

! FIRST STEM BOTTOM VGRID
  ilist = 0
  do j = nyv11, nyv21-1
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_bot_v(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_bot_v(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_bot_v(3, ilist) =  j

      list_ib_bot_v(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_bot_v(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_bot_v(6, ilist) =  j

      list_ib_bot_v(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_bot_v(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_bot_v(9, ilist) =  j

      list_ib_w_bot_v(1, ilist) =  list_ib_w_xz(5, ilist_xz) 
      list_ib_w_bot_v(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_bot_v(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do
  end do
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_bot_v(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_bot_v(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_bot_v(3, ilist) =  nyv21

      list_ib_bot_v(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_bot_v(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_bot_v(6, ilist) =  nyv21

      list_ib_bot_v(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_bot_v(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_bot_v(9, ilist) =  nyv21

      list_ib_w_bot_v(1, ilist) =  1d0 
      list_ib_w_bot_v(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_bot_v(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do

  if (ilist/=points_stem_v) then
    write(*,*) 'ERROR: ilist is not equal to points_stem 1'
    stop
  end if 
  
! FIRST STEM BOTTOM UGRID
  ilist = 0
  do j = nyu11, nyu21-1
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_bot_u(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_bot_u(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_bot_u(3, ilist) =  j

      list_ib_bot_u(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_bot_u(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_bot_u(6, ilist) =  j

      list_ib_bot_u(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_bot_u(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_bot_u(9, ilist) =  j

      list_ib_w_bot_u(1, ilist) =  list_ib_w_xz(5, ilist_xz) 
      list_ib_w_bot_u(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_bot_u(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do
  end do
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_bot_u(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_bot_u(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_bot_u(3, ilist) =  nyu21

      list_ib_bot_u(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_bot_u(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_bot_u(6, ilist) =  nyu21

      list_ib_bot_u(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_bot_u(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_bot_u(9, ilist) =  nyu21

      list_ib_w_bot_u(1, ilist) =  1d0 
      list_ib_w_bot_u(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_bot_u(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do

  if (ilist/=points_stem_u) then
    write(*,*) 'ERROR: ilist is not equal to points_stem 2'
    stop
  end if 

! FIRST STEM TOP VGRID
  ilist = 0
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_top_v(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_top_v(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_top_v(3, ilist) =  nyv12

      list_ib_top_v(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_top_v(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_top_v(6, ilist) =  nyv12

      list_ib_top_v(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_top_v(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_top_v(9, ilist) =  nyv12

      list_ib_w_top_v(1, ilist) =  1d0 
      list_ib_w_top_v(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_top_v(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do
  do j = nyv12+1, nyv22
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_top_v(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_top_v(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_top_v(3, ilist) =  j

      list_ib_top_v(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_top_v(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_top_v(6, ilist) =  j

      list_ib_top_v(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_top_v(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_top_v(9, ilist) =  j

      list_ib_w_top_v(1, ilist) =  list_ib_w_xz(5, ilist_xz) 
      list_ib_w_top_v(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_top_v(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do
  end do

  if (ilist/=points_stem_v) then
    write(*,*) 'ERROR: ilist is not equal to points_stem 3'
    stop
  end if
 
! FIRST STEM TOP UGRID
  ilist = 0
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_top_u(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_top_u(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_top_u(3, ilist) =  nyu12

      list_ib_top_u(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_top_u(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_top_u(6, ilist) =  nyu12

      list_ib_top_u(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_top_u(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_top_u(9, ilist) =  nyu12

      list_ib_w_top_u(1, ilist) =  1d0 
      list_ib_w_top_u(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_top_u(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do
  do j = nyu12+1, nyu22
    do ilist_xz = 1, nlist_xz

      ilist = ilist + 1

      list_ib_top_u(1, ilist) =  list_ib_xz(1, ilist_xz) 
      list_ib_top_u(2, ilist) =  list_ib_xz(2, ilist_xz) 
      list_ib_top_u(3, ilist) =  j

      list_ib_top_u(4, ilist) =  list_ib_xz(3, ilist_xz) 
      list_ib_top_u(5, ilist) =  list_ib_xz(4, ilist_xz) 
      list_ib_top_u(6, ilist) =  j

      list_ib_top_u(7, ilist) =  list_ib_xz(5, ilist_xz) 
      list_ib_top_u(8, ilist) =  list_ib_xz(6, ilist_xz) 
      list_ib_top_u(9, ilist) =  j

      list_ib_w_top_u(1, ilist) =  list_ib_w_xz(5, ilist_xz) 
      list_ib_w_top_u(2, ilist) =  list_ib_w_xz(1, ilist_xz) 
      list_ib_w_top_u(3, ilist) =  list_ib_w_xz(2, ilist_xz) 

    end do
  end do

  if (ilist/=points_stem_u) then
    write(*,*) 'ERROR: ilist is not equal to points_stem 4'
    stop
  end if 

! REPLICATE THE PATTERN STEM
! This section should be common for all geometries
  do ix = 1,ntilex
    do iz = 1,ntilez
      shift = points_stem_u*(iz-1) + points_stem_u*ntilez*(ix-1) 
      do ilist = 1,points_stem_u
        list_ib_bot_u(1,ilist+shift) = list_ib_bot_u(1,ilist) + dnx*(ix-1)  ! i
        list_ib_bot_u(2,ilist+shift) = 1+mod(list_ib_bot_u(2,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))
        list_ib_bot_u(3,ilist+shift) = list_ib_bot_u(3,ilist)               ! j 
        list_ib_bot_u(4,ilist+shift) = list_ib_bot_u(4,ilist) + dnx*(ix-1)  ! i
        list_ib_bot_u(5,ilist+shift) = 1+mod(list_ib_bot_u(5,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
        list_ib_bot_u(6,ilist+shift) = list_ib_bot_u(6,ilist)               ! j 
        list_ib_bot_u(7,ilist+shift) = list_ib_bot_u(7,ilist) + dnx*(ix-1)  ! i
        list_ib_bot_u(8,ilist+shift) = 1+mod(list_ib_bot_u(8,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
        list_ib_bot_u(9,ilist+shift) = list_ib_bot_u(9,ilist)               ! j 

        list_ib_w_bot_u(1,ilist+shift) = list_ib_w_bot_u(1,ilist)
        list_ib_w_bot_u(2,ilist+shift) = list_ib_w_bot_u(2,ilist)
        list_ib_w_bot_u(3,ilist+shift) = list_ib_w_bot_u(3,ilist)


        list_ib_top_u(1,ilist+shift) = list_ib_top_u(1,ilist) + dnx*(ix-1)  ! i
        list_ib_top_u(2,ilist+shift) = 1+mod(list_ib_top_u(2,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
        list_ib_top_u(3,ilist+shift) = list_ib_top_u(3,ilist)               ! j 
        list_ib_top_u(4,ilist+shift) = list_ib_top_u(4,ilist) + dnx*(ix-1)  ! i
        list_ib_top_u(5,ilist+shift) = 1+mod(list_ib_top_u(5,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))   
        list_ib_top_u(6,ilist+shift) = list_ib_top_u(6,ilist)               ! j 
        list_ib_top_u(7,ilist+shift) = list_ib_top_u(7,ilist) + dnx*(ix-1)  ! i
        list_ib_top_u(8,ilist+shift) = 1+mod(list_ib_top_u(8,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))   
        list_ib_top_u(9,ilist+shift) = list_ib_top_u(9,ilist)               ! j 

        list_ib_w_top_u(1,ilist+shift) = list_ib_w_top_u(1,ilist)                
        list_ib_w_top_u(2,ilist+shift) = list_ib_w_top_u(2,ilist)                
        list_ib_w_top_u(3,ilist+shift) = list_ib_w_top_u(3,ilist)                
      end do

      shift = points_stem_v*(iz-1) + points_stem_v*ntilez*(ix-1) 
      do ilist = 1,points_stem_v
        list_ib_bot_v(1,ilist+shift) = list_ib_bot_v(1,ilist) + dnx*(ix-1)  ! i
        list_ib_bot_v(2,ilist+shift) = 1+mod(list_ib_bot_v(2,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3)) 
        list_ib_bot_v(3,ilist+shift) = list_ib_bot_v(3,ilist)               ! j 
        list_ib_bot_v(4,ilist+shift) = list_ib_bot_v(4,ilist) + dnx*(ix-1)  ! i
        list_ib_bot_v(5,ilist+shift) = 1+mod(list_ib_bot_v(5,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
        list_ib_bot_v(6,ilist+shift) = list_ib_bot_v(6,ilist)               ! j 
        list_ib_bot_v(7,ilist+shift) = list_ib_bot_v(7,ilist) + dnx*(ix-1)  ! i
        list_ib_bot_v(8,ilist+shift) = 1+mod(list_ib_bot_v(8,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))
        list_ib_bot_v(9,ilist+shift) = list_ib_bot_v(9,ilist)               ! j 

        list_ib_w_bot_v(1,ilist+shift) = list_ib_w_bot_v(1,ilist)
        list_ib_w_bot_v(2,ilist+shift) = list_ib_w_bot_v(2,ilist)
        list_ib_w_bot_v(3,ilist+shift) = list_ib_w_bot_v(3,ilist)


        list_ib_top_v(1,ilist+shift) = list_ib_top_v(1,ilist) + dnx*(ix-1)  ! i
        list_ib_top_v(2,ilist+shift) = 1+mod(list_ib_top_v(2,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))   
        list_ib_top_v(3,ilist+shift) = list_ib_top_v(3,ilist)               ! j 
        list_ib_top_v(4,ilist+shift) = list_ib_top_v(4,ilist) + dnx*(ix-1)  ! i 
        list_ib_top_v(5,ilist+shift) = 1+mod(list_ib_top_v(5,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
        list_ib_top_v(6,ilist+shift) = list_ib_top_v(6,ilist)               ! j 
        list_ib_top_v(7,ilist+shift) = list_ib_top_v(7,ilist) + dnx*(ix-1)  ! i 
        list_ib_top_v(8,ilist+shift) = 1+mod(list_ib_top_v(8,ilist) + dnz*(iz-1) &! k
&                                            + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3)) 
        list_ib_top_v(9,ilist+shift) = list_ib_top_v(9,ilist)               ! j 

        list_ib_w_top_v(1,ilist+shift) = list_ib_w_top_v(1,ilist)                
        list_ib_w_top_v(2,ilist+shift) = list_ib_w_top_v(2,ilist)                
        list_ib_w_top_v(3,ilist+shift) = list_ib_w_top_v(3,ilist)                
      end do
    end do
  end do

! Sort the lists 
  allocate(sortInd(nlist_ib_bot_v))
  call MergeSort3D(nlist_ib_bot_v, list_ib_bot_v(1:3,:), sortInd, (/.true.,.true.,.true./))
  do ind = 1,3
    list_ib_w_bot_v(ind,:) = list_ib_w_bot_v(ind,sortInd)
    list_ib_bot_v(ind+3,:) = list_ib_bot_v(ind+3,sortInd)
    list_ib_bot_v(ind+6,:) = list_ib_bot_v(ind+6,sortInd)
  end do 
  call MergeSort3D(nlist_ib_bot_v, list_ib_top_v(1:3,:), sortInd, (/.true.,.true.,.true./))
  do ind = 1,3
    list_ib_w_top_v(ind,:) = list_ib_w_top_v(ind,sortInd)
    list_ib_top_v(ind+3,:) = list_ib_top_v(ind+3,sortInd)
    list_ib_top_v(ind+6,:) = list_ib_top_v(ind+6,sortInd)
  end do 
  deallocate(sortInd)
  allocate(sortInd(nlist_ib_bot_u))
  call MergeSort3D(nlist_ib_bot_u, list_ib_bot_u(1:3,:), sortInd, (/.true.,.true.,.true./))
  do ind = 1,3
    list_ib_w_bot_u(ind,:) = list_ib_w_bot_u(ind,sortInd)
    list_ib_bot_u(ind+3,:) = list_ib_bot_u(ind+3,sortInd)
    list_ib_bot_u(ind+6,:) = list_ib_bot_u(ind+6,sortInd)
  end do 
  call MergeSort3D(nlist_ib_bot_u, list_ib_top_u(1:3,:), sortInd, (/.true.,.true.,.true./))
  do ind = 1,3
    list_ib_w_top_u(ind,:) = list_ib_w_top_u(ind,sortInd)
    list_ib_top_u(ind+3,:) = list_ib_top_u(ind+3,sortInd)
    list_ib_top_u(ind+6,:) = list_ib_top_u(ind+6,sortInd)
  end do 
  deallocate(sortInd)
  
! Save the lists into a file
  open(10,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3//'.dat',form='unformatted',access='stream')
  write(10) Lx,Ly,Lz
  write(10) Ngal,nlist_ib_bot_v,nlist_ib_bot_u,nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
  write(10) list_ib_bot_v, list_ib_top_v, list_ib_bot_u, list_ib_top_u
  write(10) list_ib_w_bot_v, list_ib_w_top_v, list_ib_w_bot_u, list_ib_w_top_u
  write(10) r, c_i, c_k, list_ib_xz, list_ib_w_xz
  close(10)
print*, "Geom complete"

  deallocate(list_ib_xz, list_ib_w_xz)
  deallocate(list_ib_bot_v, list_ib_top_v, list_ib_bot_u, list_ib_top_u)
  deallocate(list_ib_w_bot_v, list_ib_w_top_v, list_ib_w_bot_u, list_ib_w_top_u)
 
end subroutine 

real(8) function circ(x, z, r, c_i, c_k)
  
  implicit none 

  real(8) :: x, z, r, c_i, c_k

  circ = (x - c_i)**2 + (z - c_k)**2 - r**2

end function

subroutine lin_interp(i_ele, k_ele, interpp, weicoef, zbound, xbound, radius, c_i, c_k)

  use declaration
  implicit none

  real(8) :: i_ele, k_ele, aa
  integer :: z2, x2, z3, x3, interpp(4)
  integer :: n_z_sign, n_x_sign 
  real(8) :: radius, c_i, c_k
  real(8) :: zbound, xbound, n_z, n_x
  real(8) :: weicoef(2), denominator 


  if (i_ele .ge. c_i) then
    aa = 1d0
  else
    aa = -1d0
  end if

  zbound = c_k + (k_ele-c_k)*radius/((k_ele-c_k)**2 + (i_ele-c_i)**2)**0.5d0
  xbound = aa*sqrt(radius**2-(zbound-c_k)**2) + c_i


  n_z = (zbound-c_k)/sqrt((zbound-c_k)**2+(xbound-c_i)**2)
  n_x = (xbound-c_i)/sqrt((zbound-c_k)**2+(xbound-c_i)**2)


  if (i_ele .eq. c_i) then 
  ! Back and front of circle  
    xbound     = c_i
    n_z = (zbound-c_k)/sqrt((zbound-c_k)**2+(xbound-c_i)**2)
    n_z_sign = n_z / abs(n_z)

    z2 = k_ele + n_z_sign
    x2 = i_ele
    z3 = k_ele + n_z_sign
    x3 = i_ele

    weicoef(1) = (k_ele-zbound)/(z2-zbound) 
    weicoef(2) =  0d0

  elseif (k_ele .eq. c_k) then
  ! Top and bot of circle
    n_x_sign = n_x / abs(n_x)

    z2 = k_ele 
    x2 = i_ele + n_x_sign
    z3 = k_ele 
    x3 = i_ele + n_x_sign

    weicoef(1) = (i_ele-xbound)/(x2-xbound)
    weicoef(2) =  0d0

  else  

    n_z_sign = n_z / abs(n_z)
    n_x_sign = n_x / abs(n_x)

    if (abs(i_ele-c_i) .eq. abs(k_ele-c_k)) then
  !  if (abs(n_x) .eq. abs(n_z)) then !only accurate to 3dp so doesnt work
    ! Diagonal point 
      z2 = k_ele + n_z_sign
      x2 = i_ele  
      z3 = k_ele 
      x3 = i_ele + n_x_sign 

    elseif (abs(n_x) .lt. abs(n_z)) then 
    ! Between 10.5 o clock and 1.5 o clock  
    ! Between 4.5 o clock and 7.5 o clock 
      z2 = k_ele + n_z_sign
      x2 = i_ele   
      z3 = k_ele + n_z_sign
      x3 = i_ele + n_x_sign

    elseif (abs(n_x) .gt. abs(n_z)) then 
    ! Between 1.5 o clock and 4.5 o clock 
    ! Between 1.5 o clock and 4.5 o clock 
      z2 = k_ele 
      x2 = i_ele + n_x_sign 
      z3 = k_ele + n_z_sign
      x3 = i_ele + n_x_sign

    else

    write(*,*) 'Mel fucked up'
    stop
 
    end if

    denominator=  (xbound*z2 - x2*zbound - xbound*z3 + x3*zbound + x2*z3 - x3*z2)


    weicoef(1) = -(xbound*z3 - x3*zbound - xbound*k_ele + i_ele*zbound + x3*k_ele - i_ele*z3)/denominator

    weicoef(2) =  (xbound*z2 - x2*zbound - xbound*k_ele + i_ele*zbound + x2*k_ele - i_ele*z2)/denominator 

  end if 

  interpp(1) = z2 
  interpp(2) = x2 
  interpp(3) = z3 
  interpp(4) = x3 


end subroutine


subroutine add_forcing_point(list_ib_xz,list_ib_w_xz,i_ele,k_ele,interpp,weicoef,xbound,zbound)

  use declaration
  implicit none
  integer :: i_ele, k_ele
  integer :: interpp(4)
  real(8) :: weicoef(2), zbound, xbound
  integer :: list_ib_xz(6)
  real(8) :: list_ib_w_xz(5)

  list_ib_xz(1) = 1+mod(i_ele-1d0+dnx,1d0*dnx)
  list_ib_xz(2) = 1+mod(k_ele-1d0+dnz,1d0*dnz)
  list_ib_xz(3) = 1+mod(interpp(2)-1d0+dnx,1d0*dnx)   ! i of fluid point 2
  list_ib_xz(4) = 1+mod(interpp(1)-1d0+dnz,1d0*dnz)   ! k of fluid point 2
  list_ib_xz(5) = 1+mod(interpp(4)-1d0+dnx,1d0*dnx)   ! i of fluid point 3
  list_ib_xz(6) = 1+mod(interpp(3)-1d0+dnz,1d0*dnz)   ! k of fluid point 3

  list_ib_w_xz(1) = weicoef(1) ! weighting of fluid point 2
  list_ib_w_xz(2) = weicoef(2) ! weighting of fluid point 3
  list_ib_w_xz(3) = 1+mod(xbound-1d0+dnx,1d0*dnx) ! i of circular boundary point
  list_ib_w_xz(4) = 1+mod(zbound-1d0+dnz,1d0*dnz) ! k of circular boundary point
  list_ib_w_xz(5) = 1d0 ! Laplacian coefficient 

end subroutine


subroutine add_solid_point(list_ib_xz,list_ib_w_xz,i_ele,k_ele)

  use declaration
  implicit none
  integer :: i_ele, k_ele
  integer :: list_ib_xz(6)
  real(8) :: list_ib_w_xz(5)

  list_ib_xz(1) = 1+mod(i_ele-1d0+dnx,1d0*dnx) 
  list_ib_xz(2) = 1+mod(k_ele-1d0+dnz,1d0*dnz)
  list_ib_xz(3) = 1+mod(i_ele-1d0+dnx,1d0*dnx)   ! i of fluid point 2
  list_ib_xz(4) = 1+mod(k_ele-1d0+dnz,1d0*dnz)   ! k of fluid point 2
  list_ib_xz(5) = 1+mod(i_ele-1d0+dnx,1d0*dnx)   ! i of fluid point 3
  list_ib_xz(6) = 1+mod(k_ele-1d0+dnz,1d0*dnz)   ! k of fluid point 3

  list_ib_w_xz(1) = 0d0 ! weighting of fluid point 2
  list_ib_w_xz(2) = 0d0 ! weighting of fluid point 3
  list_ib_w_xz(3) = 0d0 ! i of circular boundary point
  list_ib_w_xz(4) = 0d0 ! k of circular boundary point
  list_ib_w_xz(5) = 0d0 ! Laplacian coefficient 

end subroutine
