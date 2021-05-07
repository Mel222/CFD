Module Sort
!Merge sort algorithm
!Reference: http://www-h.eng.cam.ac.uk/help/languages/fortran/f90/mergsort.f90

    Implicit None
    Private Merging

Contains

    recursive subroutine RcrsDivide( a, ind )

        implicit none
        integer, intent(inout) :: a(:)
        integer, intent(inout) :: ind(:)

        integer :: a0(size(a))
        integer :: low, high, mid

        low  = lbound(a,1)
        high = ubound(a,1)
        if ( low < high ) then
            mid = (low+high)/2
            call RcrsDivide( a(low  :mid ), ind(low  :mid ) )
            call RcrsDivide( a(mid+1:high), ind(mid+1:high) )
            a0 = a(low:high)
            call Merging( a0(low:mid), a0(mid+1:high), a(low:high), ind(low:high) )
        end if

    end subroutine RcrsDivide

    subroutine Merging( a, b, c, ind )

        implicit none
        integer, intent(in)    :: a(:), b(:)
        integer, intent(out)   :: c(size(a)+size(b))
        integer, intent(inout) :: ind(:)

        integer :: a_ptr, a_high, a_low
        integer :: b_ptr, b_high, b_low
        integer :: c_ptr
        integer :: ind1(size(a)), ind2(size(b))

        a_low  = lbound(a,1)
        a_high = ubound(a,1)
        b_low  = lbound(b,1)
        b_high = ubound(b,1)

        a_ptr = a_low
        b_ptr = b_low
        c_ptr = 1

        ind1(:) = ind(1        :size(a)        )
        ind2(:) = ind(1+size(a):size(a)+size(b))

        do while ( a_ptr <= a_high .AND. b_ptr <= b_high )
            if ( a(a_ptr) <= b(b_ptr) ) then
                c(  c_ptr) = a(   a_ptr)
                ind(c_ptr) = ind1(a_ptr)
                a_ptr = a_ptr + 1
            else
                c(  c_ptr) = b(   b_ptr)
                ind(c_ptr) = ind2(b_ptr)
                b_ptr = b_ptr + 1
            end if
            c_ptr = c_ptr + 1
        end do

        if ( a_ptr > a_high ) then
            c(  c_ptr:) = b(   b_ptr:b_high)
            ind(c_ptr:) = ind2(b_ptr:b_high)
        else
            c(  c_ptr:) = a(   a_ptr:a_high)
            ind(c_ptr:) = ind1(a_ptr:a_high)
        end if

    end subroutine Merging

End Module Sort


Subroutine MergeSort3D( nList, ptList, sortInd, ascend )
!Sorting a set of points with (i,k,j) indices according to the order of j, k and i, successively
!Using merge sort algorithm
Use Sort

    implicit none

    integer, intent(in)            :: nList
    !"nList": Number of points participating in the sorting
    integer, intent(inout)         :: ptList(3,nList)
    !"ptList(1:3,ilist)": i/k/j of the "ilist"-th point
    integer, intent(out)           :: sortInd(nList)
    !"sortInd(ilist)": The pre-sorting index of the post-sorting "ilist"-th point
    logical, intent(in), optional  :: ascend(3)
    !"ascend(1:3)": Whther or not in ascending order of i/k/j; all ture by default

    integer :: metric1D(nList)
    integer :: i_min, k_min, j_min
    integer :: i_max, k_max, j_max
    integer :: indSgn(3)
    integer :: ilist

    indSgn(:) = 1    !Default: sorting in ascending order of i/k/j
    if ( present(ascend) ) then
        if ( .NOT. ascend(1) ) indSgn(1) = -1
        if ( .NOT. ascend(2) ) indSgn(2) = -1
        if ( .NOT. ascend(3) ) indSgn(3) = -1
    end if

    i_min = minval( ptList(1,:) )
    i_max = maxval( ptList(1,:) )
    k_min = minval( ptList(2,:) )
    k_max = maxval( ptList(2,:) )
    j_min = minval( ptList(3,:) )
    j_max = maxval( ptList(3,:) )
    metric1D(:) = &    !Monotonically mapping (i,k,j) to a number for each point
              &   indSgn(3) * ((i_max-i_min+1)*(k_max-k_min+1)) * ptList(3,:) &
              & + indSgn(2) *  (i_max-i_min+1)                  * ptList(2,:) &
              & + indSgn(1) *                                     ptList(1,:)
    
    do ilist = 1, nList
        sortInd(ilist) = ilist    !Initializing "sortInd"
    end do
    
    call RcrsDivide( metric1D, sortInd )
    !Ascendingly sort "metric1D" while record the change of point indexes in "sortInd"

    ptList(:,:) = ptList(:,sortInd)

End Subroutine MergeSort3D