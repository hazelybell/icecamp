program fem
  integer,parameter :: nr_v = 2 ! Number of variables (U,V)
  integer,parameter :: dp = kind(1.0d0)
  integer,parameter :: m = 11, n = 11 ! Grid size
  integer :: mi, ni, ti, tj, bi, e, pi, ki, kj, iter, ierr
  integer, parameter :: nr_p = m * n ! Number of points
  real(dp), dimension(nr_p,2) :: p ! List of points
  integer, parameter :: nr_t = (m-1) * (n-1) * 2 ! Number of triangles
  integer, dimension(nr_t,3) :: t ! List of triangles
  integer, parameter :: max_neighbors = 6 ! Maximum number of points adjacent to a single point
  integer, parameter :: nr_b = 2*n+2*m - 4 ! Number of points on the boundary
  integer, dimension(nr_b) :: b
  integer :: b_nw, b_ne, b_sw, b_se
  integer, dimension(m-2) :: b_n, b_s
  integer, dimension(n-2) :: b_e, b_w
  real(dp), dimension(3,3) :: Pe, C, Ke
  real(dp) :: detPe, Area, invDet, converr
  integer, parameter :: nr_k_nz = (nr_p + nr_t * 3 + nr_b) ! * 6
  real(dp), dimension(nr_k_nz) :: vK, vKb
  integer, dimension(nr_k_nz) :: iK
  integer, dimension(nr_p+1) :: jK
  real(dp), dimension(nr_p) :: F, Fb, U
  real(dp), dimension(nr_p*9) :: rwork
  integer, dimension(100) :: iwork
! For Stokes, the BMP is:
! 1 0 1
! 0 1 1
! 1 1 0
  integer, dimension(nr_v,nr_v), parameter :: bmp = reshape((/1, 0, 1, 0, 1, 1, 1, 1, 0/), shape(bmp)) ! Block Matrix Pattern
! This must be multiplied by the number of non-zero elements in the BMP
  integer :: vi, vj ! For loops over variables
  
  print *, 'm=', m, ' n=', n

! Make grid
  do ni = 0, (m-1)
    do mi = 0, (n-1)
      p(ni*(m)+mi+1,:) = (/(dble(mi)/dble(m-1)),(dble(ni)/dble(n-1))/)
    end do
  end do
!   do e = 1, nr_p
!     print *, p(e,:)
!   end do

! Make triangle list
  do ni = 1, n-1
    do mi = 1, m-1
      t((ni-1)*(m-1)+mi,:) = (/((mi-1)*n)+ni, ((mi-1)*n)+ni+1, (mi*n)+ni+1/)
      t((ni-1)*(m-1)+mi+((m-1)*(n-1)),:) = (/((mi-1)*n)+ni, (mi*n)+ni+1, (mi*n)+ni/)
    end do
  end do

!   do ti = 1, nr_t
!     print *, t(ti,:)
!   end do

! Make boundary list
  b_nw = 1
  b_ne = m
  b_n = (/(mi, mi = 2, m-1)/)
  b_s = (/(mi+(n-1)*m, mi = 2, m-1)/)
  b_sw = (n-1)*m+1
  b_se = (n-1)*m+m
  b_w = (/((ni-1)*m+1, ni = 2, n-1)/)
  b_e = (/((ni)*m, ni = 2, n-1)/)
!   b=(/(mi, mi = 1, m), ((ni-1)*m+1, ni = 2, n), ((ni-1)*m, ni = 3, n+1), (mi+(n-1)*m, mi = 2, m-1) /)
  b = (/b_nw, b_n, b_ne, b_w, b_sw, b_e, b_se, b_s/)
!   print *, b
!   do bi = 1, nr_b
!     print *, t(bi,:) 
!   end do

  call init_k_sparsity

  do e = 1, nr_t
    Pe(1,:) = (/ 1.0d0, p(t(e,1),:) /)
    Pe(2,:) = (/ 1.0d0, p(t(e,2),:) /)
    Pe(3,:) = (/ 1.0d0, p(t(e,3),:) /)
    detPe = Pe(1,1)*(Pe(2,2)*Pe(3,3) - Pe(2,3)*Pe(3,2)) &
                + Pe(1,2)*(Pe(2,3)*Pe(3,1) - Pe(2,1)*Pe(3,3)) &
                + Pe(1,3)*(Pe(2,1)*Pe(3,2) - Pe(2,2)*Pe(3,1))
    Area = dabs(detPe)/2.0d0
    ! Compute inverse of Pe
    invDet = 1.0d0/detPe
    C(1,1) = invDet * (Pe(2,2)*Pe(3,3) - Pe(2,3)*Pe(3,2))
    C(1,2) = invDet * (Pe(1,3)*Pe(3,2) - Pe(1,2)*Pe(3,3))
    C(1,3) = invDet * (Pe(1,2)*Pe(2,3) - Pe(1,3)*Pe(2,2))
    C(2,1) = invDet * (Pe(2,3)*Pe(3,1) - Pe(2,1)*Pe(3,3))
    C(2,2) = invDet * (Pe(1,1)*Pe(3,3) - Pe(1,3)*Pe(3,1))
    C(2,3) = invDet * (Pe(1,3)*Pe(2,1) - Pe(1,1)*Pe(2,3))
    C(3,1) = invDet * (Pe(2,1)*Pe(3,2) - Pe(2,2)*Pe(3,1))
    C(3,2) = invDet * (Pe(3,1)*Pe(1,2) - Pe(1,1)*Pe(3,2))
    C(3,3) = invDet * (Pe(1,1)*Pe(2,2) - Pe(1,2)*Pe(2,1))
    ! C[2,1] is the partial phi1 / partial x
    ! C[3,1] is the partial phi1 / partial x 
    ! etc. these are the same for V instead of phi!!!
    Ke(1,1) = Area * (C(2,1)*C(2,1) + C(3,1)*C(3,1))
    Ke(1,2) = Area * (C(2,1)*C(2,2) + C(3,1)*C(3,2))
    Ke(1,3) = Area * (C(2,1)*C(2,3) + C(3,1)*C(3,3))
    Ke(2,1) = Area * (C(2,2)*C(2,1) + C(3,2)*C(3,1))
    Ke(2,2) = Area * (C(2,2)*C(2,2) + C(3,2)*C(3,2))
    Ke(2,3) = Area * (C(2,2)*C(2,3) + C(3,2)*C(3,3))
    Ke(3,1) = Area * (C(2,3)*C(2,1) + C(3,3)*C(3,1))
    Ke(3,2) = Area * (C(2,3)*C(2,2) + C(3,3)*C(3,2))
    Ke(3,3) = Area * (C(2,3)*C(2,3) + C(3,3)*C(3,3))
    ! So the non-zero entries are each line twice and each point. There are nr_p points and nr_t * 3 lines so there are nr_p + nr_t * 6 non-sparse entries
    ! For a node pi, the relevant entries are (in slap column order) the diagonal (pi, pi), nw, n, w, e, s, sw
    ! Boundary nodes only have 4 or 3 points each
    do ti = 1, 3
      do tj = 1, 3
        ki = k_index(t(e,ti),t(e,tj))
        if (ki .eq. nr_k_nz+1) then
          print *, e, ti, tj, t(e,ti), t(e,tj)
        end if
        vK(ki) = vK(ki) + Ke(ti,tj)
      end do
      F(t(e,ti)) = F(t(e,ti)) + (Area / 3)
    end do
  end do
!   print *, Pe(1,:)
!   print *, Pe(2,:)
!   print *, Pe(3,:)
!   print *, Area
!   print *, C(1,:)
!   print *, C(2,:)
!   print *, C(3,:)
!   print *, Ke(1,:)
!   print *, Ke(2,:)
!   print *, Ke(3,:)
!   kj = 1
!   do ki = 1, nr_k_nz
!     if (ki .eq. jK(kj+1)) kj = kj + 1
!     print *, iK(ki), kj, vK(ki)
!   end do
  print *, 'done assy'
  F(b) = 0.0d0
   do ki = 1, nr_p
     do bi = 1, nr_b
       if (k_index(ki, b(bi)) .ne. nr_k_nz+1) then
        vK(k_index(ki, b(bi))) = 0.0d0
        vK(k_index(b(bi), ki)) = 0.0d0
       end if
     end do
   end do
  do bi = 1, nr_b
    vk(k_index(b(bi), b(bi))) = 1.0d0
  end do
  print *, vk(k_index(121,120))
  print *, vk(k_index(121,121))
  Fb = F
  vKb = vK
  U = 1.0d0
  CALL DSDBCG(nr_p, Fb, U, nr_k_nz, iK, jK, vKb, 0, 1, 1.0d-8, 10000, iter, converr, ierr, 6, rwork, nr_p * 9, iwork, 100)
  print *, iter, converr, ierr

  do ni = 1, n
    print *, U((ni-1)*m+1:ni*m-1)
  end do
  
  contains
    function k_index(i, j)
      integer, intent(in) :: i, j
      integer :: ii
      ii = jK(j)
      do while (ii .le. nr_k_nz .and. ii .lt. jK(j+1))
        if (iK(ii) .eq. i) then
          k_index = ii
          return
        end if
        ii = ii + 1
      end do
!       print *, 'error', i, j
      k_index = nr_k_nz+1
    end function

    subroutine init_k_sparsity()
      integer, dimension(nr_p) :: nr_neighbors = 0
      integer, dimension(nr_p, max_neighbors) :: neighbors = 0
      integer :: e, pi, ki, ni

      do e = 1, nr_t
        call maybe_add_neighbor(t(e,1),t(e,2),nr_neighbors(t(e,1)),neighbors(t(e,1),:))
        call maybe_add_neighbor(t(e,1),t(e,3),nr_neighbors(t(e,1)),neighbors(t(e,1),:))
        call maybe_add_neighbor(t(e,2),t(e,1),nr_neighbors(t(e,2)),neighbors(t(e,2),:))
        call maybe_add_neighbor(t(e,2),t(e,3),nr_neighbors(t(e,2)),neighbors(t(e,2),:))
        call maybe_add_neighbor(t(e,3),t(e,1),nr_neighbors(t(e,3)),neighbors(t(e,3),:))
        call maybe_add_neighbor(t(e,3),t(e,2),nr_neighbors(t(e,3)),neighbors(t(e,3),:))
      end do
!       do pi = 1, nr_p
!          print *, pi, ":", neighbors(pi, :)
!       end do
      ki = 1
      do pi = 1, nr_p
        jK(pi) = ki
        iK(ki) = pi
        vK(ki) = 0.0d0
        ki = ki + 1
        do ni = 1, max_neighbors
          if (neighbors(pi, ni) .ne. 0 .and. neighbors(pi, ni) .ne. pi) then
            print *, ki, pi, ni, neighbors(pi, ni)
            iK(ki) = neighbors(pi, ni)
            vK = 0.0d0
            ki = ki + 1
          end if
        end do
      end do
      jK(nr_p+1)=nr_k_nz+1
      do pi = 1, nr_p
        ki = jK(pi)
        do while (ki .le. nr_k_nz .and. ki .lt. jK(pi+1))
          print *, ki, iK(ki), pi, vK(ki), jK(pi+1)
          ki = ki + 1
        end do
      end do
      print *, jK(nr_p), nr_p, nr_k_nz
    end subroutine

    subroutine maybe_add_neighbor(node, neighbor, nr, list)
      integer, intent(in) :: node, neighbor
      integer, intent(inout) :: nr
      integer, dimension(max_neighbors), intent(inout) :: list
      integer :: ni
      integer :: insert_at
      insert_at = 0
      do ni = 1, nr
        if (list(ni) .eq. neighbor) then
          return
        end if
        if ((list(ni) .gt. neighbor) .and. (insert_at .eq. 0)) then
          insert_at = ni
          exit
        end if
      end do
      nr = nr + 1
      if (nr .gt. max_neighbors) then
        print *, "Too many neighbors, increase max_neighbors." 
        call exit(1)
      end if
      if (insert_at .eq. 0) then
        insert_at = nr
      end if
      if (insert_at .lt. nr) then
        do ni =  nr-1, insert_at, -1
          list(ni+1) = list(ni)
        end do
      end if
      list(insert_at) = neighbor
    end subroutine
end program


