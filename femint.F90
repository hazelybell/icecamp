!    Prototype Fortran Code
!    Copyright (C) Joshua Charles Campbell 2012

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
program fem
  integer,parameter :: dp = kind(1.0d0)
  integer,parameter :: m = 11, n = 11
  integer :: mi, ni, ti, tj, bi, e, pi, ki, kj, iter, ierr
  integer, parameter :: nr_p = m * n
  real(dp), dimension(nr_p,2) :: p
  integer, parameter :: nr_t = (m-1) * (n-1) * 2
  integer, dimension(nr_t,3) :: t
  integer, parameter :: nr_b = 2*n+2*m - 4
  integer, dimension(nr_b) :: b
  integer :: b_nw, b_ne, b_sw, b_se
  integer, dimension(m-2) :: b_n, b_s
  integer, dimension(n-2) :: b_e, b_w
  real(dp), dimension(3,3) :: Pe, C, Ke
  real(dp) :: detPe, Area, invDet, converr
  integer, parameter :: nr_k_nz = nr_p + nr_t * 3 + nr_b
  real(dp), dimension(nr_k_nz) :: vK, vKb
  integer, dimension(nr_k_nz) :: iK
  integer, dimension(nr_p+1) :: jK
  real(dp), dimension(nr_p) :: F, Fb, U
  real(dp), dimension(nr_p*9) :: rwork
  integer, dimension(100) :: iwork

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

!   do ti = 1, nr_t7
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

! Initialize K matrix
  ki = 1
  kj = 1
  ! So the non-zero entries are each line twice and each point. There are nr_p points and nr_t * 3 + nr_b lines so there are nr_p + nr_t * 3 + nr_b non-sparse entries
  ! For a node pi, the relevant entries are (in slap column order) the diagonal (pi, pi), nw, n, w, e, s, sw
  ! Boundary nodes only have 4 or 3 points each
  do pi = 1, nr_p
    iK(ki) = pi
    jK(kj) = ki
    vK(ki) = 0.0d0
    ki = ki + 1
    kj = kj + 1
    if (.not. ((pi .ge. b_nw .and. pi .le. b_ne) .or. (mod(pi-1, m) .eq. 0))) then
      iK(ki) = pi - m - 1
      vK(ki) = 0.0d0
      ki = ki + 1
    end if
    if (.not. ((pi .ge. b_nw .and. pi .le. b_ne))) then
      iK(ki) = pi - m
      vK(ki) = 0.0d0
      ki = ki + 1
    end if
    if (.not. (mod(pi-1, m) .eq. 0)) then
      iK(ki) = pi - 1
      vK(ki) = 0.0d0
      ki = ki + 1
    end if
    if (.not. (mod(pi, m) .eq. 0)) then
      iK(ki) = pi + 1
      vK(ki) = 0.0d0
      ki = ki + 1
    end if
    if (.not. ((pi .ge. b_sw .and. pi .le. b_se))) then
      iK(ki) = pi + m
      vK(ki) = 0.0d0
      ki = ki + 1
    end if
    if (.not. ((pi .ge. b_sw .and. pi .le. b_se) .or. (mod(pi, m) .eq. 0))) then
      iK(ki) = pi + m + 1
      vK(ki) = 0.0d0
      ki = ki + 1
    end if
  end do
  jK(nr_p+1)=ki
  print *, nr_k_nz, ki

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
       vK(k_index(ki, b(bi))) = 0.0d0
       vK(k_index(b(bi), ki)) = 0.0d0
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
    print *, U((ni-1)*m:ni*m-1)
  end do
  
  contains
    function k_index(i, j)
      integer :: i, j, ii, ans
      ii = jK(j)
      do while (ii .le. nr_k_nz .and. ii .le. jK(j+1))
        if (iK(ii) .eq. i) then
          k_index = ii
          return
        end if
        ii = ii + 1
      end do
!       print *, 'error', i, j
      k_index = nr_k_nz+1
    end function
end program