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
  integer,parameter :: nr_v = 2 ! Number of variables (U,V)
  integer,parameter :: dp = kind(1.0d0)
  integer,parameter :: n = 11, m = 11  ! Grid size n = rows, m = cols
  integer :: mi, ni, ti, tj, bi, bj, e, pi, ki, kj, iter, ierr
  integer, parameter :: nr_p = m * n ! Number of points
  real(dp), dimension(nr_p,2) :: p ! List of points
  integer, parameter :: nr_t = (m-1) * (n-1) * 2 ! Number of triangles
  integer, dimension(nr_t,3) :: t ! List of triangles
  integer, parameter :: max_neighbors = 7 ! Maximum number of points adjacent to a single point
  integer, parameter :: nr_b = 2*n+2*m - 4 ! Number of points on the boundary
  integer, dimension(nr_b) :: b
  integer, parameter ::  nr_b_we = 2*n
  integer :: b_nw, b_ne, b_sw, b_se
  integer, dimension(m-2) :: b_n, b_s
  integer, dimension(n-2) :: b_e, b_w
  integer, dimension(nr_b_we) :: b_we
  real(dp), dimension(3,3) :: Pe, C, Ke
  real(dp) :: detPe, Area, invDet, converr
  integer, parameter :: nr_bmp_nz = 3 ! Number of block matrix pattern elements that are non-zero
  integer, parameter :: linear_order = nr_p * nr_v
  integer, parameter :: nr_k_nz = (nr_p + nr_t * 3 + nr_b) * nr_bmp_nz ! Number of entries in the K matrix that are non-zero
  real(dp), dimension(nr_k_nz) :: vK, vKb
  integer, dimension(nr_k_nz) :: iK
  integer, dimension(linear_order+1) :: jK
  real(dp), dimension(linear_order) :: F, Fb, U
  real(dp), dimension(linear_order*9) :: rwork
  integer, dimension(100) :: iwork
! For ssa, the BMP is:
! 1 12
! 0  2
  integer, dimension(nr_v,nr_v), parameter :: bmp = reshape((/1, 0, 12, 2/), shape(bmp)) ! Block Matrix Pattern
! This must be multiplied by the number of non-zero elements in the BMP
  integer :: vi, vj ! For loops over variables

!compute row/col index from variable interleave
#define VPtoI(v, p) ((p-1)*nr_v+v)
!compute variable interleave from row/col index
#define ItoV(i) (MOD(i-1,nr_v)+1)
#define ItoP(i) (((i-1)/nr_v)+1)

  print *, 'm=', m, ' n=', n

! Make grid
  do ni = 0, (n-1) ! rows
    do mi = 0, (m-1) ! cols
      p(ni*(m)+mi+1,:) = (/(dble(mi)/dble(m-1)),(dble(ni)/dble(n-1))/)
    end do
  end do
!   do e = 1, nr_p
!     print *, p(e,:)
!   end do

! Make triangle list
  do ni = 1, n-1 ! rows
    do mi = 1, m-1 ! cols
      t((ni-1)*(m-1)+mi,:) = (/((ni-1)*m)+mi, ((ni-1)*m)+mi+1, (ni*m)+mi+1/)
      t((ni-1)*(m-1)+mi+((m-1)*(n-1)),:) = (/((ni-1)*m)+mi, (ni*m)+mi+1, (ni*m)+mi/)
    end do
  end do

  print *, "Triangles"
  do ti = 1, nr_t
    print *, t(ti,:)-1
  end do
  print *, "End Triangles"

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
  b_we = (/b_nw, b_w, b_sw, b_ne, b_e, b_se/)
!   print *, b
!   do bi = 1, nr_b
!     print *, t(bi,:) 
!   end do

  call init_k_sparsity
  
  F=0d0

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
        do vi = 1, nr_v
          ki = k_index(VPtoI(vi, t(e,ti)),VPtoI(vi, t(e,tj)))
          if (ki .eq. nr_k_nz+1) then
            print *, e, ti, tj, t(e,ti), t(e,tj)
          end if
            vK(ki) = vK(ki) - Ke(ti,tj)
        end do
      end do
!       do vi = 1, nr_v
        F(VPtoI(2, t(e,ti))) = F(VPtoI(2, t(e,ti))) + 24 * (Area / 3)
!         F(VPtoI(1, t(e,ti))) = 0
!       end do
    end do
!     print *, t(e, :)
!     print *, "(", Pe(1,2:3), "), (", Pe(2,2:3), "), (", Pe(3,2:3), ")"
    
    ! Compute the integral for the vw term analytically
    Ke(1,1) = int_lintimeslin(C(:,(/1,1/)), Pe(:,2:3))
    Ke(1,2) = int_lintimeslin(C(:,(/1,2/)), Pe(:,2:3))
    Ke(1,3) = int_lintimeslin(C(:,(/1,3/)), Pe(:,2:3))
    Ke(2,1) = int_lintimeslin(C(:,(/2,1/)), Pe(:,2:3))
    Ke(2,2) = int_lintimeslin(C(:,(/2,2/)), Pe(:,2:3))
    Ke(2,3) = int_lintimeslin(C(:,(/2,3/)), Pe(:,2:3))
    Ke(3,1) = int_lintimeslin(C(:,(/3,1/)), Pe(:,2:3))
    Ke(3,2) = int_lintimeslin(C(:,(/3,2/)), Pe(:,2:3))
    Ke(3,3) = int_lintimeslin(C(:,(/3,3/)), Pe(:,2:3))
!         print *, "Ke"
!         print *, Ke(1,:)
!         print *, Ke(2,:)
!         print *, Ke(3,:)
    do ti = 1, 3
      do tj = 1, 3
        ki = k_index(VPtoI(1, t(e,ti)),VPtoI(2, t(e,tj)))
        if (ki .eq. nr_k_nz+1) then
          print *, e, ti, tj, t(e,ti), t(e,tj), VPtoI(1, t(e,ti)),VPtoI(2, t(e,tj))
        end if
        vK(ki) = vK(ki) - Ke(ti,tj)
      end do
    end do

  end do
!   print *, F(:)
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
! BOUNDARY CONDITIONS !!!!!!
  print *, b_we
  print *, nr_b_we
  F(VPtoI(1, b_we)) = 0.0d0
  F(VPtoI(2, b_we)) = 2.0d0
  do ki = 1, nr_p
    do bi = 1, nr_b_we
      do vi = 1, nr_v
        do vj = 1, nr_v
        ! wtf strang
!           if (k_index(VPtoI(vi, ki), VPtoI(vj, b_we(bi))) .ne. nr_k_nz+1) then
!             vK(k_index(VPtoI(vi, ki), VPtoI(vj, b_we(bi)))) = 0.0d0
!           end if
          if (k_index(VPtoI(vi, b_we(bi)), VPtoI(vj, ki)) .ne. nr_k_nz+1) then
            vK(k_index(VPtoI(vi, b_we(bi)), VPtoI(vj, ki))) = 0.0d0
          end if
        end do
      end do
     end do
   end do
  do bi = 1, nr_b_we
    do vi = 1, nr_v
      vk(k_index(VPtoI(vi, b_we(bi)), VPtoI(vi, b_we(bi)))) = 1.0d0
    end do
  end do
!   print *, vk(k_index(VPtoI(1, 8),VPtoI(1, 9)))
!   print *, vk(k_index(VPtoI(1, 8),VPtoI(1, 8)))
  Fb = F
  vKb = vK
  U = 1.0d0
  CALL DSDBCG(linear_order, Fb, U, nr_k_nz, iK, jK, vKb, &
    0, 1, 1.0d-7, 10000, iter, converr, ierr, 6, rwork, linear_order * 9, iwork, 100)
  print *, iter, converr, ierr

!   do vi = 1, nr_v
!   print *, "Variable ", vi
!     do ni = 1, n
!       print *, ItoV((ni-1)*m*nr_v+vi), ItoP((ni-1)*m*nr_v+vi), &
!         ItoV(ni*m*nr_v+vi-2), ItoP(ni*m*nr_v+vi-2)
!       print *, U((ni-1)*m*nr_v+vi:ni*m*nr_v+vi-1:nr_v)
!     end do
!   end do
  
  do vi = 1, nr_v
    print *, "Variable ", vi
    do pi = 1, nr_p
      print *, p(pi,1:2), U(VPtoI(vi, pi)), F(VPtoI(vi, pi)), vK(k_index(VPtoI(vi, pi), VPtoI(vi, pi)))
    end do
    print *, "End Variable"
  end do
  
  do vi = 1, nr_v
    print *, minval(U(vi:nr_v*nr_p:nr_v)), maxval(U(vi:nr_v*nr_p:nr_v))
  end do
  contains
    function k_index(i, j)
      integer, intent(in) :: i, j
      integer :: ii
!       print *, i, j
      ii = jK(j)
      do while (ii .le. nr_k_nz .and. ii .lt. jK(j+1))
        if (iK(ii) .eq. i) then
          k_index = ii
          return
        end if
        ii = ii + 1
      end do
!      print *, 'error', i, j
      k_index = nr_k_nz+1
    end function
    
    real(dp) function int_lintimeslin(coeffs, domain)
      real(dp), dimension(3,2), intent(in) :: coeffs, domain
      real(dp), dimension(2) :: temp
      real(dp) :: accumulator, ma, oa, mb, ob, xmin, xmax, signfix
      real(dp) :: av, bv, cv, aw, bw, cw
      real(dp), dimension(3,2) :: sorted
      av = coeffs(1,1)
      bv = coeffs(2,1)
      cv = coeffs(3,1)
      aw = coeffs(1,2)
      bw = coeffs(2,2)
      cw = coeffs(3,2)
      sorted = domain
      ! sort triangle coords by x
      if (sorted(1,1) .gt. sorted(2,1)) then
        if (sorted(1,1) .gt. sorted(3,1)) then
          temp = sorted(3,:)
          sorted(3,:) = sorted(1,:)
          sorted(1,:) = temp
          if (sorted(1,1) .gt. sorted(2,1)) then
            temp = sorted(2,:)
            sorted(2,:) = sorted(1,:)
            sorted(1,:) = temp
          end if
        else
          temp = sorted(2,:)
          sorted(2,:) = sorted(1,:)
          sorted(1,:) = temp
          if (sorted(1,1) .gt. sorted(3,1)) then
            temp = sorted(3,:)
            sorted(3,:) = sorted(1,:)
            sorted(1,:) = temp
          end if
        end if
      else
        if (sorted(2,1) .gt. sorted(3,1)) then
          if (sorted(1,1) .gt. sorted(3,1)) then
            temp = sorted(1,:)
            sorted(1,:) = sorted(3,:)
            sorted(3,:) = temp
          else
            temp = sorted(2,:)
            sorted(2,:) = sorted(3,:)
            sorted(3,:) = temp
          end if
        end if
      end if
      ! Compute the integral in possibly two parts
      if (sorted(3,1) .le. sorted(1,1)) then
        call exit(2)
      end if
      ma = (sorted(3,2) - sorted(1,2))/(sorted(3,1) - sorted(1,1))
      oa = sorted(1,2) - ma * sorted(1,1)
      ! If we have a downward pointing triangle the intergral would be wrong...
      if (sorted(2,2) .gt. (oa + ma * sorted(2,1))) then
        signfix = 1.0d0
      else
        signfix = -1.0d0
      end if
      accumulator = 0.0d0
      if (sorted(2,1) .gt. sorted(1,1)) then
        mb = (sorted(2,2) - sorted(1,2))/(sorted(2,1) - sorted(1,1))
        ob = sorted(1,2) - mb * sorted(1,1)
        xmin = sorted(1,1)
        xmax = sorted(2,1)
        accumulator = accumulator &
          + xmax*(av*aw*ob - av*aw*oa + av*cw*ob**2/2 + aw*cv*ob**2/2 &
            - av*cw*oa**2/2 - aw*cv*oa**2/2 - cv*cw*oa**3/3 &
            + cv*cw*ob**3/3) &
          + xmax**2*(av*aw*mb/2 + av*bw*ob/2 + aw*bv*ob/2 &
            - av*aw*ma/2 - av*bw*oa/2 - aw*bv*oa/2 + av*cw*mb*ob/2 &
            + aw*cv*mb*ob/2 - av*cw*ma*oa/2 - aw*cv*ma*oa/2 &
            - bv*cw*oa**2/4 - bw*cv*oa**2/4 + bv*cw*ob**2/4 &
            + bw*cv*ob**2/4 + cv*cw*mb*ob**2/2 - cv*cw*ma*oa**2/2) &
          + xmax**3*(-av*bw*ma/3 - aw*bv*ma/3 - bv*bw*oa/3 &
            + av*bw*mb/3 + aw*bv*mb/3 + bv*bw*ob/3 - bv*cw*ma*oa/3 &
            - bw*cv*ma*oa/3 + bv*cw*mb*ob/3 + bw*cv*mb*ob/3 &
            - av*cw*ma**2/6 - aw*cv*ma**2/6 + av*cw*mb**2/6 &
            + aw*cv*mb**2/6 - cv*cw*oa*ma**2/3 + cv*cw*ob*mb**2/3) &
          + xmax**4*(-bv*bw*ma/4 + bv*bw*mb/4 - bv*cw*ma**2/8 &
            - bw*cv*ma**2/8 + bv*cw*mb**2/8 + bw*cv*mb**2/8 &
            - cv*cw*ma**3/12 + cv*cw*mb**3/12) &
          - xmin*(av*aw*ob - av*aw*oa + av*cw*ob**2/2 + aw*cv*ob**2/2 &
            - av*cw*oa**2/2 - aw*cv*oa**2/2 - cv*cw*oa**3/3 &
            + cv*cw*ob**3/3) &
          - xmin**2*(av*aw*mb/2 + av*bw*ob/2 + aw*bv*ob/2 &
            - av*aw*ma/2 - av*bw*oa/2 - aw*bv*oa/2 + av*cw*mb*ob/2 &
            + aw*cv*mb*ob/2 - av*cw*ma*oa/2 - aw*cv*ma*oa/2 &
            - bv*cw*oa**2/4 - bw*cv*oa**2/4 + bv*cw*ob**2/4 &
            + bw*cv*ob**2/4 + cv*cw*mb*ob**2/2 - cv*cw*ma*oa**2/2) &
          - xmin**3*(-av*bw*ma/3 - aw*bv*ma/3 - bv*bw*oa/3 &
            + av*bw*mb/3 + aw*bv*mb/3 + bv*bw*ob/3 - bv*cw*ma*oa/3 &
            - bw*cv*ma*oa/3 + bv*cw*mb*ob/3 + bw*cv*mb*ob/3 &
            - av*cw*ma**2/6 - aw*cv*ma**2/6 + av*cw*mb**2/6 &
            + aw*cv*mb**2/6 - cv*cw*oa*ma**2/3 + cv*cw*ob*mb**2/3) &
          - xmin**4*(-bv*bw*ma/4 + bv*bw*mb/4 - bv*cw*ma**2/8 &
            - bw*cv*ma**2/8 + bv*cw*mb**2/8 + bw*cv*mb**2/8 &
            - cv*cw*ma**3/12 + cv*cw*mb**3/12)
      end if
      if (sorted(3,1) .gt. sorted(2,1)) then
        mb = (sorted(3,2) - sorted(2,2))/(sorted(3,1) - sorted(2,1))
        ob = sorted(2,2) - mb * sorted(2,1)
        xmin = sorted(2,1)
        xmax = sorted(3,1)
        accumulator = accumulator &
          + xmax*(av*aw*ob - av*aw*oa + av*cw*ob**2/2 + aw*cv*ob**2/2 &
            - av*cw*oa**2/2 - aw*cv*oa**2/2 - cv*cw*oa**3/3 &
            + cv*cw*ob**3/3) &
          + xmax**2*(av*aw*mb/2 + av*bw*ob/2 + aw*bv*ob/2 &
            - av*aw*ma/2 - av*bw*oa/2 - aw*bv*oa/2 + av*cw*mb*ob/2 &
            + aw*cv*mb*ob/2 - av*cw*ma*oa/2 - aw*cv*ma*oa/2 &
            - bv*cw*oa**2/4 - bw*cv*oa**2/4 + bv*cw*ob**2/4 &
            + bw*cv*ob**2/4 + cv*cw*mb*ob**2/2 - cv*cw*ma*oa**2/2) &
          + xmax**3*(-av*bw*ma/3 - aw*bv*ma/3 - bv*bw*oa/3 &
            + av*bw*mb/3 + aw*bv*mb/3 + bv*bw*ob/3 - bv*cw*ma*oa/3 &
            - bw*cv*ma*oa/3 + bv*cw*mb*ob/3 + bw*cv*mb*ob/3 &
            - av*cw*ma**2/6 - aw*cv*ma**2/6 + av*cw*mb**2/6 &
            + aw*cv*mb**2/6 - cv*cw*oa*ma**2/3 + cv*cw*ob*mb**2/3) &
          + xmax**4*(-bv*bw*ma/4 + bv*bw*mb/4 - bv*cw*ma**2/8 &
            - bw*cv*ma**2/8 + bv*cw*mb**2/8 + bw*cv*mb**2/8 &
            - cv*cw*ma**3/12 + cv*cw*mb**3/12) &
          - xmin*(av*aw*ob - av*aw*oa + av*cw*ob**2/2 + aw*cv*ob**2/2 &
            - av*cw*oa**2/2 - aw*cv*oa**2/2 - cv*cw*oa**3/3 &
            + cv*cw*ob**3/3) &
          - xmin**2*(av*aw*mb/2 + av*bw*ob/2 + aw*bv*ob/2 &
            - av*aw*ma/2 - av*bw*oa/2 - aw*bv*oa/2 + av*cw*mb*ob/2 &
            + aw*cv*mb*ob/2 - av*cw*ma*oa/2 - aw*cv*ma*oa/2 &
            - bv*cw*oa**2/4 - bw*cv*oa**2/4 + bv*cw*ob**2/4 &
            + bw*cv*ob**2/4 + cv*cw*mb*ob**2/2 - cv*cw*ma*oa**2/2) &
          - xmin**3*(-av*bw*ma/3 - aw*bv*ma/3 - bv*bw*oa/3 &
            + av*bw*mb/3 + aw*bv*mb/3 + bv*bw*ob/3 - bv*cw*ma*oa/3 &
            - bw*cv*ma*oa/3 + bv*cw*mb*ob/3 + bw*cv*mb*ob/3 &
            - av*cw*ma**2/6 - aw*cv*ma**2/6 + av*cw*mb**2/6 &
            + aw*cv*mb**2/6 - cv*cw*oa*ma**2/3 + cv*cw*ob*mb**2/3) &
          - xmin**4*(-bv*bw*ma/4 + bv*bw*mb/4 - bv*cw*ma**2/8 &
            - bw*cv*ma**2/8 + bv*cw*mb**2/8 + bw*cv*mb**2/8 &
            - cv*cw*ma**3/12 + cv*cw*mb**3/12)
      end if
      int_lintimeslin = accumulator * signfix
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
      do pi = 1, nr_p
        call maybe_add_neighbor(pi,pi,nr_neighbors(pi),neighbors(pi,:))
      end do
!       do pi = 1, nr_p
!          print *, pi, ":", neighbors(pi, :)
!       end do
      ki = 1
      do pi = 1, nr_p
        do vj = 1, nr_v
              ! Diagonal entry, element and variable both equal
              jK(VPtoI(vj, pi)) = ki
!               print *, ki, pi, "-----------", "------------", vj, vj, VPtoI(vj, pi)
              iK(ki) = VPtoI(vj, pi)
              vK(ki) = 0.0d0
              ki = ki + 1
              ! Off-diagonal entries in order
              do ni = 1, max_neighbors
                do vi = 1, nr_v
                  if (neighbors(pi, ni) .ne. 0 .and. (neighbors(pi, ni) .ne. pi &
                    .or. vi .ne. vj) .and. bmp(vi,vj) .ne. 0) then
!                     print *, ki, pi, ni, neighbors(pi, ni), vi, vj, VPtoI(vi, neighbors(pi, ni))
                    iK(ki) = VPtoI(vi, neighbors(pi, ni))
                    vK = 0.0d0
                    ki = ki + 1
                  end if
                end do
              end do
        end do
      end do
      if (ki .ne. nr_k_nz+1) then
        print *, ki, nr_k_nz
        call exit(3)
      end if
      jK(linear_order+1)=nr_k_nz+1
!        do pi = 1, linear_order
!          ki = jK(pi)
!          if (ki .eq. 0) then
!            print *, pi, jK(pi)
!            call exit(1)
!          end if
!          do while (ki .le. nr_k_nz .and. ki .lt. jK(pi+1))
!            print *, ki, iK(ki), pi, vK(ki), jK(pi+1)
!            ki = ki + 1
!          end do
!        end do
!        print *, jK(nr_p), nr_p, nr_k_nz
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



