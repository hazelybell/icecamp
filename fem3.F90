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
! program fem3
  integer,parameter :: dp = kind(1.0d0)
  ! Number of nodes along x, y, z axis (z is height, opengl is for losers)
  integer,parameter :: m = 11, n = 11, o = 11
  integer :: mi, ni, oi, ti, tj, bi, e, pi, ki, kj, iter, ierr, actual_nz
  ! Number of nodes
  integer, parameter :: nr_p = m * n * o
  real(dp), dimension(nr_p,3) :: p
  integer, parameter :: nr_t = (m-1) * (n-1) * (o-1) * 5
  ! tetrahedrons
  integer, dimension(nr_t,4) :: t
  ! Number of boundary nodes
  integer, parameter :: nr_b = 2*m*n + 2*m*o + 2*n*o - 4*m - 4*n - 4*o + 8
  integer, dimension(nr_b) :: b
  ! North west top, north east top, ...
!   integer :: b_nwt, b_net, b_swt, b_set, b_nwb, b_neb, b_swb, b_seb
!   integer, dimension(o-2) :: b_nw, b_ne, b_sw, b_se
!   integer, dimension(m-2) :: b_nt, b_st, b_nb, b_sb
!   integer, dimension(n-2) :: b_et, b_wt, b_eb, b_wb
!   integer, dimension((m-2)*(n-2)) :: b_t, b_b
!   integer, dimension((m-2)*(o-2)) :: b_e, b_w
!   integer, dimension((o-2)*(n-2)) :: b_n, b_s
  real(dp), dimension(4,4) :: Pe, C, Ke
  real(dp) :: detPe, Volume, converr
  real(dp), dimension(12) :: subDets
  integer, parameter :: nr_k_nz = nr_p + (m-1)*n*o*2 + m*(n-1)*o*2 + m*n*(o-1)*2 &
                                  + (m-1)*(n-1)*o*4 + (m-1)*n*(o-1)*4 + m*(n-1)*(o-1)*4
  real(dp), dimension(nr_k_nz) :: vK, vKb
  integer, dimension(nr_k_nz) :: iK
  integer, dimension(nr_p+1) :: jK
  real(dp), dimension(nr_p) :: F, Fb, U
  real(dp), dimension(nr_p*9) :: rwork
  integer, dimension(100) :: iwork
  integer, dimension(5,4,3) :: tetras, tetrasn, tetrasma, tetrasmb

  print *, 'm=', m, ' n=', n, 'o=', o

  bi = 1
  ki = 1
  kj = 1
! Make grid *and* init boundaries *and* init K matrix
  do oi = 0, (o-1)
    do ni = 0, (m-1)
      do mi = 0, (n-1)
        pi = oi*m*n+ni*m+mi+1
        p(pi,:) = (/(dble(mi)/dble(m-1)),(dble(ni)/dble(n-1)),(dble(oi)/dble(o-1))/)
        ! Add node to the boundary list if its a boundary node
        if ((oi .eq. 0) .or. (ni .eq. 0) .or. (mi .eq. 0) .or. (oi .eq. (o-1)) .or. (ni .eq. (n-1)) .or. (mi .eq. (m-1))) then
          b(bi) = pi
          bi = bi + 1
        end if
        ! Start a new column of K, starting with the K(node,node) elt (per SLAP column format)
        ! (self)
        iK(ki) = pi
        jK(kj) = ki
        vK(ki) = 0.0d0
        ki = ki + 1
        kj = kj + 1
        ! The rest of the entries are up to 2 diffls in any direction (we dont need the through-the-cube 3-diffs)
        if (oi .gt. 0) then 
          ! (0,-1,-1)
          if (ni .gt. 0) then
            iK(ki) = (oi-1)*m*n+(ni-1)*m+mi+1
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (-1,0,-1)
          if (mi .gt. 0) then
            iK(ki) = (oi-1)*m*n+ni*m+mi
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (0,0,-1)
          iK(ki) = (oi-1)*m*n+ni*m+mi+1
          vK(ki) = 0.0d0
          ki = ki + 1
          ! (1,0,-1)
          if (mi .lt. (m-1)) then
            iK(ki) = (oi-1)*m*n+ni*m+mi+2
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (0,1,-1)
          if (ni .lt. (n-1)) then
            iK(ki) = (oi-1)*m*n+(ni+1)*m+mi+1
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
        end if
        ! now do all z=0
        if (ni .gt. 0) then
          ! (-1, -1, 0)
          if (mi .gt. 0) then
            iK(ki) = oi*m*n+(ni-1)*m+mi
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (0, -1, 0)
          iK(ki) = oi*m*n+(ni-1)*m+mi+1
          vK(ki) = 0.0d0
          ki = ki + 1
          ! (1, -1, 0)
          if (mi .lt. (m-1)) then
            iK(ki) = oi*m*n+(ni-1)*m+mi+2
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
        end if
        ! (-1, 0, 0)
        if (mi .gt. 0) then
          iK(ki) = oi*m*n+ni*m+mi
          vK(ki) = 0.0d0
          ki = ki + 1
        end if
        ! (1, 0, 0)
        if (mi .lt. (m-1)) then
          iK(ki) = oi*m*n+ni*m+mi+2
          vK(ki) = 0.0d0
          ki = ki + 1
        end if
        if (ni .lt. (n-1)) then
          ! (-1, 1, 0)
          if (mi .gt. 0) then
            iK(ki) = oi*m*n+(ni+1)*m+mi
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (0, 1, 0)
          iK(ki) = oi*m*n+(ni+1)*m+mi+1
          vK(ki) = 0.0d0
          ki = ki + 1
          ! (1, 1, 0)
          if (mi .lt. (m-1)) then
            iK(ki) = oi*m*n+(ni+1)*m+mi+2
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
        end if
        if (oi .lt. (o-1)) then
          ! (0,-1,+1)
          if (ni .gt. 0) then
            iK(ki) = (oi+1)*m*n+(ni-1)*m+mi+1
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (-1, 0,+1)
          if (mi .gt. 0) then
            iK(ki) = (oi+1)*m*n+ni*m+mi
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (0, 0, +1)
          iK(ki) = (oi+1)*m*n+ni*m+mi+1
          vK(ki) = 0.0d0
          ki = ki + 1
          ! (+1, 0,+1)
          if (mi .lt. (m-1)) then
            iK(ki) = (oi+1)*m*n+ni*m+mi+2
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
          ! (0,+1,+1)
          if (ni .lt. (n-1)) then
            iK(ki) = (oi+1)*m*n+(ni+1)*m+mi+1
            vK(ki) = 0.0d0
            ki = ki + 1
          end if
        end if
      end do
    end do
  end do
  jK(nr_p+1)=ki
  print *, nr_b, bi
  print *, nr_k_nz, ki
  actual_nz = ki - 1
!   do e = 1, nr_p
!     print *, p(e,:)
!   end do

  print *, 'done grid'

  ! north-west-top tetrahedron    (0,0,0) (1,0,0) (0,1,0) (0,0,1)
  tetras(1,1,:) = (/ 0, 0, 0 /)
  tetras(1,2,:) = (/ 1, 0, 0 /)
  tetras(1,3,:) = (/ 0, 1, 0 /)
  tetras(1,4,:) = (/ 0, 0, 1 /)
  ! south-east-top tetrahedron    (1,1,0) (1,0,0) (0,1,0) (1,1,1)
  tetras(2,1,:) = (/ 1, 1, 0 /)
  tetras(2,2,:) = (/ 1, 0, 0 /)
  tetras(2,3,:) = (/ 0, 1, 0 /)
  tetras(2,4,:) = (/ 1, 1, 1 /)
  ! north-east-bottom tetrahedron (1,0,1) (1,0,0) (0,0,1) (1,1,1)
  tetras(3,1,:) = (/ 1, 0, 1 /)
  tetras(3,2,:) = (/ 1, 0, 0 /)
  tetras(3,3,:) = (/ 0, 0, 1 /)
  tetras(3,4,:) = (/ 1, 1, 1 /)
  ! south-west-bottom tetrahedron (0,1,1) (0,1,0) (0,0,1) (1,1,1)
  tetras(4,1,:) = (/ 0, 1, 1 /)
  tetras(4,2,:) = (/ 0, 1, 0 /)
  tetras(4,3,:) = (/ 0, 0, 1 /)
  tetras(4,4,:) = (/ 1, 1, 1 /)
  ! internal tetrahedron          (1,0,0) (0,1,0) (0,0,1) (1,1,1)
  tetras(5,1,:) = (/ 1, 0, 0 /)
  tetras(5,2,:) = (/ 0, 1, 0 /)
  tetras(5,3,:) = (/ 0, 0, 1 /)
  tetras(5,4,:) = (/ 1, 1, 1 /)
  

! Make tetrahedron list
  ti = 1
  do oi = 1, o-1
    tetrasn = tetras
    do ni = 1, n-1
      ! flip z every x
      tetrasma = tetrasn
      tetrasmb(:,:,1:2) = tetrasn(:,:,1:2)
      tetrasmb(:,:,3) = mod(tetrasn(:,:,3) + 1, 2)
      do mi = 1, m-1
        if (mod(mi,2) .eq. 0) then
          t(ti,:) = &
            (/(oi-1+tetrasma(1,1,3))*m*n+(ni-1+tetrasma(1,1,2))*m+mi+tetrasma(1,1,1), &
              (oi-1+tetrasma(1,2,3))*m*n+(ni-1+tetrasma(1,2,2))*m+mi+tetrasma(1,2,1), &
              (oi-1+tetrasma(1,3,3))*m*n+(ni-1+tetrasma(1,3,2))*m+mi+tetrasma(1,3,1), &
              (oi-1+tetrasma(1,4,3))*m*n+(ni-1+tetrasma(1,4,2))*m+mi+tetrasma(1,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*1),:) = &
            (/(oi-1+tetrasma(2,1,3))*m*n+(ni-1+tetrasma(2,1,2))*m+mi+tetrasma(2,1,1), &
              (oi-1+tetrasma(2,2,3))*m*n+(ni-1+tetrasma(2,2,2))*m+mi+tetrasma(2,2,1), &
              (oi-1+tetrasma(2,3,3))*m*n+(ni-1+tetrasma(2,3,2))*m+mi+tetrasma(2,3,1), &
              (oi-1+tetrasma(2,4,3))*m*n+(ni-1+tetrasma(2,4,2))*m+mi+tetrasma(2,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*2),:) = &
            (/(oi-1+tetrasma(3,1,3))*m*n+(ni-1+tetrasma(3,1,2))*m+mi+tetrasma(3,1,1), &
              (oi-1+tetrasma(3,2,3))*m*n+(ni-1+tetrasma(3,2,2))*m+mi+tetrasma(3,2,1), &
              (oi-1+tetrasma(3,3,3))*m*n+(ni-1+tetrasma(3,3,2))*m+mi+tetrasma(3,3,1), &
              (oi-1+tetrasma(3,4,3))*m*n+(ni-1+tetrasma(3,4,2))*m+mi+tetrasma(3,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*3),:) = &
            (/(oi-1+tetrasma(4,1,3))*m*n+(ni-1+tetrasma(4,1,2))*m+mi+tetrasma(4,1,1), &
              (oi-1+tetrasma(4,2,3))*m*n+(ni-1+tetrasma(4,2,2))*m+mi+tetrasma(4,2,1), &
              (oi-1+tetrasma(4,3,3))*m*n+(ni-1+tetrasma(4,3,2))*m+mi+tetrasma(4,3,1), &
              (oi-1+tetrasma(4,4,3))*m*n+(ni-1+tetrasma(4,4,2))*m+mi+tetrasma(4,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*4),:) = &
            (/(oi-1+tetrasma(5,1,3))*m*n+(ni-1+tetrasma(5,1,2))*m+mi+tetrasma(5,1,1), &
              (oi-1+tetrasma(5,2,3))*m*n+(ni-1+tetrasma(5,2,2))*m+mi+tetrasma(5,2,1), &
              (oi-1+tetrasma(5,3,3))*m*n+(ni-1+tetrasma(5,3,2))*m+mi+tetrasma(5,3,1), &
              (oi-1+tetrasma(5,4,3))*m*n+(ni-1+tetrasma(5,4,2))*m+mi+tetrasma(5,4,1) /)
          ti = ti + 1
         else
          t(ti,:) = &
            (/(oi-1+tetrasmb(1,1,3))*m*n+(ni-1+tetrasmb(1,1,2))*m+mi+tetrasmb(1,1,1), &
              (oi-1+tetrasmb(1,2,3))*m*n+(ni-1+tetrasmb(1,2,2))*m+mi+tetrasmb(1,2,1), &
              (oi-1+tetrasmb(1,3,3))*m*n+(ni-1+tetrasmb(1,3,2))*m+mi+tetrasmb(1,3,1), &
              (oi-1+tetrasmb(1,4,3))*m*n+(ni-1+tetrasmb(1,4,2))*m+mi+tetrasmb(1,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*1),:) = &
            (/(oi-1+tetrasmb(2,1,3))*m*n+(ni-1+tetrasmb(2,1,2))*m+mi+tetrasmb(2,1,1), &
              (oi-1+tetrasmb(2,2,3))*m*n+(ni-1+tetrasmb(2,2,2))*m+mi+tetrasmb(2,2,1), &
              (oi-1+tetrasmb(2,3,3))*m*n+(ni-1+tetrasmb(2,3,2))*m+mi+tetrasmb(2,3,1), &
              (oi-1+tetrasmb(2,4,3))*m*n+(ni-1+tetrasmb(2,4,2))*m+mi+tetrasmb(2,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*2),:) = &
            (/(oi-1+tetrasmb(3,1,3))*m*n+(ni-1+tetrasmb(3,1,2))*m+mi+tetrasmb(3,1,1), &
              (oi-1+tetrasmb(3,2,3))*m*n+(ni-1+tetrasmb(3,2,2))*m+mi+tetrasmb(3,2,1), &
              (oi-1+tetrasmb(3,3,3))*m*n+(ni-1+tetrasmb(3,3,2))*m+mi+tetrasmb(3,3,1), &
              (oi-1+tetrasmb(3,4,3))*m*n+(ni-1+tetrasmb(3,4,2))*m+mi+tetrasmb(3,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*3),:) = &
            (/(oi-1+tetrasmb(4,1,3))*m*n+(ni-1+tetrasmb(4,1,2))*m+mi+tetrasmb(4,1,1), &
              (oi-1+tetrasmb(4,2,3))*m*n+(ni-1+tetrasmb(4,2,2))*m+mi+tetrasmb(4,2,1), &
              (oi-1+tetrasmb(4,3,3))*m*n+(ni-1+tetrasmb(4,3,2))*m+mi+tetrasmb(4,3,1), &
              (oi-1+tetrasmb(4,4,3))*m*n+(ni-1+tetrasmb(4,4,2))*m+mi+tetrasmb(4,4,1) /)
          t(ti+((m-1)*(n-1)*(o-1)*4),:) = &
            (/(oi-1+tetrasmb(5,1,3))*m*n+(ni-1+tetrasmb(5,1,2))*m+mi+tetrasmb(5,1,1), &
              (oi-1+tetrasmb(5,2,3))*m*n+(ni-1+tetrasmb(5,2,2))*m+mi+tetrasmb(5,2,1), &
              (oi-1+tetrasmb(5,3,3))*m*n+(ni-1+tetrasmb(5,3,2))*m+mi+tetrasmb(5,3,1), &
              (oi-1+tetrasmb(5,4,3))*m*n+(ni-1+tetrasmb(5,4,2))*m+mi+tetrasmb(5,4,1) /)
          ti = ti + 1
         end if
      end do
      ! flip x every y
      tetrasn(:,:,1) = mod(tetrasn(:,:,1) + 1, 2)
    end do
    ! flip y every z
    tetras(:,:,2) = mod(tetras(:,:,2) + 1, 2)
  end do
  print *, nr_t, ti

  print *, 'done tets'

!   do ti = 1, nr_t
!     print *, t(ti,:)
!   end do

  do e = 1, nr_t
    
    Pe(1,:) = (/ 1.0d0, p(t(e,1),:) /)
    Pe(2,:) = (/ 1.0d0, p(t(e,2),:) /)
    Pe(3,:) = (/ 1.0d0, p(t(e,3),:) /)
    Pe(4,:) = (/ 1.0d0, p(t(e,4),:) /)
    
    ! Invert using Laplace Expansion thm
    subDets = (/ Pe(1,1)*Pe(2,2)-Pe(1,2)*Pe(2,1), &
                 Pe(1,1)*Pe(2,3)-Pe(2,1)*Pe(1,3), &
                 Pe(1,1)*Pe(2,4)-Pe(2,1)*Pe(1,4), &
                 Pe(1,2)*Pe(2,3)-Pe(2,2)*Pe(1,3), &
                 Pe(1,2)*Pe(2,4)-Pe(2,2)*Pe(1,4), &
                 Pe(1,3)*Pe(2,4)-Pe(2,3)*Pe(1,4), &
                 Pe(3,1)*Pe(4,2)-Pe(4,1)*Pe(3,2), &
                 Pe(3,1)*Pe(4,3)-Pe(4,1)*Pe(3,3), &
                 Pe(3,1)*Pe(4,4)-Pe(4,1)*Pe(3,4), &
                 Pe(3,2)*Pe(4,3)-Pe(4,1)*Pe(3,3), &
                 Pe(3,2)*Pe(4,4)-Pe(3,4)*Pe(4,2), &
                 Pe(3,3)*Pe(4,4)-Pe(3,4)*Pe(4,3) /)
    
    detPe = subDets(1) * subDets(12) &
          - subDets(2) * subDets(11) &
          + subDets(3) * subDets(10) &
          + subDets(4) * subDets(9) &
          - subDets(5) * subDets(8) &
          + subDets(6) * subDets(7)
         
         
    Volume = dabs(detPe)/6.0d0

    C(1,:) = (/  Pe(2,2)*subDets(12)-Pe(2,3)*subDets(11)+Pe(2,4)*subDets(10), &
                -Pe(1,2)*subDets(12)+Pe(1,3)*subDets(11)-Pe(1,4)*subDets(10), &
                 Pe(4,2)*subDets(6) -Pe(4,3)*subDets(5) +Pe(4,4)*subDets(4) , &
                -Pe(3,2)*subDets(6) +Pe(3,3)*subDets(5) -Pe(3,4)*subDets(4) /)
    C(2,:) = (/ -Pe(2,1)*subDets(12)+Pe(2,3)*subDets(9) -Pe(2,4)*subDets(8) , &
                 Pe(1,1)*subDets(12)-Pe(1,3)*subDets(9) +Pe(1,4)*subDets(8) , &
                -Pe(4,1)*subDets(6) +Pe(4,3)*subDets(3) -Pe(4,4)*subDets(2) , &
                 Pe(3,1)*subDets(6) -Pe(3,3)*subDets(3) +Pe(3,4)*subDets(2) /)
    C(3,:) = (/  Pe(2,1)*subDets(11)-Pe(2,2)*subDets(9) +Pe(2,4)*subDets(7) , &
                -Pe(1,1)*subDets(11)+Pe(1,2)*subDets(9) -Pe(1,4)*subDets(7) , &
                 Pe(4,1)*subDets(5) -Pe(4,2)*subDets(3) +Pe(4,4)*subDets(1) , &
                -Pe(3,1)*subDets(5) +Pe(3,2)*subDets(3) -Pe(3,4)*subDets(1) /)
    C(4,:) = (/ -Pe(2,1)*subDets(10)+Pe(2,2)*subDets(8) -Pe(2,3)*subDets(7) , &
                 Pe(1,1)*subDets(10)-Pe(1,2)*subDets(8) +Pe(1,3)*subDets(7) , &
                -Pe(4,1)*subDets(4) +Pe(4,2)*subDets(2) -Pe(4,3)*subDets(1) , &
                 Pe(3,1)*subDets(4) -Pe(3,2)*subDets(2) +Pe(3,3)*subDets(1) /)
    
    C = C / detPe
    
    Ke = Volume * matmul(transpose(C(2:4,:)), C(2:4,:))
    
    do ti = 1, 4
      do tj = 1, 4
        ki = k_index(t(e,ti),t(e,tj))
        vK(ki) = vK(ki) + Ke(ti,tj)
      end do
      F(t(e,ti)) = F(t(e,ti)) + (Volume / 4)
    end do

    
    Fe = Volume / 4 ! is this right
 
  end do
    print *, Volume
    print *, Pe(1,:)
    print *, Pe(2,:)
    print *, Pe(3,:)
    print *, Pe(4,:)
    print *, Ke(1,:)
    print *, Ke(2,:)
    print *, Ke(3,:)
    print *, Ke(4,:)
     
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
  Fb = F
  vKb = vK
  U = 1.0d0
  CALL DSDBCG(nr_p, Fb, U, nr_k_nz, iK, jK, vKb, 0, 1, 1.0d-8, 10000, iter, converr, ierr, 6, rwork, nr_p * 9, iwork, 100)
  print *, iter, converr, ierr

  do ni = 1, n
    print *, U((o/2)*m*n+(ni-1)*m:(o/2)*m*n+ni*m-1)
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