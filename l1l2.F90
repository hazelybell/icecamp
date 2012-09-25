#define dp kind(1.0d0)
! #define LAGRANGEBC 1
! #define NEWTONS 1

!compute row/col index from variable interleave
#define VPtoI(v, p) ((p-1)*nr_v+v)
!compute variable interleave from row/col index
#define ItoV(i) (MOD(i-1,nr_v)+1)
#define ItoP(i) (((i-1)/nr_v)+1)

program fem
  call ismip_hom_c
  
    
  contains
    subroutine eismint_square_bay
        implicit none
        integer, parameter :: n = 41, m = 41  ! Grid size n = rows, m = cols
        real(dp), parameter :: domain_ns = 2d5 !m
        real(dp), parameter :: domain_ew = 2d5 !m
        external eismint_square_bay2
        call rect_grid_to_triangles(n, m, domain_ns, domain_ew, eismint_square_bay2, 0)
    end subroutine
    subroutine ismip_hom_c
        implicit none
        integer, parameter :: n = 40, m = 40  ! Grid size n = rows, m = cols
        real(dp), parameter :: ismip_size = 160
        real(dp), parameter :: domain_ns = ismip_size * (DFLOAT(m)/DFLOAT(m+1)) !m
        real(dp), parameter :: domain_ew = ismip_size * (DFLOAT(n)/DFLOAT(n+1)) !n
        external ismip_hom_c2
        call rect_grid_to_triangles(m, n, domain_ns, domain_ew, ismip_hom_c2, 1, ismip_size)
    end subroutine

end program
subroutine rect_grid_to_triangles(m, n, domain_ns, domain_ew, sub, periodic, ismip_size)
    implicit none
    integer, intent(in) :: m, n, periodic ! periodic must be 1 or 0 anything else will surely explode
    real(dp), intent(in) :: domain_ns !m
    real(dp), intent(in) :: domain_ew !n
    real(dp), intent(in) :: ismip_size 
    integer :: nr_p
    integer :: nr_t
    integer :: max_neighbors = 9
    real(dp), dimension(m * n,2) :: p ! List of points
    integer, dimension((m-1+periodic) * (n-1+periodic) * 2,3) :: t ! List of triangles, boundary flags
    integer :: mi, ni, ti, seed
!     seed = 100
    nr_p = m * n
    nr_t = (m-1+periodic) * (n-1+periodic) * 2
    print *, 'm=', m, ' n=', n
    ! Make grid
    do ni = 0, (n-1) ! rows
      do mi = 0, (m-1) ! cols
        p(mi*(m)+ni+1,:) = (/(DBLE(ni)/DBLE(n-1))*domain_ew,(DBLE(mi)/DBLE(m-1))*domain_ns/)
      end do
    end do
    ! Make triangle list
    ti = 0
    do ni = 0, n-2+periodic ! rows
      do mi = 0, m-2+periodic ! cols
        if (RAN(seed) .gt. .5) then
          ti = ti + 1
          t(ti,1) = ni*(m)+mi+1
          t(ti,2) = MOD(ni+1,n)*(m)+mi+1
          t(ti,3) = ni*(m)+MOD(mi+1,m)+1
          ti = ti + 1
          t(ti,1) = MOD(ni+1,n)*(m)+MOD(mi+1,m)+1
          t(ti,2) = MOD(ni+1,n)*(m)+mi+1
          t(ti,3) = ni*(m)+MOD(mi+1,m)+1
        else
          ti = ti + 1
          t(ti,1) = ni*(m)+mi+1
          t(ti,2) = MOD(ni+1,n)*(m)+mi+1
          t(ti,3) = MOD(ni+1,n)*(m)+MOD(mi+1,m)+1
          ti = ti + 1
          t(ti,1) = ni*(m)+mi+1
          t(ti,2) = ni*(m)+MOD(mi+1,m)+1
          t(ti,3) = MOD(ni+1,n)*(m)+MOD(mi+1,m)+1
        end if
      end do
    end do

    call sub(m, n, nr_p, p, nr_t, t, max_neighbors, ismip_size)
end subroutine


! ----------------------------------------------------------------------------
! THIS IS WHERE THE MAGIC HAPPENS
! ----------------------------------------------------------------------------

subroutine fem_imr(nr_p, p, nr_t, t, bed, thickness, basal_traction_coeff, grounded, nr_b, max_neighbors, nr_calving_front, calving_front, nr_dirichlet, dirichlet, dirichlet_val, ds_dx_fudge, ds_dy_fudge, flow_law_A, x_wrap, y_wrap, x_wrap_extra, y_wrap_extra, sliding_law)
    use m_mrgrnk, only: mrgrnk
    implicit none
!     Interface
!       Subroutine D_MRGRNK(XVALT, IRNGT)
!           Real(dp), Dimension (:), Intent (In) :: XVALT
!           Integer, Dimension (:), Intent (Out) :: IRNGT
!       End Subroutine D_MRGRNK
!     End Interface
    integer,parameter :: max_p_growth = 2
    integer,parameter :: max_t_growth = 3
    integer :: old_to_new_ratio  ! subdivide one out of every
    integer,parameter :: nr_v = 4 ! Number of variables (U,V,dU/dz,dV/dz)

    integer, intent(in) :: nr_p ! Number of points
    real(dp), dimension(nr_p,2), intent(in) :: p ! List of points
    integer, intent(in) :: nr_t ! Number of triangles
    integer, dimension(nr_t,3), intent(in) :: t ! List of triangles, boundary flags
    integer, intent(in) :: nr_b
    integer, intent(inout) :: max_neighbors
    integer, intent(in) :: nr_calving_front
    integer, dimension(nr_calving_front), intent(in) :: calving_front
    integer, intent(in) :: nr_dirichlet
    integer, dimension(nr_dirichlet), intent(in) :: dirichlet
    real(dp), dimension(nr_v, nr_dirichlet), intent(in) :: dirichlet_val
    real(dp), intent(in) :: ds_dx_fudge, ds_dy_fudge ! derivatives of the surface elev., m
    real(dp), intent(in) :: flow_law_A
    real(dp), intent(in) :: x_wrap, y_wrap, x_wrap_extra, y_wrap_extra
    real(dp), dimension(nr_p), intent(in) :: bed ! List of points, m
    real(dp), dimension(nr_p), intent(in) :: thickness! List of points, m
    real(dp), dimension(nr_p), intent(in) :: basal_traction_coeff ! m/(y*Pa^2)
    real(dp), dimension(nr_p), intent(in) :: grounded ! grounded flag 1=grounded 0=floating
    integer, intent(in) :: sliding_law

    real(dp), dimension(nr_p, nr_v) :: guess
    real(dp), dimension(nr_p, nr_v) :: error_estimate, resid_estimate
    real(dp) :: nl_resid
    real(dp) :: maxy,maxx,miny,minx
    integer :: use_guess
    integer :: max_nliter, max_nliter_imr
    integer :: vi, pi, ti, imri, ai, pa, pb, pc
    integer :: uv_ok = 0
    integer :: max_imriter = 0
    integer :: nr_p_imr, nr_t_imr, max_p, max_t
    integer, dimension(nr_t*max_t_growth) :: t_rank
    real(dp), dimension(nr_t*max_t_growth) :: t_err, t_area
    integer, dimension(3) :: t_old
    real(dp), dimension(2) :: wraps
    real(dp), dimension(2) :: wrap_extras

    real(dp), dimension(nr_p*max_p_growth,2) :: p_imr
    integer, dimension(nr_t*max_t_growth,3) :: t_imr
    real(dp), dimension(nr_p*max_p_growth) :: bed_imr
    real(dp), dimension(nr_p*max_p_growth) :: thickness_imr
    real(dp), dimension(nr_p*max_p_growth) :: basal_traction_coeff_imr
    real(dp), dimension(nr_p*max_p_growth) :: grounded_imr
    real(dp), dimension(nr_p*max_p_growth, nr_v) :: guess_imr
    real(dp), dimension(nr_p*max_p_growth, nr_v) :: error_estimate_imr, resid_estimate_imr

    p_imr = 0d0 ! List of points
    t_imr = 0 ! List of triangles, boundary flags
    bed_imr = 0d0 ! List of points, m
    thickness_imr = 0d0 ! List of points, m
    basal_traction_coeff_imr = 0d0 ! m/(y*Pa^2)
    grounded_imr = 0d0 ! grounded flag 1=grounded 0=floating
    guess_imr = 0d0
    error_estimate_imr = 0d0
    resid_estimate_imr = 0d0
#ifdef NEWTONS
    max_nliter = (INT(SQRT(DBLE(nr_p))))+2
#else
    max_nliter = (INT(SQRT(DBLE(nr_p))))+20
#endif
    guess = 0d0
    use_guess = 0
    max_p = nr_p * max_p_growth
    max_t = nr_t * max_t_growth
    wraps = (/x_wrap, y_wrap/)
    wrap_extras = (/x_wrap_extra, y_wrap_extra/)
    old_to_new_ratio = (INT(SQRT(DBLE(nr_p))))
    
    maxx = MAXVAL(p(:,1))
    minx = MINVAL(p(:,1))
    maxy = MAXVAL(p(:,2))
    miny = MINVAL(p(:,2))
    
    nr_p_imr = nr_p
    nr_t_imr = nr_t

    call fem_l1l2(nr_p, p, nr_t, t, bed, thickness, basal_traction_coeff, grounded, nr_b, max_neighbors, nr_calving_front, calving_front, nr_dirichlet, dirichlet, dirichlet_val, ds_dx_fudge, ds_dy_fudge, flow_law_A, x_wrap, y_wrap, x_wrap_extra, y_wrap_extra, max_nliter, guess, 0, error_estimate, resid_estimate, nl_resid, t_area(1:nr_t_imr), sliding_law)
    
    nr_p_imr = nr_p
    nr_t_imr = nr_t
    p_imr(1:nr_p,:) = p
    t_imr(1:nr_t,:) = t
    bed_imr(1:nr_p) = bed
    thickness_imr(1:nr_p) = thickness
    basal_traction_coeff_imr(1:nr_p) = basal_traction_coeff
    grounded_imr(1:nr_p) = grounded
    guess_imr(1:nr_p,:) = guess
    error_estimate_imr(1:nr_p,:) = error_estimate
    resid_estimate_imr(1:nr_p,:) = resid_estimate
    
    do imri = 1, max_imriter
      use_guess = use_guess + 1
      do ti = 1, nr_t_imr
        t_err(ti) = DSQRT(SUM(resid_estimate_imr(t_imr(ti,:),:)**2d0))*t_area(ti)**2d0
!         t_err(ti) = MAXVAL(DABS(resid_estimate_imr(t_imr(ti,:),:)))*t_area(ti)
      end do
      call MRGRNK(t_err(1:nr_t_imr), t_rank(1:nr_t_imr))
      a : do ti = 1, nr_t_imr
        if (t_rank(ti) .gt. nr_t_imr-(nr_t_imr/old_to_new_ratio)) then
          if (nr_p_imr .le. max_p - 1 .and. nr_t_imr .le. max_t - 2) then
            ! generate a new vertex at the centroid
            t_old = t_imr(ti,:)
            do vi = 1, 2
              if (wraps(vi) .gt. 0d0 .and. MAXVAL(p_imr(t_imr(ti,:),vi)) .gt. x_wrap * .75d0) then
                  do pi = 1, 3
                    if (p_imr(t_imr(ti,pi),vi) .lt. x_wrap * .25d0) then
                      p_imr(nr_p_imr+1,vi) = p_imr(nr_p_imr+1,vi) + p_imr(t_imr(ti,pi),vi) + wraps(vi) + wrap_extras(vi)
                    else
                      p_imr(nr_p_imr+1,vi) = p_imr(nr_p_imr+1,vi) + p_imr(t_imr(ti,pi),vi)
                    end if
                  end do
                  p_imr(nr_p_imr+1,vi) = p_imr(nr_p_imr+1,vi)/3d0
              else
                p_imr(nr_p_imr+1,vi) = (p_imr(t_imr(ti,1),vi)+p_imr(t_imr(ti,2),vi)+p_imr(t_imr(ti,3),vi))/3d0
              end if
            end do
            do vi = 1, nr_v
              guess_imr(nr_p_imr+1,vi) = (guess_imr(t_imr(ti,1),vi)+guess_imr(t_imr(ti,2),vi)+guess_imr(t_imr(ti,3),vi))/3d0
              error_estimate_imr(nr_p_imr+1,vi) = (error_estimate_imr(t_imr(ti,1),vi)+error_estimate_imr(t_imr(ti,2),vi)+error_estimate_imr(t_imr(ti,3),vi))/3d0
              resid_estimate_imr(nr_p_imr+1,vi) = (resid_estimate_imr(t_imr(ti,1),vi)+resid_estimate_imr(t_imr(ti,2),vi)+resid_estimate_imr(t_imr(ti,3),vi))/3d0
            end do
            bed_imr(nr_p_imr+1) = (bed_imr(t_imr(ti,1))+bed_imr(t_imr(ti,2))+bed_imr(t_imr(ti,3)))/3d0
            thickness_imr(nr_p_imr+1) = (thickness_imr(t_imr(ti,1))+thickness_imr(t_imr(ti,2))+thickness_imr(t_imr(ti,3)))/3d0
            basal_traction_coeff_imr(nr_p_imr+1) = (basal_traction_coeff_imr(t_imr(ti,1))+basal_traction_coeff_imr(t_imr(ti,2))+basal_traction_coeff_imr(t_imr(ti,3)))/3d0
            grounded_imr(nr_p_imr+1) = (grounded_imr(t_imr(ti,1))+grounded_imr(t_imr(ti,2))+grounded_imr(t_imr(ti,3)))/3d0
            ! generate 2 new triangles and replace 1
            t_imr(ti,:) = (/t_old(1), t_old(2), nr_p_imr+1/)
            t_imr(nr_t_imr+1,:) = (/t_old(2), t_old(3), nr_p_imr+1/)
            t_imr(nr_t_imr+2,:) = (/t_old(3), t_old(1), nr_p_imr+1/)
            ! update counts
            nr_p_imr = nr_p_imr + 1
            nr_t_imr = nr_t_imr + 2
          else
            exit a
          end if
        end if
      end do a
      max_nliter_imr = max_nliter/4+2
      max_neighbors = max_neighbors * 2
      call fem_l1l2(nr_p_imr, p_imr(1:nr_p_imr,:), nr_t_imr, t_imr(1:nr_t_imr,:), bed_imr(1:nr_p_imr), thickness_imr(1:nr_p_imr), basal_traction_coeff_imr(1:nr_p_imr), grounded_imr(1:nr_p_imr), nr_b, max_neighbors, nr_calving_front, calving_front, nr_dirichlet, dirichlet, dirichlet_val, ds_dx_fudge, ds_dy_fudge, flow_law_A, x_wrap, y_wrap, x_wrap_extra, y_wrap_extra, max_nliter_imr, guess_imr(1:nr_p_imr,:), 1, error_estimate_imr(1:nr_p_imr,:), resid_estimate_imr(1:nr_p_imr,:), nl_resid, t_area(1:nr_t_imr), sliding_law)
      print *, "Max Neighbors", max_neighbors
    end do
    
    OPEN(UNIT=2, STATUS='REPLACE', ACTION='WRITE', FILE='plot.txt')
    
    WRITE (2, *) "Triangles"
    do ti = 1, nr_t_imr
        if (.not.((MAXVAL(p_imr(t_imr(ti,:),1)) - minx .gt. (maxx - minx) * 0.75 .and. MINVAL(p_imr(t_imr(ti,:),1)) - minx .lt. (maxx - minx) * 0.25).or.(MAXVAL(p_imr(t_imr(ti,:),2)) - miny .gt. (maxy - miny) * 0.75 .and. MINVAL(p_imr(t_imr(ti,:),2)) - miny .lt. (maxy - miny) * 0.25))) then
          WRITE (2, *) t_imr(ti,:)-1
        end if
    end do
    print *, "End Triangles"
    
    do vi = 1, nr_v
      if ((MAXVAL(guess_imr(1:nr_p_imr,vi))-MINVAL(guess_imr(1:nr_p_imr,vi))) .gt. 0d0) then
        WRITE (2, *) "Variable ", vi
        do pi = 1, nr_p_imr
          WRITE (2, *) p_imr(pi,1:2), guess_imr(pi, vi)
        end do
        WRITE (2, *) "End Variable"
        uv_ok = uv_ok + 1
      else
        WRITE (2, *) "Flat Variable ", vi, " Range ", MINVAL(guess_imr(1:nr_p_imr,vi)), MAXVAL(guess_imr(1:nr_p_imr,vi))
      end if
    end do
    
    if (uv_ok .gt. 0) then
      WRITE (2, *) "Variable ", 9
      do pi = 1, nr_p_imr
        WRITE (2, *) p_imr(pi,1:2), sqrt(guess_imr(pi,1)**2d0+guess_imr(pi,2)**2d0)
      end do
      WRITE (2, *) "End Variable"
    end if
    
    uv_ok = 0
    do vi = 1, nr_v
      if ((MAXVAL(error_estimate_imr(1:nr_p_imr,vi))-MINVAL(error_estimate_imr(1:nr_p_imr,vi))) .gt. 0d0) then
        WRITE (2, *) "Variable ", vi+10
        do pi = 1, nr_p_imr
          WRITE (2, *) p_imr(pi,1:2), error_estimate_imr(pi, vi)
        end do
        WRITE (2, *) "End Variable"
        uv_ok = uv_ok + 1
      else
        WRITE (2, *) "Flat Variable ", vi+10, " Range ", MINVAL(error_estimate_imr(1:nr_p_imr,vi)), MAXVAL(error_estimate_imr(1:nr_p_imr,vi))
      end if
    end do

    if (uv_ok .gt. 0) then
      WRITE (2, *) "Variable ", 19
      do pi = 1, nr_p_imr
        WRITE (2, *) p_imr(pi,1:2), sqrt(error_estimate_imr(pi,1)**2d0+error_estimate_imr(pi,2)**2d0)
      end do
      WRITE (2, *) "End Variable"
    end if

    uv_ok = 0
    do vi = 1, nr_v
      if ((MAXVAL(resid_estimate_imr(1:nr_p_imr,vi))-MINVAL(resid_estimate_imr(1:nr_p_imr,vi))) .gt. 0d0) then
        WRITE (2, *) "Variable ", vi+20
        do pi = 1, nr_p_imr
          WRITE (2, *) p_imr(pi,1:2), resid_estimate_imr(pi, vi)
        end do
        WRITE (2, *) "End Variable"
        uv_ok = uv_ok + 1
      else
        WRITE (2, *) "Flat Variable ", vi+20, " Range ", MINVAL(resid_estimate_imr(1:nr_p_imr,vi)), MAXVAL(resid_estimate_imr(1:nr_p_imr,vi))
      end if
    end do

    if (uv_ok .gt. 0) then
      WRITE (2, *) "Variable ", 29
      do pi = 1, nr_p_imr
        WRITE (2, *) p_imr(pi,1:2), sqrt(resid_estimate_imr(pi,1)**2d0+resid_estimate_imr(pi,2)**2d0)
      end do
      WRITE (2, *) "End Variable"
    end if

    WRITE (2, *) "Done"


end subroutine

subroutine fem_l1l2(nr_p, p, nr_t, t, bed, thickness, basal_traction_coeff, grounded, nr_b, max_neighbors, nr_calving_front, calving_front, nr_dirichlet, dirichlet, dirichlet_val, ds_dx_fudge, ds_dy_fudge, flow_law_A, x_wrap, y_wrap, x_wrap_extra, y_wrap_extra, max_nliter, guess, use_guess, error_estimate, resid_estimate, nl_resid, t_area, sliding_law)
    implicit none
    integer,parameter :: nr_v = 4 ! Number of variables (U,V,dU/dz,dV/dz)
    integer,parameter :: nr_v2d = 2 ! Number of 2D meshin' variables
    integer,parameter :: nr_v0d = 2 ! Number of 0D (point) variables
    !   integer,parameter :: dp = kind(1.0d0)
    integer :: mi, ni, ti, tj, bi, bj, e, pi, ki, kj, iter, ierr, nliter
    integer, intent(in) :: nr_p ! Number of points
    integer, dimension(nr_p) :: all_p
    real(dp), dimension(nr_p,2), intent(in) :: p ! List of points
    integer, intent(in) :: nr_t ! Number of triangles
    integer, dimension(nr_t,3), intent(in) :: t ! List of triangles, boundary flags
    integer, intent(in) :: nr_b
    integer, dimension(nr_p) :: is_boundary_point
    real(dp), dimension(nr_p) :: boundary_point_nx, boundary_point_ny
    integer, intent(inout) :: max_neighbors
    integer, dimension(nr_p) :: nr_neighbors
    integer, dimension(nr_p, max_neighbors) :: neighbors ! adjacent vertices to a vertex
    integer, intent(in) :: nr_calving_front
    integer, dimension(nr_calving_front), intent(in) :: calving_front
    integer, intent(in) :: nr_dirichlet
    integer, dimension(nr_dirichlet), intent(in) :: dirichlet
    real(dp), dimension(nr_v, nr_dirichlet), intent(in) :: dirichlet_val
    real(dp), dimension(3,3) :: Pe, C, Ke, Je
    real(dp), dimension(3) :: element_thickness, element_vt, element_basal_coeff, element_basal_coeff_diff, element_u1vels, element_u2vels, basal_coeff_final
    real(dp) :: detPe, Area, invDet, converr
#ifdef NEWTONS
    integer, parameter :: nr_bmp_nz = 4 ! Number of block matrix pattern elements that are non-zero
    integer, parameter :: nr_bmp_nz = 2 ! Number of block matrix pattern elements that are non-zero
#else
    integer, parameter :: nr_bmp_nz2d = 4 ! Number of block matrix pattern elements that are non-zero
    integer, parameter :: nr_bmp_nz0d = 2 ! Number of block matrix pattern elements that are non-zero
#endif
#define LINEAR_ORDER (nr_p * nr_v)
#define NR_K_NZ ((nr_p + nr_t * 3 + nr_b) * nr_bmp_nz2d + (nr_p) * nr_bmp_nz0d)
    integer :: linear_order
    integer :: nr_k_nz ! Number of entries in the K matrix that are non-zero
    real(dp), dimension(NR_K_NZ) :: vKa, vKb ! K is normal matrix, J is the jacobian
    integer, dimension(NR_K_NZ) :: iK
    integer, dimension(LINEAR_ORDER+1) :: jK
#ifdef NEWTONS
    real(dp), dimension(NR_K_NZ) :: vJa, vJb ! K is normal matrix, J is the jacobian
    integer, dimension(NR_K_NZ) :: iJ
    integer, dimension(LINEAR_ORDER+1) :: jJ
#endif ! NEWTONS
    real(dp), dimension(LINEAR_ORDER) :: F, Fb, U, U_prev, delta_U, R, Rb
    real(dp), dimension(131 + LINEAR_ORDER * 18) :: rwork
    integer, dimension(100) :: iwork
    real(dp), parameter :: ice_density = 917d0 ! kg m-3
    real(dp), parameter :: accel_grav = 9.81d0 ! m s-2
    real(dp), parameter :: approx_ice_visc = 25d6 ! Pa a
    real(dp) :: ice_visc, du_dx, dv_dx, du_dy, dv_dy, du_dz, dv_dz
    real(dp), intent(out) :: nl_resid
    real(dp) :: ds_dx, ds_dy, rghdsdx, rghdsdy ! derivatives of the surface elev., m
    real(dp) :: pvecs(5,2) ! temporary vector storage
    real(dp), intent(in) :: ds_dx_fudge, ds_dy_fudge ! derivatives of the surface elev., m
    real(dp), intent(in) :: x_wrap, y_wrap, x_wrap_extra, y_wrap_extra
    real(dp) :: min_edge_length
    real(dp) :: A_factor, B_factor, C_factor, ice_visc_zeta
    real(dp), parameter :: seawater_density = 1028d0 ! kg m-3
#ifdef NEWTONS
    integer, dimension(nr_v,nr_v), parameter :: bmp = reshape((/1, 21, 0, 0,  12, 2, 0, 0,  0, 0, 3, 0,  0, 0, 0, 4/), shape(bmp)) ! Block Matrix Pattern
#else
    integer, dimension(nr_v,nr_v), parameter :: bmp = reshape((/1, 21, 0, 0,  12, 2, 0, 0,  0, 0, 3, 0,  0, 0, 0, 4/), shape(bmp)) ! Block Matrix Pattern
#endif
  ! This must be multiplied by the number of non-zero elements in the BMP
    integer :: vi, vj ! For loops over variables
    real(dp), parameter :: flow_law_exponent = 3d0
    real(dp), intent(in) :: flow_law_A
    real(dp) :: flow_law_B
    real(dp) :: linear_tol
    real(dp), parameter :: initial_tol = 1d-4
    real(dp), parameter :: final_tol = 2d-12
#ifdef NEWTONS
    real(dp), parameter :: mode_tol = 1d-9 ! m/a
#else
    real(dp), parameter :: mode_tol = final_tol*1d0 ! all updates w/in 2x10^-11 m/a
#endif
#ifdef LAGRANGEBC
    real(dp), parameter :: lagrange_multiplier = 1d12
#endif
    real(dp), parameter :: flow_law_epsilon = 1e-7
    real(dp), parameter :: basal_sliding_epsilon = 1e-2 ! From Beuller
    real(dp), dimension(nr_p), intent(in) :: bed ! List of points, m
    real(dp), dimension(nr_p), intent(in) :: thickness! List of points, m
    real(dp), dimension(nr_p), intent(in) :: basal_traction_coeff ! units vary: Pa*a/m ; m/(y*Pa^2)
    real(dp), parameter :: basal_sliding_exp_param = 2d0 ! "m" parameter
    real(dp), parameter :: basal_sliding_exp = ((1d0-basal_sliding_exp_param)/(2d0*basal_sliding_exp_param)) ! (1-m)/(2m)
    real(dp), parameter :: basal_sliding_exp_diff = basal_sliding_exp - 1d0 ! (1-3m)/(2m)
    real(dp), dimension(nr_p), intent(in) :: grounded ! grounded flag 1=grounded 0=floating
    real(dp), dimension(nr_p) :: basal_coeff, basal_coeff_diff, vel_sqr
    real(dp), parameter :: approx_vel_sqr = basal_sliding_epsilon**2d0 ! (m/a)^2
    integer, intent(in) :: max_nliter, sliding_law
    real(dp), intent(inout), dimension(nr_p, nr_v) :: guess
    real(dp), intent(inout), dimension(nr_p, nr_v) :: error_estimate, resid_estimate
    integer, intent(in) :: use_guess
    real(dp), dimension(nr_t), intent(out) :: t_area
    integer :: start_nliter
    real(dp), dimension(nr_p,2,2,2) :: p_ellipticals
    real(dp), dimension(nr_p,2,2) :: p_ds
    real(dp), dimension(nr_p,2) :: p_visc
    real(dp), dimension(nr_p,2) :: p_2dstrain
    real(dp), dimension(5,2) :: pvec
    real(dp) :: e2dstrain, e3dstrain
    real(dp), dimension(4) :: p_duv_dz
#define SURFACE(p) (bed(p) + thickness(p))


    ! ----- Triangles are set up by now, rest of code is geometry agnostic -----
    flow_law_B = flow_law_A**((1d0-flow_law_exponent)/(2d0*flow_law_exponent))
    nr_k_nz = NR_K_NZ
    linear_order = LINEAR_ORDER
    t_area = -1d0
    all_p = (/(pi, pi=1, nr_p)/)
    call analyze_mesh(max_neighbors, nr_t, t, nr_p, p, nr_neighbors, neighbors, is_boundary_point, boundary_point_nx, boundary_point_ny, min_edge_length)
    call init_k_sparsity(max_neighbors, nr_k_nz, linear_order, iK, jK, vKa, nr_v, nr_t, t, bmp, nr_p, nr_neighbors, neighbors)
#ifdef NEWTONS
    iJ=iK
    jJ=jK
#endif ! NEWTONS
    
    
    if (use_guess .gt. 0) then
      do vi = 1, nr_v
        U(vi::nr_v) = guess(:,vi)
        delta_U(vi::nr_v) = error_estimate(:,vi)
        R(vi::nr_v) = resid_estimate(:,vi)
      end do
      start_nliter = 2
    else
      U = 1d0
      delta_U = 1d0
      R = 1d0
      start_nliter = 1
    end if

    do nliter = start_nliter, max_nliter
      U_prev = U
      F = 0d0
#ifdef NEWTONS
      vJa = 0d0
#endif ! NEWTONS
      vKa = 0d0
      ! Pre-calculate various things for each point.
      do pi = 1, nr_p
        if (nliter .eq. 1) then
            vel_sqr(pi) = approx_vel_sqr
        else
            vel_sqr(pi) = U(VPtoI(1, pi))**2d0 + U(VPtoI(2, pi))**2d0 + basal_sliding_epsilon**2d0
        end if
        if (sliding_law .eq. 1) then
          basal_coeff(pi) = grounded(pi)*(1d0/(basal_traction_coeff(pi)**(1d0/basal_sliding_exp_param)))*(vel_sqr(pi)**basal_sliding_exp)
          basal_coeff_diff(pi) = grounded(pi)*(1d0/(basal_traction_coeff(pi)**(1d0/basal_sliding_exp_param)))*2d0*basal_sliding_exp*(vel_sqr(pi)**basal_sliding_exp_diff)
        else if (sliding_law .eq. 2) then
          basal_coeff(pi) = basal_traction_coeff(pi) * grounded(pi)
          basal_coeff_diff(pi) = basal_traction_coeff(pi) * grounded(pi)
        end if
      end do
      p_ellipticals = 0d0
      p_ds = 0d0
      p_visc = 0d0
      p_2dstrain = 0d0
      pvec = 0d0
      do e = 1, nr_t
        ! Calculate various things for each element.
        Pe(1,:) = (/ 1.0d0, p(t(e,1),:) /)
        Pe(2,:) = (/ 1.0d0, p(t(e,2),:) /)
        Pe(3,:) = (/ 1.0d0, p(t(e,3),:) /)
        if (x_wrap .gt. 0d0) then
          if (MAXVAL(Pe(:,2)) .gt. x_wrap * .75d0) then
            if (Pe(1,2) .lt. x_wrap * .25d0) then
              Pe(1,2) = Pe(1,2) + x_wrap + x_wrap_extra
            end if
            if (Pe(2,2) .lt. x_wrap * .25d0) then
              Pe(2,2) = Pe(2,2) + x_wrap + x_wrap_extra
            end if
            if (Pe(3,2) .lt. x_wrap * .25d0) then
              Pe(3,2) = Pe(3,2) + x_wrap + x_wrap_extra
            end if
          end if
        end if
        if (y_wrap .gt. 0d0) then
          if (MAXVAL(Pe(:,3)) .gt. y_wrap * .75d0) then
            if (Pe(1,3) .lt. y_wrap * .25d0) then
              Pe(1,3) = Pe(1,3) + y_wrap + y_wrap_extra
            end if
            if (Pe(2,3) .lt. y_wrap * .25d0) then
              Pe(2,3) = Pe(2,3) + y_wrap + y_wrap_extra
            end if
            if (Pe(3,3) .lt. y_wrap * .25d0) then
              Pe(3,3) = Pe(3,3) + y_wrap + y_wrap_extra
            end if
          end if
        end if
        detPe = Pe(1,1)*(Pe(2,2)*Pe(3,3) - Pe(2,3)*Pe(3,2)) &
                    + Pe(1,2)*(Pe(2,3)*Pe(3,1) - Pe(2,1)*Pe(3,3)) &
                    + Pe(1,3)*(Pe(2,1)*Pe(3,2) - Pe(2,2)*Pe(3,1))
        Area = dabs(detPe)/2.0d0
        t_area(e) = Area
!         print *, "Area", Area
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
        
        du_dx = U(VPtoI(1, t(e,1))) * C(2,1) + U(VPtoI(1, t(e,2))) * C(2,2) + U(VPtoI(1, t(e,3))) * C(2,3)
        dv_dx = U(VPtoI(2, t(e,1))) * C(2,1) + U(VPtoI(2, t(e,2))) * C(2,2) + U(VPtoI(2, t(e,3))) * C(2,3)
        du_dy = U(VPtoI(1, t(e,1))) * C(3,1) + U(VPtoI(1, t(e,2))) * C(3,2) + U(VPtoI(1, t(e,3))) * C(3,3)
        dv_dy = U(VPtoI(2, t(e,1))) * C(3,1) + U(VPtoI(2, t(e,2))) * C(3,2) + U(VPtoI(2, t(e,3))) * C(3,3)
        du_dz = (U(VPtoI(3, t(e,1))) + U(VPtoI(3, t(e,2))) + U(VPtoI(3, t(e,3))))/3d0
        dv_dz = (U(VPtoI(4, t(e,1))) + U(VPtoI(4, t(e,2))) + U(VPtoI(4, t(e,3))))/3d0
        A_factor = 2d0 * du_dx + dv_dy
        B_factor = 2d0 * dv_dy + du_dx
        C_factor = (dv_dx + du_dy) / 2d0

        if (nliter .eq. 1) then
          ice_visc = approx_ice_visc
          e2dstrain = (approx_ice_visc/(flow_law_B / 2d0))**((2d0*flow_law_exponent)/(1d0-flow_law_exponent))
          e3dstrain = e2dstrain
        else  ! GLEN'S FLOW LAW
  !         print *, du_dx, dv_dy, dv_dx, du_dy
          e2dstrain = (flow_law_epsilon + du_dx**2d0 + dv_dy**2d0 + ((du_dy + dv_dx)**2d0)/4d0 + du_dx*dv_dy)
          e3dstrain = e2dstrain + (du_dz**2d0)/4d0 + (dv_dz**2d0)/4d0
          ice_visc = (flow_law_B / 2d0) * (e3dstrain**((1d0-flow_law_exponent)/(2d0*flow_law_exponent)))
  !         print *, "nu, zeta:", ice_visc, ice_visc_zeta, (ice_visc/approx_ice_visc)
        end if
        ice_visc_zeta = ((1d0-flow_law_exponent)/(2d0*flow_law_exponent)) * (flow_law_B / 2d0) * (e3dstrain**((1d0-flow_law_exponent)/(2d0*flow_law_exponent) - 1d0))
        
        element_thickness = thickness(t(e,1:3))
        element_vt = element_thickness * ice_visc
        element_basal_coeff = basal_coeff(t(e,1:3))
        element_basal_coeff_diff = basal_coeff_diff(t(e,1:3))
  !       print *, element_basal_coeff
        element_u1vels = U(VPtoI(1, t(e,1:3)))
        element_u2vels = U(VPtoI(2, t(e,1:3)))
        ! C[2,1] is the partial phi1 / partial x
        ! C[3,1] is the partial phi1 / partial y
        ! etc. these are the same for V instead of phi!!!
#define SETUP_ELT_SUBMATIJ(MATI, MATJ, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2) ( Area * ( \
    ( ALL1(MATI) \
      * ( (TERM11)*C(2,MATI)*C(2,MATJ) \
        + (TERM12)*C(2,MATI)*C(3,MATJ) \
        + (TERM21)*C(3,MATI)*C(2,MATJ) \
        + (TERM22)*C(3,MATI)*C(3,MATJ) \
        ) \
    ) \
    + (int_quad_2b(C(:,(/MATI,MATJ/)), Pe(:,2:3)) * ((ALL2(1)+ALL2(2)+ALL2(3))/3d0)) \
  ) )
#define SETUP_ELT_SUBMAT(MAT, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2) MAT = reshape((/ \
    SETUP_ELT_SUBMATIJ(1, 1, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(2, 1, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(3, 1, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(1, 2, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(2, 2, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(3, 2, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(1, 3, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(2, 3, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2), \
    SETUP_ELT_SUBMATIJ(3, 3, TERM11, TERM12, TERM21, TERM22, ALL1, ALL2) /), shape(MAT))
    
  ! use the macro defn
        SETUP_ELT_SUBMAT(Ke, 4d0, 0d0, 0d0, 1d0, element_vt, element_basal_coeff)
  ! The above should produce:
  !       Ke(1,1) = Area * (4*C(2,1)*C(2,1) + C(3,1)*C(3,1)) * initial_thickness * ice_visc
  !       Ke(1,2) = Area * (4*C(2,1)*C(2,2) + C(3,1)*C(3,2)) * initial_thickness * ice_visc
  !       Ke(1,3) = Area * (4*C(2,1)*C(2,3) + C(3,1)*C(3,3)) * initial_thickness * ice_visc
  !       Ke(2,1) = Area * (4*C(2,2)*C(2,1) + C(3,2)*C(3,1)) * initial_thickness * ice_visc
  !       Ke(2,2) = Area * (4*C(2,2)*C(2,2) + C(3,2)*C(3,2)) * initial_thickness * ice_visc
  !       Ke(2,3) = Area * (4*C(2,2)*C(2,3) + C(3,2)*C(3,3)) * initial_thickness * ice_visc
  !       Ke(3,1) = Area * (4*C(2,3)*C(2,1) + C(3,3)*C(3,1)) * initial_thickness * ice_visc
  !       Ke(3,2) = Area * (4*C(2,3)*C(2,2) + C(3,3)*C(3,2)) * initial_thickness * ice_visc
  !       Ke(3,3) = Area * (4*C(2,3)*C(2,3) + C(3,3)*C(3,3)) * initial_thickness * ice_visc
        ! So the non-zero entries are each line twice and each point. There are nr_p points and nr_t * 3 lines so there are nr_p + nr_t * 6 non-sparse entries
        ! For a node pi, the relevant entries are (in slap column order) the diagonal (pi, pi), nw, n, w, e, s, sw
        ! Boundary nodes only have 4 or 3 points each
#ifdef NEWTONS
        basal_coeff_final = element_basal_coeff_diff*element_u1vels**2d0+element_basal_coeff
        SETUP_ELT_SUBMAT(Je, 4d0*ice_visc + 2d0*ice_visc_zeta*A_factor*A_factor, 2d0*ice_visc_zeta*A_factor*C_factor, 2d0*ice_visc_zeta*C_factor*A_factor, 1d0*ice_visc + 2d0*ice_visc_zeta*C_factor*C_factor, element_thickness, basal_coeff_final)
#endif
        do ti = 1, 3
          do tj = 1, 3
              ki = k_index(VPtoI(1, t(e,ti)),VPtoI(1, t(e,tj)))
              if (ki .eq. nr_k_nz+1) then
                print *, e, ti, tj, t(e,ti), t(e,tj)
              end if
              vKa(ki) = vKa(ki) + (Ke(ti,tj))
#ifdef NEWTONS
              vJa(ki) = vJa(ki) + (Je(ti,tj))
#endif ! NEWTONS
          end do
        end do
  !       if (e .eq. 100) then
  !       print *, Ke
  !       end if
        basal_coeff_final = 0
        SETUP_ELT_SUBMAT(Ke, 0d0, 2d0, 1d0, 0d0, element_vt, basal_coeff_final)
#ifdef NEWTONS
        basal_coeff_final = element_basal_coeff_diff*element_u1vels*element_u2vels
        SETUP_ELT_SUBMAT(Je, 2d0*ice_visc_zeta*A_factor*C_factor, 2d0*ice_visc + 2d0*ice_visc_zeta*A_factor*B_factor, 1d0*ice_visc + 2d0*ice_visc_zeta*C_factor*C_factor, 2d0*ice_visc_zeta*C_factor*B_factor, element_thickness, basal_coeff_final)
#endif ! NEWTONS
        do ti = 1, 3
          do tj = 1, 3
              ki = k_index(VPtoI(1, t(e,ti)),VPtoI(2, t(e,tj)))
              if (ki .eq. nr_k_nz+1) then
                print *, e, ti, tj, t(e,ti), t(e,tj)
              end if
              vKa(ki) = vKa(ki) + (Ke(ti,tj))
#ifdef NEWTONS
              vJa(ki) = vJa(ki) + (Je(ti,tj))
#endif ! NEWTONS
          end do
        end do

        basal_coeff_final = 0
        SETUP_ELT_SUBMAT(Ke, 0d0, 1d0, 2d0, 0d0, element_vt, basal_coeff_final)
#ifdef NEWTONS
        basal_coeff_final = element_basal_coeff_diff*element_u1vels*element_u2vels
        SETUP_ELT_SUBMAT(Je, 2d0*ice_visc_zeta*C_factor*A_factor, 1d0*ice_visc + 2d0*ice_visc_zeta*C_factor*C_factor, 2d0*ice_visc + 2d0*ice_visc_zeta*B_factor*A_factor, 2d0*ice_visc_zeta*B_factor*C_factor, element_thickness, basal_coeff_final)
#endif ! NEWTONS
        do ti = 1, 3
          do tj = 1, 3
              ki = k_index(VPtoI(2, t(e,ti)),VPtoI(1, t(e,tj)))
              if (ki .eq. nr_k_nz+1) then
                print *, e, ti, tj, t(e,ti), t(e,tj)
              end if
              vKa(ki) = vKa(ki) + (Ke(ti,tj))
#ifdef NEWTONS
              vJa(ki) = vJa(ki) + (Je(ti,tj))
#endif ! NEWTONS
          end do
        end do

        SETUP_ELT_SUBMAT(Ke, 1d0, 0d0, 0d0, 4d0, element_vt, element_basal_coeff)
#ifdef NEWTONS
        basal_coeff_final = element_basal_coeff_diff*element_u2vels**2d0+element_basal_coeff
        SETUP_ELT_SUBMAT(Je, 1d0*ice_visc + 2d0*ice_visc_zeta*C_factor*C_factor, 2d0*ice_visc_zeta*C_factor*B_factor, 2d0*ice_visc_zeta*B_factor*C_factor, 4d0*ice_visc + 2d0*ice_visc_zeta*B_factor*B_factor, element_thickness, basal_coeff_final)
#endif ! NEWTONS
        do ti = 1, 3
          do tj = 1, 3
              ki = k_index(VPtoI(2, t(e,ti)),VPtoI(2, t(e,tj)))
              if (ki .eq. nr_k_nz+1) then
                print *, e, ti, tj, t(e,ti), t(e,tj)
              end if
              vKa(ki) = vKa(ki) + (Ke(ti,tj))
#ifdef NEWTONS
              vJa(ki) = vJa(ki) + (Je(ti,tj))
#endif ! NEWTONS
          end do
        end do
        
        ds_dx = SURFACE(t(e,1)) * C(2,1) + SURFACE(t(e,2)) * C(2,2) + SURFACE(t(e,3)) * C(2,3) + ds_dx_fudge
        ds_dy = SURFACE(t(e,1)) * C(3,1) + SURFACE(t(e,2)) * C(3,2) + SURFACE(t(e,3)) * C(3,3) + ds_dy_fudge
  !       print *, ds_dx, ds_dy
        do ti = 1, 3
          pi = t(e, ti)
          pvec(1,1) = (p(t(e, MOD(ti-1+1, 3)+1),1)+p(t(e, MOD(ti-1+2, 3)+1),1)/2.0)-p(pi,1)
          pvec(1,2) = (p(t(e, MOD(ti-1+1, 3)+1),2)+p(t(e, MOD(ti-1+2, 3)+1),2)/2.0)-p(pi,2)
!           print *, pvec(1:2,:)
          pvec(2,:) = pvec(1,:) / DSQRT(pvec(1,1)**2d0+pvec(1,2)**2d0)
          pvec(3,:) = p(t(e, MOD(ti-1+1, 3)+1),:)-p(pi,:)
          pvec(4,:) = p(t(e, MOD(ti-1+2, 3)+1),:)-p(pi,:)
          pvec(5,1) = DACOS((pvec(3,1)*pvec(4,1)+pvec(3,2)*pvec(4,2))/(DSQRT(pvec(3,1)**2d0+pvec(3,2)**2d0)*DSQRT(pvec(4,1)**2d0+pvec(4,2)**2d0)))
          rghdsdx = (ice_density * accel_grav * element_thickness(ti) * ds_dx)
          rghdsdy = (ice_density * accel_grav * element_thickness(ti) * ds_dy)
          F(VPtoI(1, pi)) = F(VPtoI(1, pi)) + 1d0/3d0 * Area * rghdsdx
          F(VPtoI(2, pi)) = F(VPtoI(2, pi)) + 1d0/3d0 * Area * rghdsdy
          if (pvec(1,1) .gt. 1d-10) then
            p_ellipticals(pi,1,1,1) = p_ellipticals(pi,1,1,1) + ice_visc*((SUM(element_thickness)/3d0)/pvec(1,1))*2d0*A_factor*pvec(2,1)
            p_ellipticals(pi,1,1,2) = p_ellipticals(pi,1,1,2) + pvec(2,1)
            p_ellipticals(pi,2,2,1) = p_ellipticals(pi,2,2,1) + ice_visc*(SUM((element_thickness)/3d0)/pvec(1,1))*2d0*C_factor*pvec(2,1)
            p_ellipticals(pi,2,2,2) = p_ellipticals(pi,2,2,2) + pvec(2,1)
          end if
          if (pvec(1,2) .gt. 1d-10) then
            p_ellipticals(pi,1,2,1) = p_ellipticals(pi,1,2,1) + ice_visc*((SUM(element_thickness)/3d0)/pvec(1,2))*2d0*C_factor*pvec(2,2)
            p_ellipticals(pi,1,2,2) = p_ellipticals(pi,1,2,2) + pvec(2,2)
            p_ellipticals(pi,2,1,1) = p_ellipticals(pi,2,1,1) + ice_visc*((SUM(element_thickness)/3d0)/pvec(1,2))*2d0*B_factor*pvec(2,2)
            p_ellipticals(pi,2,1,2) = p_ellipticals(pi,2,1,2) + pvec(2,2)
          end if
          p_ds(pi,1,1) = p_ds(pi,1,1) + ds_dx * pvec(2,1)
          p_ds(pi,1,2) = p_ds(pi,1,2) + pvec(2,1)
          p_ds(pi,2,1) = p_ds(pi,2,1) + ds_dx * pvec(2,2)
          p_ds(pi,2,2) = p_ds(pi,2,2) + pvec(2,2)
          p_visc(pi,1) = p_visc(pi,1) + ice_visc * pvec(5,1)
          p_visc(pi,2) = p_visc(pi,2) + pvec(5,1)
          p_2dstrain(pi,1) = p_2dstrain(pi,1) + e2dstrain * pvec(5,1)
          p_2dstrain(pi,2) = p_2dstrain(pi,2) + pvec(5,1)
  !         F(VPtoI(1, t(e, ti))) = F(VPtoI(1, t(e, ti))) + 1d0/3d0 * Area * (basal_coeff(t(e, ti)) * U(VPtoI(1, t(e, ti))))
  !         F(VPtoI(2, t(e, ti))) = F(VPtoI(2, t(e, ti))) + 1d0/3d0 * Area * (basal_coeff(t(e, ti)) * U(VPtoI(2, t(e, ti))))
        end do

  ! BOUNDARY CONDITIONS !!!!!!
  !       do ti = 1, 3
  !         if (t_boundary(e, ti) .ne. 0) then
  !           JRHS
  !         end if
  !       end do

      end do
      do pi = 1, nr_p
        ki = k_index(VPtoI(3, pi),VPtoI(3, pi))
        vKa(ki) = 1d0
#ifdef NEWTONS
        vJa(ki) = 1d0
#endif ! NEWTONS
        ki = k_index(VPtoI(4, pi),VPtoI(4, pi))
        vKa(ki) = 1d0
#ifdef NEWTONS
        vJa(ki) = 1d0
#endif ! NEWTONS
        call duv_dz( &
          p_visc(pi,1)/p_visc(pi,2), &
          ice_density * accel_grav * thickness(pi) * (p_ds(pi,1,1)/p_ds(pi,1,2)), &
          ice_density * accel_grav * thickness(pi) * (p_ds(pi,2,1)/p_ds(pi,2,2)), &
          p_ellipticals(pi,1,1,1)/p_ellipticals(pi,1,1,2) + p_ellipticals(pi,1,2,1)/p_ellipticals(pi,1,2,2), &
          p_ellipticals(pi,2,1,1)/p_ellipticals(pi,2,1,2) + p_ellipticals(pi,2,2,1)/p_ellipticals(pi,2,2,2), &
          p_2dstrain(pi,1)/p_2dstrain(pi,2), &
          flow_law_B, flow_law_exponent,  p_duv_dz, basal_coeff(pi)*U(VPtoI(1, pi)), basal_coeff(pi)*U(VPtoI(2, pi)))
!         print *,           p_visc(pi,1)/p_visc(pi,2), &
!           ice_density * accel_grav * thickness(pi) * (p_ds(pi,1,1)/p_ds(pi,1,2)), &
!           ice_density * accel_grav * thickness(pi) * (p_ds(pi,2,1)/p_ds(pi,2,2)), &
!           p_ellipticals(pi,1,1,1),p_ellipticals(pi,1,1,2) , p_ellipticals(pi,1,2,1),p_ellipticals(pi,1,2,2), &
!           p_ellipticals(pi,2,1,1)/p_ellipticals(pi,2,1,2) + p_ellipticals(pi,2,2,1)/p_ellipticals(pi,2,2,2), &
!           p_2dstrain(pi,1)/p_2dstrain(pi,2), &
!           flow_law_B, flow_law_exponent, DACOS(-1d0)

!         print *, 'd[uv]/dz', p_duv_dz
        F(VPtoI(3, pi)) = p_duv_dz(3)
        F(VPtoI(4, pi)) = p_duv_dz(4)
      end do
  !     print *, F(:)
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
    !     print *, iK(ki), kj, vKa(ki)
    !   end do
      print *, 'done assy'
    ! ------------------------------------------------------------------------
    ! BOUNDARY CONDITIONS !!!!!!
    ! ------------------------------------------------------------------------
    do bi = 1, nr_calving_front
!       (sqrt((p(calving_front(bi-1),1)-p(calving_front(bi),1))**2+((p(calving_front(bi-1),2)-p(calving_front(bi),2))**2))/2.0d0 + sqrt((p(calving_front(bi+1),1)-p(calving_front(bi),1))**2+((p(calving_front(bi+1),2)-p(calving_front(bi),2))**2))/2.0d0)
      F(VPtoI(2, calving_front(bi))) = F(VPtoI(2, calving_front(bi))) - (ice_density*accel_grav*thickness(calving_front(bi))/(2.0d0))*(1.0d0-ice_density/seawater_density) * boundary_point_ny(calving_front(bi)) * thickness(calving_front(bi))
      F(VPtoI(1, calving_front(bi))) = F(VPtoI(1, calving_front(bi))) - (ice_density*accel_grav*thickness(calving_front(bi))/(2.0d0))*(1.0d0-ice_density/seawater_density) * boundary_point_nx(calving_front(bi)) * thickness(calving_front(bi))
    end do
#ifndef LAGRANGEBC
      do ki = 1, nr_p
        do bi = 1, nr_dirichlet
          do vi = 1, nr_v
            do vj = 1, nr_v
            ! wtf strang
  !             if (k_index(VPtoI(vi, ki), VPtoI(vj, b(bi))) .ne. nr_k_nz+1) then
  !               vKa(k_index(VPtoI(vi, ki), VPtoI(vj, b(bi)))) = 0.0d0
  !             end if
              if (k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, ki)) .ne. nr_k_nz+1) then
                vKa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, ki))) = 0d0
#ifdef NEWTONS
                vJa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, ki))) = 0d0
#endif ! NEWTONS
              end if
            end do
          end do
        end do
      end do
#endif
      do bi = 1, nr_dirichlet
        do vi = 1, nr_v
          do vj = 1, nr_v
  !         if (b(bi) .ne. b_s_all(2) .and. b(bi) .ne. b_s_all(m-1) .and. (b(bi) .ne. b_n_all(2) .and. b(bi) .ne. b_n_all(m-1))) then
#ifdef LAGRANGEBC
            vKa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, dirichlet(bi)))) = vKa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vi, dirichlet(bi)))) + lagrange_multiplier * 1.0d0
            F(VPtoI(vi, dirichlet(bi))) = dirichlet_val(vi, bi) * lagrange_multiplier
#ifdef NEWTONS
            vJa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, dirichlet(bi)))) = vJa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vi, dirichlet(bi)))) + lagrange_multiplier * 1.0d0
#endif ! NEWTONS
#else ! row-decimation bc
            vKa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, dirichlet(bi)))) = 1.0d0
            F(VPtoI(vi, dirichlet(bi))) = dirichlet_val(vi, bi)
#ifdef NEWTONS
            vJa(k_index(VPtoI(vi, dirichlet(bi)), VPtoI(vj, dirichlet(bi)))) = 1.0d0
#endif ! NEWTONS
#endif ! LAGRANGEBC
  !         end if
          end do
        end do
      end do
    !   print *, vk(k_index(VPtoI(1, 8),VPtoI(1, 9)))
    !   print *, vk(k_index(VPtoI(1, 8),VPtoI(1, 8)))
      vKb = vKa
      CALL DSMV(linear_order, U, R, nr_k_nz, iK, jK, vKb, 0)
  !       delta_U = 1d0
      R = F - R
      ! boundary conditions
      do bi = 1, nr_b
  !           R(VPtoI(2, dirichlet(bi))) = 0d0
  !           R(VPtoI(1, dirichlet(bi))) = 0d0
      end do
      if (nliter .eq. 1) then
        linear_tol = final_tol
      end if
      if (nliter .eq. 2) then
        linear_tol = initial_tol
      end if

#ifdef NEWTONS
      if (nliter .ge. 2) then
        print *, "Newton update norm: ", nl_resid
        nl_resid = sqrt(sum((R)**2d0))
        print *, "Newton resid norm: ", nl_resid
        vJb = vJa
        Rb = R
        CALL DSDBCG(linear_order, Rb, delta_U, nr_k_nz, iJ, jJ, vJb, &
        0, 1, linear_tol, 2*linear_order, iter, converr, ierr, 0, rwork, linear_order * 9, iwork, 100)
  !      CALL DSDGMR(linear_order, r, delta_U, nr_k_nz, iJ, jJ, vJb, &
  !        0, 10, 0, linear_tol, 2*linear_order, iter, converr, ierr, 0, rwork, 131 + linear_order * 18, iwork, 100)
        U = U_prev + delta_U
      else
#endif
        vKb = vKa
        Fb = F
        print *, minval(vKa), maxval(vKa)
        print *, "Picard update norm: ", nl_resid
        nl_resid = sqrt(sum((R)**2d0))
        print *, "Picard resid norm: ", nl_resid
        CALL DSDBCG(linear_order, Fb, U, nr_k_nz, iK, jK, vKb, &
          0, 1, linear_tol, 2*linear_order, iter, converr, ierr, 0, rwork, linear_order * 9, iwork, 100)
  !      CALL DSDGMR(linear_order, Fb, U, nr_k_nz, iK, jK, vKb, &
  !        0, 10, 0, linear_tol, 2*linear_order, iter, converr, ierr, 0, rwork, 131 + linear_order * 18, iwork, 100)
#ifdef NEWTONS
      endif
#else
#endif
      nl_resid = sqrt(sum((U_prev - U)**2.0d0))
      print *, nliter, iter, converr, ierr, nl_resid, linear_tol
      if (nl_resid .le. mode_tol .and. nliter .gt. 2) then
        if (linear_tol .le. final_tol) then
          exit
        else
          linear_tol = final_tol
        end if
      end if
      do vi = 1, nr_v
        print *, MINVAL(U(vi:nr_v*nr_p:nr_v)), MAXVAL(U(vi:nr_v*nr_p:nr_v))
      end do  !   do vi = 1, nr_v
    end do
    
    do vi = 1, nr_v
       guess(:,vi) = U(vi::nr_v)
       error_estimate(:,vi) = U(vi::nr_v)-U_prev(vi::nr_v)
       resid_estimate(:,vi) = R(vi::nr_v)
    end do
    max_neighbors = MAXVAL(nr_neighbors)
    
  contains
      integer function k_index(i, j)
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

#define LINTIMESLIN(CO, X, Y) ( (CO(1,1) + CO(2,1) * X + CO(3,1) * Y) * (CO(1,2) + CO(2,2) * X + CO(3,2) * Y) )
#define EVALAT(INTEGRAND, WT, CA, CB) (WT * INTEGRAND(coeffs, (CA*domain(1,1) + CB*domain(2,1) + (1-CA-CB)*domain(3,1)), (CA*domain(1,2) + CB*domain(2,2) + (1-CA-CB)*domain(3,2))))

      real(dp) function int_quad_2a(coeffs, domain)
        real(dp), dimension(3,2), intent(in) :: coeffs
        real(dp), dimension(3,2), intent(in) :: domain
        int_quad_2a = EVALAT(LINTIMESLIN, 0.333333333333333d0, 0.666666666666666d0, 0.166666666666666d0) \
          + EVALAT(LINTIMESLIN, 0.333333333333333d0, 0.166666666666666d0, 0.666666666666666d0) \
          + EVALAT(LINTIMESLIN, 0.333333333333333d0, 0.166666666666666d0, 0.166666666666666d0)
      end function

      real(dp) function int_quad_2b(coeffs, domain)
        real(dp), dimension(3,2), intent(in) :: coeffs
        real(dp), dimension(3,2), intent(in) :: domain
        int_quad_2b = EVALAT(LINTIMESLIN, 0.333333333333333d0, 0.5d0, 0.5d0) \
          + EVALAT(LINTIMESLIN, 0.333333333333333d0, 0.5d0, 0.0d0) \
          + EVALAT(LINTIMESLIN, 0.333333333333333d0, 0.0d0, 0.5d0)
      end function

      real(dp) function int_quad_13(coeffs, domain)
        real(dp), dimension(3,2), intent(in) :: coeffs
        real(dp), dimension(3,2), intent(in) :: domain
        int_quad_13 = \
          ( (   EVALAT(LINTIMESLIN, -0.149570044467670d0, 0.333333333333333d0, 0.333333333333333d0)       \
            +   EVALAT(LINTIMESLIN,  0.175615257433204d0, 0.479308067841923d0, 0.260345966079038d0) )     \
          + (   EVALAT(LINTIMESLIN,  0.175615257433204d0, 0.260345966079038d0, 0.479308067841923d0)       \
            +   EVALAT(LINTIMESLIN,  0.175615257433204d0, 0.260345966079038d0, 0.260345966079038d0) ) )   \
        + ( ( ( EVALAT(LINTIMESLIN,  0.053347235608839d0, 0.869739794195568d0, 0.065130102902216d0)       \
              + EVALAT(LINTIMESLIN,  0.053347235608839d0, 0.065130102902216d0, 0.869739794195568d0)       \
              + EVALAT(LINTIMESLIN,  0.053347235608839d0, 0.065130102902216d0, 0.065130102902216d0) )     \
            + ( EVALAT(LINTIMESLIN,  0.077113760890257d0, 0.638444188569809d0, 0.312865496004875d0)       \
              + EVALAT(LINTIMESLIN,  0.077113760890257d0, 0.312865496004875d0, 0.638444188569809d0) ) )   \
          + ( ( EVALAT(LINTIMESLIN,  0.077113760890257d0, 0.638444188569809d0, 0.486903154253160d0)       \
              + EVALAT(LINTIMESLIN,  0.077113760890257d0, 0.486903154253160d0, 0.638444188569809d0) )     \
            + ( EVALAT(LINTIMESLIN,  0.077113760890257d0, 0.312865496004875d0, 0.486903154253160d0)       \
              + EVALAT(LINTIMESLIN,  0.077113760890257d0, 0.486903154253160d0, 0.312865496004875d0) ) ) )
      end function

      real(dp) function int_lintimeslin(coeffs, domain)
        real(dp), dimension(3,2), intent(in) :: coeffs, domain
        real(dp), dimension(2) :: temp
        real(dp) :: accumulator, ma, oa, mb, ob, xmin, xmax, signfix, cmp
        real(dp) :: av, bv, cv, aw, bw, cw
        real(dp), dimension(3,2) :: sorted
        real(dp) :: blah
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
        accumulator = accumulator * signfix
        
        int_lintimeslin = accumulator
      end function
      
      subroutine maybe_add_neighbor(max_neighbors, node, neighbor, nr, list)
          implicit none
          integer, intent(in) :: max_neighbors
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
      
      subroutine analyze_mesh(max_neighbors, nr_t, t, nr_p, p, nr_neighbors, neighbors, is_boundary_point, boundary_point_nx, boundary_point_ny, min_edge_length)
          implicit none
          integer, intent(in) :: max_neighbors, nr_t, nr_p
          real(dp), dimension(nr_p,2), intent(in) :: p ! List of points
          integer, dimension(nr_t, 3), intent(in) :: t
          integer, dimension(nr_p), intent(out) :: nr_neighbors, is_boundary_point
          integer, dimension(nr_p) :: nr_point_triangles
          integer, dimension(nr_p, max_neighbors), intent(out) :: neighbors
          integer, dimension(nr_p, max_neighbors) :: point_triangles ! triangles that a vertex is on
          real(dp), dimension(nr_p), intent(out) :: boundary_point_nx, boundary_point_ny
          integer, dimension(nr_t, 3) :: triangle_neighbors
          integer, dimension(nr_t) :: boundary_triangles
          integer :: nr_boundary_triangles
          integer :: edgei, pa, pb, ti, triangle, pi, pint
          real(dp) :: edge_nx, edge_ny
          real(dp), intent(out) :: min_edge_length ! real big
          nr_neighbors = 0
          neighbors = 0
          nr_point_triangles = 0
          point_triangles = 0
          triangle_neighbors = 0
          boundary_triangles = 0
          nr_boundary_triangles = 0
          is_boundary_point = 0
          boundary_point_nx = 0d0
          boundary_point_ny = 0d0
          min_edge_length = HUGE(1d0) ! real big
          do e = 1, nr_t
            call maybe_add_neighbor(max_neighbors, t(e,1), t(e,2), nr_neighbors(t(e,1)),neighbors(t(e,1),:))
            call maybe_add_neighbor(max_neighbors, t(e,1), t(e,3), nr_neighbors(t(e,1)),neighbors(t(e,1),:))
            call maybe_add_neighbor(max_neighbors, t(e,2), t(e,1), nr_neighbors(t(e,2)),neighbors(t(e,2),:))
            call maybe_add_neighbor(max_neighbors, t(e,2), t(e,3), nr_neighbors(t(e,2)),neighbors(t(e,2),:))
            call maybe_add_neighbor(max_neighbors, t(e,3), t(e,1), nr_neighbors(t(e,3)),neighbors(t(e,3),:))
            call maybe_add_neighbor(max_neighbors, t(e,3), t(e,2), nr_neighbors(t(e,3)),neighbors(t(e,3),:))
            call maybe_add_neighbor(max_neighbors, t(e,1), e, nr_point_triangles(t(e,1)), point_triangles(t(e,1),:))
            call maybe_add_neighbor(max_neighbors, t(e,2), e, nr_point_triangles(t(e,2)), point_triangles(t(e,2),:))
            call maybe_add_neighbor(max_neighbors, t(e,3), e, nr_point_triangles(t(e,3)), point_triangles(t(e,3),:))
            if (DSQRT((p(t(e,1),1)-p(t(e,2),1))**2d0+(p(t(e,1),2)-p(t(e,2),2))**2d0) .lt. min_edge_length) then
              min_edge_length = DSQRT((p(t(e,1),1)-p(t(e,2),1))**2d0+(p(t(e,1),2)-p(t(e,2),2))**2d0)
            end if
            if (DSQRT((p(t(e,2),1)-p(t(e,3),1))**2d0+(p(t(e,2),2)-p(t(e,3),2))**2d0) .lt. min_edge_length) then
              min_edge_length = DSQRT((p(t(e,2),1)-p(t(e,3),1))**2d0+(p(t(e,2),2)-p(t(e,3),2))**2d0)
            end if
            if (DSQRT((p(t(e,1),1)-p(t(e,3),1))**2d0+(p(t(e,1),2)-p(t(e,3),2))**2d0) .lt. min_edge_length) then
              min_edge_length = DSQRT((p(t(e,1),1)-p(t(e,3),1))**2d0+(p(t(e,1),2)-p(t(e,3),2))**2d0)
            end if
          end do
          do pi = 1, nr_p
            call maybe_add_neighbor(max_neighbors,pi,pi,nr_neighbors(pi),neighbors(pi,:))
          end do
          do e = 1, nr_t
            do edgei = 1, 3 ! for every edge in the triangle
              pa = t(e,MOD(edgei-1,3)+1)
              pb = t(e,MOD(edgei-1+1,3)+1)
              loop_find_neighbor : do ti = 1, nr_point_triangles(pb) ! loop over triangles with vertex pb
                triangle = point_triangles(pb, ti)
                if (triangle .ne. e) then
                  do pi = 1, 3
                    if (t(triangle, pi) .eq. pa) then ! see if one of them also has vertex pa
                    triangle_neighbors(e, edgei) = triangle
  !                     exit loop_find_neighbor
                    end if
                  end do
                end if
              end do loop_find_neighbor
            end do
          end do
          do e = 1, nr_t
            do edgei = 1, 3
              if (triangle_neighbors(e, edgei) .eq. 0) then ! it's a boundary triangle
                nr_boundary_triangles = nr_boundary_triangles + 1
                boundary_triangles(nr_boundary_triangles) = e ! add it to the list
                pa = t(e,MOD(edgei-1,3)+1)
                pb = t(e,MOD(edgei-1+1,3)+1)
                pint = t(e,MOD(edgei-1+2,3)+1) ! "internal" point, the point that is NOT on the boundary edge
                edge_nx = p(pb, 2) - p(pa, 2) ! -y
                edge_ny = p(pa, 1) - p(pb, 1) ! x
                if (edge_nx * (p(pa, 1)-p(pint, 1)) + edge_ny * (p(pa, 2)-p(pint, 2)) .gt. 0d0) then ! dot a to int and a to b
                  edge_nx = -edge_nx ! oops dot product should have been negative
                  edge_ny = -edge_ny
                end if
                is_boundary_point(pa) = is_boundary_point(pa) + 1
                boundary_point_nx(pa) = boundary_point_nx(pa) + edge_nx/2d0 ! each edge contributing half of its length to each point
                boundary_point_ny(pa) = boundary_point_ny(pa) + edge_ny/2d0
                is_boundary_point(pb) = is_boundary_point(pb) + 1
                boundary_point_nx(pb) = boundary_point_nx(pb) + edge_nx/2d0
                boundary_point_ny(pb) = boundary_point_ny(pb) + edge_ny/2d0
              end if
            end do
          end do
          do pi = 1, nr_p
            if (is_boundary_point(pi) .gt. 0) then
              boundary_point_nx(pi) = boundary_point_nx(pi) / real(is_boundary_point(pi), dp)
              boundary_point_ny(pi) = boundary_point_ny(pi) / real(is_boundary_point(pi), dp)
!               print *, pi, "nx", boundary_point_nx(pi), "ny", boundary_point_ny(pi)
            end if
          end do
!           do e = 1, nr_boundary_triangles
!             print *, e, ":", triangle_neighbors(triangle, :)
!           end do
      end subroutine

      subroutine init_k_sparsity(max_neighbors, nr_k_nz, linear_order, iK, jK, vKa, nr_v, nr_t, t, bmp, nr_p, nr_neighbors, neighbors)
          implicit none
          integer, intent(in) :: max_neighbors
          integer, intent(in) :: linear_order
          integer, intent(in) :: nr_k_nz ! Number of entries in the K matrix that are non-zero
          integer, intent(in) :: nr_v, nr_t, nr_p
          integer, dimension(nr_t, 3), intent(in) :: t
          integer, dimension(nr_v,nr_v), intent(in) :: bmp
          integer, dimension(nr_k_nz), intent(out) :: iK
          integer, dimension(linear_order+1), intent(out) :: jK
          real(dp), dimension(nr_k_nz), intent(out) :: vKa
          integer, dimension(nr_p), intent(in) :: nr_neighbors
          integer, dimension(nr_p, max_neighbors), intent(in) :: neighbors
          integer :: e, pi, ki, ni, vi, vj
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
              vKa(ki) = 0.0d0
              ki = ki + 1
              if (vj .le. nr_v2d) then ! we dont need off-diagonal entries for 0d points
                ! Off-diagonal entries in order
                do ni = 1, max_neighbors
                  do vi = 1, nr_v
                    if (neighbors(pi, ni) .ne. 0 .and. (neighbors(pi, ni) .ne. pi &
                      .or. vi .ne. vj) .and. bmp(vi,vj) .ne. 0) then
    !                     print *, ki, pi, ni, neighbors(pi, ni), vi, vj, VPtoI(vi, neighbors(pi, ni))
                      iK(ki) = VPtoI(vi, neighbors(pi, ni))
                      vKa = 0.0d0
                      ki = ki + 1
                    end if
                  end do
                end do
              else
                do vi = 1, nr_v
                  if ((vi .ne. vj) .and. bmp(vi,vj) .ne. 0) then
                    iK(ki) = VPtoI(vi, pi)
                    vKa = 0.0d0
                    ki = ki + 1
                  end if
                end do
              end if
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
      
      subroutine duv_dz(mu, RX, RY, FX, FY, E, B, n, cgret, bX, bY)
        real(dp), intent(in) :: mu ! viscosity
        real(dp), intent(in) :: RX ! surface forcing x
        real(dp), intent(in) :: RY ! surface forcing y
        real(dp), intent(in) :: FX ! elliptical from eq 2
        real(dp), intent(in) :: FY ! elliptical from eq 1
        real(dp), intent(in) :: E  ! 2d strain rate
        real(dp), intent(in) :: B  ! flow law parameter B
        real(dp), intent(in) :: n  ! flow law exponent
        real(dp), intent(in) :: bX ! basal drag X
        real(dp), intent(in) :: bY ! basal drag Y
        real(dp), dimension(4), intent(out) :: cgret
        cgret(1) = 0.2D1 * B ** (-n) * (RY - FY) * (-dble(2 ** (0.1D1 + n)) * mu ** 2 * E * (mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RX ** 2 - 0.2D1 * (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RX * FX + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * FX ** 2 + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RY ** 2 - 0.2D1 * (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RY * FY + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * FY ** 2 + 0.4D1 * (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * mu ** 2 * E) / (n * FX ** 2 + n * RY ** 2 + n * RX ** 2 + n * FY ** 2 + RX ** 2 + FX ** 2 + RY ** 2 + FY ** 2 - 0.2D1 * RX * FX - 0.2D1 * RY * FY - 0.2D1 * n * RX * FX - 0.2D1 * n * RY * FY)
        cgret(2) = 0.2D1 * B ** (-n) * (RX - FX) * (-dble(2 ** (0.1D1 + n)) * mu ** 2 * E * (mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RX ** 2 - 0.2D1 * (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RX * FX + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * FX ** 2 + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RY ** 2 - 0.2D1 * (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * RY * FY + (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * FY ** 2 + 0.4D1 * (RX ** 2 - 0.2D1 * RX * FX + FX ** 2 + RY ** 2 - 0.2D1 * RY * FY + FY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * mu ** 2 * E) / (n * FX ** 2 + n * RY ** 2 + n * RX ** 2 + n * FY ** 2 + RX ** 2 + FX ** 2 + RY ** 2 + FY ** 2 - 0.2D1 * RX * FX - 0.2D1 * RY * FY - 0.2D1 * n * RX * FX - 0.2D1 * n * RY * FY)
        cgret(3) = -0.2D1 * bY * B ** (-n) * (dble(2 ** (0.1D1 + n)) * (mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * mu ** 2 * E - (bX ** 2 + bY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * bX ** 2 - (bX ** 2 + bY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * bY ** 2 - 0.4D1 * (bX ** 2 + bY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * mu ** 2 * E) / (bX ** 2 + bY ** 2 + n * bX ** 2 + n * bY ** 2)
        cgret(4) = -0.2D1 * bX * B ** (-n) * (dble(2 ** (0.1D1 + n)) * (mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * mu ** 2 * E - (bX ** 2 + bY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * bX ** 2 - (bX ** 2 + bY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * bY ** 2 - 0.4D1 * (bX ** 2 + bY ** 2 + 0.4D1 * mu ** 2 * E) ** (n / 0.2D1 - 0.1D1 / 0.2D1) * mu ** 2 * E) / (bX ** 2 + bY ** 2 + n * bX ** 2 + n * bY ** 2)
!         print *, cgret
      end subroutine
end subroutine


subroutine eismint_square_bay2(n, m, nr_p, p, nr_t, t, max_neighbors)
    implicit none
    integer, intent(in) :: n, m, nr_p, nr_t
    real(dp), dimension(nr_p, 2), intent(in) :: p ! List of points
    integer, dimension(nr_t, 3), intent(in) :: t ! List of triangles, boundary flags
    real(dp), dimension(nr_p) :: bed ! List of points, m
    real(dp), dimension(nr_p) :: thickness ! List of points, m
    real(dp), dimension(nr_p) :: basal_traction_coeff ! m/(y*Pa^2)
    real(dp), dimension(nr_p) :: grounded ! grounded flag 1=grounded 0=floating
#define NR_B (2*n+2*m - 4)
    integer :: nr_b ! Number of points on the boundary
    integer, intent(in) :: max_neighbors
    integer, dimension(NR_B) :: b
    integer ::  nr_b_wne
    integer :: b_nw, b_ne, b_sw, b_se
    integer, dimension(m-2) :: b_n, b_s
    integer, dimension(n-2) :: b_e, b_w
    integer, dimension(2*n+m-2) :: b_wne
    integer, dimension(m) :: b_s_all, b_n_all
    integer, dimension(m) :: b_e_all, b_w_all
    integer :: pi, mi, ni
    real(dp), dimension(2*n+m-2) :: dirichlet_val
    dirichlet_val = 0d0
    nr_b_wne = 2*n+m-2
    nr_b = NR_B
    basal_traction_coeff = 0.5d-10
    grounded = 0d0
    thickness = 500d0
    bed = 0d0
!     do pi = 1, nr_p
!         thickness(pi) = max(-1d-7*sqrt((p(pi,1)+1d5)**2+(p(pi,2)+1d5)**2)**2+700d0,500d0)
!         bed(pi) = 2d5 - p(pi,2)*0.01d0
!     end do

    ! Make boundary list
    b_se = 1
    b_sw = m
    b_s = (/(mi, mi = 2, m-1)/)
    b_n = (/(mi+(n-1)*m, mi = 2, m-1)/)
    b_ne = (n-1)*m+1
    b_nw = (n-1)*m+m
    b_e = (/((ni-1)*m+1, ni = 2, n-1)/)
    b_w = (/((ni)*m, ni = 2, n-1)/)
    !   b=(/(mi, mi = 1, m), ((ni-1)*m+1, ni = 2, n), ((ni-1)*m, ni = 3, n+1), (mi+(n-1)*m, mi = 2, m-1) /)
    b = (/b_nw, b_n, b_ne, b_w, b_sw, b_e, b_se, b_s/)
    b_wne = (/b_nw, b_w, b_sw, b_n, b_ne, b_e, b_se/)
    b_n_all = (/b_nw, b_n, b_ne/)
    b_s_all = (/b_sw, b_s, b_se/)
    b_w_all = (/b_nw, b_w, b_sw/)
    b_e_all = (/b_se, b_e, b_se/)
    !   print *, b
    !   do bi = 1, nr_b
    !     print *, t(bi,:) 
    !   end do
    ! This could be a lot faster using sorted whatevers or something...
    !         t_boundary = 0
    !         do bi = 1, nr_b
    !           do e = 1, nr_t
    !             do ti = 1, 3
    !               if (t(e, ti) .eq. b(bi)) then
    !                 if (t(e, mod(ti+1-1,3)+1) .eq. b(mod(bi+1-1,nr_b)+1)) then
    !                   t_boundary(e, ti) = 1
    !                 end if
    !               end if
    !             end do
    !           end do
    !         end do
      call fem_imr(nr_p, p, nr_t, t, bed, thickness, basal_traction_coeff, grounded, nr_b, max_neighbors, m, b_s_all, nr_b_wne, b_wne, dirichlet_val, 0d0, 0d0, 5.7d-18, -1d0, -1d0, 0d0, 0d0)
end subroutine

subroutine ismip_hom_c2(n, m, nr_p, p, nr_t, t, max_neighbors, ismip_size)
    implicit none
    integer, intent(in) :: n, m, nr_p, nr_t
    real(dp), dimension(nr_p, 2), intent(in) :: p ! List of points
    integer, dimension(nr_t, 3), intent(in) :: t ! List of triangles, boundary flags
    real(dp), dimension(nr_p) :: bed ! List of points, m
    real(dp), dimension(nr_p) :: thickness ! List of points, m
    real(dp), dimension(nr_p) :: basal_traction_coeff ! m/(y*Pa^2)
    real(dp), dimension(nr_p) :: grounded ! grounded flag 1=grounded 0=floating
    real(dp), intent(in) :: ismip_size
    integer, intent(in) :: max_neighbors
    integer, dimension(1) :: b
    integer :: pi, mi, ni
    real(dp) :: pic
    pic = 3.14159
    basal_traction_coeff = 0.5d-10
    grounded = 1d0
    thickness = 1000d0
    bed = 0d0
    do pi = 1, nr_p
!       bed = 500d0 * sin(p(pi,1)*2d0*pic/ismip_size) * sin(p(pi,2)*2d0*pic/ismip_size)
!       thickness = 1000d0 - bed
        basal_traction_coeff(pi) = 1000d0 + 1000d0 * dsin(p(pi,1)*2d0*pic/ismip_size) * dsin(p(pi,2)*2d0*pic/ismip_size)
    end do
    call fem_imr(nr_p, p, nr_t, t, bed, thickness, basal_traction_coeff, grounded, 0, max_neighbors, 0, b, 0, b, b, dtan(0.1d0*2d0*pic/360d0), 0d0, 1d-16, ismip_size * (DFLOAT(m)/DFLOAT(m+1)), ismip_size * (DFLOAT(n)/DFLOAT(n+1)), ismip_size/m, ismip_size/n, 2)
end subroutine