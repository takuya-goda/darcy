! SMAC法
! 角棒
! 無次元化

program fortran

    use, intrinsic :: iso_fortran_env
    implicit none

    !=============================================
    !                   変数宣言
    !
    ! ========== メッシュ作成 ==========
    integer, parameter :: nxc = 90
    integer, parameter :: nyc = 60
    integer, parameter :: nxd = nxc + 1
    integer, parameter :: nyd = nyc + 1
    integer, parameter :: nx = nxc*2
    integer, parameter :: ny = nyc*2
    real(real64), parameter :: lx = 0.72d0
    real(real64), parameter :: ly = 0.48d0
    real(real64), parameter :: dx = 0.008d0
    real(real64), parameter :: dy = 0.008d0
    real(real64), parameter :: uin = 0.625 !6.25e-4
    real(real64), parameter :: lt = 5d0

    ! ========== 流体 ==========
    real(real64), parameter :: dens = 1.0e+3
    real(real64), parameter :: viscosity = 1.0e-1
    real(real64), parameter :: kinematic_viscosity = viscosity/dens

    real(real64), parameter :: K  = 0.00025                         ! 透水係数（K小⇒Da小⇒u_aux減）
    real(real64), parameter :: epsilon = 0.45d0                     ! ポロシティ(気孔率)
    real(real64), parameter :: Re = uin*lx/kinematic_viscosity      ! レイノルズ数
    ! real(real64), parameter :: Da = K/(lx**2)                     ! ダルシー数
    real(real64), parameter :: Da = K*kinematic_viscosity/lx**2     ! ダルシー数
    real(real64), parameter :: Forchheimer_coef = 1.5e+6            ! Forchheimer
    real(real64), parameter :: Co = 0.1d0                           ! クーラン数

    real(real64), parameter :: vof_min = 1d-12

    integer, parameter :: ic_porous_start = nint(nxc*0.3d0)
    integer, parameter :: ic_porous_end   = nint(nxc*0.7d0)

    real(real64) :: dt
    real(real64) :: arufa, err_total
    real(real64) :: start_time, end_time
    integer :: i, j, time, nt, outputstep
    integer :: ic, jc
    real(real64) :: x, y
    real(real64) :: coord_x, coord_y


    ! ========== 物理量 ==========
    real(real64), dimension(0:nxd + 1, 0:nyc + 1) :: u
    real(real64), dimension(0:nxc + 1, 0:nyd + 1) :: v
    real(real64), dimension(0:nxd + 1, 0:nyc + 1) :: uold
    real(real64), dimension(0:nxc + 1, 0:nyd + 1) :: vold
    real(real64), dimension(0:nxd + 1, 0:nyc + 1) :: u_aux
    real(real64), dimension(0:nxc + 1, 0:nyd + 1) :: v_aux
    real(real64), dimension(0:nxc + 1, 0:nyc + 1) :: p
    real(real64), dimension(0:nxc + 1, 0:nyc + 1) :: eps
    real(real64), dimension(0:nxc + 1, 0:nyc + 1) :: phi
    real(real64), dimension(0:nxc + 1, 0:nyc + 1) :: dp
    real(real64), dimension(0:nxc + 1, 0:nyc + 1) :: vof
    real(real64), dimension(0:nxd + 1, 0:nyc + 1) :: vof_flux_x
    real(real64), dimension(0:nxc + 1, 0:nyd + 1) :: vof_flux_y
    real(real64), dimension(nxc, nyc) :: theta
    real(real64), dimension(0:nxc + 1, 0:nyc + 1) :: force

    real(real64), dimension(nxc, nyc) :: coordinate

    ! ========== 円柱解析 ==========
    real(8), parameter :: center(2) = [lx*0.5d0,ly*0.5d0]
    real(8), parameter :: scale_prm = dx*0.5d0
    real(8), parameter :: radius = 0.016d0
    integer, parameter :: vp = 6
    integer, parameter :: vpc = 3
    real(8), dimension(1-vp:nx+1+vp,1-vp:ny+1+vp) :: dist
    real(8), dimension(1-vp:nx+1+vp,1-vp:ny+1+vp) :: frac

    type :: velocity
        real(8)::x(-2:nxd+3, -2:nyc+3)
        real(8)::y(-2:nxc+3, -2:nyd+3)
    end type

    type(velocity)::aprt


     !=============================================
    !                   初期設定
    !

    ! ========== 計算 ==========
    ! dt = Co*dx/uin
    dt = 0.00001d0 ! 0.00001d0
    ! nt = nint(lt/dt)

    ! ========== 物理量（パラメーター） ==========
    err_total = 1d0/1e+12
    arufa = 1.99d0

    ! ========== 物理量（配列） ==========
    u(0, :)     = uin
    v(:, :)     = 0d0
    u_aux(:, :) = 0d0
    v_aux(:, :) = 0d0
    p(:, :)     = 0d0
    phi(:, :)   = 0d0
    vof(:,:)    = 0d0
    theta(:, :) = 0d0
    eps(:, :)   = 0d0
    eps(ic_porous_start:ic_porous_end,:) = epsilon
    coordinate(:, :)   = 0d0
    coord_x = 0d0 ; coord_y = 0.24d0

    ! 外力
    force(:,:)   = 0d0
    ! force(110:120,:) = -0.01d0 !ver1

    ! ========== 境界条件 ==========
    u(1, :) = uin
    u(:, 0) = u(:, 1)
    u(:, nyc + 1) = u(:, nyc)

    v(0, :) = v(1, :)
    v(:, 0) = 0d0
    v(:, nyd) = 0d0

    vof(1,:) = 1d0

    ! ========== 計算量 ==========
    nt = 900000
    outputstep = 1000


    call set_initial(dist,frac,aprt,eps)
    call setBoundary(u, v, uin)

    ! open(unit=30,file="test.txt")
    ! do j=1,nyc
    ! do i=1,nxc
    !     write(30,*) i,j,eps(i,j)
    !     ! write(30,*) i,j,dist(i*2,j*2)
    ! enddo
    ! enddo
    ! close(30)

! メインルーチン
    call output(u, v, p, phi, 0, theta, vof, coordinate)
    call cpu_time(start_time)
    print *, "Reynolds number = ",Re
    print *, "Darcy    number = ",Da
    print *, "couran   number = ",uin*dt/dx
    ! do time = 1, nt
    do time = 1, 10000
        print *, time,"/",nt,dt*time
        ! print *, vof_flux_x(1,10),vof_flux_x(2,10),vof_flux_x(3,10)
        call computeAuxiallyVelocity(u_aux, v_aux, u, v, p, dx, dy, dt, uin, vof, force, eps, aprt)
        call computeDivergenceAuxiallyVelocity(theta, u_aux, v_aux, dx, dy, dt)
        call computePressurePoisson(p, phi, dp, u, v, dx, dy, theta, err_total, arufa, dens, dt, aprt)
        call computeVelocity(u, v, uold, vold, u_aux, v_aux, dt, dx, dy, dens, phi, uin, vof, aprt)
        call integrate_vof(vof, vof_flux_x, vof_flux_y, u, v, aprt)
        call integrate_coordinate(coordinate, coord_x, coord_y, u, v, uold, vold)
        if (mod(time, outputstep) == 0) then
            call output(u, v, p, phi, time, theta, vof, coordinate)
        elseif (mod(time, 100) == 0 .and. time < 300) then
            call output(u, v, p, phi, time, theta, vof, coordinate)
        end if
    end do
    call cpu_time(end_time)
    write (*, *) "elapsed time = ", end_time - start_time

    open (unit=30, file="elapsed.txt")
    write (30, *) end_time - start_time
    close (30)

contains

    subroutine output(u, v, p, phi, time, theta, vof, coordinate)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(in) :: u(0:, 0:)
        real(real64), intent(in) :: v(0:, 0:)
        real(real64), intent(in) :: p(0:, 0:)
        real(real64), intent(in) :: phi(0:, 0:)
        integer, intent(in) :: time
        real(real64), intent(in) :: theta(:, :)
        real(real64), intent(in) :: vof(0:, 0:)
        real(real64), intent(in) :: coordinate(:,:)

        integer :: i, j
        real(real64) :: uout, vout
        character*128 filename1, filename2, filename3

        write (filename1, '("pres/pres",i6.6,".txt")') time
        write (filename2, '("vof/vof",i6.6,".txt")') time
        write (filename3, '("crd/coord",i5.5,".txt")') time
        open (unit=10, file=filename1)
        open (unit=20, file=filename2)
        open (unit=30, file=filename3)
        do j = 1, nyc
        do i = 1, nxc
            write (10, '(2i8,3f23.14)') i, j, p(i, j), u(i,j),v(i,j)
            write (20, '(2i8,3f23.14)') i, j, vof(i, j), u(i,j),v(i,j)
            write (30, '(2i8,3f23.14)') i, j, coordinate(i, j), u(i,j),v(i,j)
        end do
        end do
        close (10)
        close (20)
        close (30)

    end subroutine

    subroutine set_initial(dist,frac,aprt,eps)
        implicit none
        real(8),intent(inout)::dist(1-vp:nx+1+vp,1-vp:ny+1+vp)
        real(8),intent(inout)::frac(1-vp:nx+1+vp,1-vp:ny+1+vp)
        type(velocity),intent(inout)::aprt
        real(8),intent(inout)::eps(0:nxc+1,0:nyc+1)

        integer::i,j

        ! dist（距離）
        distination: block
            real(8)::r,x,y
            do j=1-vp,ny+1+vp
            do i=1-vp,nx+1+vp
                x = dble(i-1)*scale_prm - center(1)
                y = dble(j-1)*scale_prm - center(2)
                dist(i,j) = sqrt(x**2 + y**2) - radius
            enddo
            enddo
        end block distination

        ! frac(流体率)
        fraction: block !実際の格子点の定義
            integer::nx_,ny_
            nx_ = nx + 1
            ny_ = ny + 1
            do j=1-vp,ny_+1+vp,2
            do i=1-vp,nx_+1+vp,2
                if(dist(i,j) > 0d0) frac(i,j) = 1d0
            enddo
            enddo
        end block fraction

        h_face:block !X軸仮想格子点
            integer::ic
            integer::nx_,ny_
            real(8)::dist_im1,dist_ip1,ls_positive,l
            nx_ = nx + 1
            ny_ = ny + 1
            do j =1-vp  ,ny_+vp  ,2
            do ic=1-vp+1,nx_+vp-1,2
                dist_im1 = dist(ic-1,j)
                dist_ip1 = dist(ic+1,j)
                if(dist_im1 > 0d0 .and. dist_ip1 > 0d0)then
                    frac(ic,j) = 1d0
                elseif(dist_im1 <= 0d0 .and. dist_ip1 <= 0d0)then
                    frac(ic,j) = 0d0
                else
                    ls_positive = max(dist_im1,dist_ip1)
                    l = ls_positive/(abs(dist_im1) + abs(dist_ip1))
                    frac(ic,j) = min(1d0,l)
                endif
            enddo
            enddo
        end block h_face


        v_face:block !Y軸仮想格子点
            integer::jc
            integer::nx_,ny_
            real(8)::dist_im1,dist_ip1,ls_positive,l
            nx_ = nx + 1
            ny_ = ny + 1
            do jc=1-vp+1,ny_+vp-1,2
            do i =1-vp  ,nx_+vp  ,2
                dist_im1 = dist(i,jc-1)
                dist_ip1 = dist(i,jc+1)
                if(dist_im1 > 0d0 .and. dist_ip1 > 0d0)then
                    frac(i,jc) = 1d0
                elseif(dist_im1 <= 0d0 .and. dist_ip1 <= 0d0)then
                    frac(i,jc) = 0d0
                else
                    ls_positive = max(dist_im1,dist_ip1)
                    l = ls_positive/(abs(dist_im1) + abs(dist_ip1))
                    frac(i,jc) = min(1d0,l)
                endif
            enddo
            enddo
        end block v_face

        cell_center:block  !仮想格子点
            integer::ic,jc
            real(8)::half1,half2
            integer::nx_,ny_
            nx_ = nx + 1
            ny_ = ny + 1
            do jc=1-vp+1,ny_+vp-1,2
            do ic=1-vp+1,nx_+vp-1,2
                half1 = compute_volume_fraction_in_triangle(dist(ic+1,jc-1),dist(ic-1,jc-1),dist(ic-1,jc+1))
                half2 = compute_volume_fraction_in_triangle(dist(ic-1,jc+1),dist(ic+1,jc+1),dist(ic+1,jc-1))
                frac(ic,jc) = half1 + half2
            enddo
            enddo
        end block cell_center

        aprt_x:block
            integer::jc
            integer::i_fine,j_fine
            do jc = 1-vpc,nyc+vpc
            do i  = 1-vpc,nxc+vpc+1
                i_fine = 2*(i-1) + 1
                j_fine = 2*jc
                aprt%x(i,jc) = frac(i_fine,j_fine)
                ! print *, i,jc,aprt%x(i,jc)
                ! write(50,*) i,jc,aprt%x(i,jc)
            enddo
            enddo
        end block aprt_x

        aprt_y:block
            integer::ic
            integer::i_fine,j_fine
            do j  = 1-vpc,nyc+vpc+1
            do ic = 1-vpc,nxc+vpc
                i_fine = 2*ic
                j_fine = 2*(jc-1) + 1
                aprt%y(ic,j) = frac(i_fine,j_fine)
            enddo
            enddo
        end block aprt_y

        ! ! ε（透水パラメーター）
        ! epsilon_block: block
        !     integer::i_fine,j_fine
        !     do jc=1,nyc
        !     do ic=1,nxc
        !         i_fine = 2*ic
        !         j_fine = 2*jc
        !         if(dist(i_fine,j_fine) < 0 )then
        !             eps(ic,jc) = 0d0
        !         endif
        !     enddo
        !     enddo
        ! end block epsilon_block


    end subroutine

    subroutine setBoundary(u, v, uin)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(inout) :: u(0:, 0:)
        real(real64), intent(inout) :: v(0:, 0:)
        real(real64), intent(in) :: uin

        !&<
        ! ========== 流入境界 ==========
        u(1, :) = uin
        v(0, :) = v(1, :)

        ! ========== 壁面 ==========
        u(:, 0) = u(:, 1)     ! slip
        ! u(:, 0) = -u(:, 1)  ! non-slip
        v(:, 1) = 0d0

        u(:, nyc + 1) = u(:, nyc)    !non- slip
        ! u(:, nyc + 1) = -u(:, nyc) ! slip
        v(:, nyd) = 0d0

        ! ========== 流出境界 ==========
        u(nxd + 1, :) = u(nxd, :)
        v(nxc + 1, :) = v(nxc, :)

        !&>

    end subroutine

    subroutine computeAuxiallyVelocity(u_aux, v_aux, u, v, p, dx, dy, dt, uin, vof, force, eps, aprt)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(inout) :: u_aux(0:, 0:)
        real(real64), intent(inout) :: v_aux(0:, 0:)
        real(real64), intent(in) :: u(0:, 0:)
        real(real64), intent(in) :: v(0:, 0:)
        real(real64), intent(in) :: p(0:, 0:)
        real(real64), intent(in) :: vof(0:, 0:)
        real(real64), intent(in) :: force(0:, 0:)
        real(real64), intent(in) :: eps(0:, 0:)
        real(real64), intent(in) :: dx, dy, dt, uin
        type(velocity),intent(in)::aprt

        real(real64)::H_f
        real(real64)::test1,test2,test3

        !&<
        do jc = 1, nyc
        do i = 1, nxd
            j = jc
            ic = i
            ! H_f = aprt%x(i,jc)
            H_f = 1d0
            
            u_aux(i, jc) = u(i, jc) - dt*H_f*((u(i  , jc ) + u(i-1 , jc ))/2d0*(u(i  , jc  ) - u(i-1, jc  ))/dx &
                                            + (u(i+1, jc ) + u(i   , jc ))/2d0*(u(i+1, jc  ) - u(i  , jc  ))/dx &
                                            + (v(ic , j+1) + v(ic-1, j+1))/2d0*(u(i  , jc+1) - u(i  , jc  ))/dy &
                                            + (v(ic , j  ) + v(ic-1, j  ))/2d0*(u(i  , jc  ) - u(i  , jc-1))/dy)/2d0 &
                         + dt*H_f/Re*((u(i+1, jc  ) - 2d0*u(i, jc) + u(i-1, jc  ))/dx**2  &
                                    + (u(i  , jc+1) - 2d0*u(i, jc) + u(i  , jc-1))/dy**2) &
                         - dt*(p(ic, jc) - p(ic - 1, jc))/dx &
                         - dt/(Re*Da)*eps(ic,jc)*u(i,jc) &
                         - dt*eps(ic,jc)*abs(u(i,jc))*u(i,jc)/sqrt(Da)
            test2 = - dt/(Re*Da)*eps(ic,jc)*u(i,jc)
            test3 = - dt*eps(ic,jc)*abs(u(i,jc))*u(i,jc)/sqrt(Da)
        end do
        end do

        do j = 1, nyd
        do ic = 1, nxc
            jc = j
            i = ic
            ! H_f = aprt%y(ic,j)
            H_f = 1d0

            v_aux(ic, j) = v(ic, j) - dt*H_f*((u(i  , jc ) + u(i  , jc-1))/2d0*(v(ic  , j  ) - v(ic-1, j  ))/dx &
                                            + (u(i+1, jc ) + u(i+1, jc-1))/2d0*(v(ic+1, j  ) - v(ic  , j  ))/dx &
                                            + (v(ic , j+1) + v(ic , j   ))/2d0*(v(ic  , j+1) - v(ic  , j  ))/dy &
                                            + (v(ic , j  ) + v(ic , j-1 ))/2d0*(v(ic  , j  ) - v(ic  , j-1))/dy)/2d0 &
                         + dt*H_f/Re*((v(ic+1, j  ) - 2d0*v(ic, j) + v(ic-1, j  ))/dx**2  &
                                    + (v(ic  , j+1) - 2d0*v(ic, j) + v(ic  , j-1))/dy**2) &
                        - dt*(p(ic, jc) - p(ic, jc - 1))/dy &
                        - dt/(Re*Da)*eps(ic,jc)*v(ic,j) &
                         - dt*eps(ic,jc)*abs(v(ic,j))*v(ic,j)/sqrt(Da)
        end do
        end do
        !&>

        call setBoundary(u_aux, v_aux, uin)

    end subroutine

    subroutine computeDivergenceAuxiallyVelocity(theta, u_aux, v_aux, dx, dy, dt)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(inout) :: theta(:, :)
        real(real64), intent(in) :: u_aux(0:, 0:)
        real(real64), intent(in) :: v_aux(0:, 0:)

        real(real64), intent(in) :: dt, dx, dy

        do jc = 1, nyc
        do ic = 1, nxc
            j = jc
            i = ic
            theta(ic, jc) = (u_aux(i + 1, jc) - u_aux(i, jc))/dx + (v_aux(ic, j + 1) - v_aux(ic, j))/dy
        end do
        end do

    end subroutine

    subroutine computePressurePoisson(p, phi, dp, u, v, dx, dy, theta, err_total, arufa, dens, dt, aprt)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(inout) :: p(0:, 0:)
        real(real64), intent(inout) :: phi(0:, 0:)
        real(real64), intent(inout) :: dp(0:, 0:)
        real(real64), intent(in) :: u(0:, 0:)
        real(real64), intent(in) :: v(0:, 0:)
        real(real64), intent(in) :: theta(:, :)
        type(velocity),intent(in)::aprt

        real(real64), intent(in) :: dx, dy, err_total, arufa, dens, dt
        real(real64) :: err, err_phi, err_dp, err_p, test_sum
        real(real64)::im1,ip1,jm1,jp1
        real(real64)::original_pressure
        integer count

        dp(:, :) = 0d0
        phi(:, :) = 0d0 
        err = 1d0
        count = 0
        test_sum = 0d0

        !&<
        do while (err > err_total)
            err_p = 0d0
            err_dp = 0d0
            err_phi = 0d0
            do jc = 1, nyc
            do ic = 1, nxc
                j = jc
                i = ic
                ! im1 = aprt%x(i  ,jc )
                ! ip1 = aprt%x(i+1,jc )
                ! jm1 = aprt%y(ic ,j  )
                ! jp1 = aprt%y(ic ,j+1)                
                im1 = 1d0
                ip1 = 1d0
                jm1 = 1d0
                jp1 = 1d0                
                dp(ic, jc) = (( dy**2d0*(ip1*phi(ic + 1, jc) + im1*phi(ic - 1, jc)) &
                              + dx**2d0*(jp1*phi(ic, jc + 1) + jm1*phi(ic, jc - 1)) &
                              - dx**2d0*dy**2d0*theta(ic, jc)/dt )/(2d0*(dx**2d0 + dy**2d0))) - phi(ic, jc)
                phi(ic, jc) = phi(ic, jc) + arufa*dp(ic, jc)
                err_dp = err_dp + abs(dp(ic, jc))
                err_phi = err_phi + abs(phi(ic, jc))
            end do
            end do

            ! slip
            phi(:,0    ) = phi(:,1  )
            phi(:,nyc+1) = phi(:,nyc)

            original_pressure = minval(p)

            if (err_phi < 1.0d-20) then
                err_phi = 1d0
            end if
            err = err_dp/err_phi

        end do
        !&>
        p(:, :) = p(:, :) - original_pressure
        p(:, :) = p(:, :) + phi(:, :) 

    end subroutine

    subroutine computeVelocity(u, v, uold, vold, u_aux, v_aux, dt, dx, dy, dens, phi, uin, vof, aprt)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(inout) :: u(0:, 0:)
        real(real64), intent(inout) :: v(0:, 0:)
        real(real64), intent(inout) :: uold(0:, 0:)
        real(real64), intent(inout) :: vold(0:, 0:)
        real(real64), intent(inout) :: phi(0:, 0:)
        real(real64), intent(in) :: u_aux(0:, 0:)
        real(real64), intent(in) :: v_aux(0:, 0:)
        real(real64), intent(in) :: vof(0:, 0:)
        type(velocity),intent(in)::aprt

        real(real64), intent(in) :: dt, dx, dy, dens, uin
        real(real64)::H_f

        uold(:,:) = u(:,:)
        vold(:,:) = v(:,:)
        do jc = 1, nyc
        do i = 1, nxd
            j = jc
            ic = i
            ! H_f = aprt%x(i,jc)
            H_f = 1d0
            u(i, jc) = u_aux(i, jc) - dt*H_f*(phi(ic, jc) - phi(ic - 1, jc))/dx/dens
        end do
        end do

        do j = 1, nyd
        do ic = 1, nxc
            jc = j
            i = ic
            ! H_f = aprt%y(ic,j)
            H_f = 1d0
            v(ic, j) = v_aux(ic, j) - dt*H_f*(phi(ic, jc) - phi(ic, jc - 1))/dy/dens
        end do
        end do

        call setBoundary(u, v, uin)

    end subroutine

    subroutine integrate_vof(vof,vof_flux_x,vof_flux_y,u,v,aprt)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64),intent(inout)::vof(0:,0:)
        real(real64),intent(inout)::vof_flux_x(0:,0:)
        real(real64),intent(inout)::vof_flux_y(0:,0:)
        real(real64),intent(in)::u(0:,0:)
        real(real64),intent(in)::v(0:,0:)
        type(velocity),intent(in)::aprt

        integer::id,jd
        real(real64)::ip1,im1,jp1,jm1

        ! ========== X方向移流量計算 ==========
        do jc=1,nyc
        do id=1,nxd
            ic = id
            vof_flux_x(id,jc) = u(id  ,jc)*(vof(ic,jc) - vof(ic-1,jc))/dx
        enddo
        enddo

        ! ========== Y方向移流量計算 ==========
        do jd=1,nyd
        do ic=1,nxc
            jc = jd
            vof_flux_y(ic,jd) = v(ic  ,jd)*(vof(ic,jc) - vof(ic,jc-1))/dy
        enddo
        enddo

        ! ========== VOF値移流量計算 ==========
        do jc=1,nyc
        do ic=1,nxc
            jd = jc
            id = ic
            ! im1 = aprt%x(i  ,jc )
            ! ip1 = aprt%x(i+1,jc )
            ! jm1 = aprt%y(ic ,j  )
            ! jp1 = aprt%y(ic ,j+1)                
            im1 = 1d0
            ip1 = 1d0
            jm1 = 1d0
            jp1 = 1d0                
            ! vof(ic,jc) = vof(ic,jc) - dt*(ip1*vof_flux_x(id+1,jc  ) + im1*vof_flux_x(id,jc))/2d0 &
            !                         - dt*(jp1*vof_flux_y(ic  ,jd+1) + jm1*vof_flux_y(ic,jd))/2d0
            vof(ic,jc) = vof(ic,jc) - dt*(vof_flux_x(id+1,jc  ) + vof_flux_x(id,jc))/2d0 &
                                    - dt*(vof_flux_y(ic  ,jd+1) + vof_flux_y(ic,jd))/2d0
            ! 上限下限の補正
            if (vof(ic,jc) < vof_min) then
                vof(ic,jc) = 0d0
            elseif (vof(ic,jc) > 1d0 - vof_min) then
                vof(ic,jc) = 1d0
            endif
            ! if(ic==1.and.jc==10)then
            !     print *, vof(ic,jc),-dt*(vof_flux_x(id+1,jc)+vof_flux_x(id,jc))/2d0
            ! elseif(ic==2.and.jc==10)then
            !     print *, vof(ic,jc),-dt*(vof_flux_x(id+1,jc)+vof_flux_x(id,jc))/2d0
            ! elseif(ic==3.and.jc==10)then
            !     print *, vof(ic,jc),-dt*(vof_flux_x(id+1,jc)+vof_flux_x(id,jc))/2d0
            ! endif
        enddo
        enddo

        vof(1,:) = 1d0

    end subroutine

    subroutine integrate_coordinate(coordinate,coord_x,coord_y,u,v,uold,vold)
        use,intrinsic::iso_fortran_env
        implicit none
        real(real64),intent(inout) :: coordinate(:,:)
        real(real64),intent(inout) :: coord_x,coord_y
        real(real64),intent(in)::u(0:,0:)
        real(real64),intent(in)::v(0:,0:)
        real(real64),intent(in)::uold(0:,0:)
        real(real64),intent(in)::vold(0:,0:)
        real(real64)::coord_xold,coord_yold
        integer::coord_i,coord_j

        coord_xold = coord_x
        coord_yold = coord_y

        ! coord（座標）
        coordination: block
            real(8)::r,x,y
            real(8)::min_val
            integer::i,j
            min_val = 100d0
            do j=1,nyc
            do i=1,nxc
                x = dble(i-1)*dx - coord_x
                y = dble(j-1)*dy - coord_y
                coordinate(i,j) = sqrt(x**2 + y**2)
                if(min_val > coordinate(i,j)) then
                    min_val = coordinate(i,j)
                    coord_i = i
                    coord_j = j
                endif
                if(coordinate(i,j) > 0.02d0) coordinate(i,j) = 1d0
            enddo
            enddo
        end block coordination

        ! 座標更新
        coord_x = coord_xold + u(coord_i,coord_j)*dt + 0.5d0*(u(coord_i,coord_j) - uold(coord_i,coord_j))/dt*dt**2
        coord_y = coord_yold + v(coord_i,coord_j)*dt + 0.5d0*(v(coord_i,coord_j) - vold(coord_i,coord_j))/dt*dt**2

    end subroutine

    function compute_volume_fraction_in_triangle(dist1,dist2,dist3) result(frac)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64),intent(in) :: dist1,dist2,dist3
        real(real64)::dist_p1,dist_p2,dist_m1,dist_m2
        real(real64)::frac

        if(dist1 > 0d0 .and. dist2 > 0d0 .and. dist3 > 0d0)then
            frac = 0.5d0
        elseif(dist1 <= 0d0 .and. dist2 <= 0d0 .and. dist3 <= 0d0)then
            frac = 0d0
        elseif(     (dist1 > 0d0 .and. dist2 <= 0d0 .and. dist3 <= 0d0) &
               .or. (dist2 > 0d0 .and. dist3 <= 0d0 .and. dist1 <= 0d0) & 
               .or. (dist3 > 0d0 .and. dist1 <= 0d0 .and. dist2 <= 0d0) )then
                if(dist1 > 0d0)then
                    dist_p1 = dist1
                    dist_m1 = dist2
                    dist_m2 = dist3
                elseif(dist2 > 0d0)then
                    dist_p1 = dist2
                    dist_m1 = dist1
                    dist_m2 = dist3
                else
                    dist_p1 = dist3
                    dist_m1 = dist1
                    dist_m2 = dist2
                endif
                frac = dist_p1**2/((dist_m1 - dist_p1)*(dist_m2 - dist_p1))
                frac = frac/2d0
        else
            if(dist1 <= 0d0)then
                dist_m1 = dist1
                dist_p1 = dist2
                dist_p2 = dist3
            elseif(dist2 <= 0d0)then
                dist_m1 = dist2
                dist_p1 = dist1
                dist_p2 = dist3
            else
                dist_m1 = dist3
                dist_p1 = dist1
                dist_p2 = dist2
            endif
            frac = 1d0 - dist_m1**2/((dist_m1 - dist_p1)*(dist_m1 - dist_p2))
            frac = frac/2d0
        endif

    end function

end program
