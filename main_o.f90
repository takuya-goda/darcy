! SMAC法
! 角棒
! 無次元化
! 単相流(大島先生論文実装)

program fortran

    use, intrinsic :: iso_fortran_env
    implicit none

    !=============================================
    !                   変数宣言
    !
    ! ========== メッシュ作成 ==========
    integer, parameter :: nxc = 200
    integer, parameter :: nyc = 100
    integer, parameter :: nxd = nxc + 1
    integer, parameter :: nyd = nyc + 1
    integer, parameter :: nx = nxc*2
    integer, parameter :: ny = nyc*2
    real(real64), parameter :: lx = 0.2d0   ![m]
    real(real64), parameter :: ly = 0.1d0   ![m]
    real(real64), parameter :: dx = 0.001d0 ![m]
    real(real64), parameter :: dy = 0.001d0 ![m]
    real(real64), parameter :: uin = 0.01d0  ![m/s]
    real(real64), parameter :: lt = 5d0

    ! ========== 流体 ==========
    real(real64), parameter :: dens = 1.0e+3       ![kg/m^3]
    real(real64), parameter :: viscosity = 1.0e-3  ![Pa・s]
    real(real64), parameter :: kinematic_viscosity = viscosity/dens ![m^2/s]

    real(real64), parameter :: K  = 0.0005d0                        ! 透水係数（K小⇒Da小⇒u_aux減）
    real(real64), parameter :: epsilon = 0.8d0                      ! ポロシティ(気孔率)
    real(real64), parameter :: Re = uin*lx/kinematic_viscosity      ! レイノルズ数
    ! real(real64), parameter :: Da = K/(lx**2)                     ! ダルシー数
    real(real64), parameter :: Da = K*kinematic_viscosity/lx**2     ! ダルシー数
    real(real64), parameter :: Forchheimer_coef = 1.5e+4            ! Forchheimer
    real(real64), parameter :: semi_implicit_prm = 0.0d0

    real(real64), parameter :: vof_min = 1d-12
    real(real64), parameter :: frate = 0.3d0

    integer, parameter :: ic_porous_start = nint(nxc*0.1d0)
    integer, parameter :: ic_porous_end   = nint(nxc*0.5d0)

    real(real64) :: dt
    real(real64) :: arufa, err_total
    real(real64) :: start_time, end_time
    integer :: i, j, time, nt, outputstep
    integer :: ic, jc
    real(real64) :: x, y


    ! ========== 物理量 ==========
    real(real64), dimension(0:nxd + 1, 0:nyc + 1) :: u
    real(real64), dimension(0:nxc + 1, 0:nyd + 1) :: v
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
    dt = 0.00001d0

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
    nt = 1000
    outputstep = 100

    call setBoundary(u, v, uin)


! メインルーチン
    call output(u, v, p, phi, 0, theta, vof)
    call cpu_time(start_time)
    print *, "Reynolds number = ",Re
    print *, "Darcy    number = ",Da
    print *, "couran   number = ",uin*dt/dx
    print *, "K * Uin         = ",k*uin
    do time = 1, nt
    ! do time = 1, 10
        print *, time,"/",nt,dt*time
        call computeAuxiallyVelocity(u_aux, v_aux, u, v, p, dx, dy, dt, uin, vof, force, eps, aprt)
        call computeDivergenceAuxiallyVelocity(theta, u_aux, v_aux, dx, dy, dt)
        call computePressurePoisson(p, phi, dp, u, v, dx, dy, theta, err_total, arufa, dens, dt, aprt)
        call computeVelocity(u, v, u_aux, v_aux, dt, dx, dy, dens, phi, uin, vof, aprt)
        if (mod(time, outputstep) == 0) then
            call output(u, v, p, phi, time, theta, vof)
        elseif (mod(time, 10) == 0 .and. time < 100) then
            call output(u, v, p, phi, time, theta, vof)
        end if
    end do
    call cpu_time(end_time)
    write (*, *) "elapsed time = ", end_time - start_time

    open (unit=30, file="elapsed.txt")
    write (30, *) end_time - start_time
    close (30)

    open (unit=40, file="velocity.txt")
    ! write (40, *) (u(i,50),i=1,nxd)
    do i=1,nxd
        write(40,*) u(i,50)
    enddo
    close (40)

contains

    subroutine output(u, v, p, phi, time, theta, vof)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(in) :: u(0:, 0:)
        real(real64), intent(in) :: v(0:, 0:)
        real(real64), intent(in) :: p(0:, 0:)
        real(real64), intent(in) :: phi(0:, 0:)
        integer, intent(in) :: time
        real(real64), intent(in) :: theta(:, :)
        real(real64), intent(in) :: vof(0:, 0:)

        integer :: i, j
        real(real64) :: uout, vout
        character*128 filename1, filename2, filename3

        write (filename1, '("pres/pres",i6.6,".txt")') time
        ! write (filename2, '("vof/vof",i6.6,".txt")') time
        ! write (filename3, '("theta/theta",i5.5,".txt")') time
        open (unit=10, file=filename1)
        ! open (unit=20, file=filename2)
        ! open (unit=30, file=filename3)
        do j = 1, nyc
        do i = 1, nxc
            uout = (u(i, j) + u(i + 1, j))/2d0
            vout = (v(i, j) + v(i, j + 1))/2d0
            write (10, '(2i8,3f23.14)') i, j, p(i, j), uout, vout
            ! write (20, '(2i8,3f23.14)') i, j, vof(i, j), uout, vout
            ! write (20, '(2i8,3f23.14)') i, j, vof(i, j), u(i,j),v(i,j)
            ! write (30, '(2i8,3f23.14)') i, j, theta(i, j), uout, vout
        end do
        end do
        close (10)
        ! close (20)
        ! close (30)

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
                         - dt*(p(ic, jc) - p(ic - 1, jc))/dy &
                         - dt*(1d0 - semi_implicit_prm)/(Re*Da)*eps(ic,jc)*u(i,jc)

            u_aux(i, jc) = u_aux(i, jc)/(1d0 + semi_implicit_prm*dt/(Re*Da)*eps(ic,jc))
            if(i==30.and.j==30)then
                print *, u(i,jc),u_aux(i, jc)
            endif
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
                        - dt*(1d0 - semi_implicit_prm)/(Re*Da)*eps(ic,jc)*v(ic,j)

            v_aux(ic, j) = v_aux(ic, j)/(1d0 + semi_implicit_prm*dt/(Re*Da)*eps(ic,jc))
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
            theta(ic, jc) = (1d0 + semi_implicit_prm*dt/(Re*Da)*eps(ic,jc)) &
                           *(u_aux(i + 1, jc) - u_aux(i, jc))/dx + (v_aux(ic, j + 1) - v_aux(ic, j))/dy
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
        ! p(:, :) = p(:, :) - original_pressure
        p(:, :) = p(:, :) + phi(:, :) 

    end subroutine

    subroutine computeVelocity(u, v, u_aux, v_aux, dt, dx, dy, dens, phi, uin, vof, aprt)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64), intent(inout) :: u(0:, 0:)
        real(real64), intent(inout) :: v(0:, 0:)
        real(real64), intent(inout) :: phi(0:, 0:)
        real(real64), intent(in) :: u_aux(0:, 0:)
        real(real64), intent(in) :: v_aux(0:, 0:)
        real(real64), intent(in) :: vof(0:, 0:)
        type(velocity),intent(in)::aprt

        real(real64), intent(in) :: dt, dx, dy, dens, uin
        real(real64)::H_f

        do jc = 1, nyc
        do i = 1, nxd
            j = jc
            ic = i
            ! H_f = aprt%x(i,jc)
            H_f = 1d0
            u(i, jc) = u_aux(i, jc) - dt*H_f*(phi(ic, jc) - phi(ic - 1, jc))/dx &
                       /(1d0 + semi_implicit_prm*dt/(Re*Da)*eps(ic,jc))
        end do
        end do

        do j = 1, nyd
        do ic = 1, nxc
            jc = j
            i = ic
            ! H_f = aprt%y(ic,j)
            H_f = 1d0
            v(ic, j) = v_aux(ic, j) - dt*H_f*(phi(ic, jc) - phi(ic, jc - 1))/dy &
                       /(1d0 + semi_implicit_prm*dt/(Re*Da)*eps(ic,jc))
        end do
        end do

        call setBoundary(u, v, uin)

    end subroutine

end program
