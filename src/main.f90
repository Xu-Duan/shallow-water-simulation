program main
    use helper
    use timestep
    implicit none

    real(dp) :: dx, dy, omega, A
    integer :: nx, ny, it=0, iprogress = 0, ioutput = 0, imethod
    real(dp), dimension(:, :), allocatable :: h, u, v, d

    real(dp) :: t_current = 0d0, t_total, odt, dt, t1, t2
    character(8) :: date
    character(10) :: time

    call cpu_time(t1)
    call date_and_time(date, time)
    print *, '|---------------- START -----------------|'
    print *, '|--------- by Duan Xu from SJTU --------|'
    print *, '|----------- ', date(1:4), '.', date(5:6), '.', date(7:8), ' ----------------|'
    print *, '|------------ ', time(1:2), ':', time(3:4), ':', time(5:6), ' -----------------|'
    print *, 'input file is "input.txt" and "depth.txt" by default'
    print *, '<RETURN> to start calculation'
    read(*, *)

    call get_input(dx, dy, nx, ny, dt, t_total, odt, omega, A, imethod) ! basic input parameter
    allocate(h(0:nx+1, 0:ny+1), u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), d(0:nx+1, 0:ny+1))
    allocate(h_tild(0:nx+1, 0:ny+1), u_tild(0:nx+1, 0:ny+1), v_tild(0:nx+1, 0:ny+1))
    allocate(h_old(0:nx+1, 0:ny+1), u_old(0:nx+1, 0:ny+1), v_old(0:nx+1, 0:ny+1))
    
    call get_init_depth(nx, ny, d) ! get and initialize depth
    call init(nx, ny, h, u, v, d) ! initialization
    call bc(nx, ny, h, u, v, d, t_current, omega, A)   ! implement boundary condition
    call write_output(nx, ny, dt, t_current, odt, t_total, h, u, v, ioutput)

    do while (t_current < t_total)
        t_current = t_current + dt
        h_old = h
        u_old = u
        v_old = v
        if (imethod == 0) then ! Euler method
            call diffcalc(nx, ny, dx, dy, dt, h, u, v, d)
            h = h_old + h
            u = u_old + u
            v = v_old + v
            call bc(nx, ny, h, u, v, d, t_current+dt, omega, A)
        else if (imethod == 1) then ! TODO : RK 3 method
            ! 1st RK step
            call diffcalc(nx, ny, dx, dy, dt, h, u, v, d)
            h = h_old + h/2d0
            u = u_old + u/2d0
            v = v_old + v/2d0
            call bc(nx, ny, h, u, v, d, t_current+dt/2d0, omega, A)

            ! 2nd RK step
            call diffcalc(nx, ny, dx, dy, dt, h, u, v, d)
            h = h + h_old
            u = u + u_old
            v = v + v_old
            call bc(nx, ny, h, u, v, d, t_current+dt, omega, A)
        else  ! RK 4 method
            ! 1st RK step
            call diffcalc(nx, ny, dx, dy, dt, h, u, v, d)
            h = h_old + h / 6d0
            u = u_old + u / 6d0
            v = v_old + v / 6d0
            h_tild = h_old + h / 2d0
            u_tild = u_old + u / 2d0
            v_tild = v_old + v / 2d0
            call bc(nx, ny, h_tild, u_tild, v_tild, d, t_current+dt/2d0, omega, A)
            
            ! 2nd RK step
            call diffcalc(nx, ny, dx, dy, dt, h_tild, u_tild, v_tild, d)
            h = h + h_tild / 3d0
            u = u + u_tild / 3d0
            v = v + v_tild / 3d0
            h_tild = h_old + h_tild / 2d0
            u_tild = u_old + u_tild / 2d0
            v_tild = v_old + v_tild / 2d0
            call bc(nx, ny, h_tild, u_tild, v_tild, d, t_current+dt/2d0, omega, A)

            ! 3rd RK step
            call diffcalc(nx, ny, dx, dy, dt, h_tild, u_tild, v_tild, d)
            h = h + h_tild / 3d0
            u = u + u_tild / 3d0
            v = v + v_tild / 3d0
            h_tild = h_old + h_tild
            u_tild = u_old + u_tild
            v_tild = v_old + v_tild
            call bc(nx, ny, h_tild, u_tild, v_tild, d, t_current+dt, omega, A)

            ! 4th RK step
            call diffcalc(nx, ny, dx, dy, dt, h_tild, u_tild, v_tild, d)
            h = h + h_tild / 6d0
            u = u + u_tild / 6d0
            v = v + v_tild / 6d0
            call bc(nx, ny, h, u, v, d, t_current+dt, omega, A)
        end if

        it = it + 1
        call show_progress(t_current, t_total, iprogress)
        call write_output(nx, ny, dt, t_current, odt, t_total, h, u, v, ioutput)
    end do

    call cpu_time(t2)
    print *, '|-- calculation completed successfully --|'
    print *, 'output file is in directory "output" by default'
    print '(a, F6.2, a)', ' |------ time consumed: ', t2-t1,  's ----------|'
    print *, '<RETURN> to end the program'
    read(*, *)
end program main
