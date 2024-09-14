program main
    implicit none
    integer :: i, j

    open(1, file = 'depth.txt')
    write(1, *) 800, 40
    do i = 1, 800
        if (i .le. 400) then
            write(1, '(40ES20.5)') (10.0, j = 1, 40)
        else
            write(1, '(40ES20.5)') (10-(i-400d0)/40, j = 1, 40)
        end if
    end do
    close(1)
end program main