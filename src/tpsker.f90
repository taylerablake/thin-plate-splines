
SUBROUTINE tpsker(x, k, n, p, q, y)

    INTEGER n
    INTEGER p
    INTEGER q
    INTEGER i
    INTEGER j
    INTEGER d
    DOUBLE PRECISION x(n,p)
    DOUBLE PRECISION k(q,p)
    DOUBLE PRECISION y(n,q)
    DOUBLE PRECISION r

        d = 4 - p

        if(mod(p,2).eq.0) then
            do j = 1, q
                do i = 1, n
                    r = sqrt(sum((x(i,:)-k(j,:))**2))
                    if (r>0) y(i,j) = (r**d)*log(r)
                end do
            end do
        else
            do j = 1, q
                do i = 1, n
                    y(i,j) = sqrt(sum((x(i,:)-k(j,:))**2))**d
                end do
            end do
        end if

        end SUBROUTINE
