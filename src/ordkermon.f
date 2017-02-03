*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     k   integer vector of dimension m by 1
*     n   integer
*     m   integer
*     y   double matrix of dimension n by m (all 0)

      subroutine ordkermon(x, k, n, m, y)

      integer n, m, i, j, x(n), k(m)
      double precision y(n,m), a, b

      do j = 1, m-1
        a = (1.0 * j) / (1.0 * m)
          do i = 1, n
            if (x(i) <= k(j)) then
              b = 1.0
            else
              b = 0.0
            end if
            y(i,j) = a - b
        end do
      end do

      end
