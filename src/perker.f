*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     k   double vector of dimension m by 1
*     n   integer
*     m   integer
*     y   double matrix of dimension n by m (all 0)

      subroutine perker(x, k, n, m, y)

      integer n, m, i, j
      double precision x(n), k(m), y(n,m), a, b, c

      x = x - 0.5
      k = k - 0.5
      do j = 1, m
        do i = 1, n
          c = abs(x(i)-k(j)) - 0.5
          y(i,j) = -((c**4)-((c**2)/2)+(7./240))/24
        end do
      end do

      end
