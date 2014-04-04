*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <nhelwig2@illinois.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     k   double vector of dimension m by 1
*     n   integer
*     m   integer
*     y   double matrix of dimension n by m (all 0)

      subroutine cubkerz(x, k, n, m, y)

      integer n, m, i, j
      double precision x(n), k(m), y(n,m), a, b

      do j = 1, m
        do i = 1, n
          a = min(x(i),k(j))
          b = max(x(i),k(j))
          y(i,j) = (a**2)*((3*b)-a)/6
        end do
      end do

      end
