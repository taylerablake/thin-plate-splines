*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     k   double vector of dimension m by 1
*     n   integer
*     m   integer
*     y   double matrix of dimension n by m (all 0)

      subroutine cubker(x, k, n, m, y)

      integer n, m, i, j
      double precision x(n), k(m), y(n,m), a, b, c, d

      x = x - 0.5
      k = k - 0.5
      do j = 1, m
        a = ((k(j)**2)-(1./12))/2
          do i = 1, n
          b = ((x(i)**2)-(1./12))/2
          d = abs(x(i)-k(j)) - 0.5
          c = ((d**4)-((d**2)/2)+(7./240))/24
          y(i,j) =  a*b - c
        end do
      end do

      end
