*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     k   integer vector of dimension m by 1
*     n   integer
*     m   integer
*     g   integer (# of factor levels)
*     y   double matrix of dimension n by m (all 0)

      subroutine ordker(x, k, n, m, g, y)

      integer n, m, i, j, x(n), k(m), g
      double precision y(n,m), a, b, c

      c = (g-1.0)*(2.0*g-1.0)/(6.0*g)
      do j = 1, m
        a = k(j)*(k(j)-1.0)
          do i = 1, n
          b = x(i)*(x(i)-1.0)
          y(i,j) = 1.0 - MAX(x(i),k(j)) + c + (a+b)/(2.0*g)
        end do
      end do

      end
