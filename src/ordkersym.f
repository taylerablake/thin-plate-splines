*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     n   integer
*     g   integer (# of factor levels)
*     y   double matrix of dimension n by n (all 0)

      subroutine ordkersym(x, n, g, y)

      integer n, i, j, x(n), g
      double precision y(n,n), a, b, c

      c = (g-1.0)*(2.0*g-1.0)/(6.0*g)
      y(1,1) = 1.0 - x(1) + c + x(1)*(x(1)-1.0)/g
      do i = 2, n
        a = x(i)*(x(i)-1.0)
        y(i,i) =  1.0 - x(i) + c + a/g
        do j = 1, (i-1)
          b = x(j)*(x(j)-1.0)
          y(i,j) =  1.0 - MAX(x(i),x(j)) + c + (a+b)/(2.0*g)
          y(j,i) = y(i,j)
        end do
      end do

      end
