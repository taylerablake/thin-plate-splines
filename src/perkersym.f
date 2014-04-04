*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <nhelwig2@illinois.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     n   integer
*     y   double matrix of dimension n by n (all 0)

      subroutine perkersym(x, n, y)

      integer n, i, j
      double precision x(n), y(n,n), a, b

      x = x - 0.5
      a = -((0.5**4)-((0.5**2)/2)+(7./240))/24
      y(1,1) = a
      do i = 2, n
        y(i,i) = a
        do j = 1, (i-1)
          b = abs(x(i)-x(j)) - 0.5
          y(i,j) = -((b**4)-((b**2)/2)+(7./240))/24
          y(j,i) = y(i,j)
        end do
      end do

      end
