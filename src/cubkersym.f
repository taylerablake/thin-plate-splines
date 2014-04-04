*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <nhelwig2@illinois.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     n   integer
*     y   double matrix of dimension n by n (all 0)

      subroutine cubkersym(x, n, y)

      integer n, i, j
      double precision x(n), y(n,n), a, b, c, d, e

      x = x - 0.5
      a = ((x(1)**2)-(1./12))/2
      c = ((0.5**4)-((0.5**2)/2)+(7./240))/24
      y(1,1) =  a**2 - c
      do i = 2, n
        a = ((x(i)**2)-(1./12))/2
        y(i,i) =  a**2 - c
        do j = 1, (i-1)
          b = ((x(j)**2)-(1./12))/2
          d = abs(x(i)-x(j)) - 0.5
          e = ((d**4)-((d**2)/2)+(7./240))/24
          y(i,j) =  a*b - e
          y(j,i) = y(i,j)
        end do
      end do

      end
