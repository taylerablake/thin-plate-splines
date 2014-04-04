*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <nhelwig2@illinois.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     n   integer
*     y   double matrix of dimension n by n (all 0)

      subroutine cubkerzsym(x, n, y)

      integer n, i, j
      double precision x(n), y(n,n), a, b

      y(1,1) = (x(1)**3)/3
      do i = 2, n
        y(i,i) = (x(i)**3)/3
        do j = 1, (i-1)
          a = min(x(i),x(j))
          b = max(x(i),x(j))
          y(i,j) = (a**2)*((3*b)-a)/6
          y(j,i) = y(i,j)
        end do
      end do

      end
