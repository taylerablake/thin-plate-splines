*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   double vector of dimension n by 1
*     n   integer
*     y   double matrix of dimension n by n (all 0)

      subroutine linkersym(x, n, y)

      integer n, i, j
      double precision x(n), y(n,n)

      x = x - 0.5
      y(1,1) =  x(1)**2 + (0.25 - 1./12)/2
      do i = 2, n
        y(i,i) =  x(i)**2 + (0.25 - 1./12)/2
        do j = 1, (i-1)
          y(i,j) = x(i)*x(j) + ((abs(x(i)-x(j))-0.5)**2 - 1./12)/2
          y(j,i) = y(i,j)
        end do
      end do

      end
