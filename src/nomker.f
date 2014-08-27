*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     k   integer vector of dimension m by 1
*     n   integer
*     m   integer
*     g   double scalar (inverse of # of factor levels)
*     y   double matrix of dimension n by m (all 0)

      subroutine nomker(x, k, n, m, g, y)

      integer n, m, i, j, x(n), k(m)
      double precision y(n,m), g

      do j = 1, m
		do i = 1, n
          if (x(i).EQ.k(j)) then 
		    y(i,j) = 1 - g
		  else 
		    y(i,j) = - g
          end if
        end do
      end do

      end
