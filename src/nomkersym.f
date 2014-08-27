*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     n   integer
*     g   double scalar (inverse of # of factor levels)
*     y   double matrix of dimension n by n (all 0)

      subroutine nomkersym(x, n, g, y)

      integer n, i, j, x(n)
      double precision y(n,n), g

	  y(1,1) = 1 - g
      do i = 2, n
	    y(i,i) = 1 - g 
		do j = 1, (i-1)
          if (x(i).EQ.x(j)) then 
		    y(i,j) = 1 - g
		  else 
		    y(i,j) = - g
          end if
		  y(j,i)=y(i,j)
        end do
      end do

      end
