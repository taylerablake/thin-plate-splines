*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <nhelwig2@illinois.edu>
*     ***** inputs *****
*     x   double matrix of dimension n by p
*     n   integer
*     p   integer
*     y   double matrix of dimension n by n (all 0)


      subroutine tpskersym(x, n, p, y)

      integer n, p, i, j, d
      double precision x(n,p), y(n,n), c

      d = 4 - p
     
      if(mod(p,2).eq.0) then
        do i = 2, n
	      do j = 1, (i-1)
            c = sqrt(sum((x(i,:)-x(j,:))**2))
		    if (c>0) y(i,j) = (c**d)*log(c)
            y(j,i) = y(i,j)
		  end do
	    end do
      else
        do i = 2, n
          do j = 1, (i-1)
            y(i,j) = sqrt(sum((x(i,:)-x(j,:))**2))**d
            y(j,i) = y(i,j)
          end do
        end do
      end if

      end
