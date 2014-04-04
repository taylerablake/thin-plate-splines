*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <nhelwig2@illinois.edu>
*     ***** inputs *****
*     x   double matrix of dimension n by p
*     k   double matrix of dimension q by p
*     n   integer
*     p   integer
*     q   integer
*     y   double matrix of dimension n by q (all 0)

      subroutine tpsker(x, k, n, p, q, y)

      integer n, p, q, i, j, d
      double precision x(n,p), k(q,p), y(n,q), c

      d = 4 - p

      if(mod(p,2).eq.0) then
        do j = 1, q
   	      do i = 1, n
            c = sqrt(sum((x(i,:)-k(j,:))**2))
            if (c>0) y(i,j) = (c**d)*log(c)
		  end do
	    end do
      else
        do j = 1, q
          do i = 1, n
            y(i,j) = sqrt(sum((x(i,:)-k(j,:))**2))**d
          end do
        end do
      end if

      end
