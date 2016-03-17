*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     n   integer
*     g   integer (# of factor levels)
*     y   double matrix of dimension n by n (all 0)

      subroutine ordkersym(x, n, g, y)

      integer n, i, j, x(n), g
      double precision y(n,n), ID(g,g), D(n,g), rj

      ID = 0
      do j = 1, (g-1)
        rj = real(j)
        do i = 1, j
          ID(i,j) = (g-rj)/g
        end do
        do i = (j+1), g
          ID(i,j) = -rj/g
        end do
      end do
      do i = 1, g
        ID(i,g) = 1./g
      end do

      D = 0
      do j = 1, g
        do i = 1, n
          if (x(i).EQ.j) then
            D(i,j) = 1 - 1./g
          else
            D(i,j) = - 1./g
          end if
        end do
      end do

      D = MATMUL(D,ID)
      y = MATMUL(D, TRANSPOSE(D))

      end
