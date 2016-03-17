*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     x   integer vector of dimension n by 1
*     k   integer vector of dimension m by 1
*     n   integer
*     m   integer
*     g   integer (# of factor levels)
*     y   double matrix of dimension n by m (all 0)

      subroutine ordker(x, k, n, m, g, y)

      integer n, m, i, j, x(n), k(m), g
      double precision y(n,m), ID(g,g), D(n,g), E(g,m), rj


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

      E = 0
      do j = 1, g
        do i = 1, m
          if (k(i).EQ.j) then
            E(j,i) = 1 - 1./g
          else
            E(j,i) = - 1./g
          end if
        end do
      end do

      y = MATMUL(D,MATMUL(ID,MATMUL(TRANSPOSE(ID),E)))

      end
