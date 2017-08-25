*     bigsplines (www.r-project.org)
*     Nathaniel E. Helwig <helwig@umn.edu>
*     ***** inputs *****
*     y   double vector of length n
*     g   integer vector of length n
*     n   integer (length of y, x's, and g)
*     m   integer (number of possible unique points)
*     s   double vector of length m (all 0.0)
*     w   integer vector of length m (all 0)

       subroutine sumfreq(y, g, n, m, s, w)

       integer n, m, j, g(n), w(m)
       double precision y(n), s(m)

       do j = 1, n
         s(g(j)) = s(g(j)) + y(j)
         w(g(j)) = w(g(j)) + 1
       end do

       end
