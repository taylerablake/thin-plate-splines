C Output from Public domain Ratfor, version 1.0
      subroutine dgold (vmu, q, ldq, n, z, low, upp, nlaht, score, varht, info, twk, work)
C Parameters
      character vmu
C In
      integer ldq, n
      double precision q(ldq,*), z(*), 
C Out
      integer info
      double precision nlaht
C In/Out
      double precision low, upp, 
C TBD
      double precision score, varht, twk(2,*), work(*)
C Subroutine local variables
      double precision ratio, mlo, mup, tmpl, tmpu
      
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      info = 0
      if( upp .lt. low ) then
            mlo = low
            low = upp
            upp = mlo
      endif
C is vmu valid? -3 errorcode
      if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
            info = -3
            return
      endif
C is n too small or too large?
      if( n .lt. 1 .or. n .gt. ldq )then
            info = -1
            return
      endif
C set mlo to be the golden section on the interval [low, upp]
      mlo = upp - ratio * (upp - low)
C set up stuff (TODO Look into it)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)

      twk(1,1) = 10.d0**mlo
      
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      
C error check
      if( info .ne. 0 )then
            info = -2
            return
      endif
C set mup to be the upper golden section      
      mup = low + ratio * (upp - low)
      
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)

      twk(1,1) = 10.d0**mup

      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)

      if( info .ne. 0 )then
            info = -2
            return
      endif
C While(TRUE) do this:
23010 continue
C if mu is very close to mlow we stop since we found our minimizer
      if( mup - mlo .lt. 1.d-7 )then
C break out of the loop
            goto 23012
      endif
      if( tmpl .lt. tmpu )then
            upp = mup
            mup = mlo
            tmpu = tmpl
            mlo = upp - ratio * (upp - low)
            
            call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
            call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
            call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
            twk(1,1) = 10.d0**mlo
            call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
            
            if( info .ne. 0 ) then
                  info = -2
                  return
            endif
      else
            low = mlo
            mlo = mup
            tmpl = tmpu
            mup = low + ratio * (upp - low)
            
            call dset (n, 10.d0 ** (mup), twk(2,1), 2)
            call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
            call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
            twk(1,1) = 10.d0**mup
            call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
            
            if( info .ne. 0 )then
                  info = -2
                  return
            endif
      endif
23011 goto 23010

C we are done, clean up
23012 continue
      nlaht = ( mup + mlo ) / 2.d0
      
      call dset (n, 10.d0 ** (nlaht), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**nlaht
      call dtrev (vmu, twk, 2, n, z, score, varht, info, work)
      
      if( info .ne. 0 )then
            info = -2
            return
      endif
      return
      end
