C Output from Public domain Ratfor, version 1.0
C call from sspreg0 with .Fortran("dsidr0",
C                  as.integer(code),
C                  swk=as.double(s), as.integer(nobs),
C                  as.integer(nobs), as.integer(nnull),
C                  as.double(y),
C                  qwk=as.double(q), as.integer(nobs),
C                  as.double(0), as.integer(0), double(2),
C                  nlambda=double(1), score=double(1), varht=as.double(varht),
C                  c=double(nobs), d=double(nnull),
C                  qraux=double(nnull), jpvt=integer(nnull),
C                  double(3*nobs),
C                  info=integer(1)

    subroutine dsidr (vmu, s, lds, nobs, nnull, y, q, ldq, tol, job, l
     *imnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
        character vmu
        integer lds, nobs, nnull, ldq, job, jpvt(*), info
        double precision s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, 
        *score(*), varht, c(*), d(*), qraux(*), wk(*)
      
        info = 0
        if( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds .or. nobs
        * .gt. ldq )then
            info = -1
            return
        endif
        
        if( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' )then
            info = -3
            return
            endif
      
        C dstup returns ___ ___ ___ 
      
        call dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nobs, 1, 
            *info, wk)
        if( info .ne. 0 )then
            return
        endif
        
        C dcore returns ___ ___ ___
        
        call dcore (vmu, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht, 
        *score, varht, info, wk, wk(2*nobs+1))
        if( info .ne. 0 )then
            return
        endif
        
        
            
        call dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c,
        * d, info, wk)
        return
        end
