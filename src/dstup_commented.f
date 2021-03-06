C Output from Public domain Ratfor, version 1.0
    subroutine dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldqr, ld
                  *qc, nq, info, work)

        C In parameters
        integer lds, nobs, nnull, jpvt(*), ldqr, ldqc, nq,
        double precision y(*), qraux(*), q(ldqr,ldqc,*),

        C Out parameters
        integer info
        
        C In/Out parameters
        double precision s(lds,*)
        double precision work(*)
        
        C local subroutine variables
        double precision dum
        integer j

        
    info = 0
    if( nobs .lt. 1 .or. nobs .gt. lds .or. nobs .gt. ldqr .or. nobs .
        *gt. ldqc )then
        info = -1
        return
    endif
    
    j=1
    
    23002 
    if(.not.(j.le.nnull))goto 23004
    jpvt(j) = 0

    23003
    j=j+1
    goto 23002
    
    23004 
    continue
    
    C compute the QR decomposition of S 
    call dqrdc (s, lds, nobs, nnull, qraux, jpvt, work, 1)
    
    C given the QR decomposition of S, compute Y = Q'Y, work = minimizer of (Y - QR beta)
    call dqrsl (s, lds, nobs, nnull, qraux, y, dum, y, work, dum, dum,
            * 01100, info)
    if( info .ne. 0 )then
        return
    endif
    
    j=1

    23007
    if(.not.(j.le.nq))goto 23009
    C compute Q2' K Q2, Q1' K Q2  --> replace K with K + M lambda
    call dqrslm (s, lds, nobs, nnull, qraux, q(1,1,j), ldqr, 0, info, 
             *work)
    
    23008 
    j=j+1
    goto 23007
    
    23009 continue
    
    return
    end