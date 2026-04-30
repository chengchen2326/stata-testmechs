      SUBROUTINE XERBLA( SRNAME, INFO )
*     Stub for BLAS error handler.
*     Real implementation prints error and stops, but we never expect
*     dqrdc2 to call this since we don't trigger error paths.
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
      WRITE(*,*) 'XERBLA called from ', SRNAME, ' with info=', INFO
      RETURN
      END
