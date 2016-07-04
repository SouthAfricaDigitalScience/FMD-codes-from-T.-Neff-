      subroutine zdet(a,lda,n,ipvt,determinant)
      integer lda,n,ipvt(*)
      complex*16 a(lda,*), determinant

c     straightforward determinant calculation
c     using LU decomposition by zgetrf

      complex*16 det
      integer i

      det = 1.0d0

      do i = 1, n
        if (ipvt(i) .ne. i) det = -det
        det = a(i,i)*det
      enddo

      determinant = det

      return
      end
