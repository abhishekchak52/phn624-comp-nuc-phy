      subroutine factlog(n)
      implicit real*8(a-h,o-z)
      common/factoriall/factl(0:200),ifactl
      factl(0)=0.0d0
      do i=1,n
        factl(i)=factl(i-1)+dlog(dfloat(i))
      enddo
      ifactl=1
      end
