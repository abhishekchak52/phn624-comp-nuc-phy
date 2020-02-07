      real*8 function cgc(j1,j2,j,m1,m2)
      implicit real*8(a-h,o-z)
      common/factoriall/factl(0:200),ifactl
      if (ifactl.eq.0) call factlog(200)
      m=m1+m2
      ia=(j1+j2-j)/2
      ib=(j1-m1)/2
      ic=(j2+m2)/2
      id=(j-j2+m1)/2
      ie=(j-j1-m2)/2

      t1=dlog(dfloat(j+1))+factl((j1+j2-j)/2)+factl((j+j1-j2)/2)+
     &         factl((j+j2-j1)/2)-factl((j+j1+j2+2)/2)
      t2=factl((j1+m1)/2)+factl((j1-m1)/2)+factl((j2+m2)/2)+
     &         factl((j2-m2)/2)+factl((j+m)/2)+factl((j-m)/2)

      imin=min(id,ie)
      if(imin.lt.0)then
        imin=iabs(imin)
      else
        imin=0
      endif

      t3=0.0d0
      do i=imin,min(ia,ib,ic)
        term=factl(i)+factl(ia-i)+factl(ib-i)
     &        +factl(ic-i)+factl(id+i)+factl(ie+i)
        t3=t3+(-1)**i*dexp(-term)
      enddo
      cgc=dexp(0.5d0*(t1+t2))*t3
      end
