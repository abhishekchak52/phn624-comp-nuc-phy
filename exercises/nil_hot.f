c     nilsson model
      implicit real*8(a-h,o-z)
      parameter (km=100)
      dimension ham(km,km), eval(km), evec(km,km)
      dimension ial(km),ialam(km),iasig(km)
      dimension spe(100,364),vk2(364)
      dimension rkappa(0:11),rmu(0:11)
      data (rkappa(i),i=0,7) /8*0.05D0/
      data (rmu(i),i=0,5) /0.0D0,0.0D0,0.0D0,0.35D0,0.625D0,0.630D0 /
      data (rmu(i),i=6,7) /0.448D0,0.434D0 /
      open (4,file='hot.out',status='unknown')

      ip=40
      in=44
      ra=dfloat(ip+in)
      icol=0
      do n=0,14,2  !all quantum numbers multiplied by 2
        do iom=1,n+1,2
          nbas=0
          do l=n,0,-4
            do lam=-l,l,2
              isig=iom-lam
              if(abs(isig).ne.1)cycle
              nbas=nbas+1
              ial(nbas)=l
              ialam(nbas)=lam
              iasig(nbas)=isig
            enddo
          enddo
          delta=0.2d0
            fdel=((1.d0+(2.d0/3.d0)*delta)**2.d0
     &            *(1.d0-(4.d0/3.d0)*delta))**(-1.d0/6.d0)
            hw00=41.d0*ra**(-1.d0/3.d0)
            hw0=hw00*fdel
            c=-2.0D0*hw00*rkappa(n/2)
            d=rmu(n/2)*c/2.0D0
            do i=1,nbas
              l=ial(i)
              lam=ialam(i)
              isig=iasig(i)
              do j=i,nbas
                l1=ial(j)
                lam1=ialam(j)
                isig1=iasig(j)
  
                h00=0.d0
                hl2=0.d0
                hls=0.d0
                hr2=0.d0
                hy20=0.d0
                if(i.eq.j)then  !Diagonal matrix elements
                  h00=dfloat(n+3)/2.d0*hw0
                  hl2=dfloat(l*(l+2))/4.d0
                endif

                if(l1.eq.l)then
                  if((lam1.eq.lam+2).and.(isig1.eq.isig-2))then
                    hls=dsqrt(dfloat((l-lam)*(l+lam+2)))/4.d0
                  elseif((lam1.eq.lam).and.(isig1.eq.isig))then
                    hls=lam*isig/4.d0
                  elseif((lam1.eq.lam-2).and.(isig1.eq.isig+2))then
                    hls=dsqrt(dfloat((l+lam)*(l-lam+2)))/4.d0
                  endif
                  hr2=dfloat(n+3)/2.d0
                elseif(l1.eq.l-4)then
                  hr2=dsqrt(dfloat((n-l+4)*(n+l+2)))/2.d0
                elseif(l1.eq.l+4)then
                  hr2=dsqrt(dfloat((n-l)*(n+l+6)))/2.d0
                endif
                if((abs(hr2).gt.1e-5).and.(lam1.eq.lam))then
                  hy20=dsqrt(dfloat(l+1)/dfloat(l1+1))
     &                 *cgc(l,4,l1,lam,0)*cgc(l,4,l1,0,0)
                endif
                hdel=-delta*hw0*dsqrt(4.d0)/3.d0*hr2*hy20

                ham(i,j)=h00+hdel+c*hls+d*hl2
                ham(j,i)=ham(i,j)
              enddo
            enddo

            call diag(ham,km,nbas,eval,evec)
            idef=1
            do i=1,nbas
              spe(idef,icol+i)=eval(i)
            enddo 
          icol=icol+nbas
        enddo
      enddo
      call sort(spe(idef,:),1000,icol)
      idef=1
      ef0=spe(idef,in/2)
      do t=0.5d0,5.0d0,0.5d0
        call hot(spe(idef,:),44,t,ef0,ef,fe)
        write(4,'(2x,i8,f6.2,2f12.6,2i3)')22222,t,ef,fe,1,1
        ef0=ef
      enddo
      
      end

      subroutine hot(e,n,t,ef0,ef,fe)
      implicit real*8(a-h,o-z)
      dimension e(364)
      common /hot1/ei(364),ti,np
      external chemp
      do i=1,364
        ei(i)=e(i)
      enddo
      np=n
      ti=t
      ef=secant(chemp,ef0)
      fe=0.0d0
      do i=1,364
        if (e(i).lt.1e-2) exit
        rn=fn_nt(e(i),ef,t)
        fn_s=0.0d0
        if((rn.ge.1d-45).and.(rn.lt.1.0))
     &    sn=-(rn*dlog(rn)+(1.0d0-rn)*dlog(1.0d0-rn))
        fe=fe+rn*e(i)-t*sn
        write(4,'(2x,f6.2,i4,f12.6,2e15.6,f12.6)')t,i,e(i),rn,sn,fe
        if (rn.lt.1e-8) exit
      enddo
      end

      real*8 function chemp(ef)
      implicit real*8(a-h,o-z)
      common /hot1/ei(364),ti,np
      xn=0
      do i=1,364
        if (ei(i).lt.1e-2) exit
        x=fn_nt(ei(i),ef,ti)
        xn=xn+x
        if (x.lt.1e-8) exit
      enddo
      chemp=dfloat(np)-2.0d0*xn
      if (dabs(chemp).gt.100.d0) stop
      print *,ti,2*xn,np
      end

      real*8 function fn_nt(e,ef,t)
      implicit real*8(a-h,o-z)
      if ((e-ef)/t.gt.700.0) then
        fn_nt=0.0d0
      else
        fn_nt=1.0d0/(1.0d0+dexp((e-ef)/t))
      end if
      end
       

      real*8 function fn_s(e,ef,t)
      implicit real*8(a-h,o-z)
      real*8 nt
      fn_s=0.0d0
      nt=fn_nt(e,ef,t)
      if(nt.lt.1d-45) return
      if(nt.ge.1.0) return
      fn_s=-(nt*dlog(nt)+(1.0d0-nt)*dlog(1.0d0-nt))
      end
