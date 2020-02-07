c     nilsson model
      implicit real*8(a-h,o-z)
      parameter (km=100,kw=50)
      dimension ham(km,km,kw), eval(km,kw), evec(km,km,kw)
      dimension ial(km),ialam(km),iasig(km)
      dimension spe(100,1000)
      dimension rkappa(0:11),rmu(0:11)
      data (rkappa(i),i=0,7) /8*0.05D0/
      data (rmu(i),i=0,5) /0.0D0,0.0D0,0.0D0,0.35D0,0.625D0,0.630D0 /
      data (rmu(i),i=6,7) /0.448D0,0.434D0 /

      read(*,*)iz,in
      rmass=dfloat(iz+in)
      icol=0
      do n=12,12,2  !all quantum numbers multiplied by 2
        do iom=-n-1,n+1,2
          nbas=0
          do l=n,n,-4
            do lam=-l,l,2
              isig=iom-lam
              if(abs(isig).ne.1)cycle
              nbas=nbas+1
              ial(nbas)=l
              ialam(nbas)=lam
              iasig(nbas)=isig
            enddo
          enddo
          idef=0
          do delta=0.3d0,0.301d0,0.01d0
            fdel=((1.d0+(2.d0/3.d0)*delta)**2.d0
     &            *(1.d0-(4.d0/3.d0)*delta))**(-1.d0/6.d0)
            hw00=41.d0*rmass**(-1.d0/3.d0)
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

                omc=0.d0
                do k=1,13
                  ham(i,j,k)=h00+hdel+c*hls+d*hl2-omc*hw00*iom/2.d0
                  ham(j,i,k)=ham(i,j,k)
                  omc=omc+0.1d0
                enddo
              enddo
            enddo
            omc=0.d0
            do k=1,13
              call diag(ham(:,:,k),km,nbas,eval(:,k),evec(:,:,k))
              spe(k,icol+1)=eval(1,k)
              omc=omc+0.1d0
            enddo
          enddo !end of deformation loop
          icol=icol+1
        enddo
      enddo
      omc=0.d0
      do k=1,13
        write(*,'(2x,1000f12.6)')delta,omc,(spe(k,j),j=1,icol)
        omc=omc+0.1d0
      enddo
      end
