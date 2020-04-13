c     nilsson model
      implicit real*8(a-h,o-z)
      parameter (km=100)
      dimension ham(km,km), eval(km), evec(km,km)
      dimension ial(km),ialam(km),iasig(km)
      dimension delta_list(60), Epc_list(60, 2)
      dimension pes_list(60), pgp_list(60), pgn_list(60)
      dimension spe(100,364),vk2(364)
      dimension rkappa(0:11),rmu(0:11)
      data (rkappa(i),i=0,7) /8*0.05D0/
      data (rmu(i),i=0,5) /0.0D0,0.0D0,0.0D0,0.35D0,0.625D0,0.630D0 /
      data (rmu(i),i=6,7) /0.448D0,0.434D0 /
      open (4,file='bcs_data.txt',status='unknown')

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
          idef = 0
          do delta=-0.3d0, 0.301d0, 0.01d0

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
            idef=idef+1
            delta_list(idef) = delta
            do i=1,nbas
              spe(idef,icol+i)=eval(i)
            enddo
          enddo
          icol=icol+nbas
        enddo
      enddo
      open(unit=45, file='bcs_spe.txt')
      do idef=1,60
        call sort(spe(idef,:),1000,icol)
        sum=sumene(iz,spe(idef,:),1000)
        sum=sum+sumene(in,spe(idef,:),1000)
        sum=sum-sum0
        pes_list(idef)=sum
        if (idef .eq. 51)
     &   write(45, *) delta_list(idef), spe(idef,:)
      enddo
      close(45)
      idef=1
c      write(*,'(2x,1000f12.6)')delta,(spe(idef,j),j=1,icol)
      g0=19.2d0
      g1= 7.4d0
      gp=(g0+g1*dfloat(in-ip)/ra)/ra
      gn=(g0-g1*dfloat(in-ip)/ra)/ra
      write(*,*) gn, gp
c$$$      Delta vs G
      open(unit=25, file="delta_vs_g.txt")
      do g_mul=1.5d0, 0.0d0, -0.001d0
         g_mod = g_mul*gn
         call bcs(in,g_mod,spe(idef,:),bcsd,bcsl,epc,vk2)
         if (bcsd .eq. 0.0d0) exit
         write(25,*) g_mul, g_mod, bcsd
      end do
      close(25)
      open(unit=99, file="bcs_delta.txt")
      do idef=1,60
         call bcs(in,gn,spe(idef,:),bcsd,bcsl,epc,vk2)
         pgn_list(idef) = bcsd
         Epc_list(idef, 1) = epc
         call bcs(ip,gp,spe(idef,:),bcsd,bcsl,epc,vk2)
         pgp_list(idef) = bcsd
         Epc_list(idef, 2) = epc
      enddo
      do idef=1,60
         write(99, '(6f12.6)') delta_list(idef),
     &      pes_list(idef)-pes_list(31),
     &      Epc_list(idef,1)-Epc_list(31, 1),
     &      Epc_list(idef,2)-Epc_list(31, 2),
     &    pgp_list(idef), pgn_list(idef)
      enddo
      close(99)
c$$$      call bcs(in,gn,spe(idef,:),bcsd,bcsl,epc,vk2)
c$$$      write(*,'(2x,4f12.6)')gn,bcsd,bcsl,epc
c$$$      do i=1,icol
c$$$        esp=spe(idef,i)
c$$$        eqp=dsqrt((esp-bcsl)**2+bcsd**2)
c$$$        v2=vk2(i)
c$$$        uv=dsqrt(1-v2)*dsqrt(v2)
c$$$        write(4,'(2x,4f12.6)')esp,eqp,v2,uv
c$$$      enddo
      close(4)
      end

      function sumene(ip,spe,ndim)
      implicit real*8(a-h,o-z)
      dimension spe(ndim)
      sumene=0.d0
      do i=1,ip/2
        sumene=sumene+spe(i)
      end do
      if(mod(ip,2).ne.0)sumene=sumene+spe(ip/2+1)
      end
