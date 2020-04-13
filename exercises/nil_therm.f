c     nilsson model
      implicit real*8(a-h,o-z)
      parameter (km=100)
      dimension ham(km,km), eval(km), evec(km,km)
      dimension ial(km),ialam(km),iasig(km)
      dimension delta_list(60), fes_list(60, 6), temp_list(6)
      dimension spe(100,364),vk2(364)
      external secant, chem_pota, fes_calc
      dimension rkappa(0:11),rmu(0:11)
      common /param/ ei(364), temp, np
      data (rkappa(i),i=0,7) /8*0.05D0/
      data (rmu(i),i=0,5) /0.0D0,0.0D0,0.0D0,0.35D0,0.625D0,0.630D0 /
      data (rmu(i),i=6,7) /0.448D0,0.434D0 /
      data (temp_list(i), i=1,6) /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/
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

      open(unit=21, file='ef_delta.txt')
      open(unit=75, file="s_data.txt")
      do idef=1,60
         do itemp=1,6
          call sort(spe(idef,:),1000,icol)
          temp = (itemp-1)
          fes_list(idef, itemp) = fes_calc(ip, in, spe(idef,:),temp)
        enddo
      enddo
      close(21)
      open(unit=99, file="pes_data.txt")
      do itemp=1,6
         fe0 = fes_list(31, itemp)
         do idef=2,60
            fes_list(idef,itemp)=fes_list(idef,itemp)-fe0
c$$$            write(99, '(2f12.6)') delta_list(idef), fes_list(idef, itemp)
         enddo
      enddo
      write(99, '(A,6f5.2)') "Delta", temp_list(:)
      do idef=2,60
         write(99, '(7f12.6)') delta_list(idef), fes_list(idef, :)
      end do
      close(99)
      close(75)

c$$$      For neutrons

c$$$      do i =1,364
c$$$         ei(i) = spe(idef,i)
c$$$      end do
c$$$
c$$$      np=44
c$$$      ef0 = ei(np/2)
c$$$      write(*,*) ef0
c$$$      open(unit=45, file='hot_nuc.txt')
c$$$      do i=1,364
c$$$         if (ei(i) .eq. 0) exit
c$$$         fi = fermi_dist(ei(i), ef_root, temp)
c$$$         write(45,*) ei(i), fi
c$$$      end do
c$$$      close(45)
c$$$      open(unit=54, file='chempot_temp.txt')
c$$$      do temp = 0.05d0, 5.0d0, 0.05d0
c$$$         ef_root = secant(chem_pot, ef0)
c$$$         write(54, *) temp, ef_root
c$$$      end do
c$$$      close(54)
      end

      real*8 function fermi_dist(e, ef, T)
        implicit real*8(a-h,o-z)
        exp_arg = (e-ef)/T

        if (exp_arg .gt. 700.0d0) then
          fermi_dist =0.0d0
        else
          fermi_dist = 1.0d0/(1.0d0+dexp(exp_arg))
        end if
      end function

      real*8 function chem_pot(ef)
        implicit real*8(a-h,o-z)
        common /param/ ei(364), temp, np
        sumni=0.0d0
        do i=1,364
           if (ei(i) .lt. 1e-5) exit
           fni = fermi_dist(ei(i), ef, temp)
           sumni = sumni+fni
           if(fni .lt. 1e-8) exit
        end do
        chem_pot = dfloat(np) - 2.0d0*sumni
c$$$        write(*, *) temp, np, ef, 2.0d0*sumni
      end function

      real*8 function fes_calc(inp, inn, spec, t)
        implicit real*8(a-h,o-z)
        common /param/ ei(364), temp, np
        dimension spec(1000)
        external sumene, chem_pot, fermi_dist, secant
        sum = 0.0d0
        temp=t
        if (temp .eq. 0.0d0) then
          sum=sumene(inp,spec,1000)
          sum=sum+sumene(inn,spec,1000)
        else
c$$$           protons
           do i=1,364
              ei(i)=spec(i)
           end do
           np = inp
           ef0 = ei(np/2)
           efp = secant(chem_pot, ef0)
c$$$           Neutrons
           np = inn
           efn = secant(chem_pot, ef0)
           write(21, *) efp, efn
           do i=1,364
              fp = fermi_dist(ei(i), efp, temp)
              fn = fermi_dist(ei(i), efn, temp)
              if( (fp .lt. 1e-8) .or. (fn .lt. 1e-8)) exit
              if ((fp .ge. 1d-45) .and. (fp .lt. 1.0d0))
     &             sp = -(fp*dlog(fp) + (1-fp)*dlog(1-fp))
              if ((fn .ge. 1d-45) .and. (fn .lt. 1.0d0))
     &             sn = -(fn*dlog(fn) + (1-fn)*dlog(1-fn))
c$$$              write(75, '(i4,5f12.6)') i,ei(i),fp,fn,sp,sn
              sum=sum+ (fp+fn)*ei(i)-temp*(sp+sn)
           end do
        end if
        fes_calc = sum
        return
      end function

      function sumene(ip,spe,ndim)
      implicit real*8(a-h,o-z)
      dimension spe(ndim)
      sumene=0.d0
      do i=1,ip/2
        sumene=sumene+2*spe(i)
      end do
      if(mod(ip,2).ne.0)sumene=sumene+spe(ip/2+1)
      end
