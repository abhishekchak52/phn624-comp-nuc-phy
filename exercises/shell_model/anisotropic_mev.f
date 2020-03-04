        program main
        real*8 E(10000),en,f,A
        integer i,j,k,N,nz,rho,t
        A=90.
        open(1,file="anisotropic_ho2.out")

        do i=-3,3,1
           delta=real(i)/10.0
        t=0
            do N=0,5,1
                    do nz=0,N,1
                           nrho=N-nz

c                              x=0.0
c                             E0=N+1.5
                                t=t+1
        f=(((1.0+(2.0/3.0)*delta)**2)*(1.0-(4.0/3.0)*delta))**(-1.0/6.0)
       E(t)=(41.*(A**(-1./3.)))*(N+(3.0/2.0)-
     &   ((delta*(2*nz-nrho))/3.0))*f
c                           en=(N+1.5-((delta*(2*nz-nrho))/3.0))

                      end do
              end do
              write(1,*) delta,(E(j),j=1,t)
         end do
         end

