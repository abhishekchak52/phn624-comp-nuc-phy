      real*8 function secant(f,x00)
      implicit real*8(a-h,o-z)
      logical flag
      x0=x00
      f0=f(x0)
      x1=x0+x0/100.0d0
      f1=f(x1)
      flag=.true.
      do i=1,100
        x=x1-f1*(x1-x0)/(f1-f0)
        if(abs(x-x1).le.1e-4) then
          flag=.false.
          exit
        end if
        x0=x1
        f0=f1
        x1=x
        f1=f(x)
      end do
      if (flag) print *, 'No convergence'
      secant=x
      end
