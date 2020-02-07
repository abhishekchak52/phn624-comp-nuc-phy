      implicit real*8(a-h,o-z)
      external f1
      write(*,*) "enter the intial guess and step value"
      read (*,*) x0,xinc
      call bracket(f1,x0,xinc,a,b)
      write(*,*)"The brackets are ",a,b
      x=bisec(f1,a,b)
      write(*,*)"The solution is ",x,a,b
      end

      real*8 function f1(x)
      implicit real*8(a-h,o-z)
      f1=sin(x)
      end


      subroutine bracket(f,x0,xinc,xpos,xneg)
* sign of xinc should be +/- for increasing/decreasing function
      implicit real*8(a-h,o-z)
      logical flag
      flag=.true.
      x=x0
      if(f(x).gt.0.0)then
        xpos=x
        ival=1
        x=x-xinc
      else
        xneg=x
        ival=-1
        x=x+xinc
      end if
      do iterate=1,50
        if(f(x).gt.0.0)then
          xpos=x
          if(ival.eq.-1)then
            flag=.false.
            exit
          end if
          x=x-xinc
        else
          xneg=x
          if(ival.eq.1)then
            flag=.false.
            exit
          end if
          x=x+xinc
        end if
      end do
      if(flag)then
        write(*,*)'Havent got change in sign for bracketing'
        stop
      end if
      end


      real*8 function bisec(f,xpos,xneg)
      implicit real*8(a-h,o-z)
      logical flag
      flag=.true.
      if(f(xpos).lt.0)then
        t=xpos
        xpos=xneg
        xneg=t
      endif
      x0=(xpos+xneg)/2.0d0
      do iterate=1,300
        x=x0
        fx=f(x)
        if(fx.gt.0.0)then
          xpos=x
        else
          xneg=x
        end if
        x0=(xpos+xneg)/2.0d0
c        print *,iterate,x,fx
        if((dabs(x0-x).lt.1.0e-8).and.(dabs(fx).lt.1e-3))then
          flag=.false.
          bisec=x0
          exit
        end if
      end do
      if(flag)then
        write(*,*)'Bisection method did not converge'
        stop
      end if
      end
