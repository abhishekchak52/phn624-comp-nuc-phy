      subroutine sort(a,ndim,n)
      implicit double precision(a-h,o-z)
      logical iflag
      dimension a(ndim)
      do i=1,n
        iflag=.true.
        do j=1,n-1
          if (a(j).GT.a(j+1)) then
            iflag=.false.
            call swap(a(j),a(j+1))
          end if
        end do
        if (iflag) exit
      end do
      end

      subroutine swap(a,b)
      implicit double precision(a-h,o-z)
      t=a
      a=b
      b=t
      end
