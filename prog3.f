*
      subroutine prog3(a,b,c,d,x,m)
      implicit real*8 (a-h,o-z)
      dimension a(m),b(m),c(m),d(m),x(m)
      if(m.le.0) return
      do j=2,m
        e=a(j)/b(j-1)
        b(j)=b(j)-e*c(j-1)
        d(j)=d(j)-e*d(j-1)
      end do
      x(m)=d(m)/b(m)
      do j=m-1,1,-1
        x(j)=(d(j)-c(j)*x(j+1))/b(j)
      end do
      return
      end

