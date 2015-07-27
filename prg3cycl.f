*
      subroutine prg3cycl(a,b,c,y,x,m)
      implicit real*8 (a-h,o-z)
      dimension a(m),b(m),c(m),y(m),x(m)
      do j=2,m-2
        e=a(j)/b(j-1)
        b(j)=b(j)-e*c(j-1)
        y(j)=y(j)-e*y(j-1)
        a(j)=-e*a(j-1)
*
        e=c(m)/b(j-1)
        c(m)=-e*c(j-1)
        b(m)=b(m)-e*a(j-1)
        y(m)=y(m)-e*y(j-1)
      end do
      j=m-1
        e=a(j)/b(j-1)
        b(j)=b(j)-e*c(j-1)
        c(j)=c(j)-e*a(j-1)
        y(j)=y(j)-e*y(j-1)
*
        e=c(m)/b(j-1)
        a(m)=a(m)-e*c(j-1)
        b(m)=b(m)-e*a(j-1)
        y(m)=y(m)-e*y(j-1)
      j=m
        e=a(j)/b(j-1)
        b(j)=b(j)-e*c(j-1)
        y(j)=y(j)-e*y(j-1)
*
*
      x(m)=y(m)/b(m)
      x(m-1)=(y(m-1)-c(m-1)*x(m))/b(m-1)
      do j=m-2,1,-1
        x(j)=(y(j)-c(j)*x(j+1)-a(j)*x(m))/b(j)
      end do
      return
      end

