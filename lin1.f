*
      subroutine lin(tau,u,v,w,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ap(256),bp(256),cp(256),dp(256),ep(256)
      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/Re/Re
      ct=Re/tau
*
* At
*   U,V,W
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
            dp(k)=-ct*u(i,j,k)*yt(j)**2
          end do
          call prg3cycl(ap,bp,cp,dp,ep,Km)
          do k=1,Km
*            u(i,j,k)=ep(k)
          end do
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*rt(j)**2
*            dp(k)=-ct*v(i,j,k)*rt(j)**2
          end do
          call prg3cycl(ap,bp,cp,dp,ep,Km)
          do k=1,Km
            v(i,j,k)=ep(k)
          end do
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
            dp(k)=-ct*w(i,j,k)*yt(j)**2
          end do
          call prg3cycl(ap,bp,cp,dp,ep,Km)
          do k=1,Km
*            w(i,j,k)=ep(k)
          end do
        end do
      end do
      return
      end

