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
     >/pry/apy(128),bpy(128),cpy(128)
     >/prr/apr(128),bpr(128),cpr(128)
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
            u(i,j,k)=ep(k)
          end do
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*rt(j)**2
            dp(k)=-ct*v(i,j,k)*rt(j)**2
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
            w(i,j,k)=ep(k)
          end do
            w(i,j,0)=w(i,j,Km)
        end do
      end do
*
* Ar
*   U,V
      do i=1,Im
        do k=1,Km
          do j=1,Jm
            bp(j)=bpy(j)-ct
            dp(j)=-ct*u(i,j,k)
          end do
          bp(Jm)=bp(Jm)-cpy(Jm)
          call prog3(apy,bp,cpy,dp,ep,Jm)
          do j=1,Jm
            u(i,j,k)=ep(j)
          end do
          do j=1,Jm-1
            bp(j)=bpr(j)-ct
            dp(j)=-ct*rt(j)*v(i,j,k)
          end do
          call prog3(apr,bp,cpr,dp,ep,Jm-1)
          do j=1,Jm-1
            v(i,j,k)=ep(j)/rt(j)
          end do
          do j=1,Jm
            bp(j)=bpy(j)-ct
            dp(j)=-ct*w(i,j,k)
          end do
          bp(Jm)=bp(Jm)-cpy(Jm)*yt(Jm)/yt(Jm+1)
          call prog3(apy,bp,cpy,dp,ep,Jm)
          do j=1,Jm
            w(i,j,k)=ep(j)
          end do
        end do
      end do
      return
      end

