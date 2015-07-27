*
      subroutine step(t,dt,tol,u,v,w,u1,v1,w1,u2,v2,w2,u3,v3,w3
     >     ,ox,oy,oz,p,p0,Imax,Jmax,N,NZ)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,u1(0:Imax,0:Jmax,0:*)
     >,v1(0:Imax,0:Jmax,0:*)
     >,w1(0:Imax,0:Jmax,0:*)
     >,u2(0:Imax,0:Jmax,0:*)
     >,v2(0:Imax,0:Jmax,0:*)
     >,w2(0:Imax,0:Jmax,0:*)
     >,u3(0:Imax,0:Jmax,0:*)
     >,v3(0:Imax,0:Jmax,0:*)
     >,w3(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,p0(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
*
      c12=1.d0/2.d0
      c13=1.d0/3.d0
      c23=2.d0/3.d0
      c14=1.d0/4.d0
      c34=3.d0/4.d0
      c32=3.d0/2.d0
      c38=3.d0/8.d0
      c58=5.d0/8.d0
      facmin=0.5d0
      facmax=1.5d0
*
10    continue
*
      tau=c13*dt
      t1=t+c23*dt
      t2=t+dt
      dt23=c23*dt
      dt13=c13*dt
      dt34=c34*dt
* 1-st step
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u1(i,j,k)=dt23*u1(i,j,k)
            v1(i,j,k)=dt23*v1(i,j,k)
            w1(i,j,k)=dt23*w1(i,j,k)
          end do
        end do
      end do
      call lin(tau,u1,v1,w1,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u1(i,j,k)=u1(i,j,k)+u(i,j,k)
            v1(i,j,k)=v1(i,j,k)+v(i,j,k)
            w1(i,j,k)=w1(i,j,k)+w(i,j,k)
          end do
        end do
      end do
* 2-nd step
      call rp(t1,u1,v1,w1,u2,v2,w2,ox,oy,oz,Imax,Jmax)
      p(0,0,1)=0.
*      call pres(u2,v2,w2,p,Imax,Jmax)
      call presLame(u2,v2,w2,p,p0,Imax,Jmax,N,NZ)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u2(i,j,k)=dt13*u2(i,j,k)
     >               -(u1(i,j,k)-u(i,j,k))
            v2(i,j,k)=dt13*v2(i,j,k)
     >               -(v1(i,j,k)-v(i,j,k))
            w2(i,j,k)=dt13*w2(i,j,k)
     >               -(w1(i,j,k)-w(i,j,k))
          end do
        end do
      end do
      call lin(tau,u2,v2,w2,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u2(i,j,k)=u2(i,j,k)+c32*u1(i,j,k)-c12*u(i,j,k)
            v2(i,j,k)=v2(i,j,k)+c32*v1(i,j,k)-c12*v(i,j,k)
            w2(i,j,k)=w2(i,j,k)+c32*w1(i,j,k)-c12*w(i,j,k)
          end do
        end do
      end do
* 3-rd step
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u3(i,j,k)=c34*(u2(i,j,k)-u(i,j,k))
            v3(i,j,k)=c34*(v2(i,j,k)-v(i,j,k))
            w3(i,j,k)=c34*(w2(i,j,k)-w(i,j,k))
          end do
        end do
      end do
      call lin(tau,u3,v3,w3,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u3(i,j,k)=u3(i,j,k)+c32*u2(i,j,k)-c34*u1(i,j,k)+c14*u(i,j,k)
            v3(i,j,k)=v3(i,j,k)+c32*v2(i,j,k)-c34*v1(i,j,k)+c14*v(i,j,k)
            w3(i,j,k)=w3(i,j,k)+c32*w2(i,j,k)-c34*w1(i,j,k)+c14*w(i,j,k)
          end do
        end do
      end do
* 4-th step
      call rp(t1,u2,v2,w2,u1,v1,w1,ox,oy,oz,Imax,Jmax)
      call gradp(u1,v1,w1,p,p0,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u1(i,j,k)=dt34*u1(i,j,k)
     >          -(u3(i,j,k)-c38*u2(i,j,k)-c58*u(i,j,k))
            v1(i,j,k)=dt34*v1(i,j,k)
     >          -(v3(i,j,k)-c38*v2(i,j,k)-c58*v(i,j,k))
            w1(i,j,k)=dt34*w1(i,j,k)
     >          -(w3(i,j,k)-c38*w2(i,j,k)-c58*w(i,j,k))
          end do
        end do
      end do
      call lin(tau,u1,v1,w1,Imax,Jmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u1(i,j,k)=u1(i,j,k)+c12*u3(i,j,k)+c34*u2(i,j,k)-c14*u(i,j,k)
            v1(i,j,k)=v1(i,j,k)+c12*v3(i,j,k)+c34*v2(i,j,k)-c14*v(i,j,k)
            w1(i,j,k)=w1(i,j,k)+c12*w3(i,j,k)+c34*w2(i,j,k)-c14*w(i,j,k)
          end do
        end do
      end do
      p(0,0,1)=c12
*      call pres(u1,v1,w1,p,Imax,Jmax)
      call presLame(u1,v1,w1,p,p0,Imax,Jmax,N,NZ)
* Accuracy estimation
      err=0.
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            uu=u1(i,j,k)-u3(i,j,k)
            vv=v1(i,j,k)-v3(i,j,k)
            ww=w1(i,j,k)-w3(i,j,k)
            err=max(err,abs(uu),abs(vv),abs(ww))
          end do
        end do
      end do
      fac=(tol/err)**c13
      if(fac.lt.facmin) then
        dt=dt*fac
        write(*,*)'  STEP:  fac=',fac,'  dt=',dt
        call rp(t,u,v,w,u1,v1,w1,ox,oy,oz,Imax,Jmax)
        p(0,0,1)=0.
*        call pres(u1,v1,w1,p,Imax,Jmax)
      call presLame(u1,v1,w1,p,p0,Imax,Jmax,N,NZ)
        goto 10
      end if
      fac=min(fac,facmax)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u1(i,j,k)
            v(i,j,k)=v1(i,j,k)
            w(i,j,k)=w1(i,j,k)
          end do
        end do
      end do
      t=t+dt
      dt=fac*dt   
      call rp(t,u,v,w,u1,v1,w1,ox,oy,oz,Imax,Jmax)
      p(0,0,1)=0.
*      call pres(u1,v1,w1,p,Imax,Jmax)   
      call presLame(u1,v1,w1,p,p0,Imax,Jmax,N,NZ)   
      return
      end