*
      subroutine wrt2file(t,dt,u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
     >/dim/Xmax,rkap,tors,epsr
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
     >/Re/Re
     >/Dp/Dp
*
      ubulk=0.
      Ss=0.
      amp=0.
      dd=0.
      do j=1,Jm
        u0=0.
        do k=1,Km
          do i=1,Im
            hss=1.+yt(j)*smm(i,k)
            ubulk=ubulk+u(1,j,k)*hss*yt(j)*yt1(j)
            Ss=Ss+hss*yt(j)*yt1(j) 
            u0=u0+u(i,j,k)
          end do
        end do
        u0=u0/(Im*Km)
        do k=1,Km
          do i=1,Im
            hss=1.+yt(j)*smm(i,k)
            amp=amp
     >   +((u(i,j,k)-u0)**2+w(i,j,k)**2+v(i,j,k)**2)*hss*yt(j)*yt1(j)
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      ubulk=ubulk/Ss
      amp=sqrt(amp/Ss)
      ucl=0.
      do k=1,Km
        do i=1,Im
          ucl=ucl+u(i,1,k)
        end do
      end do
      ucl=ucl/(Im*Km)

	  pu = u(Im/2,Jm/2,Km/2)
      pw = w(Im/2,Jm/2,Km/2)
      pv = v(Im/2,Jm/2,Km/2)
*
 
      write(8,120)t,dt,Dp,amp,ucl,dd,pu,pv,pw,ubulk
120   format(e12.4,' ',e12.4,' ',e12.4,' ',e12.4,' ',e12.4,' ',e12.4,' ',
     > e12.4,' ',e12.4,' ',e12.4,' ',e12.4) 

      return
      end

