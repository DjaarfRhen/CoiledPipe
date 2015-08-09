*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,or,ot,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ut(0:Imax,0:Jmax,0:*)
     >,vt(0:Imax,0:Jmax,0:*)
     >,wt(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,or(0:Imax,0:Jmax,0:*)
     >,ot(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
     >/Re/Re
*
* Boundary conditions
      do k=1,Km
        k1=mod(k+kd-1,Km)+1
        do j=1,Jm
          u(0,j,k)=u(Im,j,k1)
        end do
      end do
      do i=0,Im
        do j=1,Jm
          u(i,j,0)=u(i,j,Km)
          u(i,j,Km+1)=u(i,j,1)
        end do
        do k=1,Km
          hs0=1.+yt(Jm)*snm(i,k)
          hs1=1.+yt(Jm+1)*snm(i,k)
          u(i,Jm+1,k)=-u(i,Jm,k) * hs0/hs1
        end do
      end do
	  
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=0.d0
        end do
      end do
      do j=1,Jm
        do k=1,Km
          k1=mod(k+kd-1,Km)+1
          v(0,j,k)=v(Im,j,k1)
          v(Im+1,j,k1)=v(1,j,k)
        end do
        do i=1,Im
          v(i,j,0)=v(i,j,Km)
          v(i,j,Km+1)=v(i,j,1)
        end do
      end do
	  
      do j=1,Jm
        do i=1,Im
          w(i,j,0)=w(i,j,Km)
        end do
      end do
      do k=0,Km
        k1=mod(k+kd-1,Km)+1
        do i=1,Im
          w(i,Jm+1,k)=-w(i,Jm,k)*yt(Jm)/yt(Jm+1)
        end do
        do j=1,Jm
          w(0,j,k)=w(Im,j,k1)
          w(Im+1,j,k1)=w(1,j,k)
        end do
      end do
      do j=1,Jm
      w(Im+1,j,0)=w(Im+1,j,Km)
      end do
*
* Vorticities
      do i=1,Im
        do j=1,Jm
          do k=0,Km
            w0=w(i,j,k)
            w1=w(i,j+1,k)
            v0=v(i,j,k)
            v1=v(i,j,k+1)
            ox(i,j,k)=((yt(j+1)*w1-yt(j)*w0)/rt1(j)
     >                -(v1-v0)/ht)/rt(j)
          end do
        end do
        j=0
          sw=0.d0
          do k=1,Km
            sw=sw+w(i,j+1,k)
          end do
          sw=sw/Km
          do k=0,Km
            ox(i,j,k)=2.d0*sw/yt(1)
          end do
      end do
      do k=0,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j,k+1)
          hs0=1.+yt(j)*snm(i,k)
          hs1=1.+yt(j)*snm(i,k+1)
            w0=w(i,j,k)
            w1=w(i+1,j,k)
          hss=1.+yt(j)*snn(i,k)
            or(i,j,k)=((hs1*u1-hs0*u0)/(yt(j)*ht)
     >               -(w1-w0)/hx)/hss
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=0,Im
            u0=u(i,j,k)
            u1=u(i,j+1,k)
          hs0=1.+yt(j)*snm(i,k)
          hs1=1.+yt(j+1)*snm(i,k)
            v0=v(i,j,k)
            v1=v(i+1,j,k)
          hss=1.+rt(j)*snm(i,k)
            ot(i,j,k)=((v1-v0)/hx
     >               -(hs1*u1-hs0*u0)/rt1(j))/hss
          end do
		  
*          write(*,*)'(i, j, k)', i, j, k
*          write(*,*)'(u(i,j,k), u(i,j+1,k), delta)',u0, u1, rt1(j)
*          pause
		  
        end do
      end do
	  
* C CHECK: save x to file
*     open(unit=1, file='log/vort3theta.txt') 
*     do k=1,Km
*       do j=1,Jm
*         do i=0,Im
*           write(1,*) ot(i,j,k)
*         end do
*       end do
*     end do
*     close(1)
* C CHECK: save x to file
*     open(unit=1, file='log/vort2r.txt') 
*     do k=0,Km
*       do j=1,Jm
*         do i=0,Im
*           write(1,*) or(i,j,k)
*         end do
*       end do
*     end do
*     close(1) 
*     write(*,*)'right part: vort'
*     pause;

* Nonlinear terms
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            v0=0.5d0*(v(i,j-1,k)+v(i+1,j-1,k))
            v1=0.5d0*(v(i,j,k)+v(i+1,j,k))
            ot0=rt(j-1)*rt1(j-1)*ot(i,j-1,k)
            ot1=rt(j)*rt1(j)*ot(i,j,k)
            w0=0.5d0*(w(i,j,k-1)+w(i+1,j,k-1))
            w1=0.5d0*(w(i,j,k)+w(i+1,j,k))
            or0=or(i,j,k-1)
            or1=or(i,j,k)
            ut(i,j,k)=
     >           0.5d0*((v0*ot0+v1*ot1)/(yt(j)*yt1(j))
     >               -(w0*or0+w1*or1))
          end do
        end do
      end do
      do k=1,Km
        do i=1,Im
          do j=1,Jm-1
            w0=0.5d0*(w(i,j,k-1)+w(i,j+1,k-1))
            w1=0.5d0*(w(i,j,k)+w(i,j+1,k))
          hs0=1.+rt(j)*smn(i,k-1)
          hs1=1.+rt(j)*smn(i,k)
            ox0=hs0*ox(i,j,k-1)
            ox1=hs1*ox(i,j,k)
          hs0=1.+yt(j)*snm(i-1,k)
          hs1=1.+yt(j+1)*snm(i-1,k)
            u0=0.5d0*(hs0*u(i-1,j,k)+hs1*u(i-1,j+1,k))
          hs0=1.+yt(j)*snm(i,k)
          hs1=1.+yt(j+1)*snm(i,k)
            u1=0.5d0*(hs0*u(i,j,k)+hs1*u(i,j+1,k))
            ot0=ot(i-1,j,k)
            ot1=ot(i,j,k)
          hss=1.+rt(j)*smm(i,k)
            vt(i,j,k)=
     >           0.5d0*((w0*ox0+w1*ox1)
     >                 -(u0*ot0+u1*ot1))/hss
          end do
          vt(i,Jm,k)=0.d0
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=1,Im
          hs0=1.+yt(j)*snm(i-1,k)
          hs1=1.+yt(j)*snm(i-1,k+1)
            u0=0.5d0*(hs0*u(i-1,j,k)+hs1*u(i-1,j,k+1))
          hs0=1.+yt(j)*snm(i,k)
          hs1=1.+yt(j)*snm(i,k+1)
            u1=0.5d0*(hs0*u(i,j,k)+hs1*u(i,j,k+1))
            or0=or(i-1,j,k)
            or1=or(i,j,k)
            v0=0.5d0*(v(i,j-1,k)+v(i,j-1,k+1))
            v1=0.5d0*(v(i,j,k)+v(i,j,k+1))
          hs0=1.+rt(j-1)*smn(i,k)
          hs1=1.+rt(j)*smn(i,k)
            ox0=rt(j-1)*rt1(j-1)*hs0*ox(i,j-1,k)
            ox1=rt(j)*rt1(j)*hs1*ox(i,j,k)
          hss=1.+yt(j)*smn(i,k)
            wt(i,j,k)=
     >           0.5d0*((u0*or0+u1*or1)
     >                 -(v0*ox0+v1*ox1)/(yt(j)*yt1(j)))/hss
          end do
        end do
      end do
*
* Viscous terms
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            ot0=ot(i,j-1,k)
            ot1=ot(i,j,k)
            or0=or(i,j,k-1)
            or1=or(i,j,k)
            ut(i,j,k)=          ut(i,j,k)
     >               -((rt(j)*ot1-rt(j-1)*ot0)/yt1(j)
     >               -(or1-or0)/ht)/yt(j)/Re
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            ox0=ox(i,j,k-1)
            ox1=ox(i,j,k)
            ot0=ot(i-1,j,k)
            ot1=ot(i,j,k)
            hs0=1.+rt(j)*smn(i,k-1)
          hs1=1.+rt(j)*smn(i,k)
            hss=1.+rt(j)*smm(i,k)
            vt(i,j,k)=          vt(i,j,k)
     >               -((hs1*ox1-hs0*ox0)/(rt(j)*ht)
     >                -(ot1-ot0)/hx)/Re/hss
          end do
        end do
      end do
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            ox0=ox(i,j-1,k)
            ox1=ox(i,j,k)
            or0=or(i-1,j,k)
            or1=or(i,j,k)
            hs0=1.+rt(j-1)*smn(i,k)
          hs1=1.+rt(j)*smn(i,k)
            hss=1.+yt(j)*smn(i,k)
            wt(i,j,k)=          wt(i,j,k)
     >               -((or1-or0)/hx
     >                -(hs1*ox1-hs0*ox0)/yt1(j))/Re/hss
          end do
        end do
      end do
      return
      end

