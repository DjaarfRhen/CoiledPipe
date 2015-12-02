*
*     program an
      implicit real*8 (a-h,o-z)
      parameter (Imax=257, Jmax=65, Kmax=257, N=1e5, NZ=1e7)
      character*1 getcharqq
      logical peekcharqq
      character*12 fncp,fndat
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,u1(0:Imax,0:Jmax,0:Kmax)
     >,v1(0:Imax,0:Jmax,0:Kmax)
     >,w1(0:Imax,0:Jmax,0:Kmax)
     >,u2(0:Imax,0:Jmax,0:Kmax)
     >,v2(0:Imax,0:Jmax,0:Kmax)
     >,w2(0:Imax,0:Jmax,0:Kmax)
     >,u3(0:Imax,0:Jmax,0:Kmax)
     >,v3(0:Imax,0:Jmax,0:Kmax)
     >,w3(0:Imax,0:Jmax,0:Kmax)
     >,ox(0:Imax,0:Jmax,0:Kmax)
     >,or(0:Imax,0:Jmax,0:Kmax)
     >,ot(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
     >,p0(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,rkap,tors,epsr
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
     >/Re/Re
     >/Dp/Dp
     >/pressMatr/irowInd(1e7),icolInd(1e7),irowPtr(1e5),values(1e7)
     >,irowsNum
     >/linSolver/iparam(6),rparam(5),IRFAC(1e7),JCFAC(1e7)
     >,IPVT(1e5),JPVT(1e5),FACT(1e7),ipath,NL,NFAC
*

*********************************************************************
**** init start *****************************************************
*********************************************************************

    
      
      
      
* corresponding END DO is in line ~444     
*      do ind_rad=0,0
*      do ind_pitch=0,0
      do ind_rad=1,10
      do ind_pitch=0,0
      do ind_Re=1,10
*      do ind_rad=1,40
*      do ind_pitch=0,40
      
* some of those parameters are modified during the simulation,
* so it has to be inside this loop to overwrite it every time   
      open(5,file='globalStart.car')
      read(5,*) Re
      read(5,*) rkap
      read(5,*) tors
      read(5,*) Xmax
      read(5,*) Im
      read(5,*) Jm,epsr
      read(5,*) Km
      read(5,*) dt
      read(5,*) amp
      read(5,*) fncp
      read(5,*) tol
      read(5,*) nprt
      read(5,*) nwrt
      read(5,*) tmax
      read(5,*) dtmax
      read(5,*) fndat
      close(5)
*
      one=1.d0
      pi=4.d0*atan(one)
      
      write(*,*)' ***************************************************'
      write(*,*)' ***************************************************'
      write(*,*)' PARAMETERS CHANGE'
      write(*,*)' ***************************************************'
            
      t = 0.0
      
* I wanted to do 5* and 5* for 40x40 map     
*      Rad = 5      
      Rad       = 10.0*ind_rad
      HelixStep = 2.0*ind_pitch      
* helix pitch      
      Pitch = HelixStep/(2*pi)   
      
* curvature      
      rkap = abs(Rad)/(Rad**2+Pitch**2)
*      rkap = 0
      
* torsion      
      tors = Pitch/(Rad**2+Pitch**2)
*      tors = 0      
      
      Re = 200*ind_Re
      
      write(*,*)' ***************************************************'
      write(*,198) Re,Xmax,ind_rad,ind_pitch,Rad,HelixStep,rkap,tors
      write(*,*)' ***************************************************'
198   format(' Re=',e9.3' Xmax=',e9.3, /,
     >' ind_rad=',i3' ind_pitch=',i3, /,
     >' Rad=',e9.3,' Pitch=',e9.3, /,
     >' rkap=',e9.3,' tors=',e9.3)
      
 
        open(unit=1, file='log/parameters.txt') 
        write(1,*) Re,Xmax,ind_rad,ind_pitch,Rad,HelixStep,rkap,tors
        close(1)
 
*      rkap = 0
*      tors = 0
      
      call com(u,v,w,p,p0,Imax,Jmax,N,NZ)

*
      write(*,*)' ***************************************************'
      write(*,200) t,dt,Re,Xmax,rkap,tors,epsr,kd,Im,Jm,Km
      write(*,*)' ***************************************************'
200   format(' *  t=',1pe10.3,' dt=',e9.2,' Re=',e9.2
     >'            *',/,
     >' *  Xmax=',0pf6.2,'  rkap=',f6.2,'  tors=',f6.2,' epsr=',f6.2
     >' *  kd=',i3
     >,'  Im=',i3,'  Jm=',i3,'  Km=',i3,'   *')
*
      if(Im.gt.Imax-1.or.Im.gt.256) then
        write(*,*)'  Im=',Im,'  is greater than   Imax-1=',Imax-1
        stop
      end if
      if(Jm.gt.Jmax-1.or.Jm.gt.128) then
        write(*,*)'  Jm=',Jm,'  is greater than   Jmax-1=',Jmax-1
        stop
      end if
      if(Km.gt.Kmax-1.or.Km.gt.256) then
        write(*,*)'  Km=',Km,'  is greater than   Kmax=',Kmax
        stop
      end if
*


*      pause


*
      do k=1,Km
        do j=1,Jm
            do i=1,Im
            u(i,j,k)=(1.-yt(j)**2)*sin(2.*pi*i/Im)*cos(2.*pi*k/Km)
            v(i,j,k)=0.
            w(i,j,k)=0.
          end do
        end do
      end do
      p(0,0,1)=0.
      call presLame(u,v,w,p,p0,Imax,Jmax,N,NZ)
      ss=0.
      aa=0.
      do k=1,Km
      do j=1,Jm
        do i=1,Im
          hss=1.+yt(j)*smm(i,k)
          aa=aa+hss*yt(j)*yt1(j)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
          ss=ss+hss*yt(j)*yt1(j)
        end do
      end do
      end do
      c=amp/sqrt(aa/ss)
      do k=1,Km
        do j=1,Jm
            do i=1,Im
            u(i,j,k)=c*u(i,j,k)+1.-yt(j)**2
            v(i,j,k)=c*v(i,j,k)
            w(i,j,k)=c*w(i,j,k)
          end do
        end do
      end do

*      open(9,file=fncp,form='unformatted')
*      write(9)t,dt,Re,Xmax,rkap,tors,epsr,Im,Jm,Km
*        do k=1,Km
*          do j=1,Jm
*            write(9)(u(i,j,k),i=1,Im)
*            write(9)(v(i,j,k),i=1,Im)
*            write(9)(w(i,j,k),i=1,Im)
*          end do
*      end do
*      close(9)
*      close(8)
*      write(*,*)'*************************************************'
*      write(*,*)'*               Control Point                   *'
*      write(*,*)'*************************************************'
*      pause




*********************************************************************
*** init end ********************************************************
*********************************************************************




*********************************************************************
*** pipe start ******************************************************
*********************************************************************





*
* C CHECK: save x to file
      open(unit=1, file='log/u0.txt') 
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            write(1,*) u(i,j,k)
          end do
        end do
      end do
      close(1)
	  
	  
	  
      dt=min(dt,dtmax)

*      call initLinSolver(Imax,Jmax,N,NZ)    
      
      p(0,0,1) = 0.5  
      call presLame(u,v,w,p,p0,Imax,Jmax,N,NZ)

*
* C CHECK: save x to file
      open(unit=1, file='log/u0div.txt') 
      do k=1,Km
        do j=1,Jm
          do i=0,Im
            write(1,*) u(i,j,k)
          end do
        end do
      end do
      close(1)
 
* C CHECK: save x to file
      open(unit=1, file='log/v0div.txt') 
      do k=1,Km
        do j=0,Jm
          do i=1,Im
            write(1,*) v(i,j,k)
          end do
        end do
      end do
      close(1)	  
 
* C CHECK: save x to file
      open(unit=1, file='log/w0div.txt') 
      do k=0,Km
        do j=1,Jm
          do i=1,Im
            write(1,*) w(i,j,k)
          end do
        end do
      end do
      close(1)	  



      
      call rp(t,u,v,w,u1,v1,w1,ox,or,ot,Imax,Jmax)

        open(unit=1, file='log/ut.txt') 
        do i=0,Im
          do j=1,Jm
            do k=1,Km
              write(1,*) u1(i,j,k)
            end do
          end do
        end do
        close(1)
        
        open(unit=1, file='log/vt.txt') 
        do i=1,Im
          do j=0,Jm
            do k=1,Km
              write(1,*) v1(i,j,k)
            end do
          end do
        end do
        close(1)
        
        open(unit=1, file='log/wt.txt') 
        do i=1,Im
          do j=1,Jm
            do k=0,Km
              write(1,*) w1(i,j,k)
            end do
          end do
        end do
        close(1)

*        write(*,*) 'velocity time step'
*        pause




      
*	  pause
	  
      dd=0.
      do k=1,Km-1
        do j=1,Jm-1
          do i=1,Im-1
          ox0=ox(i,j,k)
          ox1=ox(i+1,j,k)
          hnn0=1.+yt(j)*snn(i,k)
          hnn1=1.+yt(j+1)*snn(i,k)
          oy0=yt(j)*hnn0*or(i,j,k)
          oy1=yt(j+1)*hnn1*or(i,j+1,k)
          hnm0=1.+rt(j)*snm(i,k)
          hnm1=1.+rt(j)*snm(i,k+1)
          ot0=hnm0*ot(i,j,k)
          ot1=hnm1*ot(i,j,k+1)
          hss=1.+rt(j)*snn(i,k)
          d=((ox1-ox0)/hx+(oy1-oy0)/rt(j)/rt1(j)+(ot1-ot0)/ht/rt(j))/hss
          dd=max(dd,abs(d))
*          write(*,*)i,j,k,d
         end do
         end do
         end do
*         write(*,*)'  div(omega)=',dd
        d=0.
        do k=1,Km
         do j=1,Jm
          do i=1,Im
           d=d+u(i,j,k)*u1(i,j,k)*yt(j)*yt1(j)*(1.+yt(j)*snm(i,k))
     >              +v(i,j,k)*v1(i,j,k)*rt(j)*rt1(j)*(1.+rt(j)*smm(i,k))
     >              +w(i,j,k)*w1(i,j,k)*yt(j)*yt1(j)*(1.+yt(j)*smn(i,k))
          end do
          end do
          end do
           d=hx*ht*d
*        write(*,*)'  Ee=',d
*
        dd=0.
        do k=2,Km
        do j=1,Jm
        do i=2,Im
*              ut0=u1(i-1,j,k)
*              ut1=u1(i,j,k)
*              hmm0=1.+rt(j-1)*smm(i,k)
*              hmm1=1.+rt(j)*smm(i,k)
*              vt0=rt(j-1)*hmm0*v1(i,j-1,k)
*              vt1=rt(j)*hmm1*v1(i,j,k)
*              hmn0=1.+yt(j)*smn(i,k-1)
*              hmn1=1.+yt(j)*smn(i,k)
*              wt0=hmn0*w1(i,j,k-1)
*              wt1=hmn1*w1(i,j,k)
*              hss=1.+yt(j)*smm(i,k)
*        d=((ut1-ut0)/hx+(vt1-vt0)/yt(j)/yt1(j)+(wt1-wt0)/ht/yt(j))/hss
         call div(i,j,k,u1,v1,w1,d,Imax,Jmax)
          dd=max(dd,abs(d))
*               write(*,*)i,j,k,d
         end do
         end do
         end do
*         write(*,*)'  div(rp)=',dd
*             pause
      p(0,0,1)=0.d0
*                stop
*      call pres(u1,v1,w1,p,Imax,Jmax)
      p(0,0,1) = 0.5  
      call presLame(u1,v1,w1,p,p0,Imax,Jmax,N,NZ)

*      pause
*      stop
*
      open(8,file=fndat,access='append')

      call wrt2file(t,dt,u,v,w,p,Imax,Jmax)
      call prt2screen(t,dt,u,v,w,p,Imax,Jmax)
      lprt=0
      lwrt=0
      
*     MAIN LOOP      
10    continue
*
      call step(t,dt,tol,u,v,w,u1,v1,w1,u2,v2,w2,u3,v3,w3
     >    ,ox,or,ot,p,p0,Imax,Jmax,N,NZ)
     
*     stop
      lprt=lprt+1
      lwrt=lwrt+1
      dt=min(dt,dtmax)
*      if(peekcharqq()) then
*        write(*,*)'  Keeboard symbol typed'
*        if(getcharqq().eq.'w') tmax=t
*      end if
      call servis(t,u,v,w,ox,or,ot,p,0,Imax,Jmax)
*      write(*,*)'  t=',t,'  dt=',dt
      if(lprt.ge.nprt.or.t.ge.tmax) then
        lprt=0
        call servis(t,u,v,w,ox,or,ot,p,1,Imax,Jmax)
        call prt2screen(t,dt,u,v,w,p,Imax,Jmax)
      end if
      if(lwrt.ge.nwrt.or.t+0.01d0*dt.ge.tmax) then
        lwrt=0
        call wrt2file(t,dt,u,v,w,p,Imax,Jmax)
*        call servis(t,u,v,w,ox,or,ot,p,2,Imax,Jmax)
*        open(9,file=fncp,form='unformatted')
*        write(9)t,dt,Re,Xmax,rkap,tors,epsr,Im,Jm,Km
*        do k=1,Km
*          do j=1,Jm
*            write(9)(u(i,j,k),i=1,Im)
*            write(9)(v(i,j,k),i=1,Im)
*            write(9)(w(i,j,k),i=1,Im)
*          end do
*        end do
*        close(9)
*        close(8)
*        
*        open(unit=1, file='log/parameters.txt') 
*        write(1,*) Re,Xmax,rkap,tors,epsr,Im,Jm,Km
*        close(1)
*        
*        open(unit=1, file='log/u.txt') 
*        do i=0,Im
*          do j=1,Jm
*            do k=1,Km
*              write(1,*) u(i,j,k)
*            end do
*          end do
*        end do
*        close(1)
*        
*        open(unit=1, file='log/v.txt') 
*        do i=1,Im
*          do j=0,Jm
*            do k=1,Km
*              write(1,*) v(i,j,k)
*            end do
*          end do
*        end do
*        close(1)
*        
*        open(unit=1, file='log/w.txt') 
*        do i=1,Im
*          do j=1,Jm
*            do k=0,Km
*              write(1,*) w(i,j,k)
*            end do
*          end do
*        end do
*        close(1)
*        
*        open(unit=1, file='log/p.txt') 
*        do i=1,Im
*          do j=1,Jm
*            do k=1,Km
*              write(1,*) p(i,j,k)
*            end do
*          end do
*        end do
*        close(1)
        
*        write(*,*)'*************************************************'
*        write(*,*)'*               Control Point                   *'
*        write(*,*)'*************************************************'
        
*        pause
      
*        open(8,file=fndat,access='append')
      end if
      if(t+0.01d0*dt.lt.tmax) goto 10



*********************************************************************
**** pipe end *******************************************************
*********************************************************************

* END DO for parameter change loop
      end do  
      end do
      end do


      stop
1     write(*,*)'  File ',fncp,' was not found'
100   format(10(1pe12.4))

      stop
      end
      