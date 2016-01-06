*
      subroutine com(u,v,w,p,p0,Imax,Jmax,N,NZ)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,p0(0:Imax,0:Jmax,0:*)
      common
     >/dim/Xmax,rkap,tors,epsr
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
     >/pry/apy(128),bpy(128),cpy(128)
     >/prr/apr(128),bpr(128),cpr(128)
     >/rl/rlx(0:128),rlt(0:128)
     >/servst/iserv
     >/pi/pi
     >/set/set0,set1,aset,bset,iset
     >/pressMatr/irowInd(1e7),icolInd(1e7),irowPtr(1e5),values(1e7)
     >,irowsNum
     >/linSolver/iparam(6),rparam(5),IRFAC(1e7),JCFAC(1e7)
     >,IPVT(1e5),JPVT(1e5),FACT(1e7),ipath,NL,NFAC
     
      dimension rhs(1e5), x(1e5)
*

      
      iserv=0
      one=1.d0
      pi=4.d0*atan(one)
*
* dimt
      ht=2.d0*pi/Km
* dimx
      kd = 0
      if(tors.ne.0) then
        kd=tors*Xmax/ht
        Xm0=kd*ht/tors
        Xm1=(kd+1)*ht/tors
* this is round-off to the closest value
* but some times ZERO happens 
* so it is safer to enlarge Xmax in any case       
*        if (Xm1-Xmax.lt.Xmax-Xm0) then
*           kd=kd+1
*        end if
        kd=kd+1
        Xmax=kd*ht/tors
      end if
      kd=mod(kd,Km)
      hx=Xmax/Im
      
      write(*,33) Xmax
33    format(' Xmax=',e9.3)
      
      open(unit=1, file='log/parameters.txt', access='append') 
      write(1,*) Xmax,kd
      close(1)
      
* dimr
      hr=1./Jm
      if(epsr.eq.1.) then
        do j=0,Jm
          ro=j*hr
          rt(j)=ro
          rt1(j)=hr
        end do
        do j=1,Jm+1
          ro=(j-0.5)*hr
          yt(j)=ro
          yt1(j)=hr
        end do
      else
        iset=0
        set0=1.d0
        set1=epsr
        do j=0,Jm
          ro=j*hr
          rt(j)=rrt(ro,0)
          rt1(j)=rrt(ro,1)*hr
        end do
        do j=1,Jm+1
          ro=(j-0.5)*hr
          yt(j)=rrt(ro,0)
          yt1(j)=rrt(ro,1)*hr
        end do
      end if

	  
C CHECK: save rth to file
      open(unit=1, file='log/rth.txt')
      do i = 1, Jm
         write(1,*) yt(i)
      enddo
      close(1)
*
      do j=1,Jm
        c=1.d0/(yt(j)*yt1(j))
        apy(j)=c*rt(j-1)/rt1(j-1)
        cpy(j)=c*rt(j)/rt1(j)
        bpy(j)=-apy(j)-cpy(j)
      end do
      do j=1,Jm-1
        c=1.d0/(rt(j)*rt1(j))
        apr(j)=c*yt(j)/yt1(j)
        cpr(j)=c*yt(j+1)/yt1(j+1)
        bpr(j)=-apr(j)-cpr(j)
      end do
* sst
      do k=0,Km
      tt=k*ht
      do i=0,Im
        s=i*hx
        snn(i,k)=rkap*sin(tt-tors*s)
c        snn(i,k)=-rkap*cos(tt-tors*s)
      end do
      end do
      do k=0,Km+1
      tt=(k-0.5)*ht
      do i=0,Im
        s=i*hx
        snm(i,k)=rkap*sin(tt-tors*s)
c        snm(i,k)=-rkap*cos(tt-tors*s)
      end do
      end do
      do k=0,Km
      tt=k*ht
      do i=0,Im+1
        s=(i-0.5)*hx
        smn(i,k)=rkap*sin(tt-tors*s)
c        smn(i,k)=-rkap*cos(tt-tors*s)
      end do
      end do
      do k=0,Km+1
      tt=(k-0.5)*ht
      do i=0,Im+1
        s=(i-0.5)*hx
        smm(i,k)=rkap*sin(tt-tors*s)
c        smn(i,k)=-rkap*cos(tt-tors*s)
      end do
      end do     
* rl
      do i=0,Im/2
        rlx(i)=(2.d0/hx*sin(i*pi/Im))**2
      end do
      do k=0,Km/2
        rlt(k)=(2.d0/ht*sin(k*pi/Km))**2
      end do


      
*      write(*,*)' ***************************************************'
*      write(*,200) Xmax,rkap,tors,epsr,kd,Im,Jm,Km
*      write(*,*)' ***************************************************'
*200   format(
*     >'            *',/,
*     >' *  Xmax=',0pf6.2,'  rkap=',f6.2,'  tors=',f6.2,' epsr=',f6.2
*     >' *  kd=',i3
*     >,'  Im=',i3,'  Jm=',i3,'  Km=',i3,'   *')
     
     
     
     
      call initLinSolver(Imax,Jmax,N,NZ)

      
     
      
* pressure gradient    initialization  
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            hss=1.+yt(j)*snm(i,k)
            u(i,j,k) = 1.0/hss
            v(i,j,k) =  0.0
            w(i,j,k) =  0.0
          end do
        end do
      end do

      open(unit=123, file='log/u_hss.txt') 
      do k=1,Km
        do j=1,Jm
            do i=1,Im
              write(123,*) u(i,j,k)
            end do
        end do
      end do
      close(123)
      
      do k=1,Km
        do i=1,Im
          do j=1,Jm
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            rhs(i+(j-1)*Im+(k-1)*Im*Jm)=d
          end do
        end do
      end do
             
      rhs(1)=0.0
      
C CHECK: save rhs to file
      open(unit=1, file='log/com_rhs.txt')
      do i = 1, irowsNum
         write(1,*) rhs(i)
      enddo
      close(1)
***********************************************************
***********************************************************
***********************************************************
          
C Solve A * X(i) = B(i)
       CALL DLFSXG(irowsNum, NFAC, NL, FACT, IRFAC, JCFAC, IPVT, JPVT,
     >                 rhs, ipath, x)

           
      do k=1,Km
        do i=1,Im
          do j=1,Jm
            p0(i,j,k)=x(i+(j-1)*Im+(k-1)*Im*Jm)
            p(i,j,k)=0.0
          end do
        end do
      end do

***********************************************************
***********************************************************
***********************************************************

C CHECK: save rhs to file
      open(unit=1, file='log/com_x.txt')
      do i = 1, irowsNum
         write(1,*) x(i)
      enddo
      close(1)
      
      dd = 0.0;
      do k=1,Km  
      do j=1,Jm
      do i=1,Im
         v(i,Jm,k)=0.d0
         call div(i,j,k,u,v,w,d,Imax,Jmax)         
         dd=max(dd,abs(d))
      end do
      end do
      end do
      write(*,*)'  div(u_hss)=',dd
      
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k) = 0.0
            v(i,j,k) = 0.0
            w(i,j,k) = 0.0
          end do
        end do
      end do
      
      p(0,0,0) = 0.0
      call gradp(u,v,w,p,p0,Imax,Jmax)      
      
      dd = 0.0;
      do k=1,Km  
      do j=1,Jm
      do i=1,Im
         v(i,Jm,k)=0.d0
         call div(i,j,k,u,v,w,d,Imax,Jmax)         
         dd=max(dd,abs(d))
      end do
      end do
      end do
      write(*,*)'  div(u_p0)=',dd


      open(unit=123, file='log/p0.txt') 
      do k=1,Km
        do j=1,Jm
            do i=1,Im
              write(123,*) p0(i,j,k)
            end do
        end do
      end do
      close(123)

      open(unit=123, file='log/p.txt') 
      do k=1,Km
        do j=1,Jm
            do i=1,Im
              write(123,*) p(i,j,k)
            end do
        end do
      end do
      close(123)




      

        
*  Mean pressure gradient
*     ss=0.d0
*     su=0.d0
*     do j=1,Jm
*       ss=ss+yt(j)*yt1(j)
*       ssu=0.d0
*       do k=1,Km
*         do i=1,Im
*           ssu=ssu+u(i,j,k)
*         end do
*       end do
*       su=su+ssu*yt(j)*yt1(j)
*     end do
*     Dp=su/(Im*Km*ss)

      ss=0.d0
      su=0.d0
      do j=1,Jm
        ss=ss+yt(j)*yt1(j)
        do k=1,Km
          do i=1,Im
            hss=1.+yt(j)*smm(i,k)
            ss=ss+hss*yt(j)*yt1(j)
            su=su+u(i,j,k)*hss*yt(j)*yt1(j)
          end do
        end do
      end do
      Dp=su/ss

      
      write(*,*)'  Dp=',Dp
      
*      pause
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      return
      end

