      subroutine presLame(u,v,w,p,p0,Imax,Jmax,N,NZ)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,p0(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/rl/rlx(0:128),rlt(0:128)
     >/pry/apy(128),bpy(128),cpy(128)
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
     >/Dp/Dp
     >/pressMatr/irowInd(1e7),icolInd(1e7),irowPtr(1e5),values(1e7)
     >,irowsNum
     >/linSolver/iparam(6),rparam(5),IRFAC(1e7),JCFAC(1e7)
     >,IPVT(1e5),JPVT(1e5),FACT(1e7),ipath,NL,NFAC
     
      dimension rhs(1e5), x(1e5)
*
      Im2=Im/2
      Km2=Km/2
      cik=4.d0/(Im2*Km2)



*
*      do k=1,Km
*        do i=1,Im
*          v(i,Jm,k)=0.d0
*          do j=1,Jm
*            call div(i,j,k,u,v,w,d,Imax,Jmax)
*            rhs(i+(j-1)*Im+(k-1)*Im*Jm)=d
*          end do
*        end do
*      end do
      
      do k=1,Km
        do i=1,Im
          do j=1,(Jm-1)
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            rhs(i+(j-1)*Im+(k-1)*Im*Jm)=d
          end do
          call div(i,Jm,k,u,v,w,d,Imax,Jmax)
          hss = 1. + yt(Jm)  * smm(i,k)
          hsr1= 1. + rt(Jm)  * smm(i,k)
          rhs(i+(Jm-1)*Im+(k-1)*Im*Jm) = d - 
     >                  rt(Jm)*hsr1*v(i,Jm,k)/(yt(Jm)*yt1(Jm)*hss)
        end do
      end do      
             
      rhs(1)=0.0

***********************************************************
***********************************************************
***********************************************************
          
C Solve A * X(i) = B(i)
       CALL DLFSXG(irowsNum, NFAC, NL, FACT, IRFAC, JCFAC, IPVT, JPVT,
     >                 rhs, ipath, x)

           
      do k=1,Km
        do i=1,Im
          do j=1,Jm
            p(i,j,k)=x(i+(j-1)*Im+(k-1)*Im*Jm)
          end do
        end do
      end do

***********************************************************
***********************************************************
***********************************************************

*
*  Mean pressure gradient
*      Ub=p(0,0,1)
*      ss=0.d0
*      su=0.d0
*      do j=1,Jm
*        ss=ss+yt(j)*yt1(j)
*        ssu=0.d0
*        do k=1,Km
*          do i=1,Im
*            ssu=ssu+u(i,j,k)
*          end do
*        end do
*        su=su+ssu*yt(j)*yt1(j)
*      end do
*      Dp=Ub-su/(Im*Km*ss)
*      p(0,0,0)=Dp
      
	  
	  Ub=p(0,0,1)
      ss=0.d0
      su=0.d0
      do j=1,Jm
        do k=1,Km
          do i=1,Im
            hss=1.+yt(j)*smm(i,k)
            ss=ss+hss*yt(j)*yt1(j)
            su=su+u(i,j,k)*hss*yt(j)*yt1(j)
          end do
        end do
      end do
      Dp=Ub-su/ss
*      write(*,*)'  ***** s=', su/ss
      p(0,0,0)=Dp
	  
	  
*	  ss=0.
*      aa=0.
*      do k=1,Km
*      do j=1,Jm
*        do i=1,Im
*          hss=1.+yt(j)*smm(i,k)
*          aa=aa+hss*yt(j)*yt1(j)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
*          ss=ss+hss*yt(j)*yt1(j)
*        end do
*      end do
*      end do
*      c=amp/sqrt(aa/ss)
	  
	  
	  
*      write(*,*)'  Ub=',Ub
*      write(*,*)'  Ub_cur=',su/ss
*      write(*,*)'  Dp=',Dp
      
  
      
*
       

      dd1 = 0.0;
      do k=1,Km  
      do j=1,Jm
      do i=1,Im
         v(i,Jm,k)=0.d0
         call div(i,j,k,u,v,w,d,Imax,Jmax)         
c         write(*,*)'  div(',i,',',j,',',k,')=',d
         dd1=max(dd1,abs(d))
      end do
      end do
      end do
*      write(*,*)'  div(u_rp)=',dd1




      call gradp(u,v,w,p,p0,Imax,Jmax)
      
      
      
      
      dd2 = 0.0;
      do k=1,Km  
      do j=1,Jm
      do i=1,Im
         v(i,Jm,k)=0.d0
         call div(i,j,k,u,v,w,d,Imax,Jmax)         
c         write(*,*)'  div(',i,',',j,',',k,')=',d
         dd2=max(dd2,abs(d))
      end do
      end do
      end do
*      write(*,*)'  div(u_pres)=',dd2
*      write(*,*)'  div(u_pres)/div(u_rp)=',dd2/dd1
              
            
            
            
            
            
            
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
*     Dp2=Ub-su/(Im*Km*ss)
 

	  
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
      Dp2=Ub-su/ss

 
*      write(*,*)'  Ub2=',Ub
*      write(*,*)'  Ub_cur2=',su/ss
*      write(*,*)'  Dp2=',Dp2            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
      
      return
      end





*
      subroutine gradp(u,v,w,p,p0,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,p0(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
      Dp=p(0,0,0)
      do k=1,Km
        do j=1,Jm
          do i=1,Im-1
            hss=1.+yt(j)*snm(i,k)
*            u(i,j,k)=u(i,j,k)-((p(i+1,j,k)-p(i,j,k))/hx-Dp)/hss
*            u(i,j,k)=u(i,j,k) + Dp/hss
            u(i,j,k)=u(i,j,k)-(p(i+1,j,k)-p(i,j,k))/(hx*hss) + Dp/hss
     >            -Dp*(p0(i+1,j,k)-p0(i,j,k))/(hx*hss)
          end do
          i=Im
          hss=1.+yt(j)*snm(i,k)
          k10=mod(Km+k-kd-1,Km)+1
*          u(i,j,k)=u(i,j,k)-((p(1,j,k10)-p(i,j,k))/hx-Dp)/hss
*          u(i,j,k)=u(i,j,k) + Dp/hss
          u(i,j,k)=u(i,j,k)-(p(1,j,k10)-p(i,j,k))/(hx*hss) + Dp/hss
     >            -Dp*(p0(1,j,k10)-p0(i,j,k))/(hx*hss)
        end do
      end do
*
      do k=1,Km
        do j=1,Jm-1
          do i=1,Im
            v(i,j,k)=v(i,j,k)-(p(i,j+1,k)-p(i,j,k))/rt1(j)
     >            -Dp*(p0(i,j+1,k)-p0(i,j,k))/rt1(j)
          end do
        end do
      end do
*
      do j=1,Jm
        do i=1,Im
          do k=1,Km-1
            w(i,j,k)=w(i,j,k)-(p(i,j,k+1)-p(i,j,k))/(yt(j)*ht)
     >            -Dp*(p0(i,j,k+1)-p0(i,j,k))/(yt(j)*ht)
          end do
          k=Km
            w(i,j,k)=w(i,j,k)-(p(i,j,1)-p(i,j,k))/(yt(j)*ht)
     >            -Dp*(p0(i,j,1)-p0(i,j,k))/(yt(j)*ht)
        end do
      end do
    
      return
      end
*
      subroutine div(i,j,k,u,v,w,d,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
      i1=mod(Im+i-2,Im)+1
      k1=mod(Km+k-2,Km)+1
      hss=1.+yt(j)*smm(i,k)
      hsr0=1.+rt(j-1)*smm(i,k)
      hsr1=1.+rt(j)*smm(i,k)
      hst0=1.+yt(j)*smn(i,k-1)
      hst1=1.+yt(j)*smn(i,k)
      d=((u(i,j,k)-u(i1,j,k))/hx
     > +(rt(j)*hsr1*v(i,j,k)-rt(j-1)*hsr0*v(i,j-1,k))/(yt(j)*yt1(j))
     > +(hst1*w(i,j,k)-hst0*w(i,j,k1))/(yt(j)*ht))/hss
      if (i.EQ.1) then
      k10=mod(Km+k+kd-1,Km)+1
      d=((u(i,j,k)-u(Im,j,k10))/hx
     > +(rt(j)*hsr1*v(i,j,k)-rt(j-1)*hsr0*v(i,j-1,k))/(yt(j)*yt1(j))
     > +(hst1*w(i,j,k)-hst0*w(i,j,k1))/(yt(j)*ht))/hss
      end if
      return
      end
      