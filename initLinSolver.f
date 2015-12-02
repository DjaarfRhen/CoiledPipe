      subroutine initLinSolver(Imax,Jmax,N,NZ)
      implicit real*8 (a-h,o-z)
      common
     > a1(256),a2(256),a3(256),a4(256)
     >,b1(256),b2(256),b3(256),b4(256)
     >,c1(256),c2(256),c3(256),c4(256)
     >,ap(256),bp(256),cp(256)
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/rl/rlx(0:128),rlt(0:128)
     >/pry/apy(128),bpy(128),cpy(128)
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)
     >/lll/lx,lt
     >/pressMatr/irowInd(1e7),icolInd(1e7),irowPtr(1e5),values(1e7)
     >,irowsNum
     >/linSolver/iparam(6),rparam(5),IRFAC(1e7),JCFAC(1e7)
     >,IPVT(1e5),JPVT(1e5),FACT(1e7),ipath,NL,NFAC
*
      Im2=Im/2
      Km2=Km/2
      cik=4.d0/(Im2*Km2)


***********************************************************
***********************************************************
***********************************************************

C Poisson matrix generation 
      call matrGen(values, irowPtr, icolInd, innzNum, irowsNum)


C CHECK: save to file
      open(unit=1, file='log/values.txt')
      do i = 1, innzNum
         write(1,*) values(i)
      enddo
      close(1)
C CHECK: save to file
      open(unit=1, file='log/irowPtr.txt')
      do i = 1, irowsNum+1
         write(1,*) irowPtr(i)
      enddo
      close(1)
C CHECK: save to file
      open(unit=1, file='log/icolInd.txt')
      do i = 1, innzNum
         write(1,*) icolInd(i)
      enddo
      close(1)





      
      
C Poisson matrix conversion
      do iRow = 1,irowsNum
         do l = irowPtr(iRow),(irowPtr(iRow+1)-1)
            irowInd(l) = iRow
         end do
      end do

      
      
C CHECK: save to file
      open(unit=1, file='log/irowInd.txt')
      do i = 1, innzNum
         write(1,*) irowInd(i)
      enddo
      close(1)      
      
 
 
 
      
      
      ipath = 1
      CALL DL4LXG(iparam, rparam) 
      iparam(5) = 1e8           
      NFAC = 1e7
      
C Perform LU factorization
      CALL DLFTXG (irowsNum, innzNum, values, irowInd, icolInd, iparam,
     > rparam, NFAC, NL, FACT, IRFAC, JCFAC, IPVT, JPVT)

*      write(*,*)'stop'
*      stop
      
***********************************************************
***********************************************************
***********************************************************
      
      return
      end