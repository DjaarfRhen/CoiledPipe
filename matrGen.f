      subroutine matrGen(values, irowPtr, icolInd, innzNum, irowsNum)
      implicit real*8 (a-h,o-z) 
        
      integer irowPtr(*), icolInd(*)
      real*8 values(*)

      common
     >/dimx/hx,Im,kd
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km
     >/sst/snn(0:257,0:257),snm(0:257,0:257)
     >     ,smn(0:257,0:257),smm(0:257,0:257)

      integer iRow

      iRow = 1
      innzNum = 1
      irowsNum = Im*Jm*Km
C innzNum = nnzNum + 1   due to C++/Fortran conversion reasons      
      irowPtr(1) = innzNum
     
      do k = 1,Km
         do j = 1,Jm
            do i = 1,Im

               aP = 0.0
               aW = 0.0
               aE = 0.0
               aS = 0.0
               aN = 0.0
               aB = 0.0
               aT = 0.0

               iP = iRow

*              =================================================================
*              Lame
               hss  = 1. + yt(j)  * smm(i,k)
               hss0 = 1. + yt(j)  * snm(i-1,k)
               hss1 = 1. + yt(j)  * snm(i,k)
               hsr0 = 1. + rt(j-1)* smm(i,k)
               hsr1 = 1. + rt(j)  * smm(i,k)
               hst0 = 1. + yt(j)  * smn(i,k-1)
               hst1 = 1. + yt(j)  * smn(i,k)
               
*              cycle BC
               if (i.GT.1) then
                  iW = iRow - 1
                  tmp = 1.0 / (hx*hss0 * hx*hss)
                  aW = aW + tmp
                  aP = aP - tmp
               else
c                  iW = iRow + (Im - 1)
                 iW = mod(irowsNum + iRow + (Im - 1) - 1 + kd*Im*Jm,
     >irowsNum) + 1
                  tmp = 1.0 / (hx*hss0 * hx*hss)
                  aW = aW + tmp
                  aP = aP - tmp
               endif

               if (i.LT.Im) then
                  iE = iRow + 1
                  tmp = 1.0 / (hx*hss1 * hx*hss)
                  aE = aE + tmp
                  aP = aP - tmp
               else
c                  iE = iRow - (Im - 1)
                iE = mod(irowsNum + iRow - (Im - 1) - 1 - kd*Im*Jm,
     >irowsNum) + 1
                  tmp = 1.0 / (hx*hss1 * hx*hss)
                  aE = aE + tmp
                  aP = aP - tmp
               endif

*              =================================================================
*              wall BC
               if (j.GT.1) then
                  iS = iRow - Im
                  tmp = rt(j-1)*hsr0 / (rt1(j-1) * yt(j)*yt1(j)*hss)
                  aS = aS + tmp
                  aP = aP - tmp
               else
                  iS = -1
               endif 

               if (j.LT.Jm) then
                  iN_ = iRow + Im
                  tmp = rt(j)*hsr1 / (rt1(j) * yt(j)*yt1(j)*hss)
                  aN = aN + tmp
                  aP = aP - tmp
               else
                  iN_ = -1
               endif

*              =================================================================
*              cycle BC
               if (k.GT.1) then
                  iB = iRow - Im*Jm    
                  tmp = hst0 / (ht * ht*(yt(j)**2)*hss)                  
                  aB = aB + tmp
                  aP = aP - tmp
               else
                  iB = (Km - 1)*Im*Jm + (j-1)*Im + i                    
                  tmp = hst0 / (ht * ht*(yt(j)**2)*hss)
                  aB = aB + tmp
                  aP = aP - tmp
               endif

               if (k.LT.Km) then
                  iT = iRow + Im*Jm  
                  tmp = hst1 / (ht * ht*(yt(j)**2)*hss)
                  aT = aT + tmp
                  aP = aP - tmp
               else  
                  iT = (j-1)*Im + i                  
                  tmp = hst1 / (ht * ht*(yt(j)**2)*hss)
                  aT = aT + tmp
                  aP = aP - tmp
               endif


*              =================================================================

		       if (iW.EQ.iP) then
			      aP = aP + aW
			      iW = -1
               endif

		       if (iE.EQ.iP) then
			      aP = aP + aE
			      iE = -1
               endif

		       if (iW.EQ.iE) then
			      aW = aW + aE
			      iE = -1
               endif

*              =================================================================
*              data should be sorted: extra care is needed for cyclic BC

*              order
               if ((iT.GE.1) .AND. (iT.NE.(iRow+Im*Jm))) then
                  icolInd(innzNum) = iT
                  values(innzNum) = aT
                  innzNum = innzNum + 1
               endif

               if ((iB.GE.1) .AND. (iB.EQ.(iRow-Im*Jm))) then
                  icolInd(innzNum) = iB
                  values(innzNum) = aB
                  innzNum = innzNum + 1
               endif

               if (iS.GE.1) then
                  icolInd(innzNum) = iS
                  values(innzNum) = aS
                  innzNum = innzNum + 1
               endif

*              order
               if ((iE.GE.1) .AND. (iE.NE.(iRow+1))) then
                  icolInd(innzNum) = iE
                  values(innzNum) = aE
                  innzNum = innzNum + 1
               endif

               if ((iW.GE.1) .AND. (iW.EQ.(iRow-1))) then
                  icolInd(innzNum) = iW
                  values(innzNum) = aW
                  innzNum = innzNum + 1
               endif

               icolInd(innzNum) = iP
               values(innzNum) = aP
               innzNum = innzNum + 1

               if ((iE.GE.1) .AND. (iE.EQ.(iRow+1))) then
                  icolInd(innzNum) = iE
                  values(innzNum) = aE
                  innzNum = innzNum + 1
               endif

*              order
               if ((iW.GE.1) .AND. (iW.NE.(iRow-1))) then
                  icolInd(innzNum) = iW
                  values(innzNum) = aW
                  innzNum = innzNum + 1
               endif

               if (iN_.GE.1) then
                  icolInd(innzNum) = iN_
                  values(innzNum) = aN
                  innzNum = innzNum + 1
               endif
                
               if ((iT.GE.1) .AND. (iT.EQ.(iRow+Im*Jm))) then
                  icolInd(innzNum) = iT
                  values(innzNum) = aT
                  innzNum = innzNum + 1
               endif

*              order
               if ((iB.GE.1) .AND. (iB.NE.(iRow-Im*Jm))) then
                  icolInd(innzNum) = iB
                  values(innzNum) = aB
                  innzNum = innzNum + 1
               endif

               irowPtr(iRow + 1) = innzNum

               iRow = iRow + 1

            end do
         end do
      end do

      do i = irowPtr(1),(irowPtr(2)-1)    
         if (icolInd(i) .NE. 1) then
            values(i) = 0.0
         else
            values(i) = 1.0/(hx*hr*ht)
            icolInd(i) = Im*Jm*Km
         endif
      enddo

C     for fortran IMSL only !!!
      innzNum = innzNum - 1

      return
      end

