*
      function rrt(x,i)
      implicit real*8 (a-h,o-z)
      common/set/y0,sety,al,tha,iset
      if(iset.ne.0) goto 2
1     continue
      t=1./sqrt(sety)
      al=dlog(t+sqrt(t**2-1.))
      tha=tanh(al)
      iset=1
2     continue
      if(i.eq.0)rrt=tanh(al*x)/tha
      if(i.eq.1)rrt=al/cosh(al*x)**2/tha
      return
      end