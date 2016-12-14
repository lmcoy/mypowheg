* looptools_nfac.F contains the definitions of the loop integrals B0c, C0c, D0c, E0c which are introduced by "CreateESPAAmplitudes" and contribute to the non-factorizable corrections. 



      double complex function b0c(ki,null,m2bar)
      implicit none
      double complex m2bar,logc,ki,ic
      double precision m2,null,RRe
      double precision getlambda

#include "params.h"
#include "amps/ampcommon.h"

      ic=(0D0,1D0)
      m2=RRe(m2bar)
    
      b0c=ic*pi**2*ki/m2*(logc(dsqrt(getlambda())*dsqrt(m2)/-ki)+1D0)

      end

      double complex function c0c(switch,p1s,p2s,p1p2s,null,m1s,m2s)
      implicit none

      integer switch
      double complex logc,ic,LLi2
      double complex m1s,m2s,M1bar,M2bar
      double complex ki,kj,rij,p1s,p2s,p1p2s
      double precision that,eps,RRe
      double precision m,null,M1,M2     
      double precision wij,sijbar,lambdakin
      double precision gammamass,getlambda

#include "params.h"
#include "amps/ampcommon.h"
 
      eps=1d-20
      ic=(0D0,1D0)
      gammamass=dsqrt(getlambda())
      if(switch.eq.1)then
         M1bar=m1s
         M1=dsqrt(RRe(M1bar))
         m=dsqrt(RRe(m2s))
         ki=p1s 
         c0c=-ic*pi**2/M1**2*(2D0*dlog(m/M1)*logc(-ki/gammamass/M1)
     &               +dlog(m/M1)**2+pi**2/6D0)
         return
      elseif(switch.eq.2)then
         M2bar=m2s
         M2=dsqrt(RRe(M2bar))
         m=dsqrt(RRe(m1s)) 
         that=p2s-M2**2
         ki=p1p2s
         
         c0c=ic*pi**2/that*(logc(-m*M2/that+ic*eps)
     &     *(logc(ki/that)+logc(-ki/gammamass**2)
     &        +dlog(m/M2))+pi**2/6D0)
         return
      elseif(switch.eq.3)then
         M1bar=m1s
         M1=dsqrt(RRe(M1bar))
         M2bar=m2s
         M2=dsqrt(RRe(M2bar))
         
         wij=dsqrt(lambdakin(p2s,M1**2,M2**2))
         sijbar=p2s-M1**2-M2**2
         rij=1D0/2D0/M1/M2*(-sijbar+wij)*(1D0-ic*eps/wij)
         ki=p1s
         kj=p1p2s

         c0c=ic*pi**2/wij*(LLi2(ki*M2/kj/M1,1D0/rij)
     &               -LLi2(ki*M2/kj/M1,rij)+LLi2(rij,rij)+logc(rij)**2
     &               +2D0*logc(rij)*logc(-kj/M2/gammamass))

         return
      endif
      end

      double complex function d0c(switch,p1s,p2s,p3s,p4s,p1p2s,p2p3s,
     &                             null,m1s,m2s,m3s)
      implicit none
      integer switch,tau
      double complex ic
      double complex m1s,m2s,m3s,M1bar,M2bar
      double complex ki,kj,LLi2,Li2sum,rij,logc
      double complex p1s,p2s,p3s,p4s,p1p2s,p2p3s
      double precision t,that,RRe,null
      double precision m,mp,M1,M2
      double precision shatij,shatji,sbarij,wij,lambdakin,eps,Li2
      double precision gammamass,getlambda


#include "params.h"
#include "amps/ampcommon.h"

      eps=1d-20
      ic=(0D0,1D0)
      gammamass=dsqrt(getlambda())

      if(switch.eq.1) then
         m=dsqrt(RRe(m3s))
         mp=dsqrt(RRe(m1s))
         M1bar=m2s
         M1=dsqrt(RRe(M1bar))

         t=p2p3s         
         that=p2s-M1**2
         ki=p1p2s

         d0c=ic*pi**2/ki/t*(2D0*logc(-m*mp/t+ic*eps)*logc(-gammamass*M1/ki)
     &                 -logc(-mp*M1/that+ic*eps)**2
     &                 -dlog(m/M1)**2-pi**2/3D0-Li2(dcmplx(1D0-that/t))) 
         return
      elseif(switch.eq.2) then
         M1bar=m2s  
         M2bar=m1s
         m=dsqrt(RRe(m3s))
         M1=dsqrt(RRe(M1bar))
         M2=dsqrt(RRe(M2bar))
         Kj=p1s
         Ki=p1p2s
         shatji=p2p3s-M2**2
         sbarij=p2s+M1**2+M2**2
         wij=dsqrt(lambdakin(p2s,M1**2,M2**2))
         rij=1D0/2D0/M1/M2*(-sbarij+wij)*(1D0-ic*eps/wij)

         Li2sum=0D0
         do tau=1,2
            Li2sum=Li2sum+LLi2(Ki*M2/Kj*M1,rij**tau)
     &                   -LLi2(-M1*M2/shatji+ic*eps,rij**tau)
         enddo
         d0c=ic*pi**2/(ki*shatji+kj*M1**2)
     &       *(Li2sum-2D0*LLi2(Ki*M2/Kj/M1,-shatji/M1/M2-ic*eps)
     &        -dlog(m**2/M1**2)
     &            *(logc(ki*M2/kj/M1)+logc(-shatji/M1/M2-ic*eps)))
         return

      endif 
      end

      double complex function e0c(p1s,p2s,p3s,p4s,p5s,
     &                            p1p2s,p2p3s,p3p4s,p4p5s,p5p1s,
     &                            null,m1s,M1bar,M2bar,m2s)
      implicit none
      double complex M1bar,M2bar,ki,kj
      double complex d0x(0:4),d0c,ic
      double complex p1s,p2s,p3s,p4s,p5s
      double complex p1p2s,p2p3s,p3p4s,p4p5s,p5p1s
      double precision Mi,Mj,m,mp,m1s,m2s,null,RRe
      double precision sbarij,shatij,shatji,sij
      double precision detY,detY0,detY2,detY3,lambdakin,Li2
      
#include "looptools.h"
#include "params.h"
#include "amps/ampcommon.h"


      Mi=dsqrt(RRe(M1bar))
      Mj=dsqrt(RRe(M2bar))
      m=dsqrt(RRe(m1s))
      mp=dsqrt(RRe(m2s))
      
      sij=p5p1s
      sbarij=p3s-Mi**2-Mj**2
      shatij=p3p4s-Mi**2
      shatji=p2p3s-Mj**2
      ki=p4p5s
      kj=p1p2s
   
      detY0=lambdakin(sij*sbarij,shatij*shatji,Mi**2*Mj**2)
     &        +4D0*sij*Mi**2*Mj**2*(shatij+shatji-sij)
      detY2=sij*(ki*(shatij*shatji-sij*sbarij+Mi**2*Mj**2)
     &           +2D0*Kj*Mi**2*(shatij-sij))
      detY3=sij*(kj*(shatji*shatij-sij*sbarij+Mj**2*Mi**2)
     &           +2D0*Ki*Mj**2*(shatji-sij))
      detY=2D0*sij*(Ki*Kj*(sij*sbarij-shatij*shatji-Mi**2*Mj**2)
     &           +Ki**2*Mj**2*(sij-shatji)+Kj**2*Mi**2*(sij-shatij))

      D0x(0)=D0i(dd0,p2s,p3s,p4s,p5p1s,p2p3s,p3p4s,m1s,Mj**2,Mi**2,m2s)
      D0x(1)=D0c(2,p4p5s,p3s,p2s,p1s,p1p2s,p2p3s,
     &                             null,M2bar,M1bar,m1s)
      D0x(2)=D0c(1,p1s,p2p3s,p4s,p5s,p4p5s,p5p1s,
     &                             lambda,m1s,M2bar,m2s)
      D0x(3)=D0c(1,p5s,p3p4s,p2s,p1s,p1p2s,p5p1s,
     &                             lambda,m2s,M1bar,m1s)
      D0x(4)=D0c(2,p1p2s,p3s,p4s,p5s,p4p5s,p3p4s,
     &                             null,M1bar,M2bar,m2s)


      e0c=ic*pi**2*(
     &    -detY0/detY*D0x(0)
     &    -detY3/detY/kj/sij
     &       *((ki*shatji+kj*Mi**2)*D0x(1)+kj*sij*D0x(3))
     &    -detY2/detY/ki/sij
     &       *((kj*shatij+ki*Mj**2)*D0x(4)+ki*sij*D0x(2))
     &    -shatij/ki/sij*D0x(4)-shatji/kj/sij*D0x(1))
      end
      
      function RRe(z)
      implicit none
      double complex z
      double precision RRe 

      RRe=(z+Conjg(z))/2D0
      end 

      double complex function LLi2(x,y)
      implicit none
      double complex x,y,logc
      double precision Li2

      LLi2=Li2(1-x*y)+(logc(x*y)-logc(x)-logc(y))*logc(1D0-x*y)

      end

      double precision function lambdakin(x,y,z)
      implicit none
      double precision x,y,z

      lambdakin=x**2+y**2+z**2-2D0*x*y-2D0*x*z-2D0*y*z

      end

      double complex function logc(z)
      implicit none
      integer sig
      double complex z,ic
      double precision rez,imz,argz,absz

#include "params.h"
#include "amps/ampcommon.h"


      ic=(0D0,1D0)
      rez=1D0/2D0*(z+dConjg(z))
      imz=-ic/2D0*(z-dConjg(z))

      absz=dsqrt(rez**2+imz**2)
      if(rez.gt.0D0) then
        argz=datan(imz/rez)
      else
        argz=sig(imz)*pi-datan(imz/abs(rez))
      endif
        logc=dlog(absz)+ic*argz
      end

      integer function sig(x)
      implicit none
      double precision x

       if(x.gt.0D0) then
        sig=1
       else
        sig=-1
       endif

      end
