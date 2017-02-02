
        subroutine setprods(K)
        implicit none 

*calculates all spinor products from the momentum matrix "k" 
       
#include "params.h"

        integer i,j,r,t
        double precision K(np,0:3)
        double precision absk(np)
        double complex ic
        double complex spinor(np,2,2),sum,lambda1(np,2)
        integer sign,pro

#include "amps/ampcommon.h"

        ic=(0D0,1D0)
        
          do t=1,2
             do j=1,np
        	do r=1,2           
                   do i=1,np  
                      prod(i,r,j,t)=0D0
                      prodc(i,r,j,t)=0D0
                      prod(j,t,i,r)=0D0
                      prodc(j,t,i,r)=0D0
                   enddo
        	enddo
             enddo
          enddo

          do i=1,np
             do j=1,2
        	do r=1,2
                   spinor(i,j,r)=0D0
        	enddo
             enddo
          enddo 

          do i=1,nexternal(1)
	     absk(i)=dsqrt(K(i,1)**2+K(i,2)**2+K(i,3)**2)
             lambda1(i,1)=dsqrt(k(i,0)+absk(i))
             If(spinorind(i,1) .eq.2) then
        	If(procmass(1,1,i) .lt. 1D-12)then
                   lambda1(i,2)=lambda1(i,1)
        	else
                   lambda1(i,2)=dsqrt(dabs(k(i,0)-absk(i)))
        	endif
             endif
          enddo
          do i=1,nexternal(1)
             If(k(i,1)**2+k(i,2)**2.EQ.0D0) Then
		if(k(i,3).lt.0D0)then
        	   spinor(i,1,1)=(0D0,0D0)
        	   spinor(i,1,2)=(1D0,0D0)
        	   If(spinorind(i,1).eq.2) then
                      spinor(i,2,1)=(1D0,0D0)
                      spinor(i,2,2)=(0D0,0D0)
        	   endif
        	elseif(k(i,3).gt.0D0)then
        	   spinor(i,1,1)=(1D0,0D0)
        	   spinor(i,1,2)=(0D0,0D0)
        	   If(spinorind(i,1).eq.2) then
                      spinor(i,2,1)=(0D0,0D0)
                      spinor(i,2,2)=(-1D0,0D0)
        	   endif
        	endif
             Else
        	spinor(i,1,1)=(k(i,1)-ic*k(i,2))
     &  	  *dsqrt((absk(i)+k(i,3))/2D0/absk(i)/(k(i,1)**2+k(i,2)**2))
        	spinor(i,1,2)=dsqrt((k(i,1)**2+k(i,2)**2)/2D0/absk(i)/(absk(i)+k(i,3)))
        	If(spinorind(i,1).eq.2) then
        	   spinor(i,2,1)=spinor(i,1,2)
        	   spinor(i,2,2)=-dconjg(spinor(i,1,1))
        	endif
             Endif
          enddo 


          do i=1,nexternal(1)
             do j=i+1,nexternal(1)
        	do r=1,spinorind(i,1)
        	   do t=1,spinorind(j,1)
        	      prod(i,r,j,t)=lambda1(i,r)*lambda1(j,t)
     &  			   *(spinor(i,r,1)*spinor(j,t,2)
     &  			    -spinor(i,r,2)*spinor(j,t,1))
        	      prodc(i,r,j,t)=dconjg(prod(i,r,j,t))
        	      prod(j,t,i,r)=-prod(i,r,j,t)
        	      prodc(j,t,i,r)=-prodc(i,r,j,t)
        	  enddo       
        	enddo
             enddo
          enddo

          do i=1,nexternal(1)
             if(spinorind(i,1).eq.2)then
        	if(procmass(1,1,i).eq.0d0) then
        	   prod(i,1,i,2)=-2D0*k(i,0)
        	   prod(i,2,i,1)=2D0*k(i,0)
        	   prodc(i,1,i,2)=-2D0*k(i,0)
        	   prodc(i,2,i,1)=2D0*k(i,0) 
        	else
        	   prod(i,1,i,2)=-procmass(1,1,i)
        	   prod(i,2,i,1)=procmass(1,1,i)
        	   prodc(i,1,i,2)=-procmass(1,1,i)
        	   prodc(i,2,i,1)=procmass(1,1,i)
               endif
             endif
          enddo
              
      end
       

      function Prp(p2,M2,W)
      implicit none 
*calculates propagator with momentum "p2", squared mass m2 and width W

      integer i1,i2 
       
      double complex Prp,ic
      double precision M2,W,p2

      ic=(0D0,1D0)
      Prp=1D0/(p2-M2+ic*dsqrt(M2)*W)
 
      end  

 
