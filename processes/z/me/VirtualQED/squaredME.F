#define COLILib
C#define DDLib 

******************************************************************************************+

      Subroutine init_all(alfa_in,mw_in,gw_in,mz_in,gz_in,
     & me_in,mm_in,ml_in,mu_in,mc_in,mt_in,md_in,ms_in,mb_in,mh_in,gh_in,inveps)
      Implicit None     
*
* Call before first call to Msquared to provide necessary input 
* for common blocks and initialize loop routines
*
#include "params.h"
#include "amps/ampcommon.h"

      double precision alfa_in
      double precision mw_in,gw_in,mz_in,gz_in
      double precision me_in,mm_in,ml_in
      double precision mu_in,mc_in,mt_in
      double precision md_in,ms_in,mb_in
      double precision mh_in,gh_in
      double precision inveps

      double precision nflavour_out,mureg_out,Qren2
      common /ct_common/nflavour_out,mureg_out,Qren2
      
      integer i1,gen,procnr
      character*5 name(np)

      PI=3.1415926535897932384D0
      sqrt2=sqrt(2D0)

      nexternal(1)=np-1

*     renormalization scale
      mureg=1D0

*     deltareg=2/epsilon-gamma_E+log(4*pi)
      deltareg=0D0

      epsfac=inveps
      epsfac2=inveps*inveps

*     photon mass
c     This is the dimReg mu in the virtual amplitudes, they occur in Log(mu/s)
c      lambda=sqrt(Qren2)
      lambda = 1d0
c     lambda=MZ_in

*     number of light flavors
      nflavour=5D0
      
      nflavour_out=nflavour
      mureg_out=mureg

      MWr=MW_in
      MZr=MZ_in
      MH=MH_in
      GW=GW_in
      GZ=GZ_in
      GH=GH_in

      ME=ME_in
      MM=MM_in
      ML=ML_in

      MU=MU_in
      MC=MC_in
      MT=MT_in

      MD=MD_in
      MS=MS_in
      MB=MB_in

      alfa=alfa_in
      
      cwr = mwr/mzr
      cw2r = cwr**2
      sw2r = 1d0-cw2r
      swr=dsqrt(sw2r) 

      MW=MWr
      MZ=MZr

      MW2r=MWr**2
      MZ2r=MZr**2
      MH2=MH**2

      ME2=ME**2
      MM2=MM**2
      ML2=ML**2

      MU2=MU**2
      MC2=MC**2
      MT2=MT**2

      MD2=MD**2
      MS2=MS**2
      MB2=MB**2

      MLE(1)=ME
      MLE(2)=MM
      MLE(3)=ML

      MLE2(1)=ME2
      MLE2(2)=MM2
      MLE2(3)=ML2

      MQU(1)=MU
      MQU(2)=MC
      MQU(3)=MT

      MQU2(1)=MU2
      MQU2(2)=MC2
      MQU2(3)=MT2

      MQD(1)=MD
      MQD(2)=MS
      MQD(3)=MB

      MQD2(1)=MD2
      MQD2(2)=MS2
      MQD2(3)=MB2
      
      alfa2=alfa**2

      el=dsqrt(4D0*pi*alfa)

      MW2c=MW2r-(0D0,1D0)*GW*MWr
      MZ2c=MZ2r-(0D0,1D0)*GZ*MZr

      MWc=sqrt(MW2c)
      MZc=sqrt(MZ2c)

      CW2c=MW2c/MZ2c
      SW2c=1D0-CW2c
      CWc=sqrt(CW2c)
      SWc=sqrt(SW2c)

      do gen=1,1
       do procnr=1,nproctree
        do i1=1,np						   
           name(i1)=subprocs(i1,procnr,gen)		   
           procmass(gen,procnr,i1)=0D0 			   
           if(name(i1).eq." up  ".or.name(i1).eq."~up  ")then     
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." do  ".or.name(i1).eq."~do  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ch  ".or.name(i1).eq."~ch  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." st  ".or.name(i1).eq."~st  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." bo  ".or.name(i1).eq."~bo  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." to  ".or.name(i1).eq."~to  ")then 
              procmass(gen,procnr,i1)=MT
           elseif(name(i1).eq." el  ".or.name(i1).eq."~el  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ve  ".or.name(i1).eq."~ve  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." mu  ".or.name(i1).eq."~mu  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ta  ".or.name(i1).eq."~ta  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." W+  ".or.name(i1).eq."~W+  ")then 
              procmass(gen,procnr,i1)=MWr
           elseif(name(i1).eq." Z0  ")then			   
              procmass(gen,procnr,i1)=MZr
           elseif(name(i1).eq." H0  ")then			   
              procmass(gen,procnr,i1)=MH
           elseif(name(i1).eq." gl  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ga  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq."     ")then			   
              procmass(gen,procnr,i1)=0D0 		   
           else
	      Print *,"Particle identification failed in mass assignment"
	      Stop
	   endif						   
        enddo  						   
       enddo
      enddo

*initialize amplitudes
      call init_amps

      End
******************************************************************************************+
      Subroutine getdZR(CT)
      Implicit None
#include "params.h"
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h"
      double precision CT

c      call abbr_const()
      CT = dZR1
c      write(*,*) dZR1
      End

******************************************************************************************+

      Subroutine set_alphas(Alfas_in,Qren2_in)
*
* Call to provide Alphas and the squared renormalization scale
* (at each phase-space point)
*
      double precision Alfas,Alfas2,GS,Alfas_in,Pi
      double precision nflavour,mureg,Qren2,Qren2_in

      common /my_alfas_common/Alfas,Alfas2,GS
      common /ct_common/nflavour,mureg,Qren2 

      PI=3.1415926535897932384D0
      
      Alfas=Alfas_in
      Alfas2=Alfas_in**2
      gs=dsqrt(4D0*pi*Alfas)
      
      Qren2=Qren2_in

      call abbr_const_mu(Qren2_in)

      End


      subroutine abbr_const_mu(Qren2_in)
#include "params.h" 
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 
      double precision Qren2_in
      double precision dZfL1photonic(4,1,1),dZfR1photonic(4,1,1),dr
      lambda=sqrt(Qren2_in)
      Cpl180 = 1/Pi*Log(lambda**2)
        dZfL1dimreg(3,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183
        dZfL1dimreg(4,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183
        dZfR1dimreg(3,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183

        dZfR1dimreg(4,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183

      dZfL1photonic(1,1,1)=0D0

	dr=1
	If(dr.Eq.0) Then
	dZfL1photonic(2,1,1)=2*Cpl15*LInt363 + Cpl15*LInt364 + 
     -   0.25D0*Cpl15 + 0.5D0*(Cpl15*LInt361)

	dZfL1photonic(3,1,1)=0.1111111111111111D0*Cpl15 + 
     -   0.2222222222222222D0*(Cpl15*LInt366) +  
     -   0.8888888888888888D0*(Cpl15*LInt368) + 
     -   0.4444444444444444D0*(Cpl15*LInt369)

	dZfL1photonic(4,1,1)=0.027777777777777776D0*Cpl15 + 
     -   0.05555555555555555D0*(Cpl15*LInt370) + 
     -   0.2222222222222222D0*(Cpl15*LInt373) + 
     -   0.1111111111111111D0*(Cpl15*LInt374)
        Else
        dZfL1photonic(2,1,1)=
     -  Alfa*(0.25D0*Cpl180 - 0.25D0*Cpl181 + 
     -        0.25D0*Cpl182 - 0.25D0*Cpl183)

        dZfL1photonic(3,1,1)=
     -  4D0/9D0*Alfa*(0.25D0*Cpl180 - 0.25D0*Cpl181 + 
     -                0.25D0*Cpl182 - 0.25D0*Cpl183)

        dZfL1photonic(4,1,1)=
     -  1D0/9D0*Alfa*(0.25D0*Cpl180 - 0.25D0*Cpl181 + 
     -                0.25D0*Cpl182 - 0.25D0*Cpl183)
        Endif

        dZfR1photonic(1,1,1)=dZfL1photonic(1,1,1)	
	dZfR1photonic(2,1,1)=dZfL1photonic(2,1,1)     
        dZfR1photonic(3,1,1)=dZfL1photonic(3,1,1)
        dZfR1photonic(4,1,1)=dZfL1photonic(4,1,1)
     
	Cpl134 = Alfa*Pi*dZfR1photonic(3,1,1)
	Cpl135 = Alfa*Pi*(SW2c*dZfR1photonic(3,1,1))/CW2c
	Cpl136 = Alfa*Pi*dZfR1photonic(3,1,1)/CW2c
	Cpl137 = Alfa*Pi*dZfL1photonic(3,1,1)
	Cpl138 = Alfa*Pi*dZfL1photonic(3,1,1)/CW2c
	Cpl139 = Alfa*Pi*(SW2c*dZfL1photonic(3,1,1))/CW2c
	Cpl140 = Alfa*Pi*dZfL1photonic(3,1,1)/(CW2c*SW2c)
	Cpl141 = Alfa*Pi*dZfR1photonic(4,1,1)
	Cpl142 = Alfa*Pi*(SW2c*dZfR1photonic(4,1,1))/CW2c
	Cpl143 = Alfa*Pi*dZfR1photonic(4,1,1)/CW2c
	Cpl144 = Alfa*Pi*dZfL1photonic(4,1,1)
	Cpl145 = Alfa*Pi*dZfL1photonic(4,1,1)/CW2c
	Cpl146 = Alfa*Pi*(SW2c*dZfL1photonic(4,1,1))/CW2c
	Cpl147 = Alfa*Pi*dZfL1photonic(4,1,1)/(CW2c*SW2c)
	Cpl148 = Alfa*Pi*dZfR1photonic(2,1,1)
	Cpl149 = Alfa*Pi*(SW2c*dZfR1photonic(2,1,1))/CW2c
	Cpl150 = Alfa*Pi*dZfL1photonic(2,1,1)
	Cpl151 = Alfa*Pi*dZfL1photonic(2,1,1)/CW2c
	Cpl152 = Alfa*Pi*(SW2c*dZfL1photonic(2,1,1))/CW2c
	Cpl153 = Alfa*Pi*dZfR1photonic(2,1,1)/CW2c
	Cpl154 = Alfa*Pi*dZfL1photonic(2,1,1)/(CW2c*SW2c)

#ifdef COLILib
      call setmphoton2(lambda**2)       
#endif
#ifdef DDLib
*Initialize Stefan's libraries
      Call DDsetparam(deltareg,mureg,epsfac2,epsfac,
     &			      lambda**2,smallm2)
      cacc = 1d-10
      dacc = 1d-10
      mode34 = 0
      mode5  = 0
      mode6  = 0
      outlevel = 0
      call DDsetmode(cacc,dacc,mode34,mode5,mode6,outlevel)
#endif
      End

******************************************************************************************+

      subroutine Msquared(msqu0,msqu1,k,flavour,virtual_type)
      implicit none                   
* Calculates:
* msqu0(nproc): Born matrix element squared
* msqu1(nproc): virtual 1-loop QCD corrected matrix element squared
*       nproc=1: ubar u
*       nproc=2: dbar d
* for momentum configuration k(np,0:3) 
* where np is the number of particles at born + 1
* np=1: anti-quark
* np=2: quark
* np=3: anti-lepton
* np=4: lepton
* np=5: empty
* You have to call "init_all" first!
* You have to provide alphas and the renormalization scale 
* with "set_alphas"!
*
* virtual_type = 0 don't calculate virtual corrections
* virtual_type = 1 QCD corrections
* virtual_type = 2 EW corrections

c flavour = 1 down
c flavour = 2 up
c calcluate amplitude seperatly

#include "params.h"

      integer procnr,helnr,nc1
      integer resnr,nrproc,flavour      
      double precision k(np,0:3),msqu0(nproc),msqu1(nproc)
      double complex M0(nproc,ncol,nhel),M1(nproc,ncol,nhel)
      double complex M1virtQCD1(nproc,ncol,nhel)
      double complex M1virtQCD2(nproc,ncol,nhel)
      double complex M1ampsfacctQCD1(nproc,ncol,nhel)
      double complex M1ampsfacctQCD2(nproc,ncol,nhel)
      double complex M1virtEWphotonicISR1(nproc,ncol,nhel)
      double complex M1virtEWphotonicISR2(nproc,ncol,nhel)
      double complex M1virtEWphotonicFSR1(nproc,ncol,nhel)
      double complex M1virtEWphotonicFSR2(nproc,ncol,nhel)
      double complex M1virtEWphotonicINT1(nproc,ncol,nhel)
      double complex M1virtEWphotonicINT2(nproc,ncol,nhel)
      double complex M1virtEWnonphotonic1(nproc,ncol,nhel)
      double complex M1virtEWnonphotonic2(nproc,ncol,nhel)
      double complex M1ampsfacctEWphotonicISR1(nproc,ncol,nhel)
      double complex M1ampsfacctEWphotonicISR2(nproc,ncol,nhel)
      double complex M1ampsfacctEWphotonicFSR1(nproc,ncol,nhel)
      double complex M1ampsfacctEWphotonicFSR2(nproc,ncol,nhel)
      double complex M1ampsfacctEWnonphotonic1(nproc,ncol,nhel)
      double complex M1ampsfacctEWnonphotonic2(nproc,ncol,nhel)
      integer virtual_type

      resnr=0
      
      do procnr=1,nproctree
         msqu0(procnr)=0D0
         msqu1(procnr)=0D0
      enddo
 
      Do procnr=1,nproctree
	 do helnr=1,nheltree
	    do nc1=1,ncol
               M0(procnr,nc1,helnr)=(0D0,0D0)
               M1(procnr,nc1,helnr)=(0D0,0D0)
            enddo
         enddo
      enddo
      	 
      If(virtual_type.eq.1) Then
	call setzero(M1virtQCD1)
	call setzero(M1virtQCD2)
	call setzero(M1ampsfacctQCD1)
	call setzero(M1ampsfacctQCD2)
      Else IF(virtual_type.eq.2) Then
	call setzero(M1virtEWphotonicISR1)
	call setzero(M1virtEWphotonicISR2)
	call setzero(M1virtEWphotonicFSR1)
	call setzero(M1virtEWphotonicFSR2)
	call setzero(M1virtEWphotonicINT1)
	call setzero(M1virtEWphotonicINT2)
	call setzero(M1virtEWnonphotonic1)
	call setzero(M1virtEWnonphotonic2)
	call setzero(M1ampsfacctEWphotonicISR1)
	call setzero(M1ampsfacctEWphotonicISR2)
	call setzero(M1ampsfacctEWphotonicFSR1)
	call setzero(M1ampsfacctEWphotonicFSR2)
	call setzero(M1ampsfacctEWnonphotonic1)
	call setzero(M1ampsfacctEWnonphotonic2)
      Endif
       
#ifdef COLILib
c      call cacheinit(0)
#endif

      call setprods(k)

      call abbrborn1(resnr)
      call born1(M0,resnr)

      call abbrborn2(resnr)
      call born2(M0,resnr)

      If(virtual_type.eq.1) Then
c	QCD Corrections
	call abbrvirtQCD1(resnr)
	call virtQCD1(M1virtQCD1,resnr)
	call sumamplitudes(M1,M1virtQCD1)

	call abbrvirtQCD2(resnr)
	call virtQCD2(M1virtQCD2,resnr)
	call sumamplitudes(M1,M1virtQCD2)

	call abbrampsfacctQCD1(resnr)
	call ampsfacctQCD1(M1ampsfacctQCD1,resnr)
	call sumamplitudes(M1,M1ampsfacctQCD1)

	call abbrampsfacctQCD2(resnr)
	call ampsfacctQCD2(M1ampsfacctQCD2,resnr)
	call sumamplitudes(M1,M1ampsfacctQCD2)

      Else If (virtual_type.eq.2) Then
c	Electroweak corrections

c 	down amplitude
	If( flavour.eq.1 .or. flavour.eq.3 .or. flavour.eq.5) Then
		call abbrvirtEWphotonicISR2(resnr)
		call virtEWphotonicISR2(M1virtEWphotonicISR2,resnr)
		call sumamplitudes(M1,M1virtEWphotonicISR2)

		call abbrvirtEWphotonicFSR2(resnr)
		call virtEWphotonicFSR2(M1virtEWphotonicFSR2,resnr)
		call sumamplitudes(M1,M1virtEWphotonicFSR2)

		call abbrvirtEWphotonicINT2(resnr)
		call virtEWphotonicINT2(M1virtEWphotonicINT2,resnr)
		call sumamplitudes(M1,M1virtEWphotonicINT2)

		call abbrvirtEWnonphotonic2(resnr)
		call virtEWnonphotonic2(M1virtEWnonphotonic2,resnr)
		call sumamplitudes(M1,M1virtEWnonphotonic2)

		call abbrampsfacctEWphotonicISR2(resnr)
		call ampsfacctEWphotonicISR2(M1ampsfacctEWphotonicISR2,resnr)
		call sumamplitudes(M1,M1ampsfacctEWphotonicISR2)

		call abbrampsfacctEWphotonicFSR2(resnr)
		call ampsfacctEWphotonicFSR2(M1ampsfacctEWphotonicFSR2,resnr)
		call sumamplitudes(M1,M1ampsfacctEWphotonicFSR2)

		call abbrampsfacctEWnonphotonic2(resnr)
		call ampsfacctEWnonphotonic2(M1ampsfacctEWnonphotonic2,resnr)
		call sumamplitudes(M1,M1ampsfacctEWnonphotonic2)
	Endif

c	up amplitude
	if( flavour.eq.2 .or. flavour.eq.4 .or. flavour.eq.6) Then
		call abbrvirtEWphotonicISR1(resnr)
		call virtEWphotonicISR1(M1virtEWphotonicISR1,resnr)
		call sumamplitudes(M1,M1virtEWphotonicISR1)

		call abbrvirtEWphotonicFSR1(resnr)
		call virtEWphotonicFSR1(M1virtEWphotonicFSR1,resnr)
		call sumamplitudes(M1,M1virtEWphotonicFSR1)

		call abbrvirtEWphotonicINT1(resnr)
		call virtEWphotonicINT1(M1virtEWphotonicINT1,resnr)
		call sumamplitudes(M1,M1virtEWphotonicINT1)

		call abbrvirtEWnonphotonic1(resnr)
		call virtEWnonphotonic1(M1virtEWnonphotonic1,resnr)
		call sumamplitudes(M1,M1virtEWnonphotonic1)

		call abbrampsfacctEWphotonicISR1(resnr)
		call ampsfacctEWphotonicISR1(M1ampsfacctEWphotonicISR1,resnr)
		call sumamplitudes(M1,M1ampsfacctEWphotonicISR1)

		call abbrampsfacctEWphotonicFSR1(resnr)
		call ampsfacctEWphotonicFSR1(M1ampsfacctEWphotonicFSR1,resnr)
		call sumamplitudes(M1,M1ampsfacctEWphotonicFSR1)

		call abbrampsfacctEWnonphotonic1(resnr)
		call ampsfacctEWnonphotonic1(M1ampsfacctEWnonphotonic1,resnr)
		call sumamplitudes(M1,M1ampsfacctEWnonphotonic1)
	Endif
	
      Endif

      call squareamp(msqu0,msqu1,M0,M1)

c	write(*,*) "B up:", msqu0(1)/36d0,"B down:", msqu0(2)/36d0
c	write(*,*) "V up:", msqu1(1)/36d0,"V down:", msqu1(2)/36d0
      end  

****************************************************************************************************************

      subroutine squareamp(msqu0,msqu1,amp0,amp1)
      implicit none

#include "params.h"

      integer nrproc
      double precision msqu0(nproc),msqu1(nproc)
      double complex amp0(nproc,ncol,nhel),amp1(nproc,ncol,nhel)
     
      integer procnr,helnr,nc1,nc2
      double complex res0,res1
      double precision re

#include "amps/ampcommon.h"

      Do procnr=1,nproctree							    
         do helnr=1,nheltree							    
            do nc1=1,ncol							    
               res0=0D0 							    
               res1=0D0 							    
               do nc2=1,ncol							  
        	  res0=res0+Cfac(nc2,nc1,procnr,1)   
     &  		     *amp0(procnr,nc2,helnr)				    
        	  res1=res1+Cfac(nc2,nc1,procnr,1)   
     &  		     *amp1(procnr,nc2,helnr)				    
               enddo								    
               Msqu0(procnr)=Msqu0(procnr)+res0*conjg(res0)			    
               Msqu1(procnr)=Msqu1(procnr)+2D0*Re(res0,res1)			    
            enddo								    
         enddo  								    
      enddo									    
      
      end

**************************************************************************************************** 
 
      subroutine sumamplitudes(Amp,M)
      implicit none

#include "params.h"

      double complex Amp(nproc,ncol,nhel),M(nproc,ncol,nhel)
      integer procnr,helnr,colnr

       Do procnr=1,nproctree
	  do helnr=1,nheltree
	     do colnr=1,ncol
                amp(procnr,colnr,helnr)
     &           =amp(procnr,colnr,helnr)+M(procnr,colnr,helnr)
             enddo
          enddo 
       enddo
      end

******************************************************************************************+

      subroutine setzero(M)       
      implicit none

#include "params.h"

      integer procnr,helnr,nc1
      double complex M(nproc,ncol,nhel)

       Do procnr=1,nproctree
	  do helnr=1,nheltree
	     do nc1=1,ncol
                M(procnr,nc1,helnr)=(0D0,0D0)
             enddo
          enddo
       enddo
      end  

*************************************************************************************************

      subroutine init_amps
      implicit none
*Intializes the amplitudes

#include "params.h"
#include "amps/ampcommon.h"

      integer gen,ii,mode34,mode5,mode6,outlevel,iimax
      double precision getlambda,smallm2(10),cacc,dacc
      double complex commass 

      common/smallm2_common/smallm2,iimax
      
#ifdef COLILib
*Initialize Ansgar's libraries
      call defcoli(1)
      call cacheon
      call countinit
#endif

*Initialize UV and infrared regulators

#ifdef COLILib
      call setmudim2(mureg)
      call setmphoton2(lambda**2)
      call setdeltair_coli(epsfac,epsfac2)
      call setdeltauv_coli(deltareg)            
#endif
      
*set small masses for Ansgar's and Stefan's loop integrals
      do ii=1,10
         smallm2(ii)=0D0
      enddo
      ii=0
      commass=MU2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MU2

      commass=MD2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MD2

C      commass=MS2
#ifdef COLILib
C      call setcminf2(commass)
#endif
C      ii=ii+1
C      smallm2(ii)=MS2

C      commass=MC2
#ifdef COLILib
C      call setcminf2(commass)
#endif
C      ii=ii+1
C      smallm2(ii)=MC2

      commass=MB2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MB2

      commass=ME2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=ME2

      commass=ML2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=ML2

      commass=MM2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MM2
            
      iimax=ii

      call compress_smallm2(smallm2,iimax)

#ifdef DDLib
*Initialize Stefan's libraries
      Call DDsetparam(deltareg,mureg,epsfac2,epsfac,
     &			      lambda**2,smallm2)
      cacc = 1d-10
      dacc = 1d-10
      mode34 = 0
      mode5  = 0
      mode6  = 0
      outlevel = 0
      call DDsetmode(cacc,dacc,mode34,mode5,mode6,outlevel)
#endif

      end
 
******************************************************************************************+

      Subroutine compress_smallm2(smallm2,iimax)
      Implicit None
      
      double precision smallm2(10)
      integer iimax,ii,jj,kk,keep_index

      keep_index=0
      do ii=1,10
 99     continue
        If(smallm2(ii).NE.0) Then      
	  If(ii.GT.0) Then
	    do jj=1,ii-1
	       If(smallm2(ii).EQ.smallm2(jj)) Then
	         do kk=ii,9
	           smallm2(kk)=smallm2(kk+1)
	         enddo
	         smallm2(10)=0D0
C	         Print ('(A11,I2,A42)'),"small mass ",ii+keep_index,
C     &                   " is eliminated (equal to other small mass)"
		 keep_index=keep_index+1
		 goto 99
	       Endif
	    enddo
	  Endif
	Else 
	  iimax=ii-1
	  Return
	Endif
      Enddo
      
      End      
******************************************************************************************+
      function delta(s1,s2)
      implicit none
      integer delta,s1,s2
c*Kronecker delta
      if (s1.eq.s2) then
        delta=1D0
      else
         delta=0D0
      endif
      end


      function Re(a,b)
      implicit none
c*real part
      double precision Re
      double complex a,b

      Re=(conjg(a)*b+a*conjg(b))/2D0
      end

