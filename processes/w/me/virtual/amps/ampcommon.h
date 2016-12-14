      double precision pi,sqrt2,el
      double precision alfa,alfa2
      double precision mwr,mw,mw2r,gw,mzr,mz,mz2r,gz
      double precision cwr,cw2r,swr,sw2r   
      double complex mwc,mw2c,mzc,mz2c
      double complex cwc,cw2c,swc,sw2c   
      double precision me,mm,ml,mu,mc,mt,md,ms,mb
      double precision me2,mm2,ml2,mu2,mc2,mt2,md2,ms2,mb2
      double precision mh,mh2,gh
      double precision MLE(3),MLE2(3),MQU(3),MQU2(3),MQD(3),MQD2(3)
      common /ampsmpara/
     &    MLE, MQU, MQD, MLE2, MQU2, MQD2, Pi, Sqrt2,
     &    EL, alfa, alfa2, 
     &    MWr, MW, MZr, MZ, MH, GW, GZ, GH, 
     &    ME, MM, ML, MU, MD, MC, MS, MB, MT,
     &    MW2r, MZ2r, MH2, ME2, MM2, ML2, MU2, MD2, MC2, MS2, 
     &    MB2, MT2, CWr, SWr, CW2r, SW2r, CWc, SWc, CW2c, SW2c, MWc, MZc, MW2c, MZ2c

      double complex prod(np,1:2,np,1:2),prodc(np,1:2,np,1:2)
      common/ampsetproducts/prod,prodc

      double precision mureg,deltareg,lambda,epsfac,epsfac2,nflavour
      common/ampregularizationparameters/mureg,deltareg,lambda,epsfac,epsfac2,nflavour

      integer nexternal(maxg)
      common/ampgeneral/nexternal

      double precision procmass(maxg,nproc,np)
      common/ampprocmass/procmass

      integer spinorind(np,maxg)
      data spinorind/1,1,1,1,1,
     &                 1,2,1,1,2/

      double precision Cfac(ncol,ncol,nproc,maxg)
        data cfac/1.732050807568877D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            1.732050807568877D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,2.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            1.732050807568877D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,2.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            1.732050807568877D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     &            0.D0,0.D0,2.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/

      character*5 subprocs(np,nproc,maxg)
        data subprocs/" up  ","~do  "," ve  ","~el  ","   ",
     &                "   ","   ","   ","   ","   ",
     &                "   ","   ","   ","   ","   ",
     &                "   ","   ","   ","   ","   ",
     &                "   ","   ","   ","   ","   ",
     &                "   ","   ","   ","   ","   ",
     &                " up  ","~do  "," ve  ","~el  "," ga  ",
     &                " up  ","~do  "," ve  ","~el  "," gl  ",
     &                " up  "," ga  "," ve  ","~el  "," do  ",
     &                " up  "," gl  "," ve  ","~el  "," do  ",
     &                "~do  "," ga  "," ve  ","~el  ","~up  ",
     &                "~do  "," gl  "," ve  ","~el  ","~up  "/
      
