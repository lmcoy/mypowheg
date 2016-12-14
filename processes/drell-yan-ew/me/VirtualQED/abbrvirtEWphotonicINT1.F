      subroutine abbrvirtEWphotonicINT1(resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)


#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV
      double precision dr
      
      integer myBget,myCget,myDget,myEget

#include "amps/looptools.h" 
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 
      if(resnr.eq.0)then

      
	WProd46 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd47 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd53 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd54 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)**2
	WProd49 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd59 = Prod(4,1,1,1)*ProdC(3,1,2,1)
        WProd60 = Prod(2,1,1,1)*Prod(4,1,1,1)*ProdC(2,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd61 = Prod(3,1,1,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd55 = Prod(3,1,2,1)*Prod(4,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
	WProd62 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(3,1,2,1)**2
	WProd45 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd48 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd50 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd51 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,1,1)
	WProd57 = Prod(3,1,1,1)**2*ProdC(3,1,2,1)*ProdC(4,1,1,1)
	WProd56 = Prod(3,1,1,1)*ProdC(4,1,2,1)
	WProd52 = Prod(3,1,2,1)**2*ProdC(3,1,1,1)*ProdC(4,1,2,1)
        WProd58 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,2,1)
      
      
	dr=1D0
	If(dr.eq.0D0) then      
	IInt46 = myCget(WProd46,ME2,ME2,0.D0,0.D0,ME2,0)
	LInt101 = myCval(cc0,IInt46)
	IInt47 = myCget(WProd46,ME2,ME2,0.D0,MZ2r,ME2,0)
	LInt102 = myCval(cc0,IInt47)
	IInt48 = myCget(WProd46,ME2,ME2,MZ2r,0.D0,ME2,0)
	LInt103 = myCval(cc0,IInt48)
        IInt49 = myDget(MU2,WProd46,ME2,-WProd47,MU2,ME2,MU2,0.D0,
     -   0.D0,ME2,2)
	LInt104 = myDval(dd00,IInt49)
	LInt110 = myDval(dd1,IInt49)
	LInt116 = myDval(dd11,IInt49)
	LInt122 = myDval(dd12,IInt49)
	LInt128 = myDval(dd13,IInt49)
	LInt134 = myDval(dd2,IInt49)
	LInt140 = myDval(dd23,IInt49)
	LInt146 = myDval(dd3,IInt49)
	LInt152 = myDval(dd33,IInt49)
        IInt50 = myDget(MU2,WProd46,ME2,-WProd47,MU2,ME2,MU2,0.D0,
     -   MZ2r,ME2,2)
	LInt105 = myDval(dd00,IInt50)
	LInt111 = myDval(dd1,IInt50)
	LInt117 = myDval(dd11,IInt50)
	LInt123 = myDval(dd12,IInt50)
	LInt129 = myDval(dd13,IInt50)
	LInt135 = myDval(dd2,IInt50)
	LInt141 = myDval(dd23,IInt50)
	LInt147 = myDval(dd3,IInt50)
	LInt153 = myDval(dd33,IInt50)
        IInt51 = myDget(MU2,WProd46,ME2,-WProd47,MU2,ME2,MU2,
     -   MZ2r,0.D0,ME2,2)
	LInt106 = myDval(dd00,IInt51)
	LInt112 = myDval(dd1,IInt51)
	LInt118 = myDval(dd11,IInt51)
	LInt124 = myDval(dd12,IInt51)
	LInt130 = myDval(dd13,IInt51)
	LInt136 = myDval(dd2,IInt51)
	LInt142 = myDval(dd23,IInt51)
	LInt148 = myDval(dd3,IInt51)
	LInt154 = myDval(dd33,IInt51)
        IInt52 = myDget(MU2,WProd46,ME2,-WProd48,MU2,ME2,MU2,0.D0,
     -   0.D0,ME2,2)
	LInt107 = myDval(dd00,IInt52)
	LInt113 = myDval(dd1,IInt52)
	LInt119 = myDval(dd11,IInt52)
	LInt125 = myDval(dd12,IInt52)
	LInt131 = myDval(dd13,IInt52)
	LInt137 = myDval(dd2,IInt52)
	LInt143 = myDval(dd23,IInt52)
	LInt149 = myDval(dd3,IInt52)
	LInt155 = myDval(dd33,IInt52)
        IInt53 = myDget(MU2,WProd46,ME2,-WProd48,MU2,ME2,MU2,0.D0,
     -   MZ2r,ME2,2)
	LInt108 = myDval(dd00,IInt53)
	LInt114 = myDval(dd1,IInt53)
	LInt120 = myDval(dd11,IInt53)
	LInt126 = myDval(dd12,IInt53)
	LInt132 = myDval(dd13,IInt53)
	LInt138 = myDval(dd2,IInt53)
	LInt144 = myDval(dd23,IInt53)
	LInt150 = myDval(dd3,IInt53)
	LInt156 = myDval(dd33,IInt53)
        IInt54 = myDget(MU2,WProd46,ME2,-WProd48,MU2,ME2,MU2,
     -   MZ2r,0.D0,ME2,2)
	LInt109 = myDval(dd00,IInt54)
	LInt115 = myDval(dd1,IInt54)
	LInt121 = myDval(dd11,IInt54)
	LInt127 = myDval(dd12,IInt54)
	LInt133 = myDval(dd13,IInt54)
	LInt139 = myDval(dd2,IInt54)
	LInt145 = myDval(dd23,IInt54)
	LInt151 = myDval(dd3,IInt54)
	LInt157 = myDval(dd33,IInt54)
	Else
	IInt46 = myCget(WProd46,0D0,0D0,0.D0,0.D0,0D0,0)
	LInt101 = myCval(cc0,IInt46)
	IInt47 = myCget(WProd46,0D0,0D0,0.D0,MZ2r,0D0,0)
	LInt102 = myCval(cc0,IInt47)
	IInt48 = myCget(WProd46,0D0,0D0,MZ2r,0.D0,0D0,0)
	LInt103 = myCval(cc0,IInt48)
        IInt49 = myDget(0D0,WProd46,0D0,-WProd47,0D0,0D0,0D0,0.D0,
     -   0.D0,0D0,2)
	LInt104 = myDval(dd00,IInt49)
	LInt110 = myDval(dd1,IInt49)
	LInt116 = myDval(dd11,IInt49)
	LInt122 = myDval(dd12,IInt49)
	LInt128 = myDval(dd13,IInt49)
	LInt134 = myDval(dd2,IInt49)
	LInt140 = myDval(dd23,IInt49)
	LInt146 = myDval(dd3,IInt49)
	LInt152 = myDval(dd33,IInt49)
        IInt50 = myDget(0D0,WProd46,0D0,-WProd47,0D0,0D0,0D0,0.D0,
     -   MZ2r,0D0,2)
	LInt105 = myDval(dd00,IInt50)
	LInt111 = myDval(dd1,IInt50)
	LInt117 = myDval(dd11,IInt50)
	LInt123 = myDval(dd12,IInt50)
	LInt129 = myDval(dd13,IInt50)
	LInt135 = myDval(dd2,IInt50)
	LInt141 = myDval(dd23,IInt50)
	LInt147 = myDval(dd3,IInt50)
	LInt153 = myDval(dd33,IInt50)
        IInt51 = myDget(0D0,WProd46,0D0,-WProd47,0D0,0D0,0D0,
     -   MZ2r,0.D0,0D0,2)
	LInt106 = myDval(dd00,IInt51)
	LInt112 = myDval(dd1,IInt51)
	LInt118 = myDval(dd11,IInt51)
	LInt124 = myDval(dd12,IInt51)
	LInt130 = myDval(dd13,IInt51)
	LInt136 = myDval(dd2,IInt51)
	LInt142 = myDval(dd23,IInt51)
	LInt148 = myDval(dd3,IInt51)
	LInt154 = myDval(dd33,IInt51)
        IInt52 = myDget(0D0,WProd46,0D0,-WProd48,0D0,0D0,0D0,0.D0,
     -   0.D0,0D0,2)
	LInt107 = myDval(dd00,IInt52)
	LInt113 = myDval(dd1,IInt52)
	LInt119 = myDval(dd11,IInt52)
	LInt125 = myDval(dd12,IInt52)
	LInt131 = myDval(dd13,IInt52)
	LInt137 = myDval(dd2,IInt52)
	LInt143 = myDval(dd23,IInt52)
	LInt149 = myDval(dd3,IInt52)
	LInt155 = myDval(dd33,IInt52)
        IInt53 = myDget(0D0,WProd46,0D0,-WProd48,0D0,0D0,0D0,0.D0,
     -   MZ2r,0D0,2)
	LInt108 = myDval(dd00,IInt53)
	LInt114 = myDval(dd1,IInt53)
	LInt120 = myDval(dd11,IInt53)
	LInt126 = myDval(dd12,IInt53)
	LInt132 = myDval(dd13,IInt53)
	LInt138 = myDval(dd2,IInt53)
	LInt144 = myDval(dd23,IInt53)
	LInt150 = myDval(dd3,IInt53)
	LInt156 = myDval(dd33,IInt53)
        IInt54 = myDget(0D0,WProd46,0D0,-WProd48,0D0,0D0,0D0,
     -   MZ2r,0.D0,0D0,2)
	LInt109 = myDval(dd00,IInt54)
	LInt115 = myDval(dd1,IInt54)
	LInt121 = myDval(dd11,IInt54)
	LInt127 = myDval(dd12,IInt54)
	LInt133 = myDval(dd13,IInt54)
	LInt139 = myDval(dd2,IInt54)
	LInt145 = myDval(dd23,IInt54)
	LInt151 = myDval(dd3,IInt54)
	LInt157 = myDval(dd33,IInt54)
	Endif
      
      
      endif

      end


