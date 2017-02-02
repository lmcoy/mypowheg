      subroutine abbrvirtQCD2(resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV

      integer myBget,myCget,myDget,myEget

#include "amps/looptools.h" 
#include "amps/ampvars.h" 
#include "amps/ampcommon.h" 
#include "amps/ampcom.h" 
      if(resnr.eq.0)then

      
	WProd132 = Prod(2,2,2,1)
	WProd139 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd148 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd162 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd152 = ProdC(2,2,2,1)
	WProd140 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd142 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd145 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd158 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd141 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd133 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd150 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd146 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd153 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(4,1,1,1))
        WProd159 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
        WProd160 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,1,1))
	WProd143 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd134 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd151 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd147 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd135 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd154 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
        WProd136 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd137 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
	WProd155 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd161 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd156 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd157 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd144 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd138 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd149 = Prod(5,1,2,1)*ProdC(5,1,2,1)
      
      
	Prp23 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd138)
	Prp26 = 1.D0/WProd139
	Prp24 = -(1.D0/WProd143)
	Prp25 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd144)
	Prp27 = -(1.D0/WProd149)
      
	IInt123 = myBget(-WProd138,0.D0,0.D0,0)
	LInt110 = myBval(bb0,IInt123)
	IInt129 = myBget(WProd139,0.D0,0.D0,0)
	LInt158 = myBval(bb0,IInt129)
	IInt130 = myBget(WProd144,0.D0,0.D0,0)
	LInt159 = myBval(bb0,IInt130)
	IInt133 = myBget(WProd139,0.D0,0.D0,1)
	LInt171 = myBval(bb1,IInt133)
	IInt134 = myBget(-WProd149,0.D0,0.D0,0)
	LInt172 = myBval(bb0,IInt134)
	IInt137 = myBget(-WProd149,0.D0,0.D0,1)
	LInt183 = myBval(bb1,IInt137)
	IInt125 = myCget(-WProd138,0.D0,WProd144,0.D0,0.D0,0.D0,1)
	LInt115 = myCval(cc0,IInt125)
	LInt117 = myCval(cc1,IInt125)
	LInt119 = myCval(cc2,IInt125)
	IInt126 = myCget(WProd139,WProd144,0.D0,0.D0,0.D0,0.D0,1)
	LInt116 = myCval(cc0,IInt126)
	LInt118 = myCval(cc1,IInt126)
	LInt120 = myCval(cc2,IInt126)
	IInt124 = myCget(0.D0,-WProd138,0.D0,0.D0,0.D0,0.D0,2)
	LInt111 = myCval(cc0,IInt124)
	LInt112 = myCval(cc00,IInt124)
	LInt113 = myCval(cc1,IInt124)
	LInt114 = myCval(cc2,IInt124)
	IInt136 = myCget(0.D0,-WProd149,0.D0,0.D0,0.D0,0.D0,2)
	LInt174 = myCval(cc00,IInt136)
	LInt178 = myCval(cc12,IInt136)
	LInt181 = myCval(cc22,IInt136)
	LInt184 = myCval(cc2,IInt136)
	IInt132 = myCget(0.D0,WProd139,0.D0,0.D0,0.D0,0.D0,2)
	LInt162 = myCval(cc00,IInt132)
	LInt164 = myCval(cc1,IInt132)
	LInt165 = myCval(cc11,IInt132)
	LInt167 = myCval(cc12,IInt132)
	LInt169 = myCval(cc2,IInt132)
	IInt131 = myCget(0.D0,WProd144,WProd139,0.D0,0.D0,0.D0,2)
	LInt160 = myCval(cc0,IInt131)
	LInt161 = myCval(cc00,IInt131)
	LInt163 = myCval(cc1,IInt131)
	LInt166 = myCval(cc12,IInt131)
	LInt168 = myCval(cc2,IInt131)
	LInt170 = myCval(cc22,IInt131)
	IInt135 = myCget(0.D0,WProd144,-WProd149,0.D0,0.D0,0.D0,2)
	LInt173 = myCval(cc0,IInt135)
	LInt175 = myCval(cc00,IInt135)
	LInt176 = myCval(cc1,IInt135)
	LInt177 = myCval(cc11,IInt135)
	LInt179 = myCval(cc12,IInt135)
	LInt180 = myCval(cc2,IInt135)
	LInt182 = myCval(cc22,IInt135)
        IInt128 = myDget(0.D0,WProd139,WProd144,-WProd138,0.D0,0.D0,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt123 = myDval(dd00,IInt128)
	LInt125 = myDval(dd001,IInt128)
	LInt127 = myDval(dd002,IInt128)
	LInt129 = myDval(dd003,IInt128)
	LInt131 = myDval(dd1,IInt128)
	LInt133 = myDval(dd11,IInt128)
	LInt134 = myDval(dd112,IInt128)
	LInt137 = myDval(dd12,IInt128)
	LInt138 = myDval(dd122,IInt128)
	LInt140 = myDval(dd123,IInt128)
	LInt142 = myDval(dd13,IInt128)
	LInt145 = myDval(dd2,IInt128)
	LInt148 = myDval(dd223,IInt128)
	LInt150 = myDval(dd23,IInt128)
	LInt152 = myDval(dd233,IInt128)
	LInt154 = myDval(dd3,IInt128)
	LInt156 = myDval(dd33,IInt128)
        IInt127 = myDget(0.D0,-WProd138,0.D0,-WProd149,0.D0,WProd144,
     -   0.D0,0.D0,0.D0,0.D0,3)
	LInt121 = myDval(dd0,IInt127)
	LInt122 = myDval(dd00,IInt127)
	LInt124 = myDval(dd001,IInt127)
	LInt126 = myDval(dd002,IInt127)
	LInt128 = myDval(dd003,IInt127)
	LInt130 = myDval(dd1,IInt127)
	LInt132 = myDval(dd11,IInt127)
	LInt135 = myDval(dd113,IInt127)
	LInt136 = myDval(dd12,IInt127)
	LInt139 = myDval(dd123,IInt127)
	LInt141 = myDval(dd13,IInt127)
	LInt143 = myDval(dd133,IInt127)
	LInt144 = myDval(dd2,IInt127)
	LInt146 = myDval(dd22,IInt127)
	LInt147 = myDval(dd223,IInt127)
	LInt149 = myDval(dd23,IInt127)
	LInt151 = myDval(dd233,IInt127)
	LInt153 = myDval(dd3,IInt127)
	LInt155 = myDval(dd33,IInt127)
	LInt157 = myDval(dd333,IInt127)
      
      
      endif

      end


