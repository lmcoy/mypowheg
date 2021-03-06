      subroutine abbrvirtEW1(resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV
      integer dr

      integer myBget,myCget,myDget,myEget

#include "looptools.h" 
#include "amps/ampvars.h" 
#include "amps/ampcommon.h" 
#include "amps/ampcom.h" 
#include "regcommon.h"
      if(resnr.eq.0)then

      
	WProd6 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd8 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd9 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd5 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd7 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd10 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,1,1)
      
      
	Prp3 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd6)
      
	If(dr.eq.1) Then
	
		IInt31 = myBget(WProd6,0.D0,0D0,0)
	LInt127 = myBval(bb0,IInt31)
	IInt33 = myBget(WProd6,0D0,0D0,0)
	LInt129 = myBval(bb0,IInt33)
	IInt34 = myBget(WProd6,MW2r,MZ2r,0)
	LInt130 = myBval(bb0,IInt34)
	IInt47 = myCget(WProd6,0D0,0.D0,MW2r,0.D0,0D0,0)
	LInt76 = myCval(cc0,IInt47)
	IInt48 = myCget(WProd6,0D0,0.D0,MW2r,MZ2r,0D0,0)
	LInt77 = myCval(cc0,IInt48)
	IInt49 = myCget(WProd6,0D0,0.D0,MZ2r,MW2r,0.D0,0)
	LInt78 = myCval(cc0,IInt49)
      

	do Gen1=1,3
	IInt50(Gen1) = myBget(WProd6,MQD2(Gen1),MQU2(Gen1),0)
	LInt173(Gen1) = myBval(bb0,IInt50(Gen1))
	IInt51(Gen1) = myBget(WProd6,0.D0,MLE2(Gen1),1)
	LInt174(Gen1) = myBval(bb1,IInt51(Gen1))
	IInt52(Gen1) = myBget(WProd6,MQU2(Gen1),MQD2(Gen1),1)
	LInt175(Gen1) = myBval(bb1,IInt52(Gen1))
        enddo

	IInt32 = myBget(WProd6,0.D0,MW2r,1)
	LInt128 = myBval(bb0,IInt32)
	LInt171 = myBval(bb1,IInt32)
	IInt45 = myBget(WProd6,MH2,MW2r,2)
	LInt168 = myBval(bb0,IInt45)
	LInt169 = myBval(bb00,IInt45)
	IInt46 = myBget(WProd6,MZ2r,MW2r,2)
	LInt170 = myBval(bb00,IInt46)
	LInt172 = myBval(bb1,IInt46)
	IInt38 = myCget(0.D0,WProd6,0D0,0.D0,MZ2r,MW2r,2)
	LInt134 = myCval(cc00,IInt38)
	LInt144 = myCval(cc1,IInt38)
	LInt154 = myCval(cc2,IInt38)
	IInt39 = myCget(0.D0,WProd6,0D0,0D0,MW2r,0.D0,2)
	LInt135 = myCval(cc00,IInt39)
	LInt145 = myCval(cc1,IInt39)
	LInt155 = myCval(cc2,IInt39)
	IInt40 = myCget(0.D0,WProd6,0D0,0D0,MW2r,MZ2r,2)
	LInt136 = myCval(cc00,IInt40)
	LInt146 = myCval(cc1,IInt40)
	LInt156 = myCval(cc2,IInt40)
	IInt41 = myCget(0D0,WProd6,0D0,0D0,MW2r,0.D0,2)
	LInt139 = myCval(cc00,IInt41)
	LInt149 = myCval(cc1,IInt41)
	LInt159 = myCval(cc2,IInt41)
	IInt42 = myCget(0D0,WProd6,0D0,0D0,MW2r,MZ2r,2)
	LInt140 = myCval(cc00,IInt42)
	LInt150 = myCval(cc1,IInt42)
	LInt160 = myCval(cc2,IInt42)
	IInt43 = myCget(0D0,WProd6,0D0,0D0,0.D0,MW2r,2)
	LInt141 = myCval(cc00,IInt43)
	LInt151 = myCval(cc1,IInt43)
	LInt161 = myCval(cc2,IInt43)
	IInt44 = myCget(0D0,WProd6,0D0,0D0,MZ2r,MW2r,2)
	LInt142 = myCval(cc00,IInt44)
	LInt152 = myCval(cc1,IInt44)
	LInt162 = myCval(cc2,IInt44)
	IInt35 = myCget(0.D0,WProd6,0D0,MZ2r,0.D0,0D0,2)
	LInt131 = myCval(cc0,IInt35)
	LInt137 = myCval(cc00,IInt35)
	LInt147 = myCval(cc1,IInt35)
	LInt157 = myCval(cc2,IInt35)
	IInt36 = myCget(0D0,WProd6,0D0,0.D0,0D0,0D0,2)
	LInt132 = myCval(cc0,IInt36)
	LInt138 = myCval(cc00,IInt36)
	LInt148 = myCval(cc1,IInt36)
	LInt158 = myCval(cc2,IInt36)
	IInt37 = myCget(0D0,WProd6,0D0,MZ2r,0D0,0D0,2)
	LInt133 = myCval(cc0,IInt37)
	LInt143 = myCval(cc00,IInt37)
	LInt153 = myCval(cc1,IInt37)
	LInt163 = myCval(cc2,IInt37)
	IInt25 = myDget(0D0,WProd6,0D0,-WProd8,0D0,0.D0,0D0,MW2r,0.D0,0D0,2)
	LInt100 = myDval(dd13,IInt25)
	LInt106 = myDval(dd2,IInt25)
	LInt112 = myDval(dd23,IInt25)
	LInt118 = myDval(dd3,IInt25)
	LInt124 = myDval(dd33,IInt25)
	LInt82 = myDval(dd00,IInt25)
	LInt88 = myDval(dd1,IInt25)
	LInt94 = myDval(dd12,IInt25)
        IInt26 = myDget(0D0,WProd6,0D0,-WProd8,0D0,0.D0,0D0,MW2r,
     -   MZ2r,0D0,2)
	LInt101 = myDval(dd13,IInt26)
	LInt107 = myDval(dd2,IInt26)
	LInt113 = myDval(dd23,IInt26)
	LInt119 = myDval(dd3,IInt26)
	LInt125 = myDval(dd33,IInt26)
	LInt83 = myDval(dd00,IInt26)
	LInt89 = myDval(dd1,IInt26)
	LInt95 = myDval(dd12,IInt26)
        IInt27 = myDget(0D0,WProd6,0D0,-WProd8,0D0,0.D0,0D0,MZ2r,
     -   MW2r,0.D0,2)
	LInt102 = myDval(dd13,IInt27)
	LInt108 = myDval(dd2,IInt27)
	LInt114 = myDval(dd23,IInt27)
	LInt120 = myDval(dd3,IInt27)
	LInt126 = myDval(dd33,IInt27)
	LInt84 = myDval(dd00,IInt27)
	LInt90 = myDval(dd1,IInt27)
	LInt96 = myDval(dd12,IInt27)
        IInt28 = myDget(0D0,WProd6,0.D0,-WProd7,0D0,0D0,0D0,MW2r,
     -   MZ2r,0.D0,2)
	LInt103 = myDval(dd2,IInt28)
	LInt109 = myDval(dd23,IInt28)
	LInt115 = myDval(dd3,IInt28)
	LInt121 = myDval(dd33,IInt28)
	LInt79 = myDval(dd00,IInt28)
	LInt85 = myDval(dd1,IInt28)
	LInt91 = myDval(dd12,IInt28)
	LInt97 = myDval(dd13,IInt28)
	IInt29 = myDget(0D0,WProd6,0.D0,-WProd7,0D0,0D0,0D0,0.D0,MW2r,0D0,2)
	LInt104 = myDval(dd2,IInt29)
	LInt110 = myDval(dd23,IInt29)
	LInt116 = myDval(dd3,IInt29)
	LInt122 = myDval(dd33,IInt29)
	LInt80 = myDval(dd00,IInt29)
	LInt86 = myDval(dd1,IInt29)
	LInt92 = myDval(dd12,IInt29)
	LInt98 = myDval(dd13,IInt29)
        IInt30 = myDget(0D0,WProd6,0.D0,-WProd7,0D0,0D0,0D0,MZ2r,
     -   MW2r,0D0,2)
	LInt105 = myDval(dd2,IInt30)
	LInt111 = myDval(dd23,IInt30)
	LInt117 = myDval(dd3,IInt30)
	LInt123 = myDval(dd33,IInt30)
	LInt81 = myDval(dd00,IInt30)
	LInt87 = myDval(dd1,IInt30)
	LInt93 = myDval(dd12,IInt30)
	LInt99 = myDval(dd13,IInt30)     
	
	Else
	
	IInt31 = myBget(WProd6,0.D0,ME2,0)
	LInt127 = myBval(bb0,IInt31)
	IInt33 = myBget(WProd6,MD2,MU2,0)
	LInt129 = myBval(bb0,IInt33)
	IInt34 = myBget(WProd6,MW2r,MZ2r,0)
	LInt130 = myBval(bb0,IInt34)
	IInt47 = myCget(WProd6,ME2,0.D0,MW2r,0.D0,ME2,0)
	LInt76 = myCval(cc0,IInt47)
	IInt48 = myCget(WProd6,ME2,0.D0,MW2r,MZ2r,ME2,0)
	LInt77 = myCval(cc0,IInt48)
	IInt49 = myCget(WProd6,ME2,0.D0,MZ2r,MW2r,0.D0,0)
	LInt78 = myCval(cc0,IInt49)
      

	do Gen1=1,3
	IInt50(Gen1) = myBget(WProd6,MQD2(Gen1),MQU2(Gen1),0)
	LInt173(Gen1) = myBval(bb0,IInt50(Gen1))
	IInt51(Gen1) = myBget(WProd6,0.D0,MLE2(Gen1),1)
	LInt174(Gen1) = myBval(bb1,IInt51(Gen1))
	IInt52(Gen1) = myBget(WProd6,MQU2(Gen1),MQD2(Gen1),1)
	LInt175(Gen1) = myBval(bb1,IInt52(Gen1))
        enddo

	IInt32 = myBget(WProd6,0.D0,MW2r,1)
	LInt128 = myBval(bb0,IInt32)
	LInt171 = myBval(bb1,IInt32)
	IInt45 = myBget(WProd6,MH2,MW2r,2)
	LInt168 = myBval(bb0,IInt45)
	LInt169 = myBval(bb00,IInt45)
	IInt46 = myBget(WProd6,MZ2r,MW2r,2)
	LInt170 = myBval(bb00,IInt46)
	LInt172 = myBval(bb1,IInt46)
	IInt38 = myCget(0.D0,WProd6,ME2,0.D0,MZ2r,MW2r,2)
	LInt134 = myCval(cc00,IInt38)
	LInt144 = myCval(cc1,IInt38)
	LInt154 = myCval(cc2,IInt38)
	IInt39 = myCget(0.D0,WProd6,ME2,ME2,MW2r,0.D0,2)
	LInt135 = myCval(cc00,IInt39)
	LInt145 = myCval(cc1,IInt39)
	LInt155 = myCval(cc2,IInt39)
	IInt40 = myCget(0.D0,WProd6,ME2,ME2,MW2r,MZ2r,2)
	LInt136 = myCval(cc00,IInt40)
	LInt146 = myCval(cc1,IInt40)
	LInt156 = myCval(cc2,IInt40)
	IInt41 = myCget(MU2,WProd6,MD2,MD2,MW2r,0.D0,2)
	LInt139 = myCval(cc00,IInt41)
	LInt149 = myCval(cc1,IInt41)
	LInt159 = myCval(cc2,IInt41)
	IInt42 = myCget(MU2,WProd6,MD2,MD2,MW2r,MZ2r,2)
	LInt140 = myCval(cc00,IInt42)
	LInt150 = myCval(cc1,IInt42)
	LInt160 = myCval(cc2,IInt42)
	IInt43 = myCget(MU2,WProd6,MD2,MU2,0.D0,MW2r,2)
	LInt141 = myCval(cc00,IInt43)
	LInt151 = myCval(cc1,IInt43)
	LInt161 = myCval(cc2,IInt43)
	IInt44 = myCget(MU2,WProd6,MD2,MU2,MZ2r,MW2r,2)
	LInt142 = myCval(cc00,IInt44)
	LInt152 = myCval(cc1,IInt44)
	LInt162 = myCval(cc2,IInt44)
	IInt35 = myCget(0.D0,WProd6,ME2,MZ2r,0.D0,ME2,2)
	LInt131 = myCval(cc0,IInt35)
	LInt137 = myCval(cc00,IInt35)
	LInt147 = myCval(cc1,IInt35)
	LInt157 = myCval(cc2,IInt35)
	IInt36 = myCget(MU2,WProd6,MD2,0.D0,MU2,MD2,2)
	LInt132 = myCval(cc0,IInt36)
	LInt138 = myCval(cc00,IInt36)
	LInt148 = myCval(cc1,IInt36)
	LInt158 = myCval(cc2,IInt36)
	IInt37 = myCget(MU2,WProd6,MD2,MZ2r,MU2,MD2,2)
	LInt133 = myCval(cc0,IInt37)
	LInt143 = myCval(cc00,IInt37)
	LInt153 = myCval(cc1,IInt37)
	LInt163 = myCval(cc2,IInt37)
	IInt25 = myDget(MU2,WProd6,ME2,-WProd8,MD2,0.D0,MD2,MW2r,0.D0,ME2,2)
	LInt100 = myDval(dd13,IInt25)
	LInt106 = myDval(dd2,IInt25)
	LInt112 = myDval(dd23,IInt25)
	LInt118 = myDval(dd3,IInt25)
	LInt124 = myDval(dd33,IInt25)
	LInt82 = myDval(dd00,IInt25)
	LInt88 = myDval(dd1,IInt25)
	LInt94 = myDval(dd12,IInt25)
        IInt26 = myDget(MU2,WProd6,ME2,-WProd8,MD2,0.D0,MD2,MW2r,
     -   MZ2r,ME2,2)
	LInt101 = myDval(dd13,IInt26)
	LInt107 = myDval(dd2,IInt26)
	LInt113 = myDval(dd23,IInt26)
	LInt119 = myDval(dd3,IInt26)
	LInt125 = myDval(dd33,IInt26)
	LInt83 = myDval(dd00,IInt26)
	LInt89 = myDval(dd1,IInt26)
	LInt95 = myDval(dd12,IInt26)
        IInt27 = myDget(MU2,WProd6,ME2,-WProd8,MD2,0.D0,MU2,MZ2r,
     -   MW2r,0.D0,2)
	LInt102 = myDval(dd13,IInt27)
	LInt108 = myDval(dd2,IInt27)
	LInt114 = myDval(dd23,IInt27)
	LInt120 = myDval(dd3,IInt27)
	LInt126 = myDval(dd33,IInt27)
	LInt84 = myDval(dd00,IInt27)
	LInt90 = myDval(dd1,IInt27)
	LInt96 = myDval(dd12,IInt27)
        IInt28 = myDget(MU2,WProd6,0.D0,-WProd7,MD2,ME2,MD2,MW2r,
     -   MZ2r,0.D0,2)
	LInt103 = myDval(dd2,IInt28)
	LInt109 = myDval(dd23,IInt28)
	LInt115 = myDval(dd3,IInt28)
	LInt121 = myDval(dd33,IInt28)
	LInt79 = myDval(dd00,IInt28)
	LInt85 = myDval(dd1,IInt28)
	LInt91 = myDval(dd12,IInt28)
	LInt97 = myDval(dd13,IInt28)
	IInt29 = myDget(MU2,WProd6,0.D0,-WProd7,MD2,ME2,MU2,0.D0,MW2r,ME2,2)
	LInt104 = myDval(dd2,IInt29)
	LInt110 = myDval(dd23,IInt29)
	LInt116 = myDval(dd3,IInt29)
	LInt122 = myDval(dd33,IInt29)
	LInt80 = myDval(dd00,IInt29)
	LInt86 = myDval(dd1,IInt29)
	LInt92 = myDval(dd12,IInt29)
	LInt98 = myDval(dd13,IInt29)
        IInt30 = myDget(MU2,WProd6,0.D0,-WProd7,MD2,ME2,MU2,MZ2r,
     -   MW2r,ME2,2)
	LInt105 = myDval(dd2,IInt30)
	LInt111 = myDval(dd23,IInt30)
	LInt117 = myDval(dd3,IInt30)
	LInt123 = myDval(dd33,IInt30)
	LInt81 = myDval(dd00,IInt30)
	LInt87 = myDval(dd1,IInt30)
	LInt93 = myDval(dd12,IInt30)
	LInt99 = myDval(dd13,IInt30)
      Endif      
      
      endif

      end


