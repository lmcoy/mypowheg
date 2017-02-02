      subroutine abbrvirtQCD1(resnr)
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

      
	WProd124 = Prod(5,2,5,1)
	WProd112 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd118 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd116 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd105 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd119 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd120 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,1,1)
	WProd121 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(4,1,1,1)**2
	WProd117 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd122 = Prod(3,1,2,1)**2*ProdC(3,1,1,1)*ProdC(4,1,2,1)
        WProd123 = Prod(3,1,2,1)*Prod(4,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd103 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd114 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd125 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd131 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd126 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd115 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd127 = Prod(5,2,2,1)*ProdC(5,1,2,1)
        WProd128 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd129 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd130 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd106 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd113 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd107 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
	WProd108 = Prod(5,1,2,1)*ProdC(5,2,2,1)
        WProd109 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd110 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd111 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd104 = ProdC(5,2,5,1)
      
      
	Prp20 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd103)
	Prp21 = -(1.D0/WProd114)
	Prp22 = -(1.D0/WProd115)
      
	IInt111 = myBget(-WProd115,0.D0,0.D0,1)
	LInt109 = myBval(bb1,IInt111)
	IInt116 = myBget(WProd103,0.D0,0.D0,0)
	LInt72 = myBval(bb0,IInt116)
	IInt117 = myBget(-WProd114,0.D0,0.D0,0)
	LInt73 = myBval(bb0,IInt117)
	IInt121 = myBget(-WProd114,0.D0,0.D0,1)
	LInt91 = myBval(bb1,IInt121)
	IInt122 = myBget(-WProd115,0.D0,0.D0,0)
	LInt92 = myBval(bb0,IInt122)
	IInt106 = myCget(WProd112,0.D0,WProd103,0.D0,0.D0,0.D0,1)
	LInt1 = myCval(cc0,IInt106)
	LInt4 = myCval(cc1,IInt106)
	LInt7 = myCval(cc2,IInt106)
	IInt114 = myCget(WProd112,WProd103,0.D0,0.D0,0.D0,0.D0,1)
	LInt2 = myCval(cc0,IInt114)
	LInt5 = myCval(cc1,IInt114)
	LInt8 = myCval(cc2,IInt114)
	IInt115 = myCget(-WProd114,0.D0,WProd103,0.D0,0.D0,0.D0,1)
	LInt3 = myCval(cc0,IInt115)
	LInt6 = myCval(cc1,IInt115)
	LInt9 = myCval(cc2,IInt115)
	IInt110 = myCget(0.D0,-WProd115,0.D0,0.D0,0.D0,0.D0,2)
	LInt103 = myCval(cc12,IInt110)
	LInt106 = myCval(cc2,IInt110)
	LInt96 = myCval(cc00,IInt110)
	LInt99 = myCval(cc1,IInt110)
	IInt109 = myCget(0.D0,-WProd115,0.D0,0.D0,0.D0,0.D0,2)
	LInt102 = myCval(cc12,IInt109)
	LInt105 = myCval(cc2,IInt109)
	LInt93 = myCval(cc0,IInt109)
	LInt95 = myCval(cc00,IInt109)
	LInt98 = myCval(cc1,IInt109)
	IInt120 = myCget(0.D0,-WProd114,0.D0,0.D0,0.D0,0.D0,2)
	LInt78 = myCval(cc00,IInt120)
	LInt81 = myCval(cc1,IInt120)
	LInt83 = myCval(cc11,IInt120)
	LInt86 = myCval(cc12,IInt120)
	LInt89 = myCval(cc2,IInt120)
	IInt118 = myCget(0.D0,WProd103,-WProd114,0.D0,0.D0,0.D0,2)
	LInt74 = myCval(cc0,IInt118)
	LInt76 = myCval(cc00,IInt118)
	LInt79 = myCval(cc1,IInt118)
	LInt84 = myCval(cc12,IInt118)
	LInt87 = myCval(cc2,IInt118)
	LInt90 = myCval(cc22,IInt118)
	IInt119 = myCget(0.D0,-WProd114,0.D0,0.D0,0.D0,0.D0,2)
	LInt75 = myCval(cc0,IInt119)
	LInt77 = myCval(cc00,IInt119)
	LInt80 = myCval(cc1,IInt119)
	LInt82 = myCval(cc11,IInt119)
	LInt85 = myCval(cc12,IInt119)
	LInt88 = myCval(cc2,IInt119)
	IInt108 = myCget(0.D0,WProd103,-WProd115,0.D0,0.D0,0.D0,2)
	LInt100 = myCval(cc1,IInt108)
	LInt101 = myCval(cc11,IInt108)
	LInt104 = myCval(cc12,IInt108)
	LInt107 = myCval(cc2,IInt108)
	LInt108 = myCval(cc22,IInt108)
	LInt94 = myCval(cc0,IInt108)
	LInt97 = myCval(cc00,IInt108)
        IInt107 = myDget(0.D0,WProd112,0.D0,-WProd115,0.D0,WProd103,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt10 = myDval(dd0,IInt107)
	LInt13 = myDval(dd00,IInt107)
	LInt16 = myDval(dd001,IInt107)
	LInt19 = myDval(dd002,IInt107)
	LInt22 = myDval(dd003,IInt107)
	LInt25 = myDval(dd1,IInt107)
	LInt28 = myDval(dd11,IInt107)
	LInt32 = myDval(dd113,IInt107)
	LInt35 = myDval(dd12,IInt107)
	LInt39 = myDval(dd123,IInt107)
	LInt42 = myDval(dd13,IInt107)
	LInt45 = myDval(dd133,IInt107)
	LInt48 = myDval(dd2,IInt107)
	LInt51 = myDval(dd22,IInt107)
	LInt54 = myDval(dd223,IInt107)
	LInt57 = myDval(dd23,IInt107)
	LInt60 = myDval(dd233,IInt107)
	LInt63 = myDval(dd3,IInt107)
	LInt66 = myDval(dd33,IInt107)
	LInt69 = myDval(dd333,IInt107)
        IInt112 = myDget(0.D0,WProd112,WProd103,-WProd114,0.D0,0.D0,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt11 = myDval(dd0,IInt112)
	LInt14 = myDval(dd00,IInt112)
	LInt17 = myDval(dd001,IInt112)
	LInt20 = myDval(dd002,IInt112)
	LInt23 = myDval(dd003,IInt112)
	LInt26 = myDval(dd1,IInt112)
	LInt29 = myDval(dd11,IInt112)
	LInt33 = myDval(dd113,IInt112)
	LInt36 = myDval(dd12,IInt112)
	LInt40 = myDval(dd123,IInt112)
	LInt43 = myDval(dd13,IInt112)
	LInt46 = myDval(dd133,IInt112)
	LInt49 = myDval(dd2,IInt112)
	LInt52 = myDval(dd22,IInt112)
	LInt55 = myDval(dd223,IInt112)
	LInt58 = myDval(dd23,IInt112)
	LInt61 = myDval(dd233,IInt112)
	LInt64 = myDval(dd3,IInt112)
	LInt67 = myDval(dd33,IInt112)
	LInt70 = myDval(dd333,IInt112)
        IInt113 = myDget(0.D0,-WProd114,0.D0,-WProd115,0.D0,WProd103,
     -   0.D0,0.D0,0.D0,0.D0,3)
	LInt12 = myDval(dd0,IInt113)
	LInt15 = myDval(dd00,IInt113)
	LInt18 = myDval(dd001,IInt113)
	LInt21 = myDval(dd002,IInt113)
	LInt24 = myDval(dd003,IInt113)
	LInt27 = myDval(dd1,IInt113)
	LInt30 = myDval(dd11,IInt113)
	LInt31 = myDval(dd112,IInt113)
	LInt34 = myDval(dd113,IInt113)
	LInt37 = myDval(dd12,IInt113)
	LInt38 = myDval(dd122,IInt113)
	LInt41 = myDval(dd123,IInt113)
	LInt44 = myDval(dd13,IInt113)
	LInt47 = myDval(dd133,IInt113)
	LInt50 = myDval(dd2,IInt113)
	LInt53 = myDval(dd22,IInt113)
	LInt56 = myDval(dd223,IInt113)
	LInt59 = myDval(dd23,IInt113)
	LInt62 = myDval(dd233,IInt113)
	LInt65 = myDval(dd3,IInt113)
	LInt68 = myDval(dd33,IInt113)
	LInt71 = myDval(dd333,IInt113)
      
      
      endif

      end


