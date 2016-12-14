      subroutine abbrvirtEWphotonicINT2(resnr)
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

      
	WProd64 = Prod(2,1,1,1)*ProdC(2,1,1,1)
c	Print *, "WProd64", Wprod64
	WProd65 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd71 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd72 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)**2
	WProd67 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd77 = Prod(4,1,1,1)*ProdC(3,1,2,1)
        WProd78 = Prod(2,1,1,1)*Prod(4,1,1,1)*ProdC(2,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd79 = Prod(3,1,1,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd73 = Prod(3,1,2,1)*Prod(4,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
	WProd80 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(3,1,2,1)**2
	WProd63 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd66 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd68 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd69 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,1,1)
	WProd75 = Prod(3,1,1,1)**2*ProdC(3,1,2,1)*ProdC(4,1,1,1)
	WProd74 = Prod(3,1,1,1)*ProdC(4,1,2,1)
	WProd70 = Prod(3,1,2,1)**2*ProdC(3,1,1,1)*ProdC(4,1,2,1)
        WProd76 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,2,1)
      
      
	dr=1D0
	If(dr.eq.0D0) then      
	IInt81 = myCget(WProd64,ME2,ME2,0.D0,0.D0,ME2,0)
	LInt158 = myCval(cc0,IInt81)
	IInt82 = myCget(WProd64,ME2,ME2,0.D0,MZ2r,ME2,0)
	LInt159 = myCval(cc0,IInt82)
	IInt83 = myCget(WProd64,ME2,ME2,MZ2r,0.D0,ME2,0)
	LInt160 = myCval(cc0,IInt83)
        IInt84 = myDget(MD2,WProd64,ME2,-WProd65,MD2,ME2,MD2,0.D0,
     -   0.D0,ME2,2)
	LInt161 = myDval(dd00,IInt84)
	LInt167 = myDval(dd1,IInt84)
	LInt173 = myDval(dd11,IInt84)
	LInt179 = myDval(dd12,IInt84)
	LInt185 = myDval(dd13,IInt84)
	LInt191 = myDval(dd2,IInt84)
	LInt197 = myDval(dd23,IInt84)
	LInt203 = myDval(dd3,IInt84)
	LInt209 = myDval(dd33,IInt84)
        IInt85 = myDget(MD2,WProd64,ME2,-WProd65,MD2,ME2,MD2,0.D0,
     -   MZ2r,ME2,2)
	LInt162 = myDval(dd00,IInt85)
	LInt168 = myDval(dd1,IInt85)
	LInt174 = myDval(dd11,IInt85)
	LInt180 = myDval(dd12,IInt85)
	LInt186 = myDval(dd13,IInt85)
	LInt192 = myDval(dd2,IInt85)
	LInt198 = myDval(dd23,IInt85)
	LInt204 = myDval(dd3,IInt85)
	LInt210 = myDval(dd33,IInt85)
        IInt86 = myDget(MD2,WProd64,ME2,-WProd65,MD2,ME2,MD2,
     -   MZ2r,0.D0,ME2,2)
	LInt163 = myDval(dd00,IInt86)
	LInt169 = myDval(dd1,IInt86)
	LInt175 = myDval(dd11,IInt86)
	LInt181 = myDval(dd12,IInt86)
	LInt187 = myDval(dd13,IInt86)
	LInt193 = myDval(dd2,IInt86)
	LInt199 = myDval(dd23,IInt86)
	LInt205 = myDval(dd3,IInt86)
	LInt211 = myDval(dd33,IInt86)
        IInt87 = myDget(MD2,WProd64,ME2,-WProd66,MD2,ME2,MD2,0.D0,
     -   0.D0,ME2,2)
	LInt164 = myDval(dd00,IInt87)
	LInt170 = myDval(dd1,IInt87)
	LInt176 = myDval(dd11,IInt87)
	LInt182 = myDval(dd12,IInt87)
	LInt188 = myDval(dd13,IInt87)
	LInt194 = myDval(dd2,IInt87)
	LInt200 = myDval(dd23,IInt87)
	LInt206 = myDval(dd3,IInt87)
	LInt212 = myDval(dd33,IInt87)
        IInt88 = myDget(MD2,WProd64,ME2,-WProd66,MD2,ME2,MD2,0.D0,
     -   MZ2r,ME2,2)
	LInt165 = myDval(dd00,IInt88)
	LInt171 = myDval(dd1,IInt88)
	LInt177 = myDval(dd11,IInt88)
	LInt183 = myDval(dd12,IInt88)
	LInt189 = myDval(dd13,IInt88)
	LInt195 = myDval(dd2,IInt88)
	LInt201 = myDval(dd23,IInt88)
	LInt207 = myDval(dd3,IInt88)
	LInt213 = myDval(dd33,IInt88)
        IInt89 = myDget(MD2,WProd64,ME2,-WProd66,MD2,ME2,MD2,
     -   MZ2r,0.D0,ME2,2)
	LInt166 = myDval(dd00,IInt89)
	LInt172 = myDval(dd1,IInt89)
	LInt178 = myDval(dd11,IInt89)
	LInt184 = myDval(dd12,IInt89)
	LInt190 = myDval(dd13,IInt89)
	LInt196 = myDval(dd2,IInt89)
	LInt202 = myDval(dd23,IInt89)
	LInt208 = myDval(dd3,IInt89)
	LInt214 = myDval(dd33,IInt89)
	Else
	IInt81 = myCget(WProd64,0D0,0D0,0.D0,0.D0,0D0,0)
	LInt158 = myCval(cc0,IInt81)
c	Print *, "Llint158", Lint158
	IInt82 = myCget(WProd64,0D0,0D0,0.D0,MZ2r,0D0,0)
	LInt159 = myCval(cc0,IInt82)
	IInt83 = myCget(WProd64,0D0,0D0,MZ2r,0.D0,0D0,0)
	LInt160 = myCval(cc0,IInt83)
        IInt84 = myDget(0D0,WProd64,0D0,-WProd65,0D0,0D0,0D0,0.D0,
     -   0.D0,0D0,2)
	LInt161 = myDval(dd00,IInt84)
	LInt167 = myDval(dd1,IInt84)
	LInt173 = myDval(dd11,IInt84)
	LInt179 = myDval(dd12,IInt84)
	LInt185 = myDval(dd13,IInt84)
	LInt191 = myDval(dd2,IInt84)
	LInt197 = myDval(dd23,IInt84)
	LInt203 = myDval(dd3,IInt84)
	LInt209 = myDval(dd33,IInt84)
        IInt85 = myDget(0D0,WProd64,0D0,-WProd65,0D0,0D0,0D0,0.D0,
     -   MZ2r,0D0,2)
	LInt162 = myDval(dd00,IInt85)
	LInt168 = myDval(dd1,IInt85)
	LInt174 = myDval(dd11,IInt85)
	LInt180 = myDval(dd12,IInt85)
	LInt186 = myDval(dd13,IInt85)
	LInt192 = myDval(dd2,IInt85)
	LInt198 = myDval(dd23,IInt85)
	LInt204 = myDval(dd3,IInt85)
	LInt210 = myDval(dd33,IInt85)
        IInt86 = myDget(0D0,WProd64,0D0,-WProd65,0D0,0D0,0D0,
     -   MZ2r,0.D0,0D0,2)
	LInt163 = myDval(dd00,IInt86)
	LInt169 = myDval(dd1,IInt86)
	LInt175 = myDval(dd11,IInt86)
	LInt181 = myDval(dd12,IInt86)
	LInt187 = myDval(dd13,IInt86)
	LInt193 = myDval(dd2,IInt86)
	LInt199 = myDval(dd23,IInt86)
	LInt205 = myDval(dd3,IInt86)
	LInt211 = myDval(dd33,IInt86)
        IInt87 = myDget(0D0,WProd64,0D0,-WProd66,0D0,0D0,0D0,0.D0,
     -   0.D0,0D0,2)
	LInt164 = myDval(dd00,IInt87)
	LInt170 = myDval(dd1,IInt87)
	LInt176 = myDval(dd11,IInt87)
	LInt182 = myDval(dd12,IInt87)
	LInt188 = myDval(dd13,IInt87)
	LInt194 = myDval(dd2,IInt87)
	LInt200 = myDval(dd23,IInt87)
	LInt206 = myDval(dd3,IInt87)
	LInt212 = myDval(dd33,IInt87)
        IInt88 = myDget(0D0,WProd64,0D0,-WProd66,0D0,0D0,0D0,0.D0,
     -   MZ2r,0D0,2)
	LInt165 = myDval(dd00,IInt88)
	LInt171 = myDval(dd1,IInt88)
	LInt177 = myDval(dd11,IInt88)
	LInt183 = myDval(dd12,IInt88)
	LInt189 = myDval(dd13,IInt88)
	LInt195 = myDval(dd2,IInt88)
	LInt201 = myDval(dd23,IInt88)
	LInt207 = myDval(dd3,IInt88)
	LInt213 = myDval(dd33,IInt88)
        IInt89 = myDget(0D0,WProd64,0D0,-WProd66,0D0,0D0,0D0,
     -   MZ2r,0.D0,0D0,2)
	LInt166 = myDval(dd00,IInt89)
	LInt172 = myDval(dd1,IInt89)
	LInt178 = myDval(dd11,IInt89)
	LInt184 = myDval(dd12,IInt89)
	LInt190 = myDval(dd13,IInt89)
	LInt196 = myDval(dd2,IInt89)
	LInt202 = myDval(dd23,IInt89)
	LInt208 = myDval(dd3,IInt89)
	LInt214 = myDval(dd33,IInt89)
	Endif
      
      
      endif

      end


