      subroutine abbrvirtQCD3(resnr)
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

      
	WProd164 = Prod(2,2,2,1)
	WProd172 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd165 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd182 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd181 = ProdC(2,2,2,1)
	WProd177 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd174 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd167 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd183 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd178 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd166 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd179 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd168 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd184 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
        WProd185 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,1,1))
	WProd175 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd169 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd180 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd170 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd171 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd186 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
	WProd187 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd188 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd163 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd173 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd176 = Prod(5,1,2,1)*ProdC(5,1,2,1)
      
      
	Prp28 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd163)
	Prp29 = 1.D0/WProd172
	Prp30 = -(1.D0/WProd176)
      
	IInt228 = myBget(WProd163,0.D0,0.D0,0)
	LInt253 = myBval(bb0,IInt228)
	IInt229 = myBget(WProd172,0.D0,0.D0,0)
	LInt254 = myBval(bb0,IInt229)
	IInt233 = myBget(WProd172,0.D0,0.D0,1)
	LInt272 = myBval(bb1,IInt233)
	IInt234 = myBget(-WProd176,0.D0,0.D0,0)
	LInt273 = myBval(bb0,IInt234)
	IInt238 = myBget(-WProd176,0.D0,0.D0,1)
	LInt290 = myBval(bb1,IInt238)
	IInt222 = myCget(WProd172,0.D0,WProd163,0.D0,0.D0,0.D0,1)
	LInt185 = myCval(cc0,IInt222)
	LInt188 = myCval(cc1,IInt222)
	LInt191 = myCval(cc2,IInt222)
	IInt223 = myCget(WProd172,WProd163,0.D0,0.D0,0.D0,0.D0,1)
	LInt186 = myCval(cc0,IInt223)
	LInt189 = myCval(cc1,IInt223)
	LInt192 = myCval(cc2,IInt223)
	IInt224 = myCget(-WProd173,0.D0,WProd163,0.D0,0.D0,0.D0,1)
	LInt187 = myCval(cc0,IInt224)
	LInt190 = myCval(cc1,IInt224)
	LInt193 = myCval(cc2,IInt224)
	IInt237 = myCget(0.D0,-WProd176,0.D0,0.D0,0.D0,0.D0,2)
	LInt277 = myCval(cc00,IInt237)
	LInt283 = myCval(cc12,IInt237)
	LInt288 = myCval(cc22,IInt237)
	LInt291 = myCval(cc2,IInt237)
	IInt232 = myCget(0.D0,WProd172,0.D0,0.D0,0.D0,0.D0,2)
	LInt259 = myCval(cc00,IInt232)
	LInt262 = myCval(cc1,IInt232)
	LInt264 = myCval(cc11,IInt232)
	LInt267 = myCval(cc12,IInt232)
	LInt270 = myCval(cc2,IInt232)
	IInt230 = myCget(0.D0,WProd163,WProd172,0.D0,0.D0,0.D0,2)
	LInt255 = myCval(cc0,IInt230)
	LInt257 = myCval(cc00,IInt230)
	LInt260 = myCval(cc1,IInt230)
	LInt265 = myCval(cc12,IInt230)
	LInt268 = myCval(cc2,IInt230)
	LInt271 = myCval(cc22,IInt230)
	IInt231 = myCget(0.D0,WProd172,0.D0,0.D0,0.D0,0.D0,2)
	LInt256 = myCval(cc0,IInt231)
	LInt258 = myCval(cc00,IInt231)
	LInt261 = myCval(cc1,IInt231)
	LInt263 = myCval(cc11,IInt231)
	LInt266 = myCval(cc12,IInt231)
	LInt269 = myCval(cc2,IInt231)
	IInt235 = myCget(0.D0,-WProd176,0.D0,0.D0,0.D0,0.D0,2)
	LInt274 = myCval(cc0,IInt235)
	LInt276 = myCval(cc00,IInt235)
	LInt279 = myCval(cc1,IInt235)
	LInt282 = myCval(cc12,IInt235)
	LInt285 = myCval(cc2,IInt235)
	LInt287 = myCval(cc22,IInt235)
	IInt236 = myCget(0.D0,WProd163,-WProd176,0.D0,0.D0,0.D0,2)
	LInt275 = myCval(cc0,IInt236)
	LInt278 = myCval(cc00,IInt236)
	LInt280 = myCval(cc1,IInt236)
	LInt281 = myCval(cc11,IInt236)
	LInt284 = myCval(cc12,IInt236)
	LInt286 = myCval(cc2,IInt236)
	LInt289 = myCval(cc22,IInt236)
        IInt227 = myDget(0.D0,WProd172,WProd163,-WProd173,0.D0,0.D0,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt197 = myDval(dd00,IInt227)
	LInt200 = myDval(dd001,IInt227)
	LInt203 = myDval(dd002,IInt227)
	LInt206 = myDval(dd003,IInt227)
	LInt209 = myDval(dd1,IInt227)
	LInt212 = myDval(dd11,IInt227)
	LInt215 = myDval(dd112,IInt227)
	LInt219 = myDval(dd12,IInt227)
	LInt222 = myDval(dd122,IInt227)
	LInt224 = myDval(dd123,IInt227)
	LInt227 = myDval(dd13,IInt227)
	LInt232 = myDval(dd2,IInt227)
	LInt237 = myDval(dd223,IInt227)
	LInt240 = myDval(dd23,IInt227)
	LInt243 = myDval(dd233,IInt227)
	LInt246 = myDval(dd3,IInt227)
	LInt249 = myDval(dd33,IInt227)
        IInt226 = myDget(0.D0,-WProd173,0.D0,-WProd176,0.D0,WProd163,
     -   0.D0,0.D0,0.D0,0.D0,3)
	LInt195 = myDval(dd0,IInt226)
	LInt198 = myDval(dd00,IInt226)
	LInt201 = myDval(dd001,IInt226)
	LInt204 = myDval(dd002,IInt226)
	LInt207 = myDval(dd003,IInt226)
	LInt210 = myDval(dd1,IInt226)
	LInt213 = myDval(dd11,IInt226)
	LInt217 = myDval(dd113,IInt226)
	LInt220 = myDval(dd12,IInt226)
	LInt225 = myDval(dd123,IInt226)
	LInt228 = myDval(dd13,IInt226)
	LInt230 = myDval(dd133,IInt226)
	LInt233 = myDval(dd2,IInt226)
	LInt235 = myDval(dd22,IInt226)
	LInt238 = myDval(dd223,IInt226)
	LInt241 = myDval(dd23,IInt226)
	LInt244 = myDval(dd233,IInt226)
	LInt247 = myDval(dd3,IInt226)
	LInt250 = myDval(dd33,IInt226)
	LInt252 = myDval(dd333,IInt226)
        IInt225 = myDget(0.D0,WProd172,0.D0,-WProd176,0.D0,WProd163,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt194 = myDval(dd0,IInt225)
	LInt196 = myDval(dd00,IInt225)
	LInt199 = myDval(dd001,IInt225)
	LInt202 = myDval(dd002,IInt225)
	LInt205 = myDval(dd003,IInt225)
	LInt208 = myDval(dd1,IInt225)
	LInt211 = myDval(dd11,IInt225)
	LInt214 = myDval(dd112,IInt225)
	LInt216 = myDval(dd113,IInt225)
	LInt218 = myDval(dd12,IInt225)
	LInt221 = myDval(dd122,IInt225)
	LInt223 = myDval(dd123,IInt225)
	LInt226 = myDval(dd13,IInt225)
	LInt229 = myDval(dd133,IInt225)
	LInt231 = myDval(dd2,IInt225)
	LInt234 = myDval(dd22,IInt225)
	LInt236 = myDval(dd223,IInt225)
	LInt239 = myDval(dd23,IInt225)
	LInt242 = myDval(dd233,IInt225)
	LInt245 = myDval(dd3,IInt225)
	LInt248 = myDval(dd33,IInt225)
	LInt251 = myDval(dd333,IInt225)
      
      
      endif

      end


