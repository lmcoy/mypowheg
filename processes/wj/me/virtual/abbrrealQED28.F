      subroutine abbrrealQED28(resnr)
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

      
	WProd2292 = Prod(6,2,6,1)
	WProd2236 = Prod(2,2,2,1)*Prod(6,2,6,1)
	WProd2182 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd2198 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd2265 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd2288 = Prod(6,2,6,1)*ProdC(2,2,2,1)
	WProd2188 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd2189 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd2220 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd2277 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd2190 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd2191 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd2221 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd2278 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd2183 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd2192 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd2193 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd2222 = Prod(5,1,2,2)*ProdC(5,1,2,1)
        WProd2225 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd2204 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd2304 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd2205 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd2233 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd2286 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
	WProd2279 = Prod(5,1,2,1)*ProdC(5,1,2,2)
        WProd2270 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,2))
        WProd2271 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd2194 = Prod(5,1,3,1)*ProdC(5,1,3,1)
        WProd2226 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd2223 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd2280 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd2208 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
	WProd2195 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd2274 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd2227 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2224 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2218 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2184 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2281 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2260 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2219 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2186 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2261 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,2)*
     -     ProdC(5,1,4,1))
	WProd2202 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd2203 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd2237 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd2252 = Prod(6,2,2,2)*ProdC(6,1,2,1)
        WProd2243 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd2244 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd2253 = -(Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2254 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2299 = -(Prod(3,1,2,1)**2*Prod(6,2,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2255 = -(Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd2238 = -(Prod(2,2,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd2256 = -(Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2239 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2257 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2300 = -(Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2240 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
	WProd2298 = Prod(6,2,2,1)*ProdC(6,1,2,2)
        WProd2295 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2296 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2301 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,2))
        WProd2289 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,2))
        WProd2290 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2302 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2291 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
	WProd2245 = Prod(6,2,3,1)*ProdC(6,1,3,1)
	WProd2246 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd2258 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,1,4,1))
        WProd2247 = -(Prod(2,2,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,1,4,1))
        WProd2297 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,2)*
     -     ProdC(6,1,4,1))
        WProd2248 = -(Prod(3,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,3,1)*
     -     ProdC(6,1,4,1))
        WProd2293 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,3,1)*ProdC(6,1,4,1))
        WProd2241 = Prod(4,1,3,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd2249 = -(Prod(4,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd2294 = -(Prod(2,1,1,1)*Prod(4,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,1,4,1))
	WProd2216 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd2250 = Prod(3,1,1,1)*Prod(6,2,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd2251 = Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd2259 = Prod(3,1,2,2)*Prod(6,2,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd2303 = Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(4,1,2,2)*
     -   ProdC(6,1,5,1)
        WProd2242 = Prod(5,1,3,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
	WProd2185 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd2217 = Prod(6,1,2,2)*ProdC(6,2,2,1)
        WProd2206 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd2207 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd2228 = -(Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd2229 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd2282 = -(Prod(3,1,2,1)**2*Prod(6,1,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd2230 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd2187 = -(Prod(2,2,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd2231 = -(Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2196 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2232 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2283 = -(Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2197 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
	WProd2276 = Prod(6,1,2,1)*ProdC(6,2,2,2)
        WProd2272 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd2273 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd2284 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,2))
        WProd2262 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,2))
        WProd2263 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd2285 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd2264 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
	WProd2209 = Prod(6,1,3,1)*ProdC(6,2,3,1)
	WProd2210 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd2234 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,2,4,1))
        WProd2211 = -(Prod(2,2,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,2,4,1))
        WProd2275 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,2)*
     -     ProdC(6,2,4,1))
        WProd2212 = -(Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,3,1)*
     -     ProdC(6,2,4,1))
        WProd2267 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,3,1)*ProdC(6,2,4,1))
        WProd2199 = Prod(4,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,4,1)
        WProd2213 = -(Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,4,1))
        WProd2268 = -(Prod(2,1,1,1)*Prod(4,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,2,4,1))
        WProd2214 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,1)
        WProd2215 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,5,1)
        WProd2235 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,5,1)
        WProd2287 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(6,2,5,1)
        WProd2200 = Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,5,1)
	WProd2269 = ProdC(6,2,6,1)
	WProd2201 = Prod(2,2,2,1)*ProdC(6,2,6,1)
	WProd2266 = ProdC(2,2,2,1)*ProdC(6,2,6,1)
      
      
	Prp233 = 1.D0/WProd2182
	Prp234 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd2183)
	Prp238 = 1.D0/(WProd2183 - WProd2188 - WProd2190)
	Prp239 = -(1.D0/WProd2193)
	Prp240 = -(1.D0/WProd2202)
	Prp235 = 1.D0/(WProd2182 - WProd2202 - WProd2203)
	Prp237 = -(1.D0/WProd2203)
	Prp236 = 1.D0/WProd2216
      
      
      
      endif

      end


