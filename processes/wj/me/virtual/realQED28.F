      subroutine realQED28(amp,resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer helnr,colnr,resnr,procnr

      double complex amp(nproc,ncol,nhel)

      double precision Alfas,Alfas2,GS

      common /my_alfas_common/Alfas,Alfas2,GS
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 

*The Amplitude is evaluated as amp(procnr,colnr,helnr), where
*
*   Resonances:
*     resnr=0:   non-resonant
*
*   Partonic Subprocesses:
*     procnr=1:   {~up  ,  up  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=2:   {~up  , ~do  ,  ve  , ~el  , ~up  , ~up  }
*     procnr=3:   { up  ,  up  ,  ve  , ~el  ,  up  ,  do  }
*     procnr=4:   { up  , ~ch  ,  ve  , ~el  , ~ch  ,  do  }
*     procnr=5:   { up  ,  ch  ,  ve  , ~el  ,  ch  ,  do  }
*     procnr=6:   { up  , ~do  ,  ve  , ~el  , ~up  ,  up  }
*     procnr=7:   { up  , ~do  ,  ve  , ~el  , ~ch  ,  ch  }
*     procnr=8:   { up  , ~do  ,  ve  , ~el  , ~do  ,  do  }
*     procnr=9:   { up  , ~do  ,  ve  , ~el  , ~st  ,  st  }
*     procnr=10:   { up  , ~do  ,  ve  , ~el  , ~bo  ,  bo  }
*     procnr=11:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  ga  }
*     procnr=12:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  gl  }
*     procnr=13:   { up  ,  do  ,  ve  , ~el  ,  do  ,  do  }
*     procnr=14:   { up  , ~st  ,  ve  , ~el  ,  do  , ~st  }
*     procnr=15:   { up  ,  st  ,  ve  , ~el  ,  do  ,  st  }
*     procnr=16:   { up  ,  ga  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=17:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  ga  }
*     procnr=18:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=19:   {~ch  ,  ch  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=20:   {~ch  , ~do  ,  ve  , ~el  , ~up  , ~ch  }
*     procnr=21:   { ch  , ~do  ,  ve  , ~el  , ~up  ,  ch  }
*     procnr=22:   {~do  , ~do  ,  ve  , ~el  , ~up  , ~do  }
*     procnr=23:   {~do  ,  do  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=24:   {~do  , ~st  ,  ve  , ~el  , ~up  , ~st  }
*     procnr=25:   {~do  ,  st  ,  ve  , ~el  , ~up  ,  st  }
*     procnr=26:   {~do  ,  ga  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=27:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  ga  }
*     procnr=28:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=29:   {~st  ,  st  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=30:   { ga  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=31:   { gl  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*
*   Helicities:
*     helnr=1:   {-1, -1, -1, 1, -1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1, -1}
*     helnr=3:   {1, -1, -1, 1, 1, -1}
*     helnr=4:   {-1, 1, -1, 1, -1, 1}
*     helnr=5:   {1, 1, -1, 1, 1, 1}
*     helnr=6:   {-1, 1, -1, 1, -1, -1}
*     helnr=7:   {-1, 1, -1, 1, 1, 1}
*     helnr=8:   {-1, -1, -1, 1, -1, 1}
*     helnr=9:   {1, -1, -1, 1, 1, 1}
*     helnr=10:   {1, 1, -1, 1, 1, -1}
*     helnr=11:   {-1, -1, -1, 1, 1, -1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col3] SUNT[Col2, Col4]
*     colnr=2:   SUNT[Col1, Col2] SUNT[Col3, Col4]
*     colnr=3:   SUNT[Col1, Col4] SUNT[Col2, Col3]
*     colnr=4:   SUNT[Glu1, Col2, Col1]
*     colnr=5:   SUNT[Glu1, Glu2, Col2, Col1]
*     colnr=6:   SUNT[Glu2, Glu1, Col2, Col1]
*     colnr=7:   SUNT[Glu1, Col1, Col2]
*     colnr=8:   SUNT[Glu1, Glu2, Col1, Col2]
*     colnr=9:   SUNT[Glu2, Glu1, Col1, Col2]


      if(resnr.eq.0)then

        amp(28,8,3) = Prp233*Prp234*
     -   ((32*Alfas*Cpl122*Prp235*
     -        (WProd2184*WProd2185 + WProd2185*WProd2186 + 
     -          2*WProd2182*WProd2187 + 2*WProd2183*WProd2187 - 
     -          2*WProd2187*WProd2188 - 2*WProd2187*WProd2189 - 
     -          2*WProd2187*WProd2190 - 2*WProd2187*WProd2191 - 
     -          2*WProd2187*WProd2192 - 2*WProd2187*WProd2193 + 
     -          2*WProd2187*WProd2194 + 2*WProd2187*WProd2195 - 
     -          WProd2196 - WProd2197 + WProd2198*WProd2199 + 
     -          WProd2198*WProd2200))/WProd2201 + 
     -     (32*Alfas*Cpl122*Prp236*
     -        (WProd2184*WProd2185 + 2*WProd2182*WProd2187 + 
     -          WProd2183*WProd2187 - WProd2187*WProd2188 - 
     -          WProd2187*WProd2189 - WProd2187*WProd2190 - 
     -          WProd2187*WProd2191 - WProd2187*WProd2192 - 
     -          WProd2187*WProd2193 + WProd2187*WProd2194 + 
     -          WProd2187*WProd2195 - WProd2196 + 
     -          WProd2198*WProd2199 + WProd2185*WProd2204 + 
     -          WProd2185*WProd2205 - WProd2206 - WProd2207 - 
     -          WProd2198*WProd2208*WProd2209 - 
     -          WProd2198*WProd2208*WProd2210 - 
     -          2*WProd2182*WProd2211 - WProd2183*WProd2211 + 
     -          WProd2188*WProd2211 + WProd2189*WProd2211 + 
     -          WProd2190*WProd2211 + WProd2191*WProd2211 + 
     -          WProd2192*WProd2211 + WProd2193*WProd2211 - 
     -          WProd2194*WProd2211 - WProd2195*WProd2211 + 
     -          WProd2198*WProd2212 + WProd2198*WProd2213 + 
     -          WProd2198*WProd2214 + WProd2198*WProd2215))/
     -      WProd2201)
        amp(28,8,3) = Prp234*
     -    ((16*Alfas*Cpl122*Prp235*Prp237*
     -         (2*WProd2184*WProd2185 + 2*WProd2185*WProd2186 + 
     -           2*WProd2198*WProd2199 + 2*WProd2198*WProd2200 + 
     -           4*WProd2182*WProd2208*WProd2217 + 
     -           5*WProd2183*WProd2208*WProd2217 - 
     -           4*WProd2188*WProd2208*WProd2217 - 
     -           4*WProd2189*WProd2208*WProd2217 - 
     -           4*WProd2190*WProd2208*WProd2217 - 
     -           4*WProd2191*WProd2208*WProd2217 - 
     -           4*WProd2192*WProd2208*WProd2217 - 
     -           4*WProd2193*WProd2208*WProd2217 + 
     -           5*WProd2194*WProd2208*WProd2217 + 
     -           5*WProd2195*WProd2208*WProd2217 - 
     -           2*WProd2217*WProd2218 - 2*WProd2217*WProd2219 - 
     -           2*WProd2199*WProd2220 - 2*WProd2200*WProd2220 - 
     -           2*WProd2199*WProd2221 - 2*WProd2200*WProd2221 - 
     -           2*WProd2199*WProd2222 - 2*WProd2200*WProd2222))/
     -       WProd2201 + Prp238*
     -       ((16*Alfas*Cpl122*Prp237*
     -            (2*WProd2184*WProd2185 + 2*WProd2198*WProd2199 + 
     -              2*WProd2185*WProd2204 - 
     -              2*WProd2198*WProd2208*WProd2209 - 
     -              2*WProd2198*WProd2208*WProd2210 + 
     -              2*WProd2198*WProd2212 + 
     -              2*WProd2198*WProd2213 + 
     -              2*WProd2198*WProd2214 + 
     -              3*WProd2182*WProd2208*WProd2217 + 
     -              2*WProd2183*WProd2208*WProd2217 - 
     -              2*WProd2188*WProd2208*WProd2217 - 
     -              WProd2189*WProd2208*WProd2217 - 
     -              2*WProd2190*WProd2208*WProd2217 - 
     -              WProd2191*WProd2208*WProd2217 - 
     -              3*WProd2192*WProd2208*WProd2217 - 
     -              3*WProd2193*WProd2208*WProd2217 + 
     -              3*WProd2194*WProd2208*WProd2217 + 
     -              3*WProd2195*WProd2208*WProd2217 - 
     -              2*WProd2217*WProd2218 - 
     -              2*WProd2199*WProd2220 - 
     -              2*WProd2185*WProd2208*WProd2220 + 
     -              2*WProd2208*WProd2209*WProd2220 + 
     -              2*WProd2208*WProd2210*WProd2220 - 
     -              2*WProd2212*WProd2220 - 
     -              2*WProd2213*WProd2220 - 
     -              2*WProd2214*WProd2220 - 
     -              2*WProd2199*WProd2221 - 
     -              2*WProd2185*WProd2208*WProd2221 + 
     -              2*WProd2208*WProd2209*WProd2221 + 
     -              2*WProd2208*WProd2210*WProd2221 - 
     -              2*WProd2212*WProd2221 - 
     -              2*WProd2213*WProd2221 - 
     -              2*WProd2214*WProd2221 - 
     -              2*WProd2199*WProd2222 + 
     -              2*WProd2208*WProd2209*WProd2222 + 
     -              2*WProd2208*WProd2210*WProd2222 - 
     -              2*WProd2212*WProd2222 - 
     -              2*WProd2213*WProd2222 - 
     -              2*WProd2214*WProd2222 + 
     -              2*WProd2185*WProd2223 + 
     -              2*WProd2185*WProd2224 - 
     -              2*WProd2217*WProd2225 - 
     -              2*WProd2217*WProd2226 - 2*WProd2217*WProd2227))
     -           /WProd2201 + 
     -         (32*Alfas*Cpl122*Prp236*
     -            (WProd2184*WProd2185 + 2*WProd2182*WProd2187 + 
     -              WProd2183*WProd2187 - WProd2187*WProd2188 - 
     -              WProd2187*WProd2189 - WProd2187*WProd2190 - 
     -              WProd2187*WProd2191 - 2*WProd2187*WProd2192 - 
     -              2*WProd2187*WProd2193 + 
     -              2*WProd2187*WProd2194 + 
     -              2*WProd2187*WProd2195 - WProd2196 + 
     -              WProd2198*WProd2199 + WProd2185*WProd2204 - 
     -              WProd2206 + WProd2205*WProd2209 - 
     -              WProd2198*WProd2208*WProd2209 + 
     -              WProd2205*WProd2210 - 
     -              WProd2198*WProd2208*WProd2210 - 
     -              2*WProd2182*WProd2211 - WProd2183*WProd2211 + 
     -              WProd2188*WProd2211 + WProd2189*WProd2211 + 
     -              WProd2190*WProd2211 + WProd2191*WProd2211 + 
     -              2*WProd2192*WProd2211 + 
     -              2*WProd2193*WProd2211 - 
     -              2*WProd2194*WProd2211 - 
     -              2*WProd2195*WProd2211 + WProd2198*WProd2212 + 
     -              WProd2198*WProd2213 + WProd2198*WProd2214 - 
     -              WProd2199*WProd2220 - 
     -              WProd2185*WProd2208*WProd2220 + 
     -              WProd2208*WProd2209*WProd2220 + 
     -              WProd2208*WProd2210*WProd2220 - 
     -              WProd2212*WProd2220 - WProd2213*WProd2220 - 
     -              WProd2214*WProd2220 - WProd2215*WProd2220 - 
     -              WProd2199*WProd2221 - 
     -              WProd2185*WProd2208*WProd2221 + 
     -              WProd2208*WProd2209*WProd2221 + 
     -              WProd2208*WProd2210*WProd2221 - 
     -              WProd2212*WProd2221 - WProd2213*WProd2221 - 
     -              WProd2214*WProd2221 - WProd2215*WProd2221 + 
     -              WProd2185*WProd2223 + WProd2185*WProd2224 - 
     -              WProd2228 - WProd2229 + WProd2220*WProd2230 + 
     -              WProd2221*WProd2230 - WProd2231 - WProd2232 - 
     -              WProd2209*WProd2233 - WProd2210*WProd2233 + 
     -              WProd2220*WProd2234 + WProd2221*WProd2234 + 
     -              2*WProd2182*WProd2235 + WProd2183*WProd2235 - 
     -              WProd2188*WProd2235 - WProd2189*WProd2235 - 
     -              WProd2190*WProd2235 - WProd2191*WProd2235 - 
     -              2*WProd2192*WProd2235 - 
     -              2*WProd2193*WProd2235 + 
     -              2*WProd2194*WProd2235 + 2*WProd2195*WProd2235))
     -           /WProd2201)) + amp(28,8,3)

#ifdef DEBUG
	print *, 'amp(28,8,3) =', amp(28,8,3)
#endif

        amp(28,8,9) = Prp233*Prp234*
     -   ((32*Alfas*Cpl122*Prp235*
     -        (WProd2184*WProd2237 + WProd2186*WProd2237 + 
     -          2*WProd2182*WProd2238 + 2*WProd2183*WProd2238 - 
     -          2*WProd2188*WProd2238 - 2*WProd2189*WProd2238 - 
     -          2*WProd2190*WProd2238 - 2*WProd2191*WProd2238 - 
     -          2*WProd2192*WProd2238 - 2*WProd2193*WProd2238 + 
     -          2*WProd2194*WProd2238 + 2*WProd2195*WProd2238 - 
     -          WProd2239 - WProd2240 + WProd2198*WProd2241 + 
     -          WProd2198*WProd2242))/WProd2236 + 
     -     (32*Alfas*Cpl122*Prp236*
     -        (WProd2184*WProd2237 + WProd2204*WProd2237 + 
     -          WProd2205*WProd2237 + 2*WProd2182*WProd2238 + 
     -          WProd2183*WProd2238 - WProd2188*WProd2238 - 
     -          WProd2189*WProd2238 - WProd2190*WProd2238 - 
     -          WProd2191*WProd2238 - WProd2192*WProd2238 - 
     -          WProd2193*WProd2238 + WProd2194*WProd2238 + 
     -          WProd2195*WProd2238 - WProd2239 + 
     -          WProd2198*WProd2241 - WProd2243 - WProd2244 - 
     -          WProd2198*WProd2208*WProd2245 - 
     -          WProd2198*WProd2208*WProd2246 - 
     -          2*WProd2182*WProd2247 - WProd2183*WProd2247 + 
     -          WProd2188*WProd2247 + WProd2189*WProd2247 + 
     -          WProd2190*WProd2247 + WProd2191*WProd2247 + 
     -          WProd2192*WProd2247 + WProd2193*WProd2247 - 
     -          WProd2194*WProd2247 - WProd2195*WProd2247 + 
     -          WProd2198*WProd2248 + WProd2198*WProd2249 + 
     -          WProd2198*WProd2250 + WProd2198*WProd2251))/
     -      WProd2236)
        amp(28,8,9) = Prp234*
     -    ((16*Alfas*Cpl122*Prp235*Prp237*
     -         (2*WProd2184*WProd2237 + 2*WProd2186*WProd2237 + 
     -           2*WProd2198*WProd2241 - 2*WProd2220*WProd2241 - 
     -           2*WProd2221*WProd2241 - 2*WProd2222*WProd2241 + 
     -           2*WProd2198*WProd2242 - 2*WProd2220*WProd2242 - 
     -           2*WProd2221*WProd2242 - 2*WProd2222*WProd2242 + 
     -           4*WProd2182*WProd2208*WProd2252 + 
     -           5*WProd2183*WProd2208*WProd2252 - 
     -           4*WProd2188*WProd2208*WProd2252 - 
     -           4*WProd2189*WProd2208*WProd2252 - 
     -           4*WProd2190*WProd2208*WProd2252 - 
     -           4*WProd2191*WProd2208*WProd2252 - 
     -           4*WProd2192*WProd2208*WProd2252 - 
     -           4*WProd2193*WProd2208*WProd2252 + 
     -           5*WProd2194*WProd2208*WProd2252 + 
     -           5*WProd2195*WProd2208*WProd2252 - 
     -           2*WProd2218*WProd2252 - 2*WProd2219*WProd2252))/
     -       WProd2236 + Prp238*
     -       ((16*Alfas*Cpl122*Prp237*
     -            (2*WProd2184*WProd2237 + 2*WProd2204*WProd2237 - 
     -              2*WProd2208*WProd2220*WProd2237 - 
     -              2*WProd2208*WProd2221*WProd2237 + 
     -              2*WProd2223*WProd2237 + 
     -              2*WProd2224*WProd2237 + 
     -              2*WProd2198*WProd2241 - 
     -              2*WProd2220*WProd2241 - 
     -              2*WProd2221*WProd2241 - 
     -              2*WProd2222*WProd2241 - 
     -              2*WProd2198*WProd2208*WProd2245 + 
     -              2*WProd2208*WProd2220*WProd2245 + 
     -              2*WProd2208*WProd2221*WProd2245 + 
     -              2*WProd2208*WProd2222*WProd2245 - 
     -              2*WProd2198*WProd2208*WProd2246 + 
     -              2*WProd2208*WProd2220*WProd2246 + 
     -              2*WProd2208*WProd2221*WProd2246 + 
     -              2*WProd2208*WProd2222*WProd2246 + 
     -              2*WProd2198*WProd2248 - 
     -              2*WProd2220*WProd2248 - 
     -              2*WProd2221*WProd2248 - 
     -              2*WProd2222*WProd2248 + 
     -              2*WProd2198*WProd2249 - 
     -              2*WProd2220*WProd2249 - 
     -              2*WProd2221*WProd2249 - 
     -              2*WProd2222*WProd2249 + 
     -              2*WProd2198*WProd2250 - 
     -              2*WProd2220*WProd2250 - 
     -              2*WProd2221*WProd2250 - 
     -              2*WProd2222*WProd2250 + 
     -              3*WProd2182*WProd2208*WProd2252 + 
     -              2*WProd2183*WProd2208*WProd2252 - 
     -              2*WProd2188*WProd2208*WProd2252 - 
     -              WProd2189*WProd2208*WProd2252 - 
     -              2*WProd2190*WProd2208*WProd2252 - 
     -              WProd2191*WProd2208*WProd2252 - 
     -              3*WProd2192*WProd2208*WProd2252 - 
     -              3*WProd2193*WProd2208*WProd2252 + 
     -              3*WProd2194*WProd2208*WProd2252 + 
     -              3*WProd2195*WProd2208*WProd2252 - 
     -              2*WProd2218*WProd2252 - 
     -              2*WProd2225*WProd2252 - 
     -              2*WProd2226*WProd2252 - 2*WProd2227*WProd2252))
     -           /WProd2236 + 
     -         (32*Alfas*Cpl122*Prp236*
     -            (WProd2184*WProd2237 + WProd2204*WProd2237 - 
     -              WProd2208*WProd2220*WProd2237 - 
     -              WProd2208*WProd2221*WProd2237 + 
     -              WProd2223*WProd2237 + WProd2224*WProd2237 + 
     -              2*WProd2182*WProd2238 + WProd2183*WProd2238 - 
     -              WProd2188*WProd2238 - WProd2189*WProd2238 - 
     -              WProd2190*WProd2238 - WProd2191*WProd2238 - 
     -              2*WProd2192*WProd2238 - 
     -              2*WProd2193*WProd2238 + 
     -              2*WProd2194*WProd2238 + 
     -              2*WProd2195*WProd2238 - WProd2239 + 
     -              WProd2198*WProd2241 - WProd2220*WProd2241 - 
     -              WProd2221*WProd2241 - WProd2243 + 
     -              WProd2205*WProd2245 - 
     -              WProd2198*WProd2208*WProd2245 + 
     -              WProd2208*WProd2220*WProd2245 + 
     -              WProd2208*WProd2221*WProd2245 - 
     -              WProd2233*WProd2245 + WProd2205*WProd2246 - 
     -              WProd2198*WProd2208*WProd2246 + 
     -              WProd2208*WProd2220*WProd2246 + 
     -              WProd2208*WProd2221*WProd2246 - 
     -              WProd2233*WProd2246 - 2*WProd2182*WProd2247 - 
     -              WProd2183*WProd2247 + WProd2188*WProd2247 + 
     -              WProd2189*WProd2247 + WProd2190*WProd2247 + 
     -              WProd2191*WProd2247 + 2*WProd2192*WProd2247 + 
     -              2*WProd2193*WProd2247 - 
     -              2*WProd2194*WProd2247 - 
     -              2*WProd2195*WProd2247 + WProd2198*WProd2248 - 
     -              WProd2220*WProd2248 - WProd2221*WProd2248 + 
     -              WProd2198*WProd2249 - WProd2220*WProd2249 - 
     -              WProd2221*WProd2249 + WProd2198*WProd2250 - 
     -              WProd2220*WProd2250 - WProd2221*WProd2250 - 
     -              WProd2220*WProd2251 - WProd2221*WProd2251 - 
     -              WProd2253 - WProd2254 + WProd2220*WProd2255 + 
     -              WProd2221*WProd2255 - WProd2256 - WProd2257 + 
     -              WProd2220*WProd2258 + WProd2221*WProd2258 + 
     -              2*WProd2182*WProd2259 + WProd2183*WProd2259 - 
     -              WProd2188*WProd2259 - WProd2189*WProd2259 - 
     -              WProd2190*WProd2259 - WProd2191*WProd2259 - 
     -              2*WProd2192*WProd2259 - 
     -              2*WProd2193*WProd2259 + 
     -              2*WProd2194*WProd2259 + 2*WProd2195*WProd2259))
     -           /WProd2236)) + amp(28,8,9)

#ifdef DEBUG
	print *, 'amp(28,8,9) =', amp(28,8,9)
#endif

        amp(28,8,10) = Prp233*Prp234*
     -   ((32*Alfas*Cpl122*Prp235*
     -        (WProd2185*WProd2260 + WProd2185*WProd2261 + 
     -          2*WProd2182*WProd2262 + 2*WProd2183*WProd2262 - 
     -          2*WProd2188*WProd2262 - 2*WProd2189*WProd2262 - 
     -          2*WProd2190*WProd2262 - 2*WProd2191*WProd2262 - 
     -          2*WProd2192*WProd2262 - 2*WProd2193*WProd2262 + 
     -          2*WProd2194*WProd2262 + 2*WProd2195*WProd2262 - 
     -          WProd2263 - WProd2264 + WProd2199*WProd2265 + 
     -          WProd2200*WProd2265))/WProd2266 + 
     -     Prp236*((32*Alfas*Cpl122*(WProd2267 + WProd2268))/
     -         WProd2269 + 
     -        (32*Alfas*Cpl122*
     -           (WProd2185*WProd2260 + 2*WProd2182*WProd2262 + 
     -             WProd2183*WProd2262 - WProd2188*WProd2262 - 
     -             WProd2189*WProd2262 - WProd2190*WProd2262 - 
     -             WProd2191*WProd2262 - WProd2192*WProd2262 - 
     -             WProd2193*WProd2262 + WProd2194*WProd2262 + 
     -             WProd2195*WProd2262 - WProd2263 + 
     -             WProd2199*WProd2265 - 
     -             WProd2208*WProd2209*WProd2265 - 
     -             WProd2208*WProd2210*WProd2265 + 
     -             WProd2212*WProd2265 + WProd2213*WProd2265 + 
     -             WProd2214*WProd2265 + WProd2215*WProd2265 + 
     -             WProd2185*WProd2270 + WProd2185*WProd2271 - 
     -             WProd2272 - WProd2273 - WProd2209*WProd2274 - 
     -             WProd2210*WProd2274 - 2*WProd2182*WProd2275 - 
     -             WProd2183*WProd2275 + WProd2188*WProd2275 + 
     -             WProd2189*WProd2275 + WProd2190*WProd2275 + 
     -             WProd2191*WProd2275 + WProd2192*WProd2275 + 
     -             WProd2193*WProd2275 - WProd2194*WProd2275 - 
     -             WProd2195*WProd2275))/WProd2266))
        amp(28,8,10) = Prp234*
     -    ((16*Alfas*Cpl122*Prp235*Prp237*
     -         (2*WProd2185*WProd2260 + 2*WProd2185*WProd2261 + 
     -           2*WProd2199*WProd2265 + 2*WProd2200*WProd2265 + 
     -           4*WProd2182*WProd2208*WProd2276 + 
     -           5*WProd2183*WProd2208*WProd2276 - 
     -           4*WProd2188*WProd2208*WProd2276 - 
     -           4*WProd2189*WProd2208*WProd2276 - 
     -           4*WProd2190*WProd2208*WProd2276 - 
     -           4*WProd2191*WProd2208*WProd2276 - 
     -           4*WProd2192*WProd2208*WProd2276 - 
     -           4*WProd2193*WProd2208*WProd2276 + 
     -           5*WProd2194*WProd2208*WProd2276 + 
     -           5*WProd2195*WProd2208*WProd2276 - 
     -           2*WProd2218*WProd2276 - 2*WProd2219*WProd2276 - 
     -           2*WProd2199*WProd2277 - 2*WProd2200*WProd2277 - 
     -           2*WProd2199*WProd2278 - 2*WProd2200*WProd2278 - 
     -           2*WProd2199*WProd2279 - 2*WProd2200*WProd2279))/
     -       WProd2266 + Prp238*
     -       ((16*Alfas*Cpl122*Prp237*
     -            (2*WProd2185*WProd2260 + 2*WProd2199*WProd2265 - 
     -              2*WProd2208*WProd2209*WProd2265 - 
     -              2*WProd2208*WProd2210*WProd2265 + 
     -              2*WProd2212*WProd2265 + 
     -              2*WProd2213*WProd2265 + 
     -              2*WProd2214*WProd2265 + 
     -              2*WProd2185*WProd2270 + 
     -              3*WProd2182*WProd2208*WProd2276 + 
     -              2*WProd2183*WProd2208*WProd2276 - 
     -              2*WProd2188*WProd2208*WProd2276 - 
     -              WProd2189*WProd2208*WProd2276 - 
     -              2*WProd2190*WProd2208*WProd2276 - 
     -              WProd2191*WProd2208*WProd2276 - 
     -              3*WProd2192*WProd2208*WProd2276 - 
     -              3*WProd2193*WProd2208*WProd2276 + 
     -              3*WProd2194*WProd2208*WProd2276 + 
     -              3*WProd2195*WProd2208*WProd2276 - 
     -              2*WProd2218*WProd2276 - 
     -              2*WProd2225*WProd2276 - 
     -              2*WProd2226*WProd2276 - 
     -              2*WProd2227*WProd2276 - 
     -              2*WProd2199*WProd2277 - 
     -              2*WProd2185*WProd2208*WProd2277 + 
     -              2*WProd2208*WProd2209*WProd2277 + 
     -              2*WProd2208*WProd2210*WProd2277 - 
     -              2*WProd2212*WProd2277 - 
     -              2*WProd2213*WProd2277 - 
     -              2*WProd2214*WProd2277 - 
     -              2*WProd2199*WProd2278 - 
     -              2*WProd2185*WProd2208*WProd2278 + 
     -              2*WProd2208*WProd2209*WProd2278 + 
     -              2*WProd2208*WProd2210*WProd2278 - 
     -              2*WProd2212*WProd2278 - 
     -              2*WProd2213*WProd2278 - 
     -              2*WProd2214*WProd2278 - 
     -              2*WProd2199*WProd2279 + 
     -              2*WProd2208*WProd2209*WProd2279 + 
     -              2*WProd2208*WProd2210*WProd2279 - 
     -              2*WProd2212*WProd2279 - 
     -              2*WProd2213*WProd2279 - 
     -              2*WProd2214*WProd2279 + 
     -              2*WProd2185*WProd2280 + 2*WProd2185*WProd2281))
     -           /WProd2266 + 
     -         Prp236*((32*Alfas*Cpl122*(WProd2267 + WProd2268))/
     -             WProd2269 + 
     -            (32*Alfas*Cpl122*
     -               (WProd2185*WProd2260 + 
     -                 2*WProd2182*WProd2262 + 
     -                 WProd2183*WProd2262 - WProd2188*WProd2262 - 
     -                 WProd2189*WProd2262 - WProd2190*WProd2262 - 
     -                 WProd2191*WProd2262 - 
     -                 2*WProd2192*WProd2262 - 
     -                 2*WProd2193*WProd2262 + 
     -                 2*WProd2194*WProd2262 + 
     -                 2*WProd2195*WProd2262 - WProd2263 + 
     -                 WProd2199*WProd2265 - 
     -                 WProd2208*WProd2209*WProd2265 - 
     -                 WProd2208*WProd2210*WProd2265 + 
     -                 WProd2212*WProd2265 + WProd2213*WProd2265 + 
     -                 WProd2214*WProd2265 + WProd2185*WProd2270 + 
     -                 WProd2209*WProd2271 + WProd2210*WProd2271 - 
     -                 WProd2272 - WProd2209*WProd2274 - 
     -                 WProd2210*WProd2274 - 
     -                 2*WProd2182*WProd2275 - 
     -                 WProd2183*WProd2275 + WProd2188*WProd2275 + 
     -                 WProd2189*WProd2275 + WProd2190*WProd2275 + 
     -                 WProd2191*WProd2275 + 
     -                 2*WProd2192*WProd2275 + 
     -                 2*WProd2193*WProd2275 - 
     -                 2*WProd2194*WProd2275 - 
     -                 2*WProd2195*WProd2275 - 
     -                 WProd2199*WProd2277 - 
     -                 WProd2185*WProd2208*WProd2277 + 
     -                 WProd2208*WProd2209*WProd2277 + 
     -                 WProd2208*WProd2210*WProd2277 - 
     -                 WProd2212*WProd2277 - WProd2213*WProd2277 - 
     -                 WProd2214*WProd2277 - WProd2215*WProd2277 + 
     -                 WProd2230*WProd2277 + WProd2234*WProd2277 - 
     -                 WProd2199*WProd2278 - 
     -                 WProd2185*WProd2208*WProd2278 + 
     -                 WProd2208*WProd2209*WProd2278 + 
     -                 WProd2208*WProd2210*WProd2278 - 
     -                 WProd2212*WProd2278 - WProd2213*WProd2278 - 
     -                 WProd2214*WProd2278 - WProd2215*WProd2278 + 
     -                 WProd2230*WProd2278 + WProd2234*WProd2278 + 
     -                 WProd2185*WProd2280 + WProd2185*WProd2281 - 
     -                 WProd2282 - WProd2283 - WProd2284 - 
     -                 WProd2285 - WProd2209*WProd2286 - 
     -                 WProd2210*WProd2286 + 
     -                 2*WProd2182*WProd2287 + 
     -                 WProd2183*WProd2287 - WProd2188*WProd2287 - 
     -                 WProd2189*WProd2287 - WProd2190*WProd2287 - 
     -                 WProd2191*WProd2287 - 
     -                 2*WProd2192*WProd2287 - 
     -                 2*WProd2193*WProd2287 + 
     -                 2*WProd2194*WProd2287 + 
     -                 2*WProd2195*WProd2287))/WProd2266))) + 
     -   amp(28,8,10)

#ifdef DEBUG
	print *, 'amp(28,8,10) =', amp(28,8,10)
#endif

        amp(28,8,5) = Prp233*Prp234*
     -   ((32*Alfas*Cpl122*Prp235*
     -        (WProd2237*WProd2260 + WProd2237*WProd2261 + 
     -          WProd2241*WProd2265 + WProd2242*WProd2265 + 
     -          2*WProd2182*WProd2289 + 2*WProd2183*WProd2289 - 
     -          2*WProd2188*WProd2289 - 2*WProd2189*WProd2289 - 
     -          2*WProd2190*WProd2289 - 2*WProd2191*WProd2289 - 
     -          2*WProd2192*WProd2289 - 2*WProd2193*WProd2289 + 
     -          2*WProd2194*WProd2289 + 2*WProd2195*WProd2289 - 
     -          WProd2290 - WProd2291))/WProd2288 + 
     -     Prp236*((32*Alfas*Cpl122*(WProd2293 + WProd2294))/
     -         WProd2292 + 
     -        (32*Alfas*Cpl122*
     -           (WProd2237*WProd2260 + WProd2241*WProd2265 - 
     -             WProd2208*WProd2245*WProd2265 - 
     -             WProd2208*WProd2246*WProd2265 + 
     -             WProd2248*WProd2265 + WProd2249*WProd2265 + 
     -             WProd2250*WProd2265 + WProd2251*WProd2265 + 
     -             WProd2237*WProd2270 + WProd2237*WProd2271 - 
     -             WProd2245*WProd2274 - WProd2246*WProd2274 + 
     -             2*WProd2182*WProd2289 + WProd2183*WProd2289 - 
     -             WProd2188*WProd2289 - WProd2189*WProd2289 - 
     -             WProd2190*WProd2289 - WProd2191*WProd2289 - 
     -             WProd2192*WProd2289 - WProd2193*WProd2289 + 
     -             WProd2194*WProd2289 + WProd2195*WProd2289 - 
     -             WProd2290 - WProd2295 - WProd2296 - 
     -             2*WProd2182*WProd2297 - WProd2183*WProd2297 + 
     -             WProd2188*WProd2297 + WProd2189*WProd2297 + 
     -             WProd2190*WProd2297 + WProd2191*WProd2297 + 
     -             WProd2192*WProd2297 + WProd2193*WProd2297 - 
     -             WProd2194*WProd2297 - WProd2195*WProd2297))/
     -         WProd2288))
        amp(28,8,5) = Prp234*
     -    ((16*Alfas*Cpl122*Prp235*Prp237*
     -         (2*WProd2237*WProd2260 + 2*WProd2237*WProd2261 + 
     -           2*WProd2241*WProd2265 + 2*WProd2242*WProd2265 - 
     -           2*WProd2241*WProd2277 - 2*WProd2242*WProd2277 - 
     -           2*WProd2241*WProd2278 - 2*WProd2242*WProd2278 - 
     -           2*WProd2241*WProd2279 - 2*WProd2242*WProd2279 + 
     -           4*WProd2182*WProd2208*WProd2298 + 
     -           5*WProd2183*WProd2208*WProd2298 - 
     -           4*WProd2188*WProd2208*WProd2298 - 
     -           4*WProd2189*WProd2208*WProd2298 - 
     -           4*WProd2190*WProd2208*WProd2298 - 
     -           4*WProd2191*WProd2208*WProd2298 - 
     -           4*WProd2192*WProd2208*WProd2298 - 
     -           4*WProd2193*WProd2208*WProd2298 + 
     -           5*WProd2194*WProd2208*WProd2298 + 
     -           5*WProd2195*WProd2208*WProd2298 - 
     -           2*WProd2218*WProd2298 - 2*WProd2219*WProd2298))/
     -       WProd2288 + Prp238*
     -       ((16*Alfas*Cpl122*Prp237*
     -            (2*WProd2237*WProd2260 + 2*WProd2241*WProd2265 - 
     -              2*WProd2208*WProd2245*WProd2265 - 
     -              2*WProd2208*WProd2246*WProd2265 + 
     -              2*WProd2248*WProd2265 + 
     -              2*WProd2249*WProd2265 + 
     -              2*WProd2250*WProd2265 + 
     -              2*WProd2237*WProd2270 - 
     -              2*WProd2208*WProd2237*WProd2277 - 
     -              2*WProd2241*WProd2277 + 
     -              2*WProd2208*WProd2245*WProd2277 + 
     -              2*WProd2208*WProd2246*WProd2277 - 
     -              2*WProd2248*WProd2277 - 
     -              2*WProd2249*WProd2277 - 
     -              2*WProd2250*WProd2277 - 
     -              2*WProd2208*WProd2237*WProd2278 - 
     -              2*WProd2241*WProd2278 + 
     -              2*WProd2208*WProd2245*WProd2278 + 
     -              2*WProd2208*WProd2246*WProd2278 - 
     -              2*WProd2248*WProd2278 - 
     -              2*WProd2249*WProd2278 - 
     -              2*WProd2250*WProd2278 - 
     -              2*WProd2241*WProd2279 + 
     -              2*WProd2208*WProd2245*WProd2279 + 
     -              2*WProd2208*WProd2246*WProd2279 - 
     -              2*WProd2248*WProd2279 - 
     -              2*WProd2249*WProd2279 - 
     -              2*WProd2250*WProd2279 + 
     -              2*WProd2237*WProd2280 + 
     -              2*WProd2237*WProd2281 + 
     -              3*WProd2182*WProd2208*WProd2298 + 
     -              2*WProd2183*WProd2208*WProd2298 - 
     -              2*WProd2188*WProd2208*WProd2298 - 
     -              WProd2189*WProd2208*WProd2298 - 
     -              2*WProd2190*WProd2208*WProd2298 - 
     -              WProd2191*WProd2208*WProd2298 - 
     -              3*WProd2192*WProd2208*WProd2298 - 
     -              3*WProd2193*WProd2208*WProd2298 + 
     -              3*WProd2194*WProd2208*WProd2298 + 
     -              3*WProd2195*WProd2208*WProd2298 - 
     -              2*WProd2218*WProd2298 - 
     -              2*WProd2225*WProd2298 - 
     -              2*WProd2226*WProd2298 - 2*WProd2227*WProd2298))
     -           /WProd2288 + 
     -         Prp236*((32*Alfas*Cpl122*(WProd2293 + WProd2294))/
     -             WProd2292 + 
     -            (32*Alfas*Cpl122*
     -               (WProd2237*WProd2260 + WProd2241*WProd2265 - 
     -                 WProd2208*WProd2245*WProd2265 - 
     -                 WProd2208*WProd2246*WProd2265 + 
     -                 WProd2248*WProd2265 + WProd2249*WProd2265 + 
     -                 WProd2250*WProd2265 + WProd2237*WProd2270 + 
     -                 WProd2245*WProd2271 + WProd2246*WProd2271 - 
     -                 WProd2245*WProd2274 - WProd2246*WProd2274 - 
     -                 WProd2208*WProd2237*WProd2277 - 
     -                 WProd2241*WProd2277 + 
     -                 WProd2208*WProd2245*WProd2277 + 
     -                 WProd2208*WProd2246*WProd2277 - 
     -                 WProd2248*WProd2277 - WProd2249*WProd2277 - 
     -                 WProd2250*WProd2277 - WProd2251*WProd2277 + 
     -                 WProd2255*WProd2277 + WProd2258*WProd2277 - 
     -                 WProd2208*WProd2237*WProd2278 - 
     -                 WProd2241*WProd2278 + 
     -                 WProd2208*WProd2245*WProd2278 + 
     -                 WProd2208*WProd2246*WProd2278 - 
     -                 WProd2248*WProd2278 - WProd2249*WProd2278 - 
     -                 WProd2250*WProd2278 - WProd2251*WProd2278 + 
     -                 WProd2255*WProd2278 + WProd2258*WProd2278 + 
     -                 WProd2237*WProd2280 + WProd2237*WProd2281 - 
     -                 WProd2245*WProd2286 - WProd2246*WProd2286 + 
     -                 2*WProd2182*WProd2289 + 
     -                 WProd2183*WProd2289 - WProd2188*WProd2289 - 
     -                 WProd2189*WProd2289 - WProd2190*WProd2289 - 
     -                 WProd2191*WProd2289 - 
     -                 2*WProd2192*WProd2289 - 
     -                 2*WProd2193*WProd2289 + 
     -                 2*WProd2194*WProd2289 + 
     -                 2*WProd2195*WProd2289 - WProd2290 - 
     -                 WProd2295 - 2*WProd2182*WProd2297 - 
     -                 WProd2183*WProd2297 + WProd2188*WProd2297 + 
     -                 WProd2189*WProd2297 + WProd2190*WProd2297 + 
     -                 WProd2191*WProd2297 + 
     -                 2*WProd2192*WProd2297 + 
     -                 2*WProd2193*WProd2297 - 
     -                 2*WProd2194*WProd2297 - 
     -                 2*WProd2195*WProd2297 - WProd2299 - 
     -                 WProd2300 - WProd2301 - WProd2302 + 
     -                 2*WProd2182*WProd2303 + 
     -                 WProd2183*WProd2303 - WProd2188*WProd2303 - 
     -                 WProd2189*WProd2303 - WProd2190*WProd2303 - 
     -                 WProd2191*WProd2303 - 
     -                 2*WProd2192*WProd2303 - 
     -                 2*WProd2193*WProd2303 + 
     -                 2*WProd2194*WProd2303 + 
     -                 2*WProd2195*WProd2303))/WProd2288))) + 
     -   amp(28,8,5)

#ifdef DEBUG
	print *, 'amp(28,8,5) =', amp(28,8,5)
#endif

        tmp1160 = (-16*Alfas*Cpl122*Prp235*Prp237*
     -      (2*WProd2184*WProd2185 + 2*WProd2185*WProd2186 + 
     -        2*WProd2198*WProd2199 + 2*WProd2198*WProd2200 + 
     -        4*WProd2182*WProd2208*WProd2217 + 
     -        5*WProd2183*WProd2208*WProd2217 - 
     -        4*WProd2188*WProd2208*WProd2217 - 
     -        4*WProd2189*WProd2208*WProd2217 - 
     -        4*WProd2190*WProd2208*WProd2217 - 
     -        4*WProd2191*WProd2208*WProd2217 - 
     -        4*WProd2192*WProd2208*WProd2217 - 
     -        4*WProd2193*WProd2208*WProd2217 + 
     -        5*WProd2194*WProd2208*WProd2217 + 
     -        5*WProd2195*WProd2208*WProd2217 - 
     -        2*WProd2217*WProd2218 - 2*WProd2217*WProd2219 - 
     -        2*WProd2199*WProd2220 - 2*WProd2200*WProd2220 - 
     -        2*WProd2199*WProd2221 - 2*WProd2200*WProd2221 - 
     -        2*WProd2199*WProd2222 - 2*WProd2200*WProd2222))/
     -    WProd2201 - (32*Alfas*Cpl122*Prp235*Prp240*
     -      (-2*WProd2182*WProd2187 - 3*WProd2183*WProd2187 + 
     -        2*WProd2187*WProd2188 + 2*WProd2187*WProd2189 + 
     -        2*WProd2187*WProd2190 + 2*WProd2187*WProd2191 + 
     -        2*WProd2187*WProd2192 + 2*WProd2187*WProd2193 - 
     -        3*WProd2187*WProd2194 - 3*WProd2187*WProd2195 + 
     -        WProd2196 + WProd2197 + 
     -        2*WProd2182*WProd2208*WProd2217 + 
     -        3*WProd2183*WProd2208*WProd2217 - 
     -        2*WProd2188*WProd2208*WProd2217 - 
     -        2*WProd2189*WProd2208*WProd2217 - 
     -        2*WProd2190*WProd2208*WProd2217 - 
     -        2*WProd2191*WProd2208*WProd2217 - 
     -        2*WProd2192*WProd2208*WProd2217 - 
     -        2*WProd2193*WProd2208*WProd2217 + 
     -        3*WProd2194*WProd2208*WProd2217 + 
     -        3*WProd2195*WProd2208*WProd2217 - 
     -        WProd2217*WProd2218 - WProd2217*WProd2219 - 
     -        WProd2199*WProd2220 - WProd2200*WProd2220 - 
     -        WProd2199*WProd2221 - WProd2200*WProd2221 - 
     -        WProd2199*WProd2222 - WProd2200*WProd2222))/WProd2201
     -     + Prp238*((-16*Alfas*Cpl122*Prp237*
     -         (2*WProd2184*WProd2185 + 2*WProd2198*WProd2199 + 
     -           2*WProd2185*WProd2204 - 
     -           2*WProd2198*WProd2208*WProd2209 - 
     -           2*WProd2198*WProd2208*WProd2210 + 
     -           2*WProd2198*WProd2212 + 2*WProd2198*WProd2213 + 
     -           2*WProd2198*WProd2214 + 
     -           3*WProd2182*WProd2208*WProd2217 + 
     -           2*WProd2183*WProd2208*WProd2217 - 
     -           2*WProd2188*WProd2208*WProd2217 - 
     -           WProd2189*WProd2208*WProd2217 - 
     -           2*WProd2190*WProd2208*WProd2217 - 
     -           WProd2191*WProd2208*WProd2217 - 
     -           3*WProd2192*WProd2208*WProd2217 - 
     -           3*WProd2193*WProd2208*WProd2217 + 
     -           3*WProd2194*WProd2208*WProd2217 + 
     -           3*WProd2195*WProd2208*WProd2217 - 
     -           2*WProd2217*WProd2218 - 2*WProd2199*WProd2220 - 
     -           2*WProd2185*WProd2208*WProd2220 + 
     -           2*WProd2208*WProd2209*WProd2220 + 
     -           2*WProd2208*WProd2210*WProd2220 - 
     -           2*WProd2212*WProd2220 - 2*WProd2213*WProd2220 - 
     -           2*WProd2214*WProd2220 - 2*WProd2199*WProd2221 - 
     -           2*WProd2185*WProd2208*WProd2221 + 
     -           2*WProd2208*WProd2209*WProd2221 + 
     -           2*WProd2208*WProd2210*WProd2221 - 
     -           2*WProd2212*WProd2221 - 2*WProd2213*WProd2221 - 
     -           2*WProd2214*WProd2221 - 2*WProd2199*WProd2222 + 
     -           2*WProd2208*WProd2209*WProd2222 + 
     -           2*WProd2208*WProd2210*WProd2222 - 
     -           2*WProd2212*WProd2222 - 2*WProd2213*WProd2222 - 
     -           2*WProd2214*WProd2222 + 2*WProd2185*WProd2223 + 
     -           2*WProd2185*WProd2224 - 2*WProd2217*WProd2225 - 
     -           2*WProd2217*WProd2226 - 2*WProd2217*WProd2227))/
     -       WProd2201 + (32*Alfas*Cpl122*Prp239*
     -         (WProd2182*WProd2187 + WProd2183*WProd2187 - 
     -           WProd2187*WProd2188 - WProd2187*WProd2190 - 
     -           WProd2187*WProd2192 - WProd2187*WProd2193 + 
     -           WProd2187*WProd2194 + WProd2187*WProd2195 - 
     -           WProd2196 - WProd2206 + WProd2205*WProd2209 + 
     -           WProd2205*WProd2210 - WProd2182*WProd2211 - 
     -           WProd2183*WProd2211 + WProd2188*WProd2211 + 
     -           WProd2190*WProd2211 + WProd2192*WProd2211 + 
     -           WProd2193*WProd2211 - WProd2194*WProd2211 - 
     -           WProd2195*WProd2211 - 
     -           WProd2182*WProd2208*WProd2217 - 
     -           WProd2183*WProd2208*WProd2217 + 
     -           WProd2188*WProd2208*WProd2217 + 
     -           WProd2190*WProd2208*WProd2217 + 
     -           WProd2192*WProd2208*WProd2217 + 
     -           WProd2193*WProd2208*WProd2217 - 
     -           WProd2194*WProd2208*WProd2217 - 
     -           WProd2195*WProd2208*WProd2217 + 
     -           WProd2217*WProd2218 - WProd2215*WProd2220 - 
     -           WProd2215*WProd2221 + WProd2199*WProd2222 - 
     -           WProd2208*WProd2209*WProd2222 - 
     -           WProd2208*WProd2210*WProd2222 + 
     -           WProd2212*WProd2222 + WProd2213*WProd2222 + 
     -           WProd2214*WProd2222 + WProd2217*WProd2225 + 
     -           WProd2217*WProd2226 + WProd2217*WProd2227 - 
     -           WProd2228 - WProd2229 + WProd2220*WProd2230 + 
     -           WProd2221*WProd2230 - WProd2231 - WProd2232 - 
     -           WProd2209*WProd2233 - WProd2210*WProd2233 + 
     -           WProd2220*WProd2234 + WProd2221*WProd2234 + 
     -           WProd2182*WProd2235 + WProd2183*WProd2235 - 
     -           WProd2188*WProd2235 - WProd2190*WProd2235 - 
     -           WProd2192*WProd2235 - WProd2193*WProd2235 + 
     -           WProd2194*WProd2235 + WProd2195*WProd2235))/
     -       WProd2201)
        tmp1160 = tmp1160 - 
     -   (32*Alfas*Cpl122*Prp239*Prp240*
     -      (-2*WProd2182*WProd2187 - 2*WProd2183*WProd2187 + 
     -        2*WProd2187*WProd2188 + WProd2187*WProd2189 + 
     -        2*WProd2187*WProd2190 + WProd2187*WProd2191 + 
     -        2*WProd2187*WProd2192 + WProd2187*WProd2193 - 
     -        2*WProd2187*WProd2194 - 2*WProd2187*WProd2195 + 
     -        WProd2196 + WProd2197 - WProd2207 + 
     -        2*WProd2182*WProd2208*WProd2217 + 
     -        2*WProd2183*WProd2208*WProd2217 - 
     -        2*WProd2188*WProd2208*WProd2217 - 
     -        WProd2189*WProd2208*WProd2217 - 
     -        2*WProd2190*WProd2208*WProd2217 - 
     -        WProd2191*WProd2208*WProd2217 - 
     -        2*WProd2192*WProd2208*WProd2217 - 
     -        WProd2193*WProd2208*WProd2217 + 
     -        2*WProd2194*WProd2208*WProd2217 + 
     -        2*WProd2195*WProd2208*WProd2217 - 
     -        WProd2217*WProd2218 - WProd2217*WProd2219 + 
     -        WProd2215*WProd2220 + WProd2215*WProd2221 - 
     -        WProd2199*WProd2222 - WProd2200*WProd2222 + 
     -        WProd2215*WProd2222 - WProd2217*WProd2226 - 
     -        WProd2217*WProd2227 + WProd2228 + WProd2229 - 
     -        WProd2220*WProd2230 - WProd2221*WProd2230 + 
     -        WProd2231 + WProd2232 - 2*WProd2182*WProd2235 - 
     -        2*WProd2183*WProd2235 + 2*WProd2188*WProd2235 + 
     -        WProd2189*WProd2235 + 2*WProd2190*WProd2235 + 
     -        WProd2191*WProd2235 + 2*WProd2192*WProd2235 + 
     -        WProd2193*WProd2235 - 2*WProd2194*WProd2235 - 
     -        2*WProd2195*WProd2235 + WProd2217*WProd2304))/
     -    WProd2201
	amp(28,9,3) = Prp234*tmp1160

#ifdef DEBUG
	print *, 'amp(28,9,3) =', amp(28,9,3)
#endif

        tmp1161 = (-16*Alfas*Cpl122*Prp235*Prp237*
     -      (2*WProd2184*WProd2237 + 2*WProd2186*WProd2237 + 
     -        2*WProd2198*WProd2241 - 2*WProd2220*WProd2241 - 
     -        2*WProd2221*WProd2241 - 2*WProd2222*WProd2241 + 
     -        2*WProd2198*WProd2242 - 2*WProd2220*WProd2242 - 
     -        2*WProd2221*WProd2242 - 2*WProd2222*WProd2242 + 
     -        4*WProd2182*WProd2208*WProd2252 + 
     -        5*WProd2183*WProd2208*WProd2252 - 
     -        4*WProd2188*WProd2208*WProd2252 - 
     -        4*WProd2189*WProd2208*WProd2252 - 
     -        4*WProd2190*WProd2208*WProd2252 - 
     -        4*WProd2191*WProd2208*WProd2252 - 
     -        4*WProd2192*WProd2208*WProd2252 - 
     -        4*WProd2193*WProd2208*WProd2252 + 
     -        5*WProd2194*WProd2208*WProd2252 + 
     -        5*WProd2195*WProd2208*WProd2252 - 
     -        2*WProd2218*WProd2252 - 2*WProd2219*WProd2252))/
     -    WProd2236 - (32*Alfas*Cpl122*Prp235*Prp240*
     -      (-2*WProd2182*WProd2238 - 3*WProd2183*WProd2238 + 
     -        2*WProd2188*WProd2238 + 2*WProd2189*WProd2238 + 
     -        2*WProd2190*WProd2238 + 2*WProd2191*WProd2238 + 
     -        2*WProd2192*WProd2238 + 2*WProd2193*WProd2238 - 
     -        3*WProd2194*WProd2238 - 3*WProd2195*WProd2238 + 
     -        WProd2239 + WProd2240 - WProd2220*WProd2241 - 
     -        WProd2221*WProd2241 - WProd2222*WProd2241 - 
     -        WProd2220*WProd2242 - WProd2221*WProd2242 - 
     -        WProd2222*WProd2242 + 
     -        2*WProd2182*WProd2208*WProd2252 + 
     -        3*WProd2183*WProd2208*WProd2252 - 
     -        2*WProd2188*WProd2208*WProd2252 - 
     -        2*WProd2189*WProd2208*WProd2252 - 
     -        2*WProd2190*WProd2208*WProd2252 - 
     -        2*WProd2191*WProd2208*WProd2252 - 
     -        2*WProd2192*WProd2208*WProd2252 - 
     -        2*WProd2193*WProd2208*WProd2252 + 
     -        3*WProd2194*WProd2208*WProd2252 + 
     -        3*WProd2195*WProd2208*WProd2252 - 
     -        WProd2218*WProd2252 - WProd2219*WProd2252))/WProd2236
     -     + Prp238*((-16*Alfas*Cpl122*Prp237*
     -         (2*WProd2184*WProd2237 + 2*WProd2204*WProd2237 - 
     -           2*WProd2208*WProd2220*WProd2237 - 
     -           2*WProd2208*WProd2221*WProd2237 + 
     -           2*WProd2223*WProd2237 + 2*WProd2224*WProd2237 + 
     -           2*WProd2198*WProd2241 - 2*WProd2220*WProd2241 - 
     -           2*WProd2221*WProd2241 - 2*WProd2222*WProd2241 - 
     -           2*WProd2198*WProd2208*WProd2245 + 
     -           2*WProd2208*WProd2220*WProd2245 + 
     -           2*WProd2208*WProd2221*WProd2245 + 
     -           2*WProd2208*WProd2222*WProd2245 - 
     -           2*WProd2198*WProd2208*WProd2246 + 
     -           2*WProd2208*WProd2220*WProd2246 + 
     -           2*WProd2208*WProd2221*WProd2246 + 
     -           2*WProd2208*WProd2222*WProd2246 + 
     -           2*WProd2198*WProd2248 - 2*WProd2220*WProd2248 - 
     -           2*WProd2221*WProd2248 - 2*WProd2222*WProd2248 + 
     -           2*WProd2198*WProd2249 - 2*WProd2220*WProd2249 - 
     -           2*WProd2221*WProd2249 - 2*WProd2222*WProd2249 + 
     -           2*WProd2198*WProd2250 - 2*WProd2220*WProd2250 - 
     -           2*WProd2221*WProd2250 - 2*WProd2222*WProd2250 + 
     -           3*WProd2182*WProd2208*WProd2252 + 
     -           2*WProd2183*WProd2208*WProd2252 - 
     -           2*WProd2188*WProd2208*WProd2252 - 
     -           WProd2189*WProd2208*WProd2252 - 
     -           2*WProd2190*WProd2208*WProd2252 - 
     -           WProd2191*WProd2208*WProd2252 - 
     -           3*WProd2192*WProd2208*WProd2252 - 
     -           3*WProd2193*WProd2208*WProd2252 + 
     -           3*WProd2194*WProd2208*WProd2252 + 
     -           3*WProd2195*WProd2208*WProd2252 - 
     -           2*WProd2218*WProd2252 - 2*WProd2225*WProd2252 - 
     -           2*WProd2226*WProd2252 - 2*WProd2227*WProd2252))/
     -       WProd2236 + (32*Alfas*Cpl122*Prp239*
     -         (WProd2182*WProd2238 + WProd2183*WProd2238 - 
     -           WProd2188*WProd2238 - WProd2190*WProd2238 - 
     -           WProd2192*WProd2238 - WProd2193*WProd2238 + 
     -           WProd2194*WProd2238 + WProd2195*WProd2238 - 
     -           WProd2239 + WProd2222*WProd2241 - WProd2243 + 
     -           WProd2205*WProd2245 - 
     -           WProd2208*WProd2222*WProd2245 - 
     -           WProd2233*WProd2245 + WProd2205*WProd2246 - 
     -           WProd2208*WProd2222*WProd2246 - 
     -           WProd2233*WProd2246 - WProd2182*WProd2247 - 
     -           WProd2183*WProd2247 + WProd2188*WProd2247 + 
     -           WProd2190*WProd2247 + WProd2192*WProd2247 + 
     -           WProd2193*WProd2247 - WProd2194*WProd2247 - 
     -           WProd2195*WProd2247 + WProd2222*WProd2248 + 
     -           WProd2222*WProd2249 + WProd2222*WProd2250 - 
     -           WProd2220*WProd2251 - WProd2221*WProd2251 - 
     -           WProd2182*WProd2208*WProd2252 - 
     -           WProd2183*WProd2208*WProd2252 + 
     -           WProd2188*WProd2208*WProd2252 + 
     -           WProd2190*WProd2208*WProd2252 + 
     -           WProd2192*WProd2208*WProd2252 + 
     -           WProd2193*WProd2208*WProd2252 - 
     -           WProd2194*WProd2208*WProd2252 - 
     -           WProd2195*WProd2208*WProd2252 + 
     -           WProd2218*WProd2252 + WProd2225*WProd2252 + 
     -           WProd2226*WProd2252 + WProd2227*WProd2252 - 
     -           WProd2253 - WProd2254 + WProd2220*WProd2255 + 
     -           WProd2221*WProd2255 - WProd2256 - WProd2257 + 
     -           WProd2220*WProd2258 + WProd2221*WProd2258 + 
     -           WProd2182*WProd2259 + WProd2183*WProd2259 - 
     -           WProd2188*WProd2259 - WProd2190*WProd2259 - 
     -           WProd2192*WProd2259 - WProd2193*WProd2259 + 
     -           WProd2194*WProd2259 + WProd2195*WProd2259))/
     -       WProd2236)
        tmp1161 = tmp1161 - 
     -   (32*Alfas*Cpl122*Prp239*Prp240*
     -      (-2*WProd2182*WProd2238 - 2*WProd2183*WProd2238 + 
     -        2*WProd2188*WProd2238 + WProd2189*WProd2238 + 
     -        2*WProd2190*WProd2238 + WProd2191*WProd2238 + 
     -        2*WProd2192*WProd2238 + WProd2193*WProd2238 - 
     -        2*WProd2194*WProd2238 - 2*WProd2195*WProd2238 + 
     -        WProd2239 + WProd2240 - WProd2222*WProd2241 - 
     -        WProd2222*WProd2242 - WProd2244 + 
     -        WProd2220*WProd2251 + WProd2221*WProd2251 + 
     -        WProd2222*WProd2251 + 
     -        2*WProd2182*WProd2208*WProd2252 + 
     -        2*WProd2183*WProd2208*WProd2252 - 
     -        2*WProd2188*WProd2208*WProd2252 - 
     -        WProd2189*WProd2208*WProd2252 - 
     -        2*WProd2190*WProd2208*WProd2252 - 
     -        WProd2191*WProd2208*WProd2252 - 
     -        2*WProd2192*WProd2208*WProd2252 - 
     -        WProd2193*WProd2208*WProd2252 + 
     -        2*WProd2194*WProd2208*WProd2252 + 
     -        2*WProd2195*WProd2208*WProd2252 - 
     -        WProd2218*WProd2252 - WProd2219*WProd2252 - 
     -        WProd2226*WProd2252 - WProd2227*WProd2252 + 
     -        WProd2253 + WProd2254 - WProd2220*WProd2255 - 
     -        WProd2221*WProd2255 + WProd2256 + WProd2257 - 
     -        2*WProd2182*WProd2259 - 2*WProd2183*WProd2259 + 
     -        2*WProd2188*WProd2259 + WProd2189*WProd2259 + 
     -        2*WProd2190*WProd2259 + WProd2191*WProd2259 + 
     -        2*WProd2192*WProd2259 + WProd2193*WProd2259 - 
     -        2*WProd2194*WProd2259 - 2*WProd2195*WProd2259 + 
     -        WProd2252*WProd2304))/WProd2236
	amp(28,9,9) = Prp234*tmp1161

#ifdef DEBUG
	print *, 'amp(28,9,9) =', amp(28,9,9)
#endif

        tmp1162 = (-16*Alfas*Cpl122*Prp235*Prp237*
     -      (2*WProd2185*WProd2260 + 2*WProd2185*WProd2261 + 
     -        2*WProd2199*WProd2265 + 2*WProd2200*WProd2265 + 
     -        4*WProd2182*WProd2208*WProd2276 + 
     -        5*WProd2183*WProd2208*WProd2276 - 
     -        4*WProd2188*WProd2208*WProd2276 - 
     -        4*WProd2189*WProd2208*WProd2276 - 
     -        4*WProd2190*WProd2208*WProd2276 - 
     -        4*WProd2191*WProd2208*WProd2276 - 
     -        4*WProd2192*WProd2208*WProd2276 - 
     -        4*WProd2193*WProd2208*WProd2276 + 
     -        5*WProd2194*WProd2208*WProd2276 + 
     -        5*WProd2195*WProd2208*WProd2276 - 
     -        2*WProd2218*WProd2276 - 2*WProd2219*WProd2276 - 
     -        2*WProd2199*WProd2277 - 2*WProd2200*WProd2277 - 
     -        2*WProd2199*WProd2278 - 2*WProd2200*WProd2278 - 
     -        2*WProd2199*WProd2279 - 2*WProd2200*WProd2279))/
     -    WProd2266 - (32*Alfas*Cpl122*Prp235*Prp240*
     -      (-2*WProd2182*WProd2262 - 3*WProd2183*WProd2262 + 
     -        2*WProd2188*WProd2262 + 2*WProd2189*WProd2262 + 
     -        2*WProd2190*WProd2262 + 2*WProd2191*WProd2262 + 
     -        2*WProd2192*WProd2262 + 2*WProd2193*WProd2262 - 
     -        3*WProd2194*WProd2262 - 3*WProd2195*WProd2262 + 
     -        WProd2263 + WProd2264 + 
     -        2*WProd2182*WProd2208*WProd2276 + 
     -        3*WProd2183*WProd2208*WProd2276 - 
     -        2*WProd2188*WProd2208*WProd2276 - 
     -        2*WProd2189*WProd2208*WProd2276 - 
     -        2*WProd2190*WProd2208*WProd2276 - 
     -        2*WProd2191*WProd2208*WProd2276 - 
     -        2*WProd2192*WProd2208*WProd2276 - 
     -        2*WProd2193*WProd2208*WProd2276 + 
     -        3*WProd2194*WProd2208*WProd2276 + 
     -        3*WProd2195*WProd2208*WProd2276 - 
     -        WProd2218*WProd2276 - WProd2219*WProd2276 - 
     -        WProd2199*WProd2277 - WProd2200*WProd2277 - 
     -        WProd2199*WProd2278 - WProd2200*WProd2278 - 
     -        WProd2199*WProd2279 - WProd2200*WProd2279))/WProd2266
     -     + Prp238*((-16*Alfas*Cpl122*Prp237*
     -         (2*WProd2185*WProd2260 + 2*WProd2199*WProd2265 - 
     -           2*WProd2208*WProd2209*WProd2265 - 
     -           2*WProd2208*WProd2210*WProd2265 + 
     -           2*WProd2212*WProd2265 + 2*WProd2213*WProd2265 + 
     -           2*WProd2214*WProd2265 + 2*WProd2185*WProd2270 + 
     -           3*WProd2182*WProd2208*WProd2276 + 
     -           2*WProd2183*WProd2208*WProd2276 - 
     -           2*WProd2188*WProd2208*WProd2276 - 
     -           WProd2189*WProd2208*WProd2276 - 
     -           2*WProd2190*WProd2208*WProd2276 - 
     -           WProd2191*WProd2208*WProd2276 - 
     -           3*WProd2192*WProd2208*WProd2276 - 
     -           3*WProd2193*WProd2208*WProd2276 + 
     -           3*WProd2194*WProd2208*WProd2276 + 
     -           3*WProd2195*WProd2208*WProd2276 - 
     -           2*WProd2218*WProd2276 - 2*WProd2225*WProd2276 - 
     -           2*WProd2226*WProd2276 - 2*WProd2227*WProd2276 - 
     -           2*WProd2199*WProd2277 - 
     -           2*WProd2185*WProd2208*WProd2277 + 
     -           2*WProd2208*WProd2209*WProd2277 + 
     -           2*WProd2208*WProd2210*WProd2277 - 
     -           2*WProd2212*WProd2277 - 2*WProd2213*WProd2277 - 
     -           2*WProd2214*WProd2277 - 2*WProd2199*WProd2278 - 
     -           2*WProd2185*WProd2208*WProd2278 + 
     -           2*WProd2208*WProd2209*WProd2278 + 
     -           2*WProd2208*WProd2210*WProd2278 - 
     -           2*WProd2212*WProd2278 - 2*WProd2213*WProd2278 - 
     -           2*WProd2214*WProd2278 - 2*WProd2199*WProd2279 + 
     -           2*WProd2208*WProd2209*WProd2279 + 
     -           2*WProd2208*WProd2210*WProd2279 - 
     -           2*WProd2212*WProd2279 - 2*WProd2213*WProd2279 - 
     -           2*WProd2214*WProd2279 + 2*WProd2185*WProd2280 + 
     -           2*WProd2185*WProd2281))/WProd2266 + 
     -      Prp239*((32*Alfas*Cpl122*(WProd2267 + WProd2268))/
     -          WProd2269 - 
     -         (32*Alfas*Cpl122*
     -            (-(WProd2182*WProd2262) - WProd2183*WProd2262 + 
     -              WProd2188*WProd2262 + WProd2190*WProd2262 + 
     -              WProd2192*WProd2262 + WProd2193*WProd2262 - 
     -              WProd2194*WProd2262 - WProd2195*WProd2262 + 
     -              WProd2263 - WProd2209*WProd2271 - 
     -              WProd2210*WProd2271 + WProd2272 + 
     -              WProd2209*WProd2274 + WProd2210*WProd2274 + 
     -              WProd2182*WProd2275 + WProd2183*WProd2275 - 
     -              WProd2188*WProd2275 - WProd2190*WProd2275 - 
     -              WProd2192*WProd2275 - WProd2193*WProd2275 + 
     -              WProd2194*WProd2275 + WProd2195*WProd2275 + 
     -              WProd2182*WProd2208*WProd2276 + 
     -              WProd2183*WProd2208*WProd2276 - 
     -              WProd2188*WProd2208*WProd2276 - 
     -              WProd2190*WProd2208*WProd2276 - 
     -              WProd2192*WProd2208*WProd2276 - 
     -              WProd2193*WProd2208*WProd2276 + 
     -              WProd2194*WProd2208*WProd2276 + 
     -              WProd2195*WProd2208*WProd2276 - 
     -              WProd2218*WProd2276 - WProd2225*WProd2276 - 
     -              WProd2226*WProd2276 - WProd2227*WProd2276 + 
     -              WProd2215*WProd2277 - WProd2230*WProd2277 - 
     -              WProd2234*WProd2277 + WProd2215*WProd2278 - 
     -              WProd2230*WProd2278 - WProd2234*WProd2278 - 
     -              WProd2199*WProd2279 + 
     -              WProd2208*WProd2209*WProd2279 + 
     -              WProd2208*WProd2210*WProd2279 - 
     -              WProd2212*WProd2279 - WProd2213*WProd2279 - 
     -              WProd2214*WProd2279 + WProd2282 + WProd2283 + 
     -              WProd2284 + WProd2285 + WProd2209*WProd2286 + 
     -              WProd2210*WProd2286 - WProd2182*WProd2287 - 
     -              WProd2183*WProd2287 + WProd2188*WProd2287 + 
     -              WProd2190*WProd2287 + WProd2192*WProd2287 + 
     -              WProd2193*WProd2287 - WProd2194*WProd2287 - 
     -              WProd2195*WProd2287))/WProd2266))
        tmp1162 = tmp1162 - 
     -   (32*Alfas*Cpl122*Prp239*Prp240*
     -      (-2*WProd2182*WProd2262 - 2*WProd2183*WProd2262 + 
     -        2*WProd2188*WProd2262 + WProd2189*WProd2262 + 
     -        2*WProd2190*WProd2262 + WProd2191*WProd2262 + 
     -        2*WProd2192*WProd2262 + WProd2193*WProd2262 - 
     -        2*WProd2194*WProd2262 - 2*WProd2195*WProd2262 + 
     -        WProd2263 + WProd2264 - WProd2273 + 
     -        2*WProd2182*WProd2208*WProd2276 + 
     -        2*WProd2183*WProd2208*WProd2276 - 
     -        2*WProd2188*WProd2208*WProd2276 - 
     -        WProd2189*WProd2208*WProd2276 - 
     -        2*WProd2190*WProd2208*WProd2276 - 
     -        WProd2191*WProd2208*WProd2276 - 
     -        2*WProd2192*WProd2208*WProd2276 - 
     -        WProd2193*WProd2208*WProd2276 + 
     -        2*WProd2194*WProd2208*WProd2276 + 
     -        2*WProd2195*WProd2208*WProd2276 - 
     -        WProd2218*WProd2276 - WProd2219*WProd2276 - 
     -        WProd2226*WProd2276 - WProd2227*WProd2276 + 
     -        WProd2215*WProd2277 - WProd2230*WProd2277 + 
     -        WProd2215*WProd2278 - WProd2230*WProd2278 - 
     -        WProd2199*WProd2279 - WProd2200*WProd2279 + 
     -        WProd2215*WProd2279 + WProd2282 + WProd2283 + 
     -        WProd2284 + WProd2285 - 2*WProd2182*WProd2287 - 
     -        2*WProd2183*WProd2287 + 2*WProd2188*WProd2287 + 
     -        WProd2189*WProd2287 + 2*WProd2190*WProd2287 + 
     -        WProd2191*WProd2287 + 2*WProd2192*WProd2287 + 
     -        WProd2193*WProd2287 - 2*WProd2194*WProd2287 - 
     -        2*WProd2195*WProd2287 + WProd2276*WProd2304))/
     -    WProd2266
	amp(28,9,10) = Prp234*tmp1162

#ifdef DEBUG
	print *, 'amp(28,9,10) =', amp(28,9,10)
#endif

        tmp1163 = (-16*Alfas*Cpl122*Prp235*Prp237*
     -      (2*WProd2237*WProd2260 + 2*WProd2237*WProd2261 + 
     -        2*WProd2241*WProd2265 + 2*WProd2242*WProd2265 - 
     -        2*WProd2241*WProd2277 - 2*WProd2242*WProd2277 - 
     -        2*WProd2241*WProd2278 - 2*WProd2242*WProd2278 - 
     -        2*WProd2241*WProd2279 - 2*WProd2242*WProd2279 + 
     -        4*WProd2182*WProd2208*WProd2298 + 
     -        5*WProd2183*WProd2208*WProd2298 - 
     -        4*WProd2188*WProd2208*WProd2298 - 
     -        4*WProd2189*WProd2208*WProd2298 - 
     -        4*WProd2190*WProd2208*WProd2298 - 
     -        4*WProd2191*WProd2208*WProd2298 - 
     -        4*WProd2192*WProd2208*WProd2298 - 
     -        4*WProd2193*WProd2208*WProd2298 + 
     -        5*WProd2194*WProd2208*WProd2298 + 
     -        5*WProd2195*WProd2208*WProd2298 - 
     -        2*WProd2218*WProd2298 - 2*WProd2219*WProd2298))/
     -    WProd2288 - (32*Alfas*Cpl122*Prp235*Prp240*
     -      (-(WProd2241*WProd2277) - WProd2242*WProd2277 - 
     -        WProd2241*WProd2278 - WProd2242*WProd2278 - 
     -        WProd2241*WProd2279 - WProd2242*WProd2279 - 
     -        2*WProd2182*WProd2289 - 3*WProd2183*WProd2289 + 
     -        2*WProd2188*WProd2289 + 2*WProd2189*WProd2289 + 
     -        2*WProd2190*WProd2289 + 2*WProd2191*WProd2289 + 
     -        2*WProd2192*WProd2289 + 2*WProd2193*WProd2289 - 
     -        3*WProd2194*WProd2289 - 3*WProd2195*WProd2289 + 
     -        WProd2290 + WProd2291 + 
     -        2*WProd2182*WProd2208*WProd2298 + 
     -        3*WProd2183*WProd2208*WProd2298 - 
     -        2*WProd2188*WProd2208*WProd2298 - 
     -        2*WProd2189*WProd2208*WProd2298 - 
     -        2*WProd2190*WProd2208*WProd2298 - 
     -        2*WProd2191*WProd2208*WProd2298 - 
     -        2*WProd2192*WProd2208*WProd2298 - 
     -        2*WProd2193*WProd2208*WProd2298 + 
     -        3*WProd2194*WProd2208*WProd2298 + 
     -        3*WProd2195*WProd2208*WProd2298 - 
     -        WProd2218*WProd2298 - WProd2219*WProd2298))/WProd2288
     -     + Prp238*((-16*Alfas*Cpl122*Prp237*
     -         (2*WProd2237*WProd2260 + 2*WProd2241*WProd2265 - 
     -           2*WProd2208*WProd2245*WProd2265 - 
     -           2*WProd2208*WProd2246*WProd2265 + 
     -           2*WProd2248*WProd2265 + 2*WProd2249*WProd2265 + 
     -           2*WProd2250*WProd2265 + 2*WProd2237*WProd2270 - 
     -           2*WProd2208*WProd2237*WProd2277 - 
     -           2*WProd2241*WProd2277 + 
     -           2*WProd2208*WProd2245*WProd2277 + 
     -           2*WProd2208*WProd2246*WProd2277 - 
     -           2*WProd2248*WProd2277 - 2*WProd2249*WProd2277 - 
     -           2*WProd2250*WProd2277 - 
     -           2*WProd2208*WProd2237*WProd2278 - 
     -           2*WProd2241*WProd2278 + 
     -           2*WProd2208*WProd2245*WProd2278 + 
     -           2*WProd2208*WProd2246*WProd2278 - 
     -           2*WProd2248*WProd2278 - 2*WProd2249*WProd2278 - 
     -           2*WProd2250*WProd2278 - 2*WProd2241*WProd2279 + 
     -           2*WProd2208*WProd2245*WProd2279 + 
     -           2*WProd2208*WProd2246*WProd2279 - 
     -           2*WProd2248*WProd2279 - 2*WProd2249*WProd2279 - 
     -           2*WProd2250*WProd2279 + 2*WProd2237*WProd2280 + 
     -           2*WProd2237*WProd2281 + 
     -           3*WProd2182*WProd2208*WProd2298 + 
     -           2*WProd2183*WProd2208*WProd2298 - 
     -           2*WProd2188*WProd2208*WProd2298 - 
     -           WProd2189*WProd2208*WProd2298 - 
     -           2*WProd2190*WProd2208*WProd2298 - 
     -           WProd2191*WProd2208*WProd2298 - 
     -           3*WProd2192*WProd2208*WProd2298 - 
     -           3*WProd2193*WProd2208*WProd2298 + 
     -           3*WProd2194*WProd2208*WProd2298 + 
     -           3*WProd2195*WProd2208*WProd2298 - 
     -           2*WProd2218*WProd2298 - 2*WProd2225*WProd2298 - 
     -           2*WProd2226*WProd2298 - 2*WProd2227*WProd2298))/
     -       WProd2288 + Prp239*
     -       ((32*Alfas*Cpl122*(WProd2293 + WProd2294))/
     -          WProd2292 - 
     -         (32*Alfas*Cpl122*
     -            (-(WProd2245*WProd2271) - WProd2246*WProd2271 + 
     -              WProd2245*WProd2274 + WProd2246*WProd2274 + 
     -              WProd2251*WProd2277 - WProd2255*WProd2277 - 
     -              WProd2258*WProd2277 + WProd2251*WProd2278 - 
     -              WProd2255*WProd2278 - WProd2258*WProd2278 - 
     -              WProd2241*WProd2279 + 
     -              WProd2208*WProd2245*WProd2279 + 
     -              WProd2208*WProd2246*WProd2279 - 
     -              WProd2248*WProd2279 - WProd2249*WProd2279 - 
     -              WProd2250*WProd2279 + WProd2245*WProd2286 + 
     -              WProd2246*WProd2286 - WProd2182*WProd2289 - 
     -              WProd2183*WProd2289 + WProd2188*WProd2289 + 
     -              WProd2190*WProd2289 + WProd2192*WProd2289 + 
     -              WProd2193*WProd2289 - WProd2194*WProd2289 - 
     -              WProd2195*WProd2289 + WProd2290 + WProd2295 + 
     -              WProd2182*WProd2297 + WProd2183*WProd2297 - 
     -              WProd2188*WProd2297 - WProd2190*WProd2297 - 
     -              WProd2192*WProd2297 - WProd2193*WProd2297 + 
     -              WProd2194*WProd2297 + WProd2195*WProd2297 + 
     -              WProd2182*WProd2208*WProd2298 + 
     -              WProd2183*WProd2208*WProd2298 - 
     -              WProd2188*WProd2208*WProd2298 - 
     -              WProd2190*WProd2208*WProd2298 - 
     -              WProd2192*WProd2208*WProd2298 - 
     -              WProd2193*WProd2208*WProd2298 + 
     -              WProd2194*WProd2208*WProd2298 + 
     -              WProd2195*WProd2208*WProd2298 - 
     -              WProd2218*WProd2298 - WProd2225*WProd2298 - 
     -              WProd2226*WProd2298 - WProd2227*WProd2298 + 
     -              WProd2299 + WProd2300 + WProd2301 + 
     -              WProd2302 - WProd2182*WProd2303 - 
     -              WProd2183*WProd2303 + WProd2188*WProd2303 + 
     -              WProd2190*WProd2303 + WProd2192*WProd2303 + 
     -              WProd2193*WProd2303 - WProd2194*WProd2303 - 
     -              WProd2195*WProd2303))/WProd2288))
        tmp1163 = tmp1163 - 
     -   (32*Alfas*Cpl122*Prp239*Prp240*
     -      (WProd2251*WProd2277 - WProd2255*WProd2277 + 
     -        WProd2251*WProd2278 - WProd2255*WProd2278 - 
     -        WProd2241*WProd2279 - WProd2242*WProd2279 + 
     -        WProd2251*WProd2279 - 2*WProd2182*WProd2289 - 
     -        2*WProd2183*WProd2289 + 2*WProd2188*WProd2289 + 
     -        WProd2189*WProd2289 + 2*WProd2190*WProd2289 + 
     -        WProd2191*WProd2289 + 2*WProd2192*WProd2289 + 
     -        WProd2193*WProd2289 - 2*WProd2194*WProd2289 - 
     -        2*WProd2195*WProd2289 + WProd2290 + WProd2291 - 
     -        WProd2296 + 2*WProd2182*WProd2208*WProd2298 + 
     -        2*WProd2183*WProd2208*WProd2298 - 
     -        2*WProd2188*WProd2208*WProd2298 - 
     -        WProd2189*WProd2208*WProd2298 - 
     -        2*WProd2190*WProd2208*WProd2298 - 
     -        WProd2191*WProd2208*WProd2298 - 
     -        2*WProd2192*WProd2208*WProd2298 - 
     -        WProd2193*WProd2208*WProd2298 + 
     -        2*WProd2194*WProd2208*WProd2298 + 
     -        2*WProd2195*WProd2208*WProd2298 - 
     -        WProd2218*WProd2298 - WProd2219*WProd2298 - 
     -        WProd2226*WProd2298 - WProd2227*WProd2298 + 
     -        WProd2299 + WProd2300 + WProd2301 + WProd2302 - 
     -        2*WProd2182*WProd2303 - 2*WProd2183*WProd2303 + 
     -        2*WProd2188*WProd2303 + WProd2189*WProd2303 + 
     -        2*WProd2190*WProd2303 + WProd2191*WProd2303 + 
     -        2*WProd2192*WProd2303 + WProd2193*WProd2303 - 
     -        2*WProd2194*WProd2303 - 2*WProd2195*WProd2303 + 
     -        WProd2298*WProd2304))/WProd2288
	amp(28,9,5) = Prp234*tmp1163

#ifdef DEBUG
	print *, 'amp(28,9,5) =', amp(28,9,5)
#endif

      endif

      end


