      subroutine realQED13(amp,resnr)
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

        amp(13,3,1) = Prp143*
     -   (-16*Alfas*Cpl122*Prp146*Prp147*
     -      (WProd1204 - WProd1206 + WProd1211) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp144*Prp145*
     -        (WProd1203 + WProd1204 + WProd1205 - WProd1206 - 
     -          WProd1207)) + 
     -     Prp148*(-16*Alfas*Cpl122*Prp147*
     -         (-WProd1204 + WProd1206 + WProd1216 - WProd1217 - 
     -           WProd1218 + WProd1219) - 
     -        5.333333333333333D0*
     -         (Alfas*Cpl122*Prp145*
     -           (WProd1203 - WProd1204 + WProd1206 - WProd1216 + 
     -             WProd1217 - WProd1218))))

#ifdef DEBUG
	print *, 'amp(13,3,1) =', amp(13,3,1)
#endif

        amp(13,3,4) = Prp143*
     -   (-16*Alfas*Cpl122*Prp147*Prp148*
     -      (WProd1220 - WProd1221 - WProd1222 + WProd1223 - 
     -        WProd1224 + WProd1225) - 
     -     16*Alfas*Cpl122*Prp146*Prp147*
     -      (WProd1222 - WProd1223 + WProd1226))

#ifdef DEBUG
	print *, 'amp(13,3,4) =', amp(13,3,4)
#endif

        amp(13,3,2) = Prp143*
     -   (-(5.333333333333333D0*
     -        (Alfas*Cpl122*Prp145*Prp148*
     -          (WProd1227 - WProd1228 - WProd1229 + WProd1230 - 
     -            WProd1231 + WProd1232))) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp144*Prp145*
     -        (WProd1227 + WProd1229 - WProd1230 + WProd1233 - 
     -          WProd1234)))

#ifdef DEBUG
	print *, 'amp(13,3,2) =', amp(13,3,2)
#endif

        amp(13,1,1) = Prp143*
     -   (16*Alfas*Cpl122*Prp144*Prp145*
     -      (WProd1203 + WProd1204 + WProd1205 - WProd1206 - 
     -        WProd1207) + 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp146*Prp147*
     -        (WProd1204 - WProd1206 + WProd1211)) + 
     -     Prp148*(16*Alfas*Cpl122*Prp145*
     -         (WProd1203 - WProd1204 + WProd1206 - WProd1216 + 
     -           WProd1217 - WProd1218) + 
     -        5.333333333333333D0*
     -         (Alfas*Cpl122*Prp147*
     -           (-WProd1204 + WProd1206 + WProd1216 - WProd1217 - 
     -             WProd1218 + WProd1219))))

#ifdef DEBUG
	print *, 'amp(13,1,1) =', amp(13,1,1)
#endif

        amp(13,1,4) = Prp143*
     -   (5.333333333333333D0*
     -      (Alfas*Cpl122*Prp147*Prp148*
     -        (WProd1220 - WProd1221 - WProd1222 + WProd1223 - 
     -          WProd1224 + WProd1225)) + 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp146*Prp147*
     -        (WProd1222 - WProd1223 + WProd1226)))

#ifdef DEBUG
	print *, 'amp(13,1,4) =', amp(13,1,4)
#endif

        amp(13,1,2) = Prp143*
     -   (16*Alfas*Cpl122*Prp145*Prp148*
     -      (WProd1227 - WProd1228 - WProd1229 + WProd1230 - 
     -        WProd1231 + WProd1232) + 
     -     16*Alfas*Cpl122*Prp144*Prp145*
     -      (WProd1227 + WProd1229 - WProd1230 + WProd1233 - 
     -        WProd1234))

#ifdef DEBUG
	print *, 'amp(13,1,2) =', amp(13,1,2)
#endif

      endif

      end


