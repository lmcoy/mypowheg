      subroutine realQED16(amp,resnr)
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

        tmp1142 = (-4*Cpl2*GS*Prp159*Prp160*
     -      (-2*WProd1277*WProd1278 - 2*WProd1278*WProd1279 + 
     -        2*WProd1281*WProd1282 - 4*WProd1282*WProd1283 - 
     -        2*WProd1282*WProd1284 - 4*WProd1282*WProd1285 - 
     -        2*WProd1282*WProd1286 + 4*WProd1282*WProd1287 - 
     -        3*WProd1282*WProd1288 - 3*WProd1282*WProd1289 + 
     -        3*WProd1282*WProd1290 + 3*WProd1282*WProd1291 + 
     -        2*WProd1277*WProd1297 + 2*WProd1279*WProd1297 + 
     -        2*WProd1299 - 2*WProd1301 - 2*WProd1302 - 
     -        2*WProd1304 - 2*WProd1305 - 2*WProd1277*WProd1311 - 
     -        2*WProd1279*WProd1311 + 2*WProd1277*WProd1315 + 
     -        2*WProd1279*WProd1315 + 2*WProd1278*WProd1316 - 
     -        2*WProd1318*WProd1319 + 2*WProd1319*WProd1320 + 
     -        2*WProd1319*WProd1321))/WProd1312 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp158*Prp162*
     -       (WProd1278*WProd1280 - WProd1280*WProd1292 - 
     -         WProd1294*WProd1296 + WProd1299 - 
     -         WProd1280*WProd1293*WProd1306 - 
     -         WProd1280*WProd1293*WProd1307 + 
     -         WProd1280*WProd1309 + WProd1280*WProd1310 + 
     -         WProd1280*WProd1311 - WProd1294*WProd1322 + 
     -         WProd1294*WProd1323 - 2*WProd1281*WProd1324 + 
     -         2*WProd1283*WProd1324 + 2*WProd1284*WProd1324 + 
     -         2*WProd1285*WProd1324 + 2*WProd1286*WProd1324 - 
     -         2*WProd1287*WProd1324 + 2*WProd1288*WProd1324 + 
     -         3*WProd1289*WProd1324 - 2*WProd1290*WProd1324 - 
     -         2*WProd1291*WProd1324 + WProd1325 - WProd1326 + 
     -         2*WProd1281*WProd1327 - 2*WProd1283*WProd1327 - 
     -         2*WProd1284*WProd1327 - 2*WProd1285*WProd1327 - 
     -         2*WProd1286*WProd1327 + 2*WProd1287*WProd1327 - 
     -         2*WProd1288*WProd1327 - 3*WProd1289*WProd1327 + 
     -         2*WProd1290*WProd1327 + 2*WProd1291*WProd1327))/
     -     WProd1312
        tmp1142 = tmp1142 + 
     -   Prp164*(5.333333333333333D0*
     -       (Cpl2*GS*Prp165*
     -          (WProd1277*WProd1278 + WProd1278*WProd1279 - 
     -            2*WProd1281*WProd1282 + 2*WProd1282*WProd1283 + 
     -            WProd1282*WProd1284 + 2*WProd1282*WProd1285 + 
     -            WProd1282*WProd1286 - 2*WProd1282*WProd1287 + 
     -            2*WProd1282*WProd1288 + 2*WProd1282*WProd1289 - 
     -            2*WProd1282*WProd1290 - 2*WProd1282*WProd1291 - 
     -            WProd1277*WProd1297 - WProd1279*WProd1297 + 
     -            WProd1301 + WProd1302 + WProd1304 + WProd1305 + 
     -            WProd1295*WProd1306 - WProd1296*WProd1306 + 
     -            WProd1295*WProd1307 - WProd1296*WProd1307 - 
     -            WProd1277*WProd1308 - WProd1279*WProd1308 - 
     -            WProd1292*WProd1316 - 
     -            WProd1293*WProd1306*WProd1316 - 
     -            WProd1293*WProd1307*WProd1316 + 
     -            WProd1309*WProd1316 + WProd1310*WProd1316 + 
     -            WProd1311*WProd1316 + 
     -            2*WProd1281*WProd1293*WProd1319 - 
     -            2*WProd1283*WProd1293*WProd1319 - 
     -            WProd1284*WProd1293*WProd1319 - 
     -            2*WProd1285*WProd1293*WProd1319 - 
     -            WProd1286*WProd1293*WProd1319 + 
     -            2*WProd1287*WProd1293*WProd1319 - 
     -            2*WProd1288*WProd1293*WProd1319 - 
     -            2*WProd1289*WProd1293*WProd1319 + 
     -            2*WProd1290*WProd1293*WProd1319 + 
     -            2*WProd1291*WProd1293*WProd1319 - 
     -            WProd1319*WProd1320 - WProd1319*WProd1321 - 
     -            2*WProd1281*WProd1324 + 2*WProd1283*WProd1324 + 
     -            WProd1284*WProd1324 + 2*WProd1285*WProd1324 + 
     -            WProd1286*WProd1324 - 2*WProd1287*WProd1324 + 
     -            2*WProd1288*WProd1324 + 2*WProd1289*WProd1324 - 
     -            2*WProd1290*WProd1324 - 2*WProd1291*WProd1324 + 
     -            WProd1325 - WProd1326 + 2*WProd1281*WProd1327 - 
     -            2*WProd1283*WProd1327 - WProd1284*WProd1327 - 
     -            2*WProd1285*WProd1327 - WProd1286*WProd1327 + 
     -            2*WProd1287*WProd1327 - 2*WProd1288*WProd1327 - 
     -            2*WProd1289*WProd1327 + 2*WProd1290*WProd1327 + 
     -            2*WProd1291*WProd1327 - WProd1319*WProd1328 + 
     -            WProd1319*WProd1329))/WProd1312 - 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp160*
     -          (-(WProd1277*WProd1278) - WProd1278*WProd1279 + 
     -            2*WProd1281*WProd1282 - 3*WProd1282*WProd1283 - 
     -            2*WProd1282*WProd1284 - 3*WProd1282*WProd1285 - 
     -            2*WProd1282*WProd1286 + 3*WProd1282*WProd1287 - 
     -            3*WProd1282*WProd1288 - 3*WProd1282*WProd1289 + 
     -            3*WProd1282*WProd1290 + 3*WProd1282*WProd1291 + 
     -            WProd1277*WProd1297 + WProd1279*WProd1297 + 
     -            WProd1299 - WProd1301 - WProd1302 - WProd1304 - 
     -            WProd1305 + WProd1278*WProd1316 - 
     -            WProd1311*WProd1316 + WProd1315*WProd1316 - 
     -            2*WProd1281*WProd1293*WProd1319 + 
     -            3*WProd1283*WProd1293*WProd1319 + 
     -            2*WProd1284*WProd1293*WProd1319 + 
     -            3*WProd1285*WProd1293*WProd1319 + 
     -            2*WProd1286*WProd1293*WProd1319 - 
     -            3*WProd1287*WProd1293*WProd1319 + 
     -            3*WProd1288*WProd1293*WProd1319 + 
     -            3*WProd1289*WProd1293*WProd1319 - 
     -            3*WProd1290*WProd1293*WProd1319 - 
     -            3*WProd1291*WProd1293*WProd1319 - 
     -            WProd1318*WProd1319 + WProd1319*WProd1320 + 
     -            WProd1319*WProd1321 + 2*WProd1281*WProd1324 - 
     -            3*WProd1283*WProd1324 - 2*WProd1284*WProd1324 - 
     -            3*WProd1285*WProd1324 - 2*WProd1286*WProd1324 + 
     -            3*WProd1287*WProd1324 - 3*WProd1288*WProd1324 - 
     -            3*WProd1289*WProd1324 + 3*WProd1290*WProd1324 + 
     -            3*WProd1291*WProd1324 - WProd1325 + 
     -            WProd1319*WProd1328 - WProd1319*WProd1331 + 
     -            WProd1332))/WProd1312) + 
     -   Prp163*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp160*
     -          (WProd1277*WProd1311 + WProd1279*WProd1311 - 
     -            WProd1277*WProd1315 - WProd1279*WProd1315 - 
     -            WProd1311*WProd1316 + WProd1315*WProd1316 + 
     -            WProd1283*WProd1293*WProd1319 + 
     -            WProd1285*WProd1293*WProd1319 - 
     -            WProd1287*WProd1293*WProd1319 - 
     -            WProd1283*WProd1324 - WProd1285*WProd1324 + 
     -            WProd1287*WProd1324 - WProd1325 + 
     -            WProd1319*WProd1328 - WProd1319*WProd1331 + 
     -            WProd1332))/WProd1312 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp162*
     -          (-(WProd1280*WProd1311) + WProd1280*WProd1315 + 
     -            WProd1294*WProd1322 + WProd1281*WProd1324 - 
     -            WProd1283*WProd1324 - WProd1284*WProd1324 - 
     -            WProd1285*WProd1324 - WProd1286*WProd1324 + 
     -            WProd1287*WProd1324 - WProd1288*WProd1324 - 
     -            WProd1289*WProd1324 + WProd1290*WProd1324 + 
     -            WProd1291*WProd1324 - WProd1325 + WProd1332 - 
     -            WProd1294*WProd1333))/WProd1312)
        tmp1142 = tmp1142 + 
     -   Prp158*((4*Cpl2*GS*Prp159*
     -         (2*WProd1277*WProd1278 + 2*WProd1278*WProd1279 - 
     -           2*WProd1278*WProd1280 - 4*WProd1281*WProd1282 + 
     -           4*WProd1282*WProd1283 + 3*WProd1282*WProd1284 + 
     -           4*WProd1282*WProd1285 + 3*WProd1282*WProd1286 - 
     -           3*WProd1282*WProd1287 + 4*WProd1282*WProd1288 + 
     -           5*WProd1282*WProd1289 - 3*WProd1282*WProd1290 - 
     -           3*WProd1282*WProd1291 - 2*WProd1277*WProd1292 - 
     -           2*WProd1279*WProd1292 + 
     -           2*WProd1277*WProd1293*WProd1294 + 
     -           2*WProd1279*WProd1293*WProd1294 + 
     -           2*WProd1294*WProd1296 - 2*WProd1277*WProd1297 - 
     -           2*WProd1279*WProd1297 - 2*WProd1299 + 
     -           2*WProd1301 + 2*WProd1302 + 2*WProd1304 + 
     -           2*WProd1305 - 2*WProd1277*WProd1293*WProd1306 - 
     -           2*WProd1279*WProd1293*WProd1306 + 
     -           2*WProd1295*WProd1306 - 2*WProd1296*WProd1306 - 
     -           2*WProd1277*WProd1293*WProd1307 - 
     -           2*WProd1279*WProd1293*WProd1307 + 
     -           2*WProd1295*WProd1307 - 2*WProd1296*WProd1307 - 
     -           2*WProd1277*WProd1308 - 2*WProd1279*WProd1308 + 
     -           2*WProd1277*WProd1309 + 2*WProd1279*WProd1309 + 
     -           2*WProd1277*WProd1310 + 2*WProd1279*WProd1310 + 
     -           2*WProd1277*WProd1311 + 2*WProd1279*WProd1311 - 
     -           2*WProd1294*WProd1334 - 2*WProd1294*WProd1335))/
     -       WProd1312 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp165*
     -          (WProd1277*WProd1278 + WProd1278*WProd1279 - 
     -            2*WProd1281*WProd1282 + 2*WProd1282*WProd1283 + 
     -            WProd1282*WProd1284 + 2*WProd1282*WProd1285 + 
     -            WProd1282*WProd1286 - WProd1282*WProd1287 + 
     -            2*WProd1282*WProd1288 + 2*WProd1282*WProd1289 - 
     -            WProd1282*WProd1290 - WProd1282*WProd1291 - 
     -            WProd1277*WProd1292 - WProd1279*WProd1292 - 
     -            WProd1280*WProd1292 + 
     -            WProd1277*WProd1293*WProd1294 + 
     -            WProd1279*WProd1293*WProd1294 - 
     -            WProd1277*WProd1297 - WProd1279*WProd1297 + 
     -            WProd1301 + WProd1302 + WProd1304 + WProd1305 - 
     -            WProd1277*WProd1293*WProd1306 - 
     -            WProd1279*WProd1293*WProd1306 - 
     -            WProd1280*WProd1293*WProd1306 + 
     -            WProd1295*WProd1306 - WProd1296*WProd1306 - 
     -            WProd1277*WProd1293*WProd1307 - 
     -            WProd1279*WProd1293*WProd1307 - 
     -            WProd1280*WProd1293*WProd1307 + 
     -            WProd1295*WProd1307 - WProd1296*WProd1307 - 
     -            WProd1277*WProd1308 - WProd1279*WProd1308 + 
     -            WProd1277*WProd1309 + WProd1279*WProd1309 + 
     -            WProd1280*WProd1309 + WProd1277*WProd1310 + 
     -            WProd1279*WProd1310 + WProd1280*WProd1310 + 
     -            WProd1277*WProd1311 + WProd1279*WProd1311 + 
     -            WProd1280*WProd1311 - WProd1294*WProd1322 + 
     -            WProd1294*WProd1323 - 2*WProd1281*WProd1324 + 
     -            2*WProd1283*WProd1324 + WProd1284*WProd1324 + 
     -            2*WProd1285*WProd1324 + WProd1286*WProd1324 - 
     -            WProd1287*WProd1324 + 2*WProd1288*WProd1324 + 
     -            2*WProd1289*WProd1324 - WProd1290*WProd1324 - 
     -            WProd1291*WProd1324 + WProd1325 - WProd1326 + 
     -            2*WProd1281*WProd1327 - 2*WProd1283*WProd1327 - 
     -            WProd1284*WProd1327 - 2*WProd1285*WProd1327 - 
     -            WProd1286*WProd1327 + WProd1287*WProd1327 - 
     -            2*WProd1288*WProd1327 - 2*WProd1289*WProd1327 + 
     -            WProd1290*WProd1327 + WProd1291*WProd1327 - 
     -            WProd1294*WProd1334 - WProd1294*WProd1335))/
     -        WProd1312)
        amp(16,4,1) = Prp161*tmp1142 + 
     -   (8*Cpl2*GS*Prp157*Prp158*Prp159*
     -      (WProd1277*WProd1278 + WProd1278*WProd1279 - 
     -        WProd1278*WProd1280 + 2*WProd1281*WProd1282 - 
     -        2*WProd1282*WProd1283 - 3*WProd1282*WProd1284 - 
     -        2*WProd1282*WProd1285 - 2*WProd1282*WProd1286 + 
     -        2*WProd1282*WProd1287 - 2*WProd1282*WProd1288 - 
     -        WProd1282*WProd1289 + 2*WProd1282*WProd1290 + 
     -        2*WProd1282*WProd1291 - WProd1279*WProd1292 + 
     -        WProd1279*WProd1293*WProd1294 - 
     -        WProd1294*WProd1295 + WProd1294*WProd1296 - 
     -        WProd1279*WProd1297 + WProd1298 - WProd1299 - 
     -        WProd1300 + WProd1301 + WProd1302 - WProd1303 + 
     -        WProd1304 + WProd1305 - 
     -        WProd1279*WProd1293*WProd1306 + 
     -        WProd1295*WProd1306 - WProd1296*WProd1306 - 
     -        WProd1279*WProd1293*WProd1307 + 
     -        WProd1295*WProd1307 - WProd1296*WProd1307 - 
     -        WProd1279*WProd1308 + WProd1279*WProd1309 + 
     -        WProd1279*WProd1310 + WProd1279*WProd1311))/WProd1312
     -     - (8*Cpl2*GS*Prp157*Prp159*Prp160*
     -      (-(WProd1277*WProd1278) - WProd1278*WProd1279 - 
     -        2*WProd1281*WProd1282 + 2*WProd1282*WProd1283 + 
     -        3*WProd1282*WProd1284 + 2*WProd1282*WProd1285 + 
     -        3*WProd1282*WProd1286 - 2*WProd1282*WProd1287 + 
     -        2*WProd1282*WProd1288 + 2*WProd1282*WProd1289 - 
     -        2*WProd1282*WProd1290 - 2*WProd1282*WProd1291 + 
     -        WProd1279*WProd1297 + WProd1299 - WProd1301 - 
     -        WProd1302 - WProd1304 - WProd1305 - 
     -        WProd1279*WProd1311 + WProd1279*WProd1315 + 
     -        WProd1278*WProd1316 - WProd1317 - 
     -        WProd1318*WProd1319 + WProd1319*WProd1320 + 
     -        WProd1319*WProd1321))/WProd1312

#ifdef DEBUG
	print *, 'amp(16,4,1) =', amp(16,4,1)
#endif

        tmp1143 = (-4*Cpl2*GS*Prp159*Prp160*
     -      (-2*WProd1277*WProd1337 - 2*WProd1279*WProd1337 + 
     -        2*WProd1316*WProd1337 + 2*WProd1281*WProd1338 - 
     -        4*WProd1283*WProd1338 - 2*WProd1284*WProd1338 - 
     -        4*WProd1285*WProd1338 - 2*WProd1286*WProd1338 + 
     -        4*WProd1287*WProd1338 - 3*WProd1288*WProd1338 - 
     -        3*WProd1289*WProd1338 + 3*WProd1290*WProd1338 + 
     -        3*WProd1291*WProd1338 + 2*WProd1277*WProd1341 + 
     -        2*WProd1279*WProd1341 + 2*WProd1343 - 2*WProd1345 - 
     -        2*WProd1346 - 2*WProd1348 - 2*WProd1349 - 
     -        2*WProd1277*WProd1355 - 2*WProd1279*WProd1355 + 
     -        2*WProd1277*WProd1356 + 2*WProd1279*WProd1356 - 
     -        2*WProd1318*WProd1358 + 2*WProd1320*WProd1358 + 
     -        2*WProd1321*WProd1358))/WProd1336 + 
     -   Prp158*((4*Cpl2*GS*Prp159*
     -         (2*WProd1277*WProd1337 + 2*WProd1279*WProd1337 - 
     -           2*WProd1280*WProd1337 - 4*WProd1281*WProd1338 + 
     -           4*WProd1283*WProd1338 + 3*WProd1284*WProd1338 + 
     -           4*WProd1285*WProd1338 + 3*WProd1286*WProd1338 - 
     -           3*WProd1287*WProd1338 + 4*WProd1288*WProd1338 + 
     -           5*WProd1289*WProd1338 - 3*WProd1290*WProd1338 - 
     -           3*WProd1291*WProd1338 - 2*WProd1277*WProd1339 - 
     -           2*WProd1279*WProd1339 + 
     -           2*WProd1277*WProd1293*WProd1340 + 
     -           2*WProd1279*WProd1293*WProd1340 + 
     -           2*WProd1296*WProd1340 - 2*WProd1334*WProd1340 - 
     -           2*WProd1335*WProd1340 - 2*WProd1277*WProd1341 - 
     -           2*WProd1279*WProd1341 - 2*WProd1343 + 
     -           2*WProd1345 + 2*WProd1346 + 2*WProd1348 + 
     -           2*WProd1349 - 2*WProd1277*WProd1293*WProd1350 - 
     -           2*WProd1279*WProd1293*WProd1350 + 
     -           2*WProd1295*WProd1350 - 2*WProd1296*WProd1350 - 
     -           2*WProd1277*WProd1293*WProd1351 - 
     -           2*WProd1279*WProd1293*WProd1351 + 
     -           2*WProd1295*WProd1351 - 2*WProd1296*WProd1351 - 
     -           2*WProd1277*WProd1352 - 2*WProd1279*WProd1352 + 
     -           2*WProd1277*WProd1353 + 2*WProd1279*WProd1353 + 
     -           2*WProd1277*WProd1354 + 2*WProd1279*WProd1354 + 
     -           2*WProd1277*WProd1355 + 2*WProd1279*WProd1355))/
     -       WProd1336 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp165*
     -          (WProd1277*WProd1337 + WProd1279*WProd1337 - 
     -            2*WProd1281*WProd1338 + 2*WProd1283*WProd1338 + 
     -            WProd1284*WProd1338 + 2*WProd1285*WProd1338 + 
     -            WProd1286*WProd1338 - WProd1287*WProd1338 + 
     -            2*WProd1288*WProd1338 + 2*WProd1289*WProd1338 - 
     -            WProd1290*WProd1338 - WProd1291*WProd1338 - 
     -            WProd1277*WProd1339 - WProd1279*WProd1339 - 
     -            WProd1280*WProd1339 + 
     -            WProd1277*WProd1293*WProd1340 + 
     -            WProd1279*WProd1293*WProd1340 - 
     -            WProd1322*WProd1340 + WProd1323*WProd1340 - 
     -            WProd1334*WProd1340 - WProd1335*WProd1340 - 
     -            WProd1277*WProd1341 - WProd1279*WProd1341 + 
     -            WProd1345 + WProd1346 + WProd1348 + WProd1349 - 
     -            WProd1277*WProd1293*WProd1350 - 
     -            WProd1279*WProd1293*WProd1350 - 
     -            WProd1280*WProd1293*WProd1350 + 
     -            WProd1295*WProd1350 - WProd1296*WProd1350 - 
     -            WProd1277*WProd1293*WProd1351 - 
     -            WProd1279*WProd1293*WProd1351 - 
     -            WProd1280*WProd1293*WProd1351 + 
     -            WProd1295*WProd1351 - WProd1296*WProd1351 - 
     -            WProd1277*WProd1352 - WProd1279*WProd1352 + 
     -            WProd1277*WProd1353 + WProd1279*WProd1353 + 
     -            WProd1280*WProd1353 + WProd1277*WProd1354 + 
     -            WProd1279*WProd1354 + WProd1280*WProd1354 + 
     -            WProd1277*WProd1355 + WProd1279*WProd1355 + 
     -            WProd1280*WProd1355 - 2*WProd1281*WProd1359 + 
     -            2*WProd1283*WProd1359 + WProd1284*WProd1359 + 
     -            2*WProd1285*WProd1359 + WProd1286*WProd1359 - 
     -            WProd1287*WProd1359 + 2*WProd1288*WProd1359 + 
     -            2*WProd1289*WProd1359 - WProd1290*WProd1359 - 
     -            WProd1291*WProd1359 + WProd1360 - WProd1361 + 
     -            2*WProd1281*WProd1362 - 2*WProd1283*WProd1362 - 
     -            WProd1284*WProd1362 - 2*WProd1285*WProd1362 - 
     -            WProd1286*WProd1362 + WProd1287*WProd1362 - 
     -            2*WProd1288*WProd1362 - 2*WProd1289*WProd1362 + 
     -            WProd1290*WProd1362 + WProd1291*WProd1362))/
     -        WProd1336) - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp158*Prp162*
     -       (WProd1280*WProd1337 - WProd1280*WProd1339 - 
     -         WProd1296*WProd1340 - WProd1322*WProd1340 + 
     -         WProd1323*WProd1340 + WProd1343 - 
     -         WProd1280*WProd1293*WProd1350 - 
     -         WProd1280*WProd1293*WProd1351 + 
     -         WProd1280*WProd1353 + WProd1280*WProd1354 + 
     -         WProd1280*WProd1355 - 2*WProd1281*WProd1359 + 
     -         2*WProd1283*WProd1359 + 2*WProd1284*WProd1359 + 
     -         2*WProd1285*WProd1359 + 2*WProd1286*WProd1359 - 
     -         2*WProd1287*WProd1359 + 2*WProd1288*WProd1359 + 
     -         3*WProd1289*WProd1359 - 2*WProd1290*WProd1359 - 
     -         2*WProd1291*WProd1359 + WProd1360 - WProd1361 + 
     -         2*WProd1281*WProd1362 - 2*WProd1283*WProd1362 - 
     -         2*WProd1284*WProd1362 - 2*WProd1285*WProd1362 - 
     -         2*WProd1286*WProd1362 + 2*WProd1287*WProd1362 - 
     -         2*WProd1288*WProd1362 - 3*WProd1289*WProd1362 + 
     -         2*WProd1290*WProd1362 + 2*WProd1291*WProd1362))/
     -     WProd1336
        tmp1143 = tmp1143 + 
     -   Prp163*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp160*
     -          (WProd1277*WProd1355 + WProd1279*WProd1355 - 
     -            WProd1316*WProd1355 - WProd1277*WProd1356 - 
     -            WProd1279*WProd1356 + WProd1316*WProd1356 + 
     -            WProd1283*WProd1293*WProd1358 + 
     -            WProd1285*WProd1293*WProd1358 - 
     -            WProd1287*WProd1293*WProd1358 + 
     -            WProd1328*WProd1358 - WProd1331*WProd1358 - 
     -            WProd1283*WProd1359 - WProd1285*WProd1359 + 
     -            WProd1287*WProd1359 - WProd1360 + WProd1363))/
     -        WProd1336 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp162*
     -          (WProd1322*WProd1340 - WProd1333*WProd1340 - 
     -            WProd1280*WProd1355 + WProd1280*WProd1356 + 
     -            WProd1281*WProd1359 - WProd1283*WProd1359 - 
     -            WProd1284*WProd1359 - WProd1285*WProd1359 - 
     -            WProd1286*WProd1359 + WProd1287*WProd1359 - 
     -            WProd1288*WProd1359 - WProd1289*WProd1359 + 
     -            WProd1290*WProd1359 + WProd1291*WProd1359 - 
     -            WProd1360 + WProd1363))/WProd1336) + 
     -   Prp164*(5.333333333333333D0*
     -       (Cpl2*GS*Prp165*
     -          (WProd1277*WProd1337 + WProd1279*WProd1337 - 
     -            2*WProd1281*WProd1338 + 2*WProd1283*WProd1338 + 
     -            WProd1284*WProd1338 + 2*WProd1285*WProd1338 + 
     -            WProd1286*WProd1338 - 2*WProd1287*WProd1338 + 
     -            2*WProd1288*WProd1338 + 2*WProd1289*WProd1338 - 
     -            2*WProd1290*WProd1338 - 2*WProd1291*WProd1338 - 
     -            WProd1316*WProd1339 - WProd1277*WProd1341 - 
     -            WProd1279*WProd1341 + WProd1345 + WProd1346 + 
     -            WProd1348 + WProd1349 + WProd1295*WProd1350 - 
     -            WProd1296*WProd1350 - 
     -            WProd1293*WProd1316*WProd1350 + 
     -            WProd1295*WProd1351 - WProd1296*WProd1351 - 
     -            WProd1293*WProd1316*WProd1351 - 
     -            WProd1277*WProd1352 - WProd1279*WProd1352 + 
     -            WProd1316*WProd1353 + WProd1316*WProd1354 + 
     -            WProd1316*WProd1355 + 
     -            2*WProd1281*WProd1293*WProd1358 - 
     -            2*WProd1283*WProd1293*WProd1358 - 
     -            WProd1284*WProd1293*WProd1358 - 
     -            2*WProd1285*WProd1293*WProd1358 - 
     -            WProd1286*WProd1293*WProd1358 + 
     -            2*WProd1287*WProd1293*WProd1358 - 
     -            2*WProd1288*WProd1293*WProd1358 - 
     -            2*WProd1289*WProd1293*WProd1358 + 
     -            2*WProd1290*WProd1293*WProd1358 + 
     -            2*WProd1291*WProd1293*WProd1358 - 
     -            WProd1320*WProd1358 - WProd1321*WProd1358 - 
     -            WProd1328*WProd1358 + WProd1329*WProd1358 - 
     -            2*WProd1281*WProd1359 + 2*WProd1283*WProd1359 + 
     -            WProd1284*WProd1359 + 2*WProd1285*WProd1359 + 
     -            WProd1286*WProd1359 - 2*WProd1287*WProd1359 + 
     -            2*WProd1288*WProd1359 + 2*WProd1289*WProd1359 - 
     -            2*WProd1290*WProd1359 - 2*WProd1291*WProd1359 + 
     -            WProd1360 - WProd1361 + 2*WProd1281*WProd1362 - 
     -            2*WProd1283*WProd1362 - WProd1284*WProd1362 - 
     -            2*WProd1285*WProd1362 - WProd1286*WProd1362 + 
     -            2*WProd1287*WProd1362 - 2*WProd1288*WProd1362 - 
     -            2*WProd1289*WProd1362 + 2*WProd1290*WProd1362 + 
     -            2*WProd1291*WProd1362))/WProd1336 - 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp160*
     -          (-(WProd1277*WProd1337) - WProd1279*WProd1337 + 
     -            WProd1316*WProd1337 + 2*WProd1281*WProd1338 - 
     -            3*WProd1283*WProd1338 - 2*WProd1284*WProd1338 - 
     -            3*WProd1285*WProd1338 - 2*WProd1286*WProd1338 + 
     -            3*WProd1287*WProd1338 - 3*WProd1288*WProd1338 - 
     -            3*WProd1289*WProd1338 + 3*WProd1290*WProd1338 + 
     -            3*WProd1291*WProd1338 + WProd1277*WProd1341 + 
     -            WProd1279*WProd1341 + WProd1343 - WProd1345 - 
     -            WProd1346 - WProd1348 - WProd1349 - 
     -            WProd1316*WProd1355 + WProd1316*WProd1356 - 
     -            2*WProd1281*WProd1293*WProd1358 + 
     -            3*WProd1283*WProd1293*WProd1358 + 
     -            2*WProd1284*WProd1293*WProd1358 + 
     -            3*WProd1285*WProd1293*WProd1358 + 
     -            2*WProd1286*WProd1293*WProd1358 - 
     -            3*WProd1287*WProd1293*WProd1358 + 
     -            3*WProd1288*WProd1293*WProd1358 + 
     -            3*WProd1289*WProd1293*WProd1358 - 
     -            3*WProd1290*WProd1293*WProd1358 - 
     -            3*WProd1291*WProd1293*WProd1358 - 
     -            WProd1318*WProd1358 + WProd1320*WProd1358 + 
     -            WProd1321*WProd1358 + WProd1328*WProd1358 - 
     -            WProd1331*WProd1358 + 2*WProd1281*WProd1359 - 
     -            3*WProd1283*WProd1359 - 2*WProd1284*WProd1359 - 
     -            3*WProd1285*WProd1359 - 2*WProd1286*WProd1359 + 
     -            3*WProd1287*WProd1359 - 3*WProd1288*WProd1359 - 
     -            3*WProd1289*WProd1359 + 3*WProd1290*WProd1359 + 
     -            3*WProd1291*WProd1359 - WProd1360 + WProd1363))/
     -        WProd1336)
        amp(16,4,8) = Prp161*tmp1143 + 
     -   (8*Cpl2*GS*Prp157*Prp158*Prp159*
     -      (WProd1277*WProd1337 + WProd1279*WProd1337 - 
     -        WProd1280*WProd1337 + 2*WProd1281*WProd1338 - 
     -        2*WProd1283*WProd1338 - 3*WProd1284*WProd1338 - 
     -        2*WProd1285*WProd1338 - 2*WProd1286*WProd1338 + 
     -        2*WProd1287*WProd1338 - 2*WProd1288*WProd1338 - 
     -        WProd1289*WProd1338 + 2*WProd1290*WProd1338 + 
     -        2*WProd1291*WProd1338 - WProd1279*WProd1339 + 
     -        WProd1279*WProd1293*WProd1340 - 
     -        WProd1295*WProd1340 + WProd1296*WProd1340 - 
     -        WProd1279*WProd1341 + WProd1342 - WProd1343 - 
     -        WProd1344 + WProd1345 + WProd1346 - WProd1347 + 
     -        WProd1348 + WProd1349 - 
     -        WProd1279*WProd1293*WProd1350 + 
     -        WProd1295*WProd1350 - WProd1296*WProd1350 - 
     -        WProd1279*WProd1293*WProd1351 + 
     -        WProd1295*WProd1351 - WProd1296*WProd1351 - 
     -        WProd1279*WProd1352 + WProd1279*WProd1353 + 
     -        WProd1279*WProd1354 + WProd1279*WProd1355))/WProd1336
     -     - (8*Cpl2*GS*Prp157*Prp159*Prp160*
     -      (-(WProd1277*WProd1337) - WProd1279*WProd1337 + 
     -        WProd1316*WProd1337 - 2*WProd1281*WProd1338 + 
     -        2*WProd1283*WProd1338 + 3*WProd1284*WProd1338 + 
     -        2*WProd1285*WProd1338 + 3*WProd1286*WProd1338 - 
     -        2*WProd1287*WProd1338 + 2*WProd1288*WProd1338 + 
     -        2*WProd1289*WProd1338 - 2*WProd1290*WProd1338 - 
     -        2*WProd1291*WProd1338 + WProd1279*WProd1341 + 
     -        WProd1343 - WProd1345 - WProd1346 - WProd1348 - 
     -        WProd1349 - WProd1279*WProd1355 + 
     -        WProd1279*WProd1356 - WProd1357 - 
     -        WProd1318*WProd1358 + WProd1320*WProd1358 + 
     -        WProd1321*WProd1358))/WProd1336

#ifdef DEBUG
	print *, 'amp(16,4,8) =', amp(16,4,8)
#endif

        tmp1144 = (-4*Cpl2*GS*Prp159*Prp160*
     -      (-2*WProd1278*WProd1368 + 2*WProd1297*WProd1368 - 
     -        2*WProd1311*WProd1368 + 2*WProd1315*WProd1368 - 
     -        2*WProd1278*WProd1369 + 2*WProd1297*WProd1369 - 
     -        2*WProd1311*WProd1369 + 2*WProd1315*WProd1369 + 
     -        2*WProd1281*WProd1371 - 4*WProd1283*WProd1371 - 
     -        2*WProd1284*WProd1371 - 4*WProd1285*WProd1371 - 
     -        2*WProd1286*WProd1371 + 4*WProd1287*WProd1371 - 
     -        3*WProd1288*WProd1371 - 3*WProd1289*WProd1371 + 
     -        3*WProd1290*WProd1371 + 3*WProd1291*WProd1371 - 
     -        2*WProd1375 - 2*WProd1376 - 2*WProd1378 - 
     -        2*WProd1380 + 2*WProd1278*WProd1382 - 
     -        2*WProd1318*WProd1384 + 2*WProd1320*WProd1384 + 
     -        2*WProd1321*WProd1384 + 2*WProd1385))/WProd1381 + 
     -   Prp158*Prp162*(2.6666666666666665D0*
     -       (Cpl2*GS*(WProd1386 + WProd1387))/WProd1367 - 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*(WProd1278*WProd1370 - WProd1292*WProd1370 - 
     -            WProd1293*WProd1306*WProd1370 - 
     -            WProd1293*WProd1307*WProd1370 + 
     -            WProd1309*WProd1370 + WProd1310*WProd1370 + 
     -            WProd1311*WProd1370 - WProd1294*WProd1373 + 
     -            WProd1385 - WProd1294*WProd1388 + 
     -            WProd1294*WProd1389 - 2*WProd1281*WProd1390 + 
     -            2*WProd1283*WProd1390 + 2*WProd1284*WProd1390 + 
     -            2*WProd1285*WProd1390 + 2*WProd1286*WProd1390 - 
     -            2*WProd1287*WProd1390 + 2*WProd1288*WProd1390 + 
     -            3*WProd1289*WProd1390 - 2*WProd1290*WProd1390 - 
     -            2*WProd1291*WProd1390 + WProd1391 - WProd1392 + 
     -            WProd1306*WProd1393 + WProd1307*WProd1393 + 
     -            2*WProd1281*WProd1394 - 2*WProd1283*WProd1394 - 
     -            2*WProd1284*WProd1394 - 2*WProd1285*WProd1394 - 
     -            2*WProd1286*WProd1394 + 2*WProd1287*WProd1394 - 
     -            2*WProd1288*WProd1394 - 3*WProd1289*WProd1394 + 
     -            2*WProd1290*WProd1394 + 2*WProd1291*WProd1394))/
     -        WProd1381)
        tmp1144 = tmp1144 + 
     -   Prp164*(Prp165*(-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1386 + WProd1387))/WProd1367) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1278*WProd1368 - 
     -               WProd1297*WProd1368 - WProd1308*WProd1368 + 
     -               WProd1278*WProd1369 - WProd1297*WProd1369 - 
     -               WProd1308*WProd1369 - 2*WProd1281*WProd1371 + 
     -               2*WProd1283*WProd1371 + WProd1284*WProd1371 + 
     -               2*WProd1285*WProd1371 + WProd1286*WProd1371 - 
     -               2*WProd1287*WProd1371 + 
     -               2*WProd1288*WProd1371 + 
     -               2*WProd1289*WProd1371 - 
     -               2*WProd1290*WProd1371 - 
     -               2*WProd1291*WProd1371 - WProd1306*WProd1373 - 
     -               WProd1307*WProd1373 + WProd1306*WProd1374 + 
     -               WProd1307*WProd1374 + WProd1375 + WProd1376 + 
     -               WProd1378 + WProd1380 - WProd1292*WProd1382 - 
     -               WProd1293*WProd1306*WProd1382 - 
     -               WProd1293*WProd1307*WProd1382 + 
     -               WProd1309*WProd1382 + WProd1310*WProd1382 + 
     -               WProd1311*WProd1382 + 
     -               2*WProd1281*WProd1293*WProd1384 - 
     -               2*WProd1283*WProd1293*WProd1384 - 
     -               WProd1284*WProd1293*WProd1384 - 
     -               2*WProd1285*WProd1293*WProd1384 - 
     -               WProd1286*WProd1293*WProd1384 + 
     -               2*WProd1287*WProd1293*WProd1384 - 
     -               2*WProd1288*WProd1293*WProd1384 - 
     -               2*WProd1289*WProd1293*WProd1384 + 
     -               2*WProd1290*WProd1293*WProd1384 + 
     -               2*WProd1291*WProd1293*WProd1384 - 
     -               WProd1320*WProd1384 - WProd1321*WProd1384 - 
     -               WProd1328*WProd1384 + WProd1329*WProd1384 - 
     -               2*WProd1281*WProd1390 + 
     -               2*WProd1283*WProd1390 + WProd1284*WProd1390 + 
     -               2*WProd1285*WProd1390 + WProd1286*WProd1390 - 
     -               2*WProd1287*WProd1390 + 
     -               2*WProd1288*WProd1390 + 
     -               2*WProd1289*WProd1390 - 
     -               2*WProd1290*WProd1390 - 
     -               2*WProd1291*WProd1390 + WProd1391 - 
     -               WProd1392 + WProd1306*WProd1393 + 
     -               WProd1307*WProd1393 + 2*WProd1281*WProd1394 - 
     -               2*WProd1283*WProd1394 - WProd1284*WProd1394 - 
     -               2*WProd1285*WProd1394 - WProd1286*WProd1394 + 
     -               2*WProd1287*WProd1394 - 
     -               2*WProd1288*WProd1394 - 
     -               2*WProd1289*WProd1394 + 
     -               2*WProd1290*WProd1394 + 2*WProd1291*WProd1394)
     -             )/WProd1381) - 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp160*
     -          (-(WProd1278*WProd1368) + WProd1297*WProd1368 - 
     -            WProd1278*WProd1369 + WProd1297*WProd1369 + 
     -            2*WProd1281*WProd1371 - 3*WProd1283*WProd1371 - 
     -            2*WProd1284*WProd1371 - 3*WProd1285*WProd1371 - 
     -            2*WProd1286*WProd1371 + 3*WProd1287*WProd1371 - 
     -            3*WProd1288*WProd1371 - 3*WProd1289*WProd1371 + 
     -            3*WProd1290*WProd1371 + 3*WProd1291*WProd1371 - 
     -            WProd1375 - WProd1376 - WProd1378 - WProd1380 + 
     -            WProd1278*WProd1382 - WProd1311*WProd1382 + 
     -            WProd1315*WProd1382 - 
     -            2*WProd1281*WProd1293*WProd1384 + 
     -            3*WProd1283*WProd1293*WProd1384 + 
     -            2*WProd1284*WProd1293*WProd1384 + 
     -            3*WProd1285*WProd1293*WProd1384 + 
     -            2*WProd1286*WProd1293*WProd1384 - 
     -            3*WProd1287*WProd1293*WProd1384 + 
     -            3*WProd1288*WProd1293*WProd1384 + 
     -            3*WProd1289*WProd1293*WProd1384 - 
     -            3*WProd1290*WProd1293*WProd1384 - 
     -            3*WProd1291*WProd1293*WProd1384 - 
     -            WProd1318*WProd1384 + WProd1320*WProd1384 + 
     -            WProd1321*WProd1384 + WProd1328*WProd1384 - 
     -            WProd1331*WProd1384 + WProd1385 + 
     -            2*WProd1281*WProd1390 - 3*WProd1283*WProd1390 - 
     -            2*WProd1284*WProd1390 - 3*WProd1285*WProd1390 - 
     -            2*WProd1286*WProd1390 + 3*WProd1287*WProd1390 - 
     -            3*WProd1288*WProd1390 - 3*WProd1289*WProd1390 + 
     -            3*WProd1290*WProd1390 + 3*WProd1291*WProd1390 - 
     -            WProd1391 + WProd1395))/WProd1381) + 
     -   Prp163*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp160*
     -          (WProd1311*WProd1368 - WProd1315*WProd1368 + 
     -            WProd1311*WProd1369 - WProd1315*WProd1369 - 
     -            WProd1311*WProd1382 + WProd1315*WProd1382 + 
     -            WProd1283*WProd1293*WProd1384 + 
     -            WProd1285*WProd1293*WProd1384 - 
     -            WProd1287*WProd1293*WProd1384 + 
     -            WProd1328*WProd1384 - WProd1331*WProd1384 - 
     -            WProd1283*WProd1390 - WProd1285*WProd1390 + 
     -            WProd1287*WProd1390 - WProd1391 + WProd1395))/
     -        WProd1381 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp162*
     -          (-(WProd1311*WProd1370) + WProd1315*WProd1370 + 
     -            WProd1294*WProd1388 + WProd1281*WProd1390 - 
     -            WProd1283*WProd1390 - WProd1284*WProd1390 - 
     -            WProd1285*WProd1390 - WProd1286*WProd1390 + 
     -            WProd1287*WProd1390 - WProd1288*WProd1390 - 
     -            WProd1289*WProd1390 + WProd1290*WProd1390 + 
     -            WProd1291*WProd1390 - WProd1391 + WProd1395 - 
     -            WProd1294*WProd1396))/WProd1381)
        tmp1144 = tmp1144 + 
     -   Prp158*((4*Cpl2*GS*Prp159*
     -         (2*WProd1278*WProd1368 - 2*WProd1292*WProd1368 + 
     -           2*WProd1293*WProd1294*WProd1368 - 
     -           2*WProd1297*WProd1368 - 
     -           2*WProd1293*WProd1306*WProd1368 - 
     -           2*WProd1293*WProd1307*WProd1368 - 
     -           2*WProd1308*WProd1368 + 2*WProd1309*WProd1368 + 
     -           2*WProd1310*WProd1368 + 2*WProd1311*WProd1368 + 
     -           2*WProd1278*WProd1369 - 2*WProd1292*WProd1369 + 
     -           2*WProd1293*WProd1294*WProd1369 - 
     -           2*WProd1297*WProd1369 - 
     -           2*WProd1293*WProd1306*WProd1369 - 
     -           2*WProd1293*WProd1307*WProd1369 - 
     -           2*WProd1308*WProd1369 + 2*WProd1309*WProd1369 + 
     -           2*WProd1310*WProd1369 + 2*WProd1311*WProd1369 - 
     -           2*WProd1278*WProd1370 - 4*WProd1281*WProd1371 + 
     -           4*WProd1283*WProd1371 + 3*WProd1284*WProd1371 + 
     -           4*WProd1285*WProd1371 + 3*WProd1286*WProd1371 - 
     -           3*WProd1287*WProd1371 + 4*WProd1288*WProd1371 + 
     -           5*WProd1289*WProd1371 - 3*WProd1290*WProd1371 - 
     -           3*WProd1291*WProd1371 + 2*WProd1294*WProd1373 - 
     -           2*WProd1306*WProd1373 - 2*WProd1307*WProd1373 + 
     -           2*WProd1306*WProd1374 + 2*WProd1307*WProd1374 + 
     -           2*WProd1375 + 2*WProd1376 + 2*WProd1378 + 
     -           2*WProd1380 - 2*WProd1385 - 
     -           2*WProd1294*WProd1397 - 2*WProd1294*WProd1398))/
     -       WProd1381 + Prp165*
     -       (-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1386 + WProd1387))/WProd1367) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1278*WProd1368 - 
     -               WProd1292*WProd1368 + 
     -               WProd1293*WProd1294*WProd1368 - 
     -               WProd1297*WProd1368 - 
     -               WProd1293*WProd1306*WProd1368 - 
     -               WProd1293*WProd1307*WProd1368 - 
     -               WProd1308*WProd1368 + WProd1309*WProd1368 + 
     -               WProd1310*WProd1368 + WProd1311*WProd1368 + 
     -               WProd1278*WProd1369 - WProd1292*WProd1369 + 
     -               WProd1293*WProd1294*WProd1369 - 
     -               WProd1297*WProd1369 - 
     -               WProd1293*WProd1306*WProd1369 - 
     -               WProd1293*WProd1307*WProd1369 - 
     -               WProd1308*WProd1369 + WProd1309*WProd1369 + 
     -               WProd1310*WProd1369 + WProd1311*WProd1369 - 
     -               WProd1292*WProd1370 - 
     -               WProd1293*WProd1306*WProd1370 - 
     -               WProd1293*WProd1307*WProd1370 + 
     -               WProd1309*WProd1370 + WProd1310*WProd1370 + 
     -               WProd1311*WProd1370 - 2*WProd1281*WProd1371 + 
     -               2*WProd1283*WProd1371 + WProd1284*WProd1371 + 
     -               2*WProd1285*WProd1371 + WProd1286*WProd1371 - 
     -               WProd1287*WProd1371 + 2*WProd1288*WProd1371 + 
     -               2*WProd1289*WProd1371 - WProd1290*WProd1371 - 
     -               WProd1291*WProd1371 - WProd1306*WProd1373 - 
     -               WProd1307*WProd1373 + WProd1306*WProd1374 + 
     -               WProd1307*WProd1374 + WProd1375 + WProd1376 + 
     -               WProd1378 + WProd1380 - WProd1294*WProd1388 + 
     -               WProd1294*WProd1389 - 2*WProd1281*WProd1390 + 
     -               2*WProd1283*WProd1390 + WProd1284*WProd1390 + 
     -               2*WProd1285*WProd1390 + WProd1286*WProd1390 - 
     -               WProd1287*WProd1390 + 2*WProd1288*WProd1390 + 
     -               2*WProd1289*WProd1390 - WProd1290*WProd1390 - 
     -               WProd1291*WProd1390 + WProd1391 - WProd1392 + 
     -               WProd1306*WProd1393 + WProd1307*WProd1393 + 
     -               2*WProd1281*WProd1394 - 
     -               2*WProd1283*WProd1394 - WProd1284*WProd1394 - 
     -               2*WProd1285*WProd1394 - WProd1286*WProd1394 + 
     -               WProd1287*WProd1394 - 2*WProd1288*WProd1394 - 
     -               2*WProd1289*WProd1394 + WProd1290*WProd1394 + 
     -               WProd1291*WProd1394 - WProd1294*WProd1397 - 
     -               WProd1294*WProd1398))/WProd1381))
        amp(16,4,6) = Prp161*tmp1144 + 
     -   Prp157*Prp158*Prp159*
     -    ((-8*Cpl2*GS*(WProd1364 - WProd1365 - WProd1366))/
     -       WProd1367 + (8*Cpl2*GS*
     -         (WProd1278*WProd1368 + WProd1278*WProd1369 - 
     -           WProd1292*WProd1369 + 
     -           WProd1293*WProd1294*WProd1369 - 
     -           WProd1297*WProd1369 - 
     -           WProd1293*WProd1306*WProd1369 - 
     -           WProd1293*WProd1307*WProd1369 - 
     -           WProd1308*WProd1369 + WProd1309*WProd1369 + 
     -           WProd1310*WProd1369 + WProd1311*WProd1369 - 
     -           WProd1278*WProd1370 + 2*WProd1281*WProd1371 - 
     -           2*WProd1283*WProd1371 - 3*WProd1284*WProd1371 - 
     -           2*WProd1285*WProd1371 - 2*WProd1286*WProd1371 + 
     -           2*WProd1287*WProd1371 - 2*WProd1288*WProd1371 - 
     -           WProd1289*WProd1371 + 2*WProd1290*WProd1371 + 
     -           2*WProd1291*WProd1371 + WProd1294*WProd1372 - 
     -           WProd1306*WProd1372 - WProd1307*WProd1372 + 
     -           WProd1294*WProd1373 - WProd1306*WProd1373 - 
     -           WProd1307*WProd1373 - WProd1294*WProd1374 + 
     -           WProd1306*WProd1374 + WProd1307*WProd1374 + 
     -           WProd1375 + WProd1376 - WProd1377 + WProd1378 - 
     -           WProd1379 + WProd1380))/WProd1381) + 
     -   Prp157*Prp159*Prp160*
     -    ((8*Cpl2*GS*(WProd1365 + WProd1366))/WProd1367 - 
     -      (8*Cpl2*GS*(-(WProd1278*WProd1368) - 
     -           WProd1278*WProd1369 + WProd1297*WProd1369 - 
     -           WProd1311*WProd1369 + WProd1315*WProd1369 - 
     -           2*WProd1281*WProd1371 + 2*WProd1283*WProd1371 + 
     -           3*WProd1284*WProd1371 + 2*WProd1285*WProd1371 + 
     -           3*WProd1286*WProd1371 - 2*WProd1287*WProd1371 + 
     -           2*WProd1288*WProd1371 + 2*WProd1289*WProd1371 - 
     -           2*WProd1290*WProd1371 - 2*WProd1291*WProd1371 - 
     -           WProd1375 - WProd1376 - WProd1378 - WProd1380 + 
     -           WProd1278*WProd1382 - WProd1383 - 
     -           WProd1318*WProd1384 + WProd1320*WProd1384 + 
     -           WProd1321*WProd1384 + WProd1385))/WProd1381)

#ifdef DEBUG
	print *, 'amp(16,4,6) =', amp(16,4,6)
#endif

        tmp1145 = (-4*Cpl2*GS*Prp159*Prp160*
     -      (-2*WProd1337*WProd1368 + 2*WProd1341*WProd1368 - 
     -        2*WProd1355*WProd1368 + 2*WProd1356*WProd1368 - 
     -        2*WProd1337*WProd1369 + 2*WProd1341*WProd1369 - 
     -        2*WProd1355*WProd1369 + 2*WProd1356*WProd1369 + 
     -        2*WProd1337*WProd1382 + 2*WProd1281*WProd1404 - 
     -        4*WProd1283*WProd1404 - 2*WProd1284*WProd1404 - 
     -        4*WProd1285*WProd1404 - 2*WProd1286*WProd1404 + 
     -        4*WProd1287*WProd1404 - 3*WProd1288*WProd1404 - 
     -        3*WProd1289*WProd1404 + 3*WProd1290*WProd1404 + 
     -        3*WProd1291*WProd1404 - 2*WProd1405 - 2*WProd1406 - 
     -        2*WProd1408 - 2*WProd1410 - 2*WProd1318*WProd1412 + 
     -        2*WProd1320*WProd1412 + 2*WProd1321*WProd1412 + 
     -        2*WProd1413))/WProd1403 + 
     -   Prp158*Prp162*(2.6666666666666665D0*
     -       (Cpl2*GS*(WProd1414 + WProd1415))/WProd1399 - 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*(WProd1337*WProd1370 - WProd1339*WProd1370 - 
     -            WProd1293*WProd1350*WProd1370 - 
     -            WProd1293*WProd1351*WProd1370 + 
     -            WProd1353*WProd1370 + WProd1354*WProd1370 + 
     -            WProd1355*WProd1370 - WProd1340*WProd1373 - 
     -            WProd1340*WProd1388 + WProd1340*WProd1389 + 
     -            WProd1350*WProd1393 + WProd1351*WProd1393 + 
     -            WProd1413 - 2*WProd1281*WProd1416 + 
     -            2*WProd1283*WProd1416 + 2*WProd1284*WProd1416 + 
     -            2*WProd1285*WProd1416 + 2*WProd1286*WProd1416 - 
     -            2*WProd1287*WProd1416 + 2*WProd1288*WProd1416 + 
     -            3*WProd1289*WProd1416 - 2*WProd1290*WProd1416 - 
     -            2*WProd1291*WProd1416 + WProd1417 - WProd1418 + 
     -            2*WProd1281*WProd1419 - 2*WProd1283*WProd1419 - 
     -            2*WProd1284*WProd1419 - 2*WProd1285*WProd1419 - 
     -            2*WProd1286*WProd1419 + 2*WProd1287*WProd1419 - 
     -            2*WProd1288*WProd1419 - 3*WProd1289*WProd1419 + 
     -            2*WProd1290*WProd1419 + 2*WProd1291*WProd1419))/
     -        WProd1403)
        tmp1145 = tmp1145 + 
     -   Prp158*((4*Cpl2*GS*Prp159*
     -         (2*WProd1337*WProd1368 - 2*WProd1339*WProd1368 + 
     -           2*WProd1293*WProd1340*WProd1368 - 
     -           2*WProd1341*WProd1368 - 
     -           2*WProd1293*WProd1350*WProd1368 - 
     -           2*WProd1293*WProd1351*WProd1368 - 
     -           2*WProd1352*WProd1368 + 2*WProd1353*WProd1368 + 
     -           2*WProd1354*WProd1368 + 2*WProd1355*WProd1368 + 
     -           2*WProd1337*WProd1369 - 2*WProd1339*WProd1369 + 
     -           2*WProd1293*WProd1340*WProd1369 - 
     -           2*WProd1341*WProd1369 - 
     -           2*WProd1293*WProd1350*WProd1369 - 
     -           2*WProd1293*WProd1351*WProd1369 - 
     -           2*WProd1352*WProd1369 + 2*WProd1353*WProd1369 + 
     -           2*WProd1354*WProd1369 + 2*WProd1355*WProd1369 - 
     -           2*WProd1337*WProd1370 + 2*WProd1340*WProd1373 - 
     -           2*WProd1350*WProd1373 - 2*WProd1351*WProd1373 + 
     -           2*WProd1350*WProd1374 + 2*WProd1351*WProd1374 - 
     -           2*WProd1340*WProd1397 - 2*WProd1340*WProd1398 - 
     -           4*WProd1281*WProd1404 + 4*WProd1283*WProd1404 + 
     -           3*WProd1284*WProd1404 + 4*WProd1285*WProd1404 + 
     -           3*WProd1286*WProd1404 - 3*WProd1287*WProd1404 + 
     -           4*WProd1288*WProd1404 + 5*WProd1289*WProd1404 - 
     -           3*WProd1290*WProd1404 - 3*WProd1291*WProd1404 + 
     -           2*WProd1405 + 2*WProd1406 + 2*WProd1408 + 
     -           2*WProd1410 - 2*WProd1413))/WProd1403 + 
     -      Prp165*(-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1414 + WProd1415))/WProd1399) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1337*WProd1368 - 
     -               WProd1339*WProd1368 + 
     -               WProd1293*WProd1340*WProd1368 - 
     -               WProd1341*WProd1368 - 
     -               WProd1293*WProd1350*WProd1368 - 
     -               WProd1293*WProd1351*WProd1368 - 
     -               WProd1352*WProd1368 + WProd1353*WProd1368 + 
     -               WProd1354*WProd1368 + WProd1355*WProd1368 + 
     -               WProd1337*WProd1369 - WProd1339*WProd1369 + 
     -               WProd1293*WProd1340*WProd1369 - 
     -               WProd1341*WProd1369 - 
     -               WProd1293*WProd1350*WProd1369 - 
     -               WProd1293*WProd1351*WProd1369 - 
     -               WProd1352*WProd1369 + WProd1353*WProd1369 + 
     -               WProd1354*WProd1369 + WProd1355*WProd1369 - 
     -               WProd1339*WProd1370 - 
     -               WProd1293*WProd1350*WProd1370 - 
     -               WProd1293*WProd1351*WProd1370 + 
     -               WProd1353*WProd1370 + WProd1354*WProd1370 + 
     -               WProd1355*WProd1370 - WProd1350*WProd1373 - 
     -               WProd1351*WProd1373 + WProd1350*WProd1374 + 
     -               WProd1351*WProd1374 - WProd1340*WProd1388 + 
     -               WProd1340*WProd1389 + WProd1350*WProd1393 + 
     -               WProd1351*WProd1393 - WProd1340*WProd1397 - 
     -               WProd1340*WProd1398 - 2*WProd1281*WProd1404 + 
     -               2*WProd1283*WProd1404 + WProd1284*WProd1404 + 
     -               2*WProd1285*WProd1404 + WProd1286*WProd1404 - 
     -               WProd1287*WProd1404 + 2*WProd1288*WProd1404 + 
     -               2*WProd1289*WProd1404 - WProd1290*WProd1404 - 
     -               WProd1291*WProd1404 + WProd1405 + WProd1406 + 
     -               WProd1408 + WProd1410 - 
     -               2*WProd1281*WProd1416 + 
     -               2*WProd1283*WProd1416 + WProd1284*WProd1416 + 
     -               2*WProd1285*WProd1416 + WProd1286*WProd1416 - 
     -               WProd1287*WProd1416 + 2*WProd1288*WProd1416 + 
     -               2*WProd1289*WProd1416 - WProd1290*WProd1416 - 
     -               WProd1291*WProd1416 + WProd1417 - WProd1418 + 
     -               2*WProd1281*WProd1419 - 
     -               2*WProd1283*WProd1419 - WProd1284*WProd1419 - 
     -               2*WProd1285*WProd1419 - WProd1286*WProd1419 + 
     -               WProd1287*WProd1419 - 2*WProd1288*WProd1419 - 
     -               2*WProd1289*WProd1419 + WProd1290*WProd1419 + 
     -               WProd1291*WProd1419))/WProd1403)) + 
     -   Prp163*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp160*
     -          (WProd1355*WProd1368 - WProd1356*WProd1368 + 
     -            WProd1355*WProd1369 - WProd1356*WProd1369 - 
     -            WProd1355*WProd1382 + WProd1356*WProd1382 + 
     -            WProd1283*WProd1293*WProd1412 + 
     -            WProd1285*WProd1293*WProd1412 - 
     -            WProd1287*WProd1293*WProd1412 + 
     -            WProd1328*WProd1412 - WProd1331*WProd1412 - 
     -            WProd1283*WProd1416 - WProd1285*WProd1416 + 
     -            WProd1287*WProd1416 - WProd1417 + WProd1420))/
     -        WProd1403 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp162*
     -          (-(WProd1355*WProd1370) + WProd1356*WProd1370 + 
     -            WProd1340*WProd1388 - WProd1340*WProd1396 + 
     -            WProd1281*WProd1416 - WProd1283*WProd1416 - 
     -            WProd1284*WProd1416 - WProd1285*WProd1416 - 
     -            WProd1286*WProd1416 + WProd1287*WProd1416 - 
     -            WProd1288*WProd1416 - WProd1289*WProd1416 + 
     -            WProd1290*WProd1416 + WProd1291*WProd1416 - 
     -            WProd1417 + WProd1420))/WProd1403)
        tmp1145 = tmp1145 + 
     -   Prp164*(Prp165*(-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1414 + WProd1415))/WProd1399) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1337*WProd1368 - 
     -               WProd1341*WProd1368 - WProd1352*WProd1368 + 
     -               WProd1337*WProd1369 - WProd1341*WProd1369 - 
     -               WProd1352*WProd1369 - WProd1350*WProd1373 - 
     -               WProd1351*WProd1373 + WProd1350*WProd1374 + 
     -               WProd1351*WProd1374 - WProd1339*WProd1382 - 
     -               WProd1293*WProd1350*WProd1382 - 
     -               WProd1293*WProd1351*WProd1382 + 
     -               WProd1353*WProd1382 + WProd1354*WProd1382 + 
     -               WProd1355*WProd1382 + WProd1350*WProd1393 + 
     -               WProd1351*WProd1393 - 2*WProd1281*WProd1404 + 
     -               2*WProd1283*WProd1404 + WProd1284*WProd1404 + 
     -               2*WProd1285*WProd1404 + WProd1286*WProd1404 - 
     -               2*WProd1287*WProd1404 + 
     -               2*WProd1288*WProd1404 + 
     -               2*WProd1289*WProd1404 - 
     -               2*WProd1290*WProd1404 - 
     -               2*WProd1291*WProd1404 + WProd1405 + 
     -               WProd1406 + WProd1408 + WProd1410 + 
     -               2*WProd1281*WProd1293*WProd1412 - 
     -               2*WProd1283*WProd1293*WProd1412 - 
     -               WProd1284*WProd1293*WProd1412 - 
     -               2*WProd1285*WProd1293*WProd1412 - 
     -               WProd1286*WProd1293*WProd1412 + 
     -               2*WProd1287*WProd1293*WProd1412 - 
     -               2*WProd1288*WProd1293*WProd1412 - 
     -               2*WProd1289*WProd1293*WProd1412 + 
     -               2*WProd1290*WProd1293*WProd1412 + 
     -               2*WProd1291*WProd1293*WProd1412 - 
     -               WProd1320*WProd1412 - WProd1321*WProd1412 - 
     -               WProd1328*WProd1412 + WProd1329*WProd1412 - 
     -               2*WProd1281*WProd1416 + 
     -               2*WProd1283*WProd1416 + WProd1284*WProd1416 + 
     -               2*WProd1285*WProd1416 + WProd1286*WProd1416 - 
     -               2*WProd1287*WProd1416 + 
     -               2*WProd1288*WProd1416 + 
     -               2*WProd1289*WProd1416 - 
     -               2*WProd1290*WProd1416 - 
     -               2*WProd1291*WProd1416 + WProd1417 - 
     -               WProd1418 + 2*WProd1281*WProd1419 - 
     -               2*WProd1283*WProd1419 - WProd1284*WProd1419 - 
     -               2*WProd1285*WProd1419 - WProd1286*WProd1419 + 
     -               2*WProd1287*WProd1419 - 
     -               2*WProd1288*WProd1419 - 
     -               2*WProd1289*WProd1419 + 
     -               2*WProd1290*WProd1419 + 2*WProd1291*WProd1419)
     -             )/WProd1403) - 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp160*
     -          (-(WProd1337*WProd1368) + WProd1341*WProd1368 - 
     -            WProd1337*WProd1369 + WProd1341*WProd1369 + 
     -            WProd1337*WProd1382 - WProd1355*WProd1382 + 
     -            WProd1356*WProd1382 + 2*WProd1281*WProd1404 - 
     -            3*WProd1283*WProd1404 - 2*WProd1284*WProd1404 - 
     -            3*WProd1285*WProd1404 - 2*WProd1286*WProd1404 + 
     -            3*WProd1287*WProd1404 - 3*WProd1288*WProd1404 - 
     -            3*WProd1289*WProd1404 + 3*WProd1290*WProd1404 + 
     -            3*WProd1291*WProd1404 - WProd1405 - WProd1406 - 
     -            WProd1408 - WProd1410 - 
     -            2*WProd1281*WProd1293*WProd1412 + 
     -            3*WProd1283*WProd1293*WProd1412 + 
     -            2*WProd1284*WProd1293*WProd1412 + 
     -            3*WProd1285*WProd1293*WProd1412 + 
     -            2*WProd1286*WProd1293*WProd1412 - 
     -            3*WProd1287*WProd1293*WProd1412 + 
     -            3*WProd1288*WProd1293*WProd1412 + 
     -            3*WProd1289*WProd1293*WProd1412 - 
     -            3*WProd1290*WProd1293*WProd1412 - 
     -            3*WProd1291*WProd1293*WProd1412 - 
     -            WProd1318*WProd1412 + WProd1320*WProd1412 + 
     -            WProd1321*WProd1412 + WProd1328*WProd1412 - 
     -            WProd1331*WProd1412 + WProd1413 + 
     -            2*WProd1281*WProd1416 - 3*WProd1283*WProd1416 - 
     -            2*WProd1284*WProd1416 - 3*WProd1285*WProd1416 - 
     -            2*WProd1286*WProd1416 + 3*WProd1287*WProd1416 - 
     -            3*WProd1288*WProd1416 - 3*WProd1289*WProd1416 + 
     -            3*WProd1290*WProd1416 + 3*WProd1291*WProd1416 - 
     -            WProd1417 + WProd1420))/WProd1403)
        amp(16,4,4) = Prp161*tmp1145 + 
     -   Prp157*Prp158*Prp159*
     -    ((-8*Cpl2*GS*(WProd1400 - WProd1401 - WProd1402))/
     -       WProd1399 + (8*Cpl2*GS*
     -         (WProd1337*WProd1368 + WProd1337*WProd1369 - 
     -           WProd1339*WProd1369 + 
     -           WProd1293*WProd1340*WProd1369 - 
     -           WProd1341*WProd1369 - 
     -           WProd1293*WProd1350*WProd1369 - 
     -           WProd1293*WProd1351*WProd1369 - 
     -           WProd1352*WProd1369 + WProd1353*WProd1369 + 
     -           WProd1354*WProd1369 + WProd1355*WProd1369 - 
     -           WProd1337*WProd1370 + WProd1340*WProd1372 - 
     -           WProd1350*WProd1372 - WProd1351*WProd1372 + 
     -           WProd1340*WProd1373 - WProd1350*WProd1373 - 
     -           WProd1351*WProd1373 - WProd1340*WProd1374 + 
     -           WProd1350*WProd1374 + WProd1351*WProd1374 + 
     -           2*WProd1281*WProd1404 - 2*WProd1283*WProd1404 - 
     -           3*WProd1284*WProd1404 - 2*WProd1285*WProd1404 - 
     -           2*WProd1286*WProd1404 + 2*WProd1287*WProd1404 - 
     -           2*WProd1288*WProd1404 - WProd1289*WProd1404 + 
     -           2*WProd1290*WProd1404 + 2*WProd1291*WProd1404 + 
     -           WProd1405 + WProd1406 - WProd1407 + WProd1408 - 
     -           WProd1409 + WProd1410))/WProd1403) + 
     -   Prp157*Prp159*Prp160*
     -    ((8*Cpl2*GS*(WProd1401 + WProd1402))/WProd1399 - 
     -      (8*Cpl2*GS*(-(WProd1337*WProd1368) - 
     -           WProd1337*WProd1369 + WProd1341*WProd1369 - 
     -           WProd1355*WProd1369 + WProd1356*WProd1369 + 
     -           WProd1337*WProd1382 - 2*WProd1281*WProd1404 + 
     -           2*WProd1283*WProd1404 + 3*WProd1284*WProd1404 + 
     -           2*WProd1285*WProd1404 + 3*WProd1286*WProd1404 - 
     -           2*WProd1287*WProd1404 + 2*WProd1288*WProd1404 + 
     -           2*WProd1289*WProd1404 - 2*WProd1290*WProd1404 - 
     -           2*WProd1291*WProd1404 - WProd1405 - WProd1406 - 
     -           WProd1408 - WProd1410 - WProd1411 - 
     -           WProd1318*WProd1412 + WProd1320*WProd1412 + 
     -           WProd1321*WProd1412 + WProd1413))/WProd1403)

#ifdef DEBUG
	print *, 'amp(16,4,4) =', amp(16,4,4)
#endif

      endif

      end


