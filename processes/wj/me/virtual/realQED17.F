      subroutine realQED17(amp,resnr)
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

        amp(17,4,1) = (8*Cpl2*GS*Prp166*Prp167*Prp168*
     -     (WProd1421*WProd1422 + WProd1421*WProd1423 + 
     -       WProd1421*WProd1424 - WProd1421*WProd1425 - 
     -       WProd1426*WProd1427 - WProd1427*WProd1428 + 
     -       WProd1427*WProd1429 + WProd1421*WProd1430 + 
     -       WProd1431 + WProd1432 - WProd1433 - 
     -       WProd1421*WProd1434*WProd1435 - 
     -       WProd1421*WProd1434*WProd1436 + WProd1437*WProd1438 - 
     -       WProd1438*WProd1439 - WProd1438*WProd1440 - 
     -       WProd1438*WProd1441 - WProd1438*WProd1442 + 
     -       WProd1438*WProd1443 - WProd1438*WProd1444 - 
     -       2*WProd1438*WProd1445 + WProd1438*WProd1446 + 
     -       WProd1438*WProd1447 + WProd1421*WProd1448))/WProd1449
        amp(17,4,1) = amp(17,4,1) + 
     -   Prp169*((8*Cpl2*GS*Prp166*Prp168*
     -         (WProd1431 + WProd1432 - WProd1433 - 
     -           WProd1428*WProd1435 - WProd1428*WProd1436 - 
     -           WProd1438*WProd1439 - WProd1438*WProd1441 - 
     -           WProd1438*WProd1442 + WProd1438*WProd1443 - 
     -           WProd1438*WProd1444 - 2*WProd1438*WProd1445 + 
     -           WProd1438*WProd1446 + WProd1438*WProd1447 + 
     -           WProd1422*WProd1451 + WProd1423*WProd1451 + 
     -           WProd1424*WProd1451 - WProd1425*WProd1451 + 
     -           WProd1430*WProd1451 - 
     -           WProd1434*WProd1435*WProd1451 - 
     -           WProd1434*WProd1436*WProd1451 + 
     -           WProd1448*WProd1451 + WProd1452 - WProd1453 - 
     -           WProd1454 + WProd1455 + WProd1456 + 
     -           WProd1435*WProd1457 + WProd1436*WProd1457 + 
     -           WProd1421*WProd1458 - WProd1451*WProd1458 - 
     -           WProd1458*WProd1459))/WProd1449 + 
     -      Prp170*((4*Cpl2*GS*Prp166*
     -            (2*WProd1431 + 2*WProd1432 - 2*WProd1433 - 
     -              2*WProd1428*WProd1435 - 
     -              2*WProd1428*WProd1436 + 
     -              4*WProd1437*WProd1438 - 
     -              5*WProd1438*WProd1439 - 
     -              3*WProd1438*WProd1440 - 
     -              5*WProd1438*WProd1441 - 
     -              3*WProd1438*WProd1442 + 
     -              5*WProd1438*WProd1443 - 
     -              5*WProd1438*WProd1444 - 
     -              5*WProd1438*WProd1445 + 
     -              5*WProd1438*WProd1446 + 
     -              5*WProd1438*WProd1447 + 
     -              2*WProd1423*WProd1451 + 
     -              2*WProd1424*WProd1451 - 
     -              2*WProd1425*WProd1451 - 
     -              2*WProd1434*WProd1435*WProd1451 - 
     -              2*WProd1434*WProd1436*WProd1451 + 
     -              2*WProd1448*WProd1451 + 
     -              2*WProd1435*WProd1457 + 
     -              2*WProd1436*WProd1457 - 
     -              2*WProd1458*WProd1459 - 
     -              2*WProd1460*WProd1461 - 
     -              2*WProd1461*WProd1462 + 
     -              2*WProd1461*WProd1463 - 
     -              2*WProd1458*WProd1464 + 2*WProd1451*WProd1465))
     -           /WProd1449 + 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp172*
     -             (WProd1431 + WProd1432 + WProd1423*WProd1451 + 
     -               WProd1424*WProd1451 - WProd1424*WProd1459 - 
     -               2*WProd1434*WProd1437*WProd1461 + 
     -               3*WProd1434*WProd1439*WProd1461 + 
     -               2*WProd1434*WProd1440*WProd1461 + 
     -               3*WProd1434*WProd1441*WProd1461 + 
     -               2*WProd1434*WProd1442*WProd1461 - 
     -               3*WProd1434*WProd1443*WProd1461 + 
     -               3*WProd1434*WProd1444*WProd1461 + 
     -               3*WProd1434*WProd1445*WProd1461 - 
     -               3*WProd1434*WProd1446*WProd1461 - 
     -               3*WProd1434*WProd1447*WProd1461 - 
     -               WProd1460*WProd1461 - WProd1461*WProd1462 - 
     -               WProd1424*WProd1464 + 2*WProd1437*WProd1466 - 
     -               3*WProd1439*WProd1466 - 
     -               2*WProd1440*WProd1466 - 
     -               3*WProd1441*WProd1466 - 
     -               2*WProd1442*WProd1466 + 
     -               3*WProd1443*WProd1466 - 
     -               3*WProd1444*WProd1466 - 
     -               3*WProd1445*WProd1466 + 
     -               3*WProd1446*WProd1466 + 
     -               3*WProd1447*WProd1466 + WProd1461*WProd1467 + 
     -               WProd1461*WProd1468 + WProd1461*WProd1469 + 
     -               2*WProd1437*WProd1470 - 
     -               3*WProd1439*WProd1470 - 
     -               2*WProd1440*WProd1470 - 
     -               3*WProd1441*WProd1470 - 
     -               2*WProd1442*WProd1470 + 
     -               3*WProd1443*WProd1470 - 
     -               3*WProd1444*WProd1470 - 
     -               3*WProd1445*WProd1470 + 
     -               3*WProd1446*WProd1470 + 
     -               3*WProd1447*WProd1470 + WProd1459*WProd1471 + 
     -               WProd1464*WProd1471 - WProd1472 - WProd1473 - 
     -               WProd1474 - WProd1475 - WProd1476 - 
     -               WProd1451*WProd1477))/WProd1449 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*Prp171*
     -             (-WProd1433 - WProd1428*WProd1435 - 
     -               WProd1428*WProd1436 + 2*WProd1437*WProd1438 - 
     -               2*WProd1438*WProd1439 - WProd1438*WProd1440 - 
     -               2*WProd1438*WProd1441 - WProd1438*WProd1442 + 
     -               2*WProd1438*WProd1443 - 
     -               2*WProd1438*WProd1444 - 
     -               2*WProd1438*WProd1445 + 
     -               2*WProd1438*WProd1446 + 
     -               2*WProd1438*WProd1447 - WProd1425*WProd1451 - 
     -               WProd1434*WProd1435*WProd1451 - 
     -               WProd1434*WProd1436*WProd1451 + 
     -               WProd1448*WProd1451 + WProd1435*WProd1457 + 
     -               WProd1436*WProd1457 + WProd1424*WProd1459 - 
     -               WProd1458*WProd1459 + 
     -               2*WProd1434*WProd1437*WProd1461 - 
     -               2*WProd1434*WProd1439*WProd1461 - 
     -               WProd1434*WProd1440*WProd1461 - 
     -               2*WProd1434*WProd1441*WProd1461 - 
     -               WProd1434*WProd1442*WProd1461 + 
     -               2*WProd1434*WProd1443*WProd1461 - 
     -               2*WProd1434*WProd1444*WProd1461 - 
     -               2*WProd1434*WProd1445*WProd1461 + 
     -               2*WProd1434*WProd1446*WProd1461 + 
     -               2*WProd1434*WProd1447*WProd1461 + 
     -               WProd1461*WProd1463 + WProd1424*WProd1464 - 
     -               WProd1458*WProd1464 + WProd1451*WProd1465 - 
     -               2*WProd1437*WProd1466 + 
     -               2*WProd1439*WProd1466 + WProd1440*WProd1466 + 
     -               2*WProd1441*WProd1466 + WProd1442*WProd1466 - 
     -               2*WProd1443*WProd1466 + 
     -               2*WProd1444*WProd1466 + 
     -               2*WProd1445*WProd1466 - 
     -               2*WProd1446*WProd1466 - 
     -               2*WProd1447*WProd1466 - WProd1461*WProd1467 - 
     -               WProd1461*WProd1468 - WProd1461*WProd1469 - 
     -               2*WProd1437*WProd1470 + 
     -               2*WProd1439*WProd1470 + WProd1440*WProd1470 + 
     -               2*WProd1441*WProd1470 + WProd1442*WProd1470 - 
     -               2*WProd1443*WProd1470 + 
     -               2*WProd1444*WProd1470 + 
     -               2*WProd1445*WProd1470 - 
     -               2*WProd1446*WProd1470 - 
     -               2*WProd1447*WProd1470 - WProd1459*WProd1471 - 
     -               WProd1464*WProd1471 + WProd1472 + WProd1473 + 
     -               WProd1474 + WProd1475 + WProd1476 + 
     -               WProd1451*WProd1477))/WProd1449))
        amp(17,4,1) = amp(17,4,1) + 
     -   Prp170*((4*Cpl2*GS*Prp166*Prp167*
     -         (2*WProd1421*WProd1423 + 2*WProd1421*WProd1424 - 
     -           2*WProd1421*WProd1425 - 2*WProd1426*WProd1427 - 
     -           2*WProd1427*WProd1428 + 2*WProd1427*WProd1429 + 
     -           2*WProd1431 + 2*WProd1432 - 2*WProd1433 - 
     -           2*WProd1421*WProd1434*WProd1435 - 
     -           2*WProd1421*WProd1434*WProd1436 + 
     -           3*WProd1437*WProd1438 - 3*WProd1438*WProd1439 - 
     -           3*WProd1438*WProd1440 - 3*WProd1438*WProd1441 - 
     -           3*WProd1438*WProd1442 + 3*WProd1438*WProd1443 - 
     -           3*WProd1438*WProd1444 - 4*WProd1438*WProd1445 + 
     -           3*WProd1438*WProd1446 + 3*WProd1438*WProd1447 + 
     -           2*WProd1421*WProd1448 + 2*WProd1421*WProd1465))/
     -       WProd1449 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp167*Prp173*
     -          (WProd1421*WProd1424 - WProd1421*WProd1425 - 
     -            WProd1427*WProd1428 + WProd1427*WProd1429 + 
     -            WProd1432 - WProd1433 - 
     -            WProd1421*WProd1434*WProd1435 - 
     -            WProd1421*WProd1434*WProd1436 + 
     -            2*WProd1437*WProd1438 - 2*WProd1438*WProd1439 - 
     -            2*WProd1438*WProd1440 - 2*WProd1438*WProd1441 - 
     -            2*WProd1438*WProd1442 + 2*WProd1438*WProd1443 - 
     -            2*WProd1438*WProd1444 - 3*WProd1438*WProd1445 + 
     -            2*WProd1438*WProd1446 + 2*WProd1438*WProd1447 + 
     -            WProd1421*WProd1448 + WProd1421*WProd1465 - 
     -            2*WProd1437*WProd1470 + 2*WProd1439*WProd1470 + 
     -            2*WProd1440*WProd1470 + 2*WProd1441*WProd1470 + 
     -            2*WProd1442*WProd1470 - 2*WProd1443*WProd1470 + 
     -            2*WProd1444*WProd1470 + 3*WProd1445*WProd1470 - 
     -            2*WProd1446*WProd1470 - 2*WProd1447*WProd1470 + 
     -            WProd1474 + WProd1421*WProd1477 - 
     -            WProd1427*WProd1481))/WProd1449 + 
     -      Prp174*(2.6666666666666665D0*
     -          (Cpl2*GS*Prp172*
     -             (WProd1431 + WProd1423*WProd1451 - 
     -               WProd1423*WProd1459 + 
     -               WProd1434*WProd1439*WProd1461 + 
     -               WProd1434*WProd1441*WProd1461 - 
     -               WProd1434*WProd1443*WProd1461 - 
     -               WProd1460*WProd1461 - WProd1423*WProd1464 + 
     -               WProd1461*WProd1468 - WProd1439*WProd1470 - 
     -               WProd1441*WProd1470 + WProd1443*WProd1470 - 
     -               WProd1474 - WProd1451*WProd1477 + 
     -               WProd1459*WProd1477 + WProd1464*WProd1477))/
     -           WProd1449 + 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp167*
     -             (WProd1421*WProd1423 - WProd1426*WProd1427 + 
     -               WProd1431 + WProd1437*WProd1470 - 
     -               WProd1439*WProd1470 - WProd1440*WProd1470 - 
     -               WProd1441*WProd1470 - WProd1442*WProd1470 + 
     -               WProd1443*WProd1470 - WProd1444*WProd1470 - 
     -               WProd1445*WProd1470 + WProd1446*WProd1470 + 
     -               WProd1447*WProd1470 - WProd1474 - 
     -               WProd1421*WProd1477 + WProd1427*WProd1481))/
     -           WProd1449) + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp171*Prp173*
     -          (-(WProd1421*WProd1425) + WProd1427*WProd1429 - 
     -            WProd1433 - WProd1428*WProd1435 - 
     -            WProd1421*WProd1434*WProd1435 - 
     -            WProd1428*WProd1436 - 
     -            WProd1421*WProd1434*WProd1436 + 
     -            2*WProd1437*WProd1438 - 2*WProd1438*WProd1439 - 
     -            WProd1438*WProd1440 - 2*WProd1438*WProd1441 - 
     -            WProd1438*WProd1442 + WProd1438*WProd1443 - 
     -            2*WProd1438*WProd1444 - 2*WProd1438*WProd1445 + 
     -            WProd1438*WProd1446 + WProd1438*WProd1447 + 
     -            WProd1421*WProd1448 + WProd1435*WProd1457 + 
     -            WProd1436*WProd1457 + WProd1424*WProd1459 - 
     -            WProd1425*WProd1459 + 
     -            WProd1427*WProd1434*WProd1459 - 
     -            WProd1434*WProd1435*WProd1459 - 
     -            WProd1434*WProd1436*WProd1459 + 
     -            WProd1448*WProd1459 - WProd1458*WProd1459 + 
     -            WProd1424*WProd1464 - WProd1425*WProd1464 + 
     -            WProd1427*WProd1434*WProd1464 - 
     -            WProd1434*WProd1435*WProd1464 - 
     -            WProd1434*WProd1436*WProd1464 + 
     -            WProd1448*WProd1464 - WProd1458*WProd1464 + 
     -            WProd1421*WProd1465 + WProd1459*WProd1465 + 
     -            WProd1464*WProd1465 - 2*WProd1437*WProd1466 + 
     -            2*WProd1439*WProd1466 + WProd1440*WProd1466 + 
     -            2*WProd1441*WProd1466 + WProd1442*WProd1466 - 
     -            WProd1443*WProd1466 + 2*WProd1444*WProd1466 + 
     -            2*WProd1445*WProd1466 - WProd1446*WProd1466 - 
     -            WProd1447*WProd1466 - 2*WProd1437*WProd1470 + 
     -            2*WProd1439*WProd1470 + WProd1440*WProd1470 + 
     -            2*WProd1441*WProd1470 + WProd1442*WProd1470 - 
     -            WProd1443*WProd1470 + 2*WProd1444*WProd1470 + 
     -            2*WProd1445*WProd1470 - WProd1446*WProd1470 - 
     -            WProd1447*WProd1470 - WProd1459*WProd1471 - 
     -            WProd1464*WProd1471 + WProd1472 + WProd1473 + 
     -            WProd1474 + WProd1475 + WProd1476 + 
     -            WProd1421*WProd1477 + WProd1459*WProd1477 + 
     -            WProd1464*WProd1477 - WProd1427*WProd1481 - 
     -            WProd1427*WProd1482 - WProd1427*WProd1483))/
     -        WProd1449)

#ifdef DEBUG
	print *, 'amp(17,4,1) =', amp(17,4,1)
#endif

        amp(17,4,8) = (8*Cpl2*GS*Prp166*Prp167*Prp168*
     -     (WProd1421*WProd1485 + WProd1421*WProd1486 + 
     -       WProd1421*WProd1487 - WProd1421*WProd1488 - 
     -       WProd1426*WProd1489 - WProd1428*WProd1489 + 
     -       WProd1429*WProd1489 + WProd1421*WProd1490 + 
     -       WProd1491 + WProd1492 - WProd1493 - 
     -       WProd1421*WProd1434*WProd1494 - 
     -       WProd1421*WProd1434*WProd1495 + WProd1437*WProd1496 - 
     -       WProd1439*WProd1496 - WProd1440*WProd1496 - 
     -       WProd1441*WProd1496 - WProd1442*WProd1496 + 
     -       WProd1443*WProd1496 - WProd1444*WProd1496 - 
     -       2*WProd1445*WProd1496 + WProd1446*WProd1496 + 
     -       WProd1447*WProd1496 + WProd1421*WProd1497))/WProd1484
        amp(17,4,8) = amp(17,4,8) + 
     -   Prp169*((8*Cpl2*GS*Prp166*Prp168*
     -         (WProd1451*WProd1485 + WProd1451*WProd1486 + 
     -           WProd1451*WProd1487 - WProd1451*WProd1488 + 
     -           WProd1451*WProd1490 + WProd1491 + WProd1492 - 
     -           WProd1493 - WProd1428*WProd1494 - 
     -           WProd1434*WProd1451*WProd1494 + 
     -           WProd1457*WProd1494 - WProd1428*WProd1495 - 
     -           WProd1434*WProd1451*WProd1495 + 
     -           WProd1457*WProd1495 - WProd1439*WProd1496 - 
     -           WProd1441*WProd1496 - WProd1442*WProd1496 + 
     -           WProd1443*WProd1496 - WProd1444*WProd1496 - 
     -           2*WProd1445*WProd1496 + WProd1446*WProd1496 + 
     -           WProd1447*WProd1496 + WProd1451*WProd1497 + 
     -           WProd1498 - WProd1499 - WProd1500 + WProd1501 + 
     -           WProd1502 + WProd1421*WProd1503 - 
     -           WProd1451*WProd1503 - WProd1459*WProd1503))/
     -       WProd1484 + Prp170*
     -       ((4*Cpl2*GS*Prp166*
     -            (2*WProd1451*WProd1486 + 2*WProd1451*WProd1487 - 
     -              2*WProd1451*WProd1488 + 2*WProd1491 + 
     -              2*WProd1492 - 2*WProd1493 - 
     -              2*WProd1428*WProd1494 - 
     -              2*WProd1434*WProd1451*WProd1494 + 
     -              2*WProd1457*WProd1494 - 
     -              2*WProd1428*WProd1495 - 
     -              2*WProd1434*WProd1451*WProd1495 + 
     -              2*WProd1457*WProd1495 + 
     -              4*WProd1437*WProd1496 - 
     -              5*WProd1439*WProd1496 - 
     -              3*WProd1440*WProd1496 - 
     -              5*WProd1441*WProd1496 - 
     -              3*WProd1442*WProd1496 + 
     -              5*WProd1443*WProd1496 - 
     -              5*WProd1444*WProd1496 - 
     -              5*WProd1445*WProd1496 + 
     -              5*WProd1446*WProd1496 + 
     -              5*WProd1447*WProd1496 + 
     -              2*WProd1451*WProd1497 - 
     -              2*WProd1459*WProd1503 - 
     -              2*WProd1464*WProd1503 - 
     -              2*WProd1460*WProd1504 - 
     -              2*WProd1462*WProd1504 + 
     -              2*WProd1463*WProd1504 + 2*WProd1451*WProd1505))
     -           /WProd1484 + 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp172*
     -             (WProd1451*WProd1486 + WProd1451*WProd1487 - 
     -               WProd1459*WProd1487 - WProd1464*WProd1487 + 
     -               WProd1491 + WProd1492 - 
     -               2*WProd1434*WProd1437*WProd1504 + 
     -               3*WProd1434*WProd1439*WProd1504 + 
     -               2*WProd1434*WProd1440*WProd1504 + 
     -               3*WProd1434*WProd1441*WProd1504 + 
     -               2*WProd1434*WProd1442*WProd1504 - 
     -               3*WProd1434*WProd1443*WProd1504 + 
     -               3*WProd1434*WProd1444*WProd1504 + 
     -               3*WProd1434*WProd1445*WProd1504 - 
     -               3*WProd1434*WProd1446*WProd1504 - 
     -               3*WProd1434*WProd1447*WProd1504 - 
     -               WProd1460*WProd1504 - WProd1462*WProd1504 + 
     -               WProd1467*WProd1504 + WProd1468*WProd1504 + 
     -               WProd1469*WProd1504 + 2*WProd1437*WProd1506 - 
     -               3*WProd1439*WProd1506 - 
     -               2*WProd1440*WProd1506 - 
     -               3*WProd1441*WProd1506 - 
     -               2*WProd1442*WProd1506 + 
     -               3*WProd1443*WProd1506 - 
     -               3*WProd1444*WProd1506 - 
     -               3*WProd1445*WProd1506 + 
     -               3*WProd1446*WProd1506 + 
     -               3*WProd1447*WProd1506 + 
     -               2*WProd1437*WProd1507 - 
     -               3*WProd1439*WProd1507 - 
     -               2*WProd1440*WProd1507 - 
     -               3*WProd1441*WProd1507 - 
     -               2*WProd1442*WProd1507 + 
     -               3*WProd1443*WProd1507 - 
     -               3*WProd1444*WProd1507 - 
     -               3*WProd1445*WProd1507 + 
     -               3*WProd1446*WProd1507 + 
     -               3*WProd1447*WProd1507 + WProd1459*WProd1508 + 
     -               WProd1464*WProd1508 - WProd1509 - WProd1510 - 
     -               WProd1511 - WProd1512 - WProd1513 - 
     -               WProd1451*WProd1514))/WProd1484 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*Prp171*
     -             (WProd1459*WProd1487 + WProd1464*WProd1487 - 
     -               WProd1451*WProd1488 - WProd1493 - 
     -               WProd1428*WProd1494 - 
     -               WProd1434*WProd1451*WProd1494 + 
     -               WProd1457*WProd1494 - WProd1428*WProd1495 - 
     -               WProd1434*WProd1451*WProd1495 + 
     -               WProd1457*WProd1495 + 2*WProd1437*WProd1496 - 
     -               2*WProd1439*WProd1496 - WProd1440*WProd1496 - 
     -               2*WProd1441*WProd1496 - WProd1442*WProd1496 + 
     -               2*WProd1443*WProd1496 - 
     -               2*WProd1444*WProd1496 - 
     -               2*WProd1445*WProd1496 + 
     -               2*WProd1446*WProd1496 + 
     -               2*WProd1447*WProd1496 + WProd1451*WProd1497 - 
     -               WProd1459*WProd1503 - WProd1464*WProd1503 + 
     -               2*WProd1434*WProd1437*WProd1504 - 
     -               2*WProd1434*WProd1439*WProd1504 - 
     -               WProd1434*WProd1440*WProd1504 - 
     -               2*WProd1434*WProd1441*WProd1504 - 
     -               WProd1434*WProd1442*WProd1504 + 
     -               2*WProd1434*WProd1443*WProd1504 - 
     -               2*WProd1434*WProd1444*WProd1504 - 
     -               2*WProd1434*WProd1445*WProd1504 + 
     -               2*WProd1434*WProd1446*WProd1504 + 
     -               2*WProd1434*WProd1447*WProd1504 + 
     -               WProd1463*WProd1504 - WProd1467*WProd1504 - 
     -               WProd1468*WProd1504 - WProd1469*WProd1504 + 
     -               WProd1451*WProd1505 - 2*WProd1437*WProd1506 + 
     -               2*WProd1439*WProd1506 + WProd1440*WProd1506 + 
     -               2*WProd1441*WProd1506 + WProd1442*WProd1506 - 
     -               2*WProd1443*WProd1506 + 
     -               2*WProd1444*WProd1506 + 
     -               2*WProd1445*WProd1506 - 
     -               2*WProd1446*WProd1506 - 
     -               2*WProd1447*WProd1506 - 
     -               2*WProd1437*WProd1507 + 
     -               2*WProd1439*WProd1507 + WProd1440*WProd1507 + 
     -               2*WProd1441*WProd1507 + WProd1442*WProd1507 - 
     -               2*WProd1443*WProd1507 + 
     -               2*WProd1444*WProd1507 + 
     -               2*WProd1445*WProd1507 - 
     -               2*WProd1446*WProd1507 - 
     -               2*WProd1447*WProd1507 - WProd1459*WProd1508 - 
     -               WProd1464*WProd1508 + WProd1509 + WProd1510 + 
     -               WProd1511 + WProd1512 + WProd1513 + 
     -               WProd1451*WProd1514))/WProd1484))
        amp(17,4,8) = amp(17,4,8) + 
     -   Prp170*((4*Cpl2*GS*Prp166*Prp167*
     -         (2*WProd1421*WProd1486 + 2*WProd1421*WProd1487 - 
     -           2*WProd1421*WProd1488 - 2*WProd1426*WProd1489 - 
     -           2*WProd1428*WProd1489 + 2*WProd1429*WProd1489 + 
     -           2*WProd1491 + 2*WProd1492 - 2*WProd1493 - 
     -           2*WProd1421*WProd1434*WProd1494 - 
     -           2*WProd1421*WProd1434*WProd1495 + 
     -           3*WProd1437*WProd1496 - 3*WProd1439*WProd1496 - 
     -           3*WProd1440*WProd1496 - 3*WProd1441*WProd1496 - 
     -           3*WProd1442*WProd1496 + 3*WProd1443*WProd1496 - 
     -           3*WProd1444*WProd1496 - 4*WProd1445*WProd1496 + 
     -           3*WProd1446*WProd1496 + 3*WProd1447*WProd1496 + 
     -           2*WProd1421*WProd1497 + 2*WProd1421*WProd1505))/
     -       WProd1484 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp167*Prp173*
     -          (WProd1421*WProd1487 - WProd1421*WProd1488 - 
     -            WProd1428*WProd1489 + WProd1429*WProd1489 - 
     -            WProd1481*WProd1489 + WProd1492 - WProd1493 - 
     -            WProd1421*WProd1434*WProd1494 - 
     -            WProd1421*WProd1434*WProd1495 + 
     -            2*WProd1437*WProd1496 - 2*WProd1439*WProd1496 - 
     -            2*WProd1440*WProd1496 - 2*WProd1441*WProd1496 - 
     -            2*WProd1442*WProd1496 + 2*WProd1443*WProd1496 - 
     -            2*WProd1444*WProd1496 - 3*WProd1445*WProd1496 + 
     -            2*WProd1446*WProd1496 + 2*WProd1447*WProd1496 + 
     -            WProd1421*WProd1497 + WProd1421*WProd1505 - 
     -            2*WProd1437*WProd1507 + 2*WProd1439*WProd1507 + 
     -            2*WProd1440*WProd1507 + 2*WProd1441*WProd1507 + 
     -            2*WProd1442*WProd1507 - 2*WProd1443*WProd1507 + 
     -            2*WProd1444*WProd1507 + 3*WProd1445*WProd1507 - 
     -            2*WProd1446*WProd1507 - 2*WProd1447*WProd1507 + 
     -            WProd1511 + WProd1421*WProd1514))/WProd1484 + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp171*Prp173*
     -          (WProd1459*WProd1487 + WProd1464*WProd1487 - 
     -            WProd1421*WProd1488 - WProd1459*WProd1488 - 
     -            WProd1464*WProd1488 + WProd1429*WProd1489 + 
     -            WProd1434*WProd1459*WProd1489 + 
     -            WProd1434*WProd1464*WProd1489 - 
     -            WProd1481*WProd1489 - WProd1482*WProd1489 - 
     -            WProd1483*WProd1489 - WProd1493 - 
     -            WProd1428*WProd1494 - 
     -            WProd1421*WProd1434*WProd1494 + 
     -            WProd1457*WProd1494 - 
     -            WProd1434*WProd1459*WProd1494 - 
     -            WProd1434*WProd1464*WProd1494 - 
     -            WProd1428*WProd1495 - 
     -            WProd1421*WProd1434*WProd1495 + 
     -            WProd1457*WProd1495 - 
     -            WProd1434*WProd1459*WProd1495 - 
     -            WProd1434*WProd1464*WProd1495 + 
     -            2*WProd1437*WProd1496 - 2*WProd1439*WProd1496 - 
     -            WProd1440*WProd1496 - 2*WProd1441*WProd1496 - 
     -            WProd1442*WProd1496 + WProd1443*WProd1496 - 
     -            2*WProd1444*WProd1496 - 2*WProd1445*WProd1496 + 
     -            WProd1446*WProd1496 + WProd1447*WProd1496 + 
     -            WProd1421*WProd1497 + WProd1459*WProd1497 + 
     -            WProd1464*WProd1497 - WProd1459*WProd1503 - 
     -            WProd1464*WProd1503 + WProd1421*WProd1505 + 
     -            WProd1459*WProd1505 + WProd1464*WProd1505 - 
     -            2*WProd1437*WProd1506 + 2*WProd1439*WProd1506 + 
     -            WProd1440*WProd1506 + 2*WProd1441*WProd1506 + 
     -            WProd1442*WProd1506 - WProd1443*WProd1506 + 
     -            2*WProd1444*WProd1506 + 2*WProd1445*WProd1506 - 
     -            WProd1446*WProd1506 - WProd1447*WProd1506 - 
     -            2*WProd1437*WProd1507 + 2*WProd1439*WProd1507 + 
     -            WProd1440*WProd1507 + 2*WProd1441*WProd1507 + 
     -            WProd1442*WProd1507 - WProd1443*WProd1507 + 
     -            2*WProd1444*WProd1507 + 2*WProd1445*WProd1507 - 
     -            WProd1446*WProd1507 - WProd1447*WProd1507 - 
     -            WProd1459*WProd1508 - WProd1464*WProd1508 + 
     -            WProd1509 + WProd1510 + WProd1511 + WProd1512 + 
     -            WProd1513 + WProd1421*WProd1514 + 
     -            WProd1459*WProd1514 + WProd1464*WProd1514))/
     -        WProd1484 + Prp174*
     -       (2.6666666666666665D0*
     -          (Cpl2*GS*Prp167*
     -             (WProd1421*WProd1486 - WProd1426*WProd1489 + 
     -               WProd1481*WProd1489 + WProd1491 + 
     -               WProd1437*WProd1507 - WProd1439*WProd1507 - 
     -               WProd1440*WProd1507 - WProd1441*WProd1507 - 
     -               WProd1442*WProd1507 + WProd1443*WProd1507 - 
     -               WProd1444*WProd1507 - WProd1445*WProd1507 + 
     -               WProd1446*WProd1507 + WProd1447*WProd1507 - 
     -               WProd1511 - WProd1421*WProd1514))/WProd1484 + 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp172*
     -             (WProd1451*WProd1486 - WProd1459*WProd1486 - 
     -               WProd1464*WProd1486 + WProd1491 + 
     -               WProd1434*WProd1439*WProd1504 + 
     -               WProd1434*WProd1441*WProd1504 - 
     -               WProd1434*WProd1443*WProd1504 - 
     -               WProd1460*WProd1504 + WProd1468*WProd1504 - 
     -               WProd1439*WProd1507 - WProd1441*WProd1507 + 
     -               WProd1443*WProd1507 - WProd1511 - 
     -               WProd1451*WProd1514 + WProd1459*WProd1514 + 
     -               WProd1464*WProd1514))/WProd1484))

#ifdef DEBUG
	print *, 'amp(17,4,8) =', amp(17,4,8)
#endif

        tmp1146 = Prp166*Prp168*
     -   ((-8*Cpl2*GS*WProd1515)/WProd1516 + 
     -     (8*Cpl2*GS*(WProd1458*WProd1517 - WProd1518 - 
     -          WProd1435*WProd1520 - WProd1436*WProd1520 + 
     -          WProd1522 - WProd1524 + WProd1435*WProd1525 + 
     -          WProd1436*WProd1525 - WProd1439*WProd1526 - 
     -          WProd1441*WProd1526 - WProd1442*WProd1526 + 
     -          WProd1443*WProd1526 - WProd1444*WProd1526 - 
     -          2*WProd1445*WProd1526 + WProd1446*WProd1526 + 
     -          WProd1447*WProd1526 + WProd1527 + 
     -          WProd1422*WProd1529 + WProd1423*WProd1529 + 
     -          WProd1424*WProd1529 - WProd1425*WProd1529 + 
     -          WProd1430*WProd1529 - 
     -          WProd1434*WProd1435*WProd1529 - 
     -          WProd1434*WProd1436*WProd1529 + 
     -          WProd1448*WProd1529 - WProd1458*WProd1529 + 
     -          WProd1530 + WProd1531 - WProd1532 + WProd1533 + 
     -          WProd1435*WProd1534 + WProd1436*WProd1534 - 
     -          WProd1458*WProd1535))/WProd1528)
        tmp1146 = tmp1146 + 
     -   Prp170*(Prp166*((-8*Cpl2*GS*(WProd1515 + WProd1536))/
     -          WProd1516 + 
     -         (4*Cpl2*GS*(-2*WProd1435*WProd1520 - 
     -              2*WProd1436*WProd1520 + 2*WProd1522 + 
     -              2*WProd1523 - 2*WProd1524 + 
     -              2*WProd1435*WProd1525 + 
     -              2*WProd1436*WProd1525 + 
     -              4*WProd1437*WProd1526 - 
     -              5*WProd1439*WProd1526 - 
     -              3*WProd1440*WProd1526 - 
     -              5*WProd1441*WProd1526 - 
     -              3*WProd1442*WProd1526 + 
     -              5*WProd1443*WProd1526 - 
     -              5*WProd1444*WProd1526 - 
     -              5*WProd1445*WProd1526 + 
     -              5*WProd1446*WProd1526 + 
     -              5*WProd1447*WProd1526 + 
     -              2*WProd1423*WProd1529 + 
     -              2*WProd1424*WProd1529 - 
     -              2*WProd1425*WProd1529 - 
     -              2*WProd1434*WProd1435*WProd1529 - 
     -              2*WProd1434*WProd1436*WProd1529 + 
     -              2*WProd1448*WProd1529 + 
     -              2*WProd1465*WProd1529 + 
     -              2*WProd1435*WProd1534 + 
     -              2*WProd1436*WProd1534 - 
     -              2*WProd1458*WProd1535 - 
     -              2*WProd1460*WProd1537 - 
     -              2*WProd1462*WProd1537 + 
     -              2*WProd1463*WProd1537 - 2*WProd1458*WProd1538))
     -           /WProd1528) + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp172*
     -          (WProd1522 + WProd1523 + WProd1423*WProd1529 + 
     -            WProd1424*WProd1529 - WProd1477*WProd1529 - 
     -            WProd1424*WProd1535 + WProd1471*WProd1535 - 
     -            2*WProd1434*WProd1437*WProd1537 + 
     -            3*WProd1434*WProd1439*WProd1537 + 
     -            2*WProd1434*WProd1440*WProd1537 + 
     -            3*WProd1434*WProd1441*WProd1537 + 
     -            2*WProd1434*WProd1442*WProd1537 - 
     -            3*WProd1434*WProd1443*WProd1537 + 
     -            3*WProd1434*WProd1444*WProd1537 + 
     -            3*WProd1434*WProd1445*WProd1537 - 
     -            3*WProd1434*WProd1446*WProd1537 - 
     -            3*WProd1434*WProd1447*WProd1537 - 
     -            WProd1460*WProd1537 - WProd1462*WProd1537 + 
     -            WProd1467*WProd1537 + WProd1468*WProd1537 + 
     -            WProd1469*WProd1537 - WProd1424*WProd1538 + 
     -            WProd1471*WProd1538 + 2*WProd1437*WProd1539 - 
     -            3*WProd1439*WProd1539 - 2*WProd1440*WProd1539 - 
     -            3*WProd1441*WProd1539 - 2*WProd1442*WProd1539 + 
     -            3*WProd1443*WProd1539 - 3*WProd1444*WProd1539 - 
     -            3*WProd1445*WProd1539 + 3*WProd1446*WProd1539 + 
     -            3*WProd1447*WProd1539 - WProd1540 - WProd1541 + 
     -            2*WProd1437*WProd1542 - 3*WProd1439*WProd1542 - 
     -            2*WProd1440*WProd1542 - 3*WProd1441*WProd1542 - 
     -            2*WProd1442*WProd1542 + 3*WProd1443*WProd1542 - 
     -            3*WProd1444*WProd1542 - 3*WProd1445*WProd1542 + 
     -            3*WProd1446*WProd1542 + 3*WProd1447*WProd1542 - 
     -            WProd1543 - WProd1544 - WProd1545))/WProd1528 + 
     -      Prp171*(-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1515 + WProd1536))/WProd1516) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd1435*WProd1520) - 
     -               WProd1436*WProd1520 - WProd1524 + 
     -               WProd1435*WProd1525 + WProd1436*WProd1525 + 
     -               2*WProd1437*WProd1526 - 
     -               2*WProd1439*WProd1526 - WProd1440*WProd1526 - 
     -               2*WProd1441*WProd1526 - WProd1442*WProd1526 + 
     -               2*WProd1443*WProd1526 - 
     -               2*WProd1444*WProd1526 - 
     -               2*WProd1445*WProd1526 + 
     -               2*WProd1446*WProd1526 + 
     -               2*WProd1447*WProd1526 - WProd1425*WProd1529 - 
     -               WProd1434*WProd1435*WProd1529 - 
     -               WProd1434*WProd1436*WProd1529 + 
     -               WProd1448*WProd1529 + WProd1465*WProd1529 + 
     -               WProd1477*WProd1529 + WProd1435*WProd1534 + 
     -               WProd1436*WProd1534 + WProd1424*WProd1535 - 
     -               WProd1458*WProd1535 - WProd1471*WProd1535 + 
     -               2*WProd1434*WProd1437*WProd1537 - 
     -               2*WProd1434*WProd1439*WProd1537 - 
     -               WProd1434*WProd1440*WProd1537 - 
     -               2*WProd1434*WProd1441*WProd1537 - 
     -               WProd1434*WProd1442*WProd1537 + 
     -               2*WProd1434*WProd1443*WProd1537 - 
     -               2*WProd1434*WProd1444*WProd1537 - 
     -               2*WProd1434*WProd1445*WProd1537 + 
     -               2*WProd1434*WProd1446*WProd1537 + 
     -               2*WProd1434*WProd1447*WProd1537 + 
     -               WProd1463*WProd1537 - WProd1467*WProd1537 - 
     -               WProd1468*WProd1537 - WProd1469*WProd1537 + 
     -               WProd1424*WProd1538 - WProd1458*WProd1538 - 
     -               WProd1471*WProd1538 - 2*WProd1437*WProd1539 + 
     -               2*WProd1439*WProd1539 + WProd1440*WProd1539 + 
     -               2*WProd1441*WProd1539 + WProd1442*WProd1539 - 
     -               2*WProd1443*WProd1539 + 
     -               2*WProd1444*WProd1539 + 
     -               2*WProd1445*WProd1539 - 
     -               2*WProd1446*WProd1539 - 
     -               2*WProd1447*WProd1539 + WProd1540 + 
     -               WProd1541 - 2*WProd1437*WProd1542 + 
     -               2*WProd1439*WProd1542 + WProd1440*WProd1542 + 
     -               2*WProd1441*WProd1542 + WProd1442*WProd1542 - 
     -               2*WProd1443*WProd1542 + 
     -               2*WProd1444*WProd1542 + 
     -               2*WProd1445*WProd1542 - 
     -               2*WProd1446*WProd1542 - 
     -               2*WProd1447*WProd1542 + WProd1543 + 
     -               WProd1544 + WProd1545))/WProd1528))
        amp(17,4,6) = Prp169*tmp1146 + 
     -   Prp166*Prp167*Prp168*
     -    ((-8*Cpl2*GS*WProd1515)/WProd1516 + 
     -      (8*Cpl2*GS*(WProd1422*WProd1517 + 
     -           WProd1423*WProd1517 + WProd1424*WProd1517 - 
     -           WProd1425*WProd1517 + WProd1430*WProd1517 - 
     -           WProd1434*WProd1435*WProd1517 - 
     -           WProd1434*WProd1436*WProd1517 + 
     -           WProd1448*WProd1517 - WProd1518 - 
     -           WProd1427*WProd1519 - WProd1427*WProd1520 + 
     -           WProd1427*WProd1521 + WProd1522 + WProd1523 - 
     -           WProd1524 + WProd1435*WProd1525 + 
     -           WProd1436*WProd1525 + WProd1437*WProd1526 - 
     -           WProd1439*WProd1526 - WProd1440*WProd1526 - 
     -           WProd1441*WProd1526 - WProd1442*WProd1526 + 
     -           WProd1443*WProd1526 - WProd1444*WProd1526 - 
     -           2*WProd1445*WProd1526 + WProd1446*WProd1526 + 
     -           WProd1447*WProd1526 + WProd1527))/WProd1528)
        amp(17,4,6) = amp(17,4,6) + 
     -   Prp170*(Prp166*Prp167*
     -       ((4*Cpl2*GS*(2*WProd1423*WProd1517 + 
     -              2*WProd1424*WProd1517 - 
     -              2*WProd1425*WProd1517 - 
     -              2*WProd1434*WProd1435*WProd1517 - 
     -              2*WProd1434*WProd1436*WProd1517 + 
     -              2*WProd1448*WProd1517 + 
     -              2*WProd1465*WProd1517 - 
     -              2*WProd1427*WProd1519 - 
     -              2*WProd1427*WProd1520 + 
     -              2*WProd1427*WProd1521 + 2*WProd1522 + 
     -              2*WProd1523 - 2*WProd1524 + 
     -              2*WProd1435*WProd1525 + 
     -              2*WProd1436*WProd1525 + 
     -              3*WProd1437*WProd1526 - 
     -              3*WProd1439*WProd1526 - 
     -              3*WProd1440*WProd1526 - 
     -              3*WProd1441*WProd1526 - 
     -              3*WProd1442*WProd1526 + 
     -              3*WProd1443*WProd1526 - 
     -              3*WProd1444*WProd1526 - 
     -              4*WProd1445*WProd1526 + 
     -              3*WProd1446*WProd1526 + 3*WProd1447*WProd1526))
     -           /WProd1528 - 
     -         (8*Cpl2*GS*(WProd1515 + WProd1536))/WProd1516) + 
     -      Prp167*Prp173*(-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1515 + WProd1536))/WProd1516) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1424*WProd1517 - 
     -               WProd1425*WProd1517 - 
     -               WProd1434*WProd1435*WProd1517 - 
     -               WProd1434*WProd1436*WProd1517 + 
     -               WProd1448*WProd1517 + WProd1465*WProd1517 + 
     -               WProd1477*WProd1517 - WProd1427*WProd1520 + 
     -               WProd1427*WProd1521 + WProd1523 - WProd1524 + 
     -               WProd1435*WProd1525 + WProd1436*WProd1525 + 
     -               2*WProd1437*WProd1526 - 
     -               2*WProd1439*WProd1526 - 
     -               2*WProd1440*WProd1526 - 
     -               2*WProd1441*WProd1526 - 
     -               2*WProd1442*WProd1526 + 
     -               2*WProd1443*WProd1526 - 
     -               2*WProd1444*WProd1526 - 
     -               3*WProd1445*WProd1526 + 
     -               2*WProd1446*WProd1526 + 
     -               2*WProd1447*WProd1526 - 
     -               2*WProd1437*WProd1542 + 
     -               2*WProd1439*WProd1542 + 
     -               2*WProd1440*WProd1542 + 
     -               2*WProd1441*WProd1542 + 
     -               2*WProd1442*WProd1542 - 
     -               2*WProd1443*WProd1542 + 
     -               2*WProd1444*WProd1542 + 
     -               3*WProd1445*WProd1542 - 
     -               2*WProd1446*WProd1542 - 
     -               2*WProd1447*WProd1542 + WProd1544 - 
     -               WProd1427*WProd1546))/WProd1528) + 
     -      Prp174*(2.6666666666666665D0*
     -          (Cpl2*GS*Prp172*
     -             (WProd1522 + WProd1423*WProd1529 - 
     -               WProd1477*WProd1529 - WProd1423*WProd1535 + 
     -               WProd1477*WProd1535 + 
     -               WProd1434*WProd1439*WProd1537 + 
     -               WProd1434*WProd1441*WProd1537 - 
     -               WProd1434*WProd1443*WProd1537 - 
     -               WProd1460*WProd1537 + WProd1468*WProd1537 - 
     -               WProd1423*WProd1538 + WProd1477*WProd1538 - 
     -               WProd1439*WProd1542 - WProd1441*WProd1542 + 
     -               WProd1443*WProd1542 - WProd1544))/WProd1528 + 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp167*
     -             (WProd1423*WProd1517 - WProd1477*WProd1517 - 
     -               WProd1427*WProd1519 + WProd1522 + 
     -               WProd1437*WProd1542 - WProd1439*WProd1542 - 
     -               WProd1440*WProd1542 - WProd1441*WProd1542 - 
     -               WProd1442*WProd1542 + WProd1443*WProd1542 - 
     -               WProd1444*WProd1542 - WProd1445*WProd1542 + 
     -               WProd1446*WProd1542 + WProd1447*WProd1542 - 
     -               WProd1544 + WProd1427*WProd1546))/WProd1528)+
     -        Prp171*Prp173*
     -       (-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1515 + WProd1536))/WProd1516) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd1425*WProd1517) - 
     -               WProd1434*WProd1435*WProd1517 - 
     -               WProd1434*WProd1436*WProd1517 + 
     -               WProd1448*WProd1517 + WProd1465*WProd1517 + 
     -               WProd1477*WProd1517 - WProd1435*WProd1520 - 
     -               WProd1436*WProd1520 + WProd1427*WProd1521 - 
     -               WProd1524 + WProd1435*WProd1525 + 
     -               WProd1436*WProd1525 + 2*WProd1437*WProd1526 - 
     -               2*WProd1439*WProd1526 - WProd1440*WProd1526 - 
     -               2*WProd1441*WProd1526 - WProd1442*WProd1526 + 
     -               WProd1443*WProd1526 - 2*WProd1444*WProd1526 - 
     -               2*WProd1445*WProd1526 + WProd1446*WProd1526 + 
     -               WProd1447*WProd1526 + WProd1435*WProd1534 + 
     -               WProd1436*WProd1534 + WProd1424*WProd1535 - 
     -               WProd1425*WProd1535 + 
     -               WProd1427*WProd1434*WProd1535 - 
     -               WProd1434*WProd1435*WProd1535 - 
     -               WProd1434*WProd1436*WProd1535 + 
     -               WProd1448*WProd1535 - WProd1458*WProd1535 + 
     -               WProd1465*WProd1535 - WProd1471*WProd1535 + 
     -               WProd1477*WProd1535 + WProd1424*WProd1538 - 
     -               WProd1425*WProd1538 + 
     -               WProd1427*WProd1434*WProd1538 - 
     -               WProd1434*WProd1435*WProd1538 - 
     -               WProd1434*WProd1436*WProd1538 + 
     -               WProd1448*WProd1538 - WProd1458*WProd1538 + 
     -               WProd1465*WProd1538 - WProd1471*WProd1538 + 
     -               WProd1477*WProd1538 - 2*WProd1437*WProd1539 + 
     -               2*WProd1439*WProd1539 + WProd1440*WProd1539 + 
     -               2*WProd1441*WProd1539 + WProd1442*WProd1539 - 
     -               WProd1443*WProd1539 + 2*WProd1444*WProd1539 + 
     -               2*WProd1445*WProd1539 - WProd1446*WProd1539 - 
     -               WProd1447*WProd1539 + WProd1540 + WProd1541 - 
     -               2*WProd1437*WProd1542 + 
     -               2*WProd1439*WProd1542 + WProd1440*WProd1542 + 
     -               2*WProd1441*WProd1542 + WProd1442*WProd1542 - 
     -               WProd1443*WProd1542 + 2*WProd1444*WProd1542 + 
     -               2*WProd1445*WProd1542 - WProd1446*WProd1542 - 
     -               WProd1447*WProd1542 + WProd1543 + WProd1544 + 
     -               WProd1545 - WProd1427*WProd1546 - 
     -               WProd1427*WProd1547 - WProd1427*WProd1548))/
     -           WProd1528))

#ifdef DEBUG
	print *, 'amp(17,4,6) =', amp(17,4,6)
#endif

        tmp1147 = Prp166*Prp168*
     -   ((-8*Cpl2*GS*WProd1550)/WProd1549 + 
     -     (8*Cpl2*GS*(WProd1503*WProd1517 - WProd1494*WProd1520 - 
     -          WProd1495*WProd1520 + WProd1494*WProd1525 + 
     -          WProd1495*WProd1525 + WProd1485*WProd1529 + 
     -          WProd1486*WProd1529 + WProd1487*WProd1529 - 
     -          WProd1488*WProd1529 + WProd1490*WProd1529 - 
     -          WProd1434*WProd1494*WProd1529 - 
     -          WProd1434*WProd1495*WProd1529 + 
     -          WProd1497*WProd1529 - WProd1503*WProd1529 + 
     -          WProd1494*WProd1534 + WProd1495*WProd1534 - 
     -          WProd1503*WProd1535 - WProd1552 + WProd1553 - 
     -          WProd1555 - WProd1439*WProd1556 - 
     -          WProd1441*WProd1556 - WProd1442*WProd1556 + 
     -          WProd1443*WProd1556 - WProd1444*WProd1556 - 
     -          2*WProd1445*WProd1556 + WProd1446*WProd1556 + 
     -          WProd1447*WProd1556 + WProd1557 + WProd1558 + 
     -          WProd1559 - WProd1560 + WProd1561))/WProd1551)
        tmp1147 = tmp1147 + 
     -   Prp170*(Prp166*((-8*Cpl2*GS*(WProd1550 + WProd1562))/
     -          WProd1549 + 
     -         (4*Cpl2*GS*(-2*WProd1494*WProd1520 - 
     -              2*WProd1495*WProd1520 + 
     -              2*WProd1494*WProd1525 + 
     -              2*WProd1495*WProd1525 + 
     -              2*WProd1486*WProd1529 + 
     -              2*WProd1487*WProd1529 - 
     -              2*WProd1488*WProd1529 - 
     -              2*WProd1434*WProd1494*WProd1529 - 
     -              2*WProd1434*WProd1495*WProd1529 + 
     -              2*WProd1497*WProd1529 + 
     -              2*WProd1505*WProd1529 + 
     -              2*WProd1494*WProd1534 + 
     -              2*WProd1495*WProd1534 - 
     -              2*WProd1503*WProd1535 - 
     -              2*WProd1503*WProd1538 + 2*WProd1553 + 
     -              2*WProd1554 - 2*WProd1555 + 
     -              4*WProd1437*WProd1556 - 
     -              5*WProd1439*WProd1556 - 
     -              3*WProd1440*WProd1556 - 
     -              5*WProd1441*WProd1556 - 
     -              3*WProd1442*WProd1556 + 
     -              5*WProd1443*WProd1556 - 
     -              5*WProd1444*WProd1556 - 
     -              5*WProd1445*WProd1556 + 
     -              5*WProd1446*WProd1556 + 
     -              5*WProd1447*WProd1556 - 
     -              2*WProd1460*WProd1563 - 
     -              2*WProd1462*WProd1563 + 2*WProd1463*WProd1563))
     -           /WProd1551) + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp172*
     -          (WProd1486*WProd1529 + WProd1487*WProd1529 - 
     -            WProd1514*WProd1529 - WProd1487*WProd1535 + 
     -            WProd1508*WProd1535 - WProd1487*WProd1538 + 
     -            WProd1508*WProd1538 + WProd1553 + WProd1554 - 
     -            2*WProd1434*WProd1437*WProd1563 + 
     -            3*WProd1434*WProd1439*WProd1563 + 
     -            2*WProd1434*WProd1440*WProd1563 + 
     -            3*WProd1434*WProd1441*WProd1563 + 
     -            2*WProd1434*WProd1442*WProd1563 - 
     -            3*WProd1434*WProd1443*WProd1563 + 
     -            3*WProd1434*WProd1444*WProd1563 + 
     -            3*WProd1434*WProd1445*WProd1563 - 
     -            3*WProd1434*WProd1446*WProd1563 - 
     -            3*WProd1434*WProd1447*WProd1563 - 
     -            WProd1460*WProd1563 - WProd1462*WProd1563 + 
     -            WProd1467*WProd1563 + WProd1468*WProd1563 + 
     -            WProd1469*WProd1563 + 2*WProd1437*WProd1564 - 
     -            3*WProd1439*WProd1564 - 2*WProd1440*WProd1564 - 
     -            3*WProd1441*WProd1564 - 2*WProd1442*WProd1564 + 
     -            3*WProd1443*WProd1564 - 3*WProd1444*WProd1564 - 
     -            3*WProd1445*WProd1564 + 3*WProd1446*WProd1564 + 
     -            3*WProd1447*WProd1564 - WProd1565 - WProd1566 + 
     -            2*WProd1437*WProd1567 - 3*WProd1439*WProd1567 - 
     -            2*WProd1440*WProd1567 - 3*WProd1441*WProd1567 - 
     -            2*WProd1442*WProd1567 + 3*WProd1443*WProd1567 - 
     -            3*WProd1444*WProd1567 - 3*WProd1445*WProd1567 + 
     -            3*WProd1446*WProd1567 + 3*WProd1447*WProd1567 - 
     -            WProd1568 - WProd1569 - WProd1570))/WProd1551 + 
     -      Prp171*(-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1550 + WProd1562))/WProd1549) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd1494*WProd1520) - 
     -               WProd1495*WProd1520 + WProd1494*WProd1525 + 
     -               WProd1495*WProd1525 - WProd1488*WProd1529 - 
     -               WProd1434*WProd1494*WProd1529 - 
     -               WProd1434*WProd1495*WProd1529 + 
     -               WProd1497*WProd1529 + WProd1505*WProd1529 + 
     -               WProd1514*WProd1529 + WProd1494*WProd1534 + 
     -               WProd1495*WProd1534 + WProd1487*WProd1535 - 
     -               WProd1503*WProd1535 - WProd1508*WProd1535 + 
     -               WProd1487*WProd1538 - WProd1503*WProd1538 - 
     -               WProd1508*WProd1538 - WProd1555 + 
     -               2*WProd1437*WProd1556 - 
     -               2*WProd1439*WProd1556 - WProd1440*WProd1556 - 
     -               2*WProd1441*WProd1556 - WProd1442*WProd1556 + 
     -               2*WProd1443*WProd1556 - 
     -               2*WProd1444*WProd1556 - 
     -               2*WProd1445*WProd1556 + 
     -               2*WProd1446*WProd1556 + 
     -               2*WProd1447*WProd1556 + 
     -               2*WProd1434*WProd1437*WProd1563 - 
     -               2*WProd1434*WProd1439*WProd1563 - 
     -               WProd1434*WProd1440*WProd1563 - 
     -               2*WProd1434*WProd1441*WProd1563 - 
     -               WProd1434*WProd1442*WProd1563 + 
     -               2*WProd1434*WProd1443*WProd1563 - 
     -               2*WProd1434*WProd1444*WProd1563 - 
     -               2*WProd1434*WProd1445*WProd1563 + 
     -               2*WProd1434*WProd1446*WProd1563 + 
     -               2*WProd1434*WProd1447*WProd1563 + 
     -               WProd1463*WProd1563 - WProd1467*WProd1563 - 
     -               WProd1468*WProd1563 - WProd1469*WProd1563 - 
     -               2*WProd1437*WProd1564 + 
     -               2*WProd1439*WProd1564 + WProd1440*WProd1564 + 
     -               2*WProd1441*WProd1564 + WProd1442*WProd1564 - 
     -               2*WProd1443*WProd1564 + 
     -               2*WProd1444*WProd1564 + 
     -               2*WProd1445*WProd1564 - 
     -               2*WProd1446*WProd1564 - 
     -               2*WProd1447*WProd1564 + WProd1565 + 
     -               WProd1566 - 2*WProd1437*WProd1567 + 
     -               2*WProd1439*WProd1567 + WProd1440*WProd1567 + 
     -               2*WProd1441*WProd1567 + WProd1442*WProd1567 - 
     -               2*WProd1443*WProd1567 + 
     -               2*WProd1444*WProd1567 + 
     -               2*WProd1445*WProd1567 - 
     -               2*WProd1446*WProd1567 - 
     -               2*WProd1447*WProd1567 + WProd1568 + 
     -               WProd1569 + WProd1570))/WProd1551))
        amp(17,4,4) = Prp169*tmp1147 + 
     -   Prp166*Prp167*Prp168*
     -    ((-8*Cpl2*GS*WProd1550)/WProd1549 + 
     -      (8*Cpl2*GS*(WProd1485*WProd1517 + 
     -           WProd1486*WProd1517 + WProd1487*WProd1517 - 
     -           WProd1488*WProd1517 + WProd1490*WProd1517 - 
     -           WProd1434*WProd1494*WProd1517 - 
     -           WProd1434*WProd1495*WProd1517 + 
     -           WProd1497*WProd1517 - WProd1489*WProd1519 - 
     -           WProd1489*WProd1520 + WProd1489*WProd1521 + 
     -           WProd1494*WProd1525 + WProd1495*WProd1525 - 
     -           WProd1552 + WProd1553 + WProd1554 - WProd1555 + 
     -           WProd1437*WProd1556 - WProd1439*WProd1556 - 
     -           WProd1440*WProd1556 - WProd1441*WProd1556 - 
     -           WProd1442*WProd1556 + WProd1443*WProd1556 - 
     -           WProd1444*WProd1556 - 2*WProd1445*WProd1556 + 
     -           WProd1446*WProd1556 + WProd1447*WProd1556 + 
     -           WProd1557))/WProd1551)
        amp(17,4,4) = amp(17,4,4) + 
     -   Prp170*(Prp166*Prp167*
     -       ((4*Cpl2*GS*(2*WProd1486*WProd1517 + 
     -              2*WProd1487*WProd1517 - 
     -              2*WProd1488*WProd1517 - 
     -              2*WProd1434*WProd1494*WProd1517 - 
     -              2*WProd1434*WProd1495*WProd1517 + 
     -              2*WProd1497*WProd1517 + 
     -              2*WProd1505*WProd1517 - 
     -              2*WProd1489*WProd1519 - 
     -              2*WProd1489*WProd1520 + 
     -              2*WProd1489*WProd1521 + 
     -              2*WProd1494*WProd1525 + 
     -              2*WProd1495*WProd1525 + 2*WProd1553 + 
     -              2*WProd1554 - 2*WProd1555 + 
     -              3*WProd1437*WProd1556 - 
     -              3*WProd1439*WProd1556 - 
     -              3*WProd1440*WProd1556 - 
     -              3*WProd1441*WProd1556 - 
     -              3*WProd1442*WProd1556 + 
     -              3*WProd1443*WProd1556 - 
     -              3*WProd1444*WProd1556 - 
     -              4*WProd1445*WProd1556 + 
     -              3*WProd1446*WProd1556 + 3*WProd1447*WProd1556))
     -           /WProd1551 - 
     -         (8*Cpl2*GS*(WProd1550 + WProd1562))/WProd1549) + 
     -      Prp174*(2.6666666666666665D0*
     -          (Cpl2*GS*Prp172*
     -             (WProd1486*WProd1529 - WProd1514*WProd1529 - 
     -               WProd1486*WProd1535 + WProd1514*WProd1535 - 
     -               WProd1486*WProd1538 + WProd1514*WProd1538 + 
     -               WProd1553 + WProd1434*WProd1439*WProd1563 + 
     -               WProd1434*WProd1441*WProd1563 - 
     -               WProd1434*WProd1443*WProd1563 - 
     -               WProd1460*WProd1563 + WProd1468*WProd1563 - 
     -               WProd1439*WProd1567 - WProd1441*WProd1567 + 
     -               WProd1443*WProd1567 - WProd1569))/WProd1551 + 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp167*
     -             (WProd1486*WProd1517 - WProd1514*WProd1517 - 
     -               WProd1489*WProd1519 + WProd1489*WProd1546 + 
     -               WProd1553 + WProd1437*WProd1567 - 
     -               WProd1439*WProd1567 - WProd1440*WProd1567 - 
     -               WProd1441*WProd1567 - WProd1442*WProd1567 + 
     -               WProd1443*WProd1567 - WProd1444*WProd1567 - 
     -               WProd1445*WProd1567 + WProd1446*WProd1567 + 
     -               WProd1447*WProd1567 - WProd1569))/WProd1551)+
     -        Prp167*Prp173*
     -       (-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1550 + WProd1562))/WProd1549) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1487*WProd1517 - 
     -               WProd1488*WProd1517 - 
     -               WProd1434*WProd1494*WProd1517 - 
     -               WProd1434*WProd1495*WProd1517 + 
     -               WProd1497*WProd1517 + WProd1505*WProd1517 + 
     -               WProd1514*WProd1517 - WProd1489*WProd1520 + 
     -               WProd1489*WProd1521 + WProd1494*WProd1525 + 
     -               WProd1495*WProd1525 - WProd1489*WProd1546 + 
     -               WProd1554 - WProd1555 + 
     -               2*WProd1437*WProd1556 - 
     -               2*WProd1439*WProd1556 - 
     -               2*WProd1440*WProd1556 - 
     -               2*WProd1441*WProd1556 - 
     -               2*WProd1442*WProd1556 + 
     -               2*WProd1443*WProd1556 - 
     -               2*WProd1444*WProd1556 - 
     -               3*WProd1445*WProd1556 + 
     -               2*WProd1446*WProd1556 + 
     -               2*WProd1447*WProd1556 - 
     -               2*WProd1437*WProd1567 + 
     -               2*WProd1439*WProd1567 + 
     -               2*WProd1440*WProd1567 + 
     -               2*WProd1441*WProd1567 + 
     -               2*WProd1442*WProd1567 - 
     -               2*WProd1443*WProd1567 + 
     -               2*WProd1444*WProd1567 + 
     -               3*WProd1445*WProd1567 - 
     -               2*WProd1446*WProd1567 - 
     -               2*WProd1447*WProd1567 + WProd1569))/WProd1551)
     -        + Prp171*Prp173*
     -       (-(5.333333333333333D0*
     -            (Cpl2*GS*(WProd1550 + WProd1562))/WProd1549) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd1488*WProd1517) - 
     -               WProd1434*WProd1494*WProd1517 - 
     -               WProd1434*WProd1495*WProd1517 + 
     -               WProd1497*WProd1517 + WProd1505*WProd1517 + 
     -               WProd1514*WProd1517 - WProd1494*WProd1520 - 
     -               WProd1495*WProd1520 + WProd1489*WProd1521 + 
     -               WProd1494*WProd1525 + WProd1495*WProd1525 + 
     -               WProd1494*WProd1534 + WProd1495*WProd1534 + 
     -               WProd1487*WProd1535 - WProd1488*WProd1535 + 
     -               WProd1434*WProd1489*WProd1535 - 
     -               WProd1434*WProd1494*WProd1535 - 
     -               WProd1434*WProd1495*WProd1535 + 
     -               WProd1497*WProd1535 - WProd1503*WProd1535 + 
     -               WProd1505*WProd1535 - WProd1508*WProd1535 + 
     -               WProd1514*WProd1535 + WProd1487*WProd1538 - 
     -               WProd1488*WProd1538 + 
     -               WProd1434*WProd1489*WProd1538 - 
     -               WProd1434*WProd1494*WProd1538 - 
     -               WProd1434*WProd1495*WProd1538 + 
     -               WProd1497*WProd1538 - WProd1503*WProd1538 + 
     -               WProd1505*WProd1538 - WProd1508*WProd1538 + 
     -               WProd1514*WProd1538 - WProd1489*WProd1546 - 
     -               WProd1489*WProd1547 - WProd1489*WProd1548 - 
     -               WProd1555 + 2*WProd1437*WProd1556 - 
     -               2*WProd1439*WProd1556 - WProd1440*WProd1556 - 
     -               2*WProd1441*WProd1556 - WProd1442*WProd1556 + 
     -               WProd1443*WProd1556 - 2*WProd1444*WProd1556 - 
     -               2*WProd1445*WProd1556 + WProd1446*WProd1556 + 
     -               WProd1447*WProd1556 - 2*WProd1437*WProd1564 + 
     -               2*WProd1439*WProd1564 + WProd1440*WProd1564 + 
     -               2*WProd1441*WProd1564 + WProd1442*WProd1564 - 
     -               WProd1443*WProd1564 + 2*WProd1444*WProd1564 + 
     -               2*WProd1445*WProd1564 - WProd1446*WProd1564 - 
     -               WProd1447*WProd1564 + WProd1565 + WProd1566 - 
     -               2*WProd1437*WProd1567 + 
     -               2*WProd1439*WProd1567 + WProd1440*WProd1567 + 
     -               2*WProd1441*WProd1567 + WProd1442*WProd1567 - 
     -               WProd1443*WProd1567 + 2*WProd1444*WProd1567 + 
     -               2*WProd1445*WProd1567 - WProd1446*WProd1567 - 
     -               WProd1447*WProd1567 + WProd1568 + WProd1569 + 
     -               WProd1570))/WProd1551))

#ifdef DEBUG
	print *, 'amp(17,4,4) =', amp(17,4,4)
#endif

      endif

      end


