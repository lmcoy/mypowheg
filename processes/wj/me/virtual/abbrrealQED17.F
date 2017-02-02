      subroutine abbrrealQED17(resnr)
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

      
	WProd1549 = Prod(6,2,6,1)
	WProd1484 = Prod(2,2,2,1)*Prod(6,2,6,1)
	WProd1437 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1451 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd1529 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd1551 = Prod(6,2,6,1)*ProdC(2,2,2,1)
	WProd1439 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1440 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd1459 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd1535 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd1441 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1434 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd1460 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1426 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1525 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
        WProd1519 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,1,1))
	WProd1442 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd1464 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd1462 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1457 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1428 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1520 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
        WProd1467 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1482 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1468 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1481 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1469 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1483 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
	WProd1538 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd1534 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1547 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1546 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1548 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd1443 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1444 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd1445 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd1421 = Prod(5,1,2,2)*ProdC(5,1,2,1)
	WProd1517 = Prod(5,1,2,1)*ProdC(5,1,2,2)
	WProd1446 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd1447 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd1463 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1429 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1521 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(2,2,1,1)*
     -     ProdC(5,1,4,1))
	WProd1478 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1485 = -(Prod(5,1,1,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,1,1))
        WProd1486 = Prod(3,1,1,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd1552 = -(Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(2,1,2,2)*ProdC(4,1,1,1)*ProdC(6,1,1,1))
        WProd1487 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1506 = Prod(3,1,2,2)*Prod(6,2,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1498 = Prod(2,1,1,1)*Prod(5,1,2,2)*Prod(6,2,3,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,1,1)
        WProd1558 = Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(6,1,1,1)
        WProd1564 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,2,2)*
     -   ProdC(6,1,1,1)
        WProd1488 = Prod(5,1,3,1)*Prod(6,2,5,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
	WProd1479 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd1489 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd1504 = Prod(6,2,2,2)*ProdC(6,1,2,1)
        WProd1490 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,1))
        WProd1507 = -(Prod(5,1,2,2)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,1))
        WProd1508 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd1499 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,1))
        WProd1491 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,1))
        WProd1500 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1492 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1509 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1510 = Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,2,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1)
        WProd1511 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1512 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1513 = Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1)
        WProd1565 = Prod(3,1,2,1)**2*Prod(6,2,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,2)*ProdC(6,1,2,1)
        WProd1566 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,2)*ProdC(6,1,2,1)
        WProd1493 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd1501 = -(Prod(5,1,2,1)*Prod(5,1,2,2)*Prod(6,2,3,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd1502 = Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1)
        WProd1559 = -(Prod(5,1,2,1)**2*Prod(6,2,3,1)*
     -     ProdC(2,2,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
	WProd1563 = Prod(6,2,2,1)*ProdC(6,1,2,2)
        WProd1567 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,2))
        WProd1560 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,2))
        WProd1553 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,2))
        WProd1554 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1568 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1569 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1570 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1555 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd1561 = Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2)
	WProd1494 = Prod(6,2,3,1)*ProdC(6,1,3,1)
	WProd1450 = Prod(6,1,4,1)*ProdC(6,1,4,1)
	WProd1495 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd1503 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1496 = -(Prod(5,1,2,2)*Prod(6,2,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1556 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(2,2,1,1)*
     -     ProdC(6,1,4,1))
        WProd1550 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(6,1,4,1))
        WProd1497 = -(Prod(5,1,3,1)*Prod(6,2,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1562 = -(Prod(4,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,1,4,1))
        WProd1505 = -(Prod(5,1,4,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1514 = -(Prod(4,1,3,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
	WProd1480 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1557 = Prod(5,1,2,1)**2*Prod(6,2,3,1)*
     -   ProdC(2,1,2,2)*ProdC(4,1,1,1)*ProdC(6,1,5,1)
        WProd1422 = -(Prod(5,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,1,1))
        WProd1423 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,1,1)
        WProd1518 = -(Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(2,1,2,2)*ProdC(4,1,1,1)*ProdC(6,2,1,1))
        WProd1424 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1466 = Prod(3,1,2,2)*Prod(6,1,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1452 = Prod(2,1,1,1)*Prod(5,1,2,2)*Prod(6,1,3,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,1,1)
        WProd1530 = Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(6,2,1,1)
        WProd1539 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,2,2)*
     -   ProdC(6,2,1,1)
        WProd1425 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
	WProd1427 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd1461 = Prod(6,1,2,2)*ProdC(6,2,2,1)
        WProd1430 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,1))
        WProd1470 = -(Prod(5,1,2,2)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,1))
        WProd1471 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,2,1)
        WProd1453 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,1))
        WProd1431 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,1))
        WProd1454 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1432 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1472 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1473 = Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,1,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1)
        WProd1474 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1475 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1476 = Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1)
        WProd1540 = Prod(3,1,2,1)**2*Prod(6,1,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,2)*ProdC(6,2,2,1)
        WProd1541 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,2)*ProdC(6,2,2,1)
        WProd1433 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1455 = -(Prod(5,1,2,1)*Prod(5,1,2,2)*Prod(6,1,3,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1456 = Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1)
        WProd1531 = -(Prod(5,1,2,1)**2*Prod(6,1,3,1)*
     -     ProdC(2,2,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
	WProd1537 = Prod(6,1,2,1)*ProdC(6,2,2,2)
        WProd1542 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,2))
        WProd1532 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,2))
        WProd1522 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,2))
        WProd1523 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1543 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1544 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1545 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1524 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd1533 = Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2)
	WProd1435 = Prod(6,1,3,1)*ProdC(6,2,3,1)
	WProd1436 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd1458 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1438 = -(Prod(5,1,2,2)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1526 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(2,2,1,1)*
     -     ProdC(6,2,4,1))
        WProd1515 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(6,2,4,1))
        WProd1448 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1536 = -(Prod(4,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,2,4,1))
        WProd1465 = -(Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1477 = -(Prod(4,1,3,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1527 = Prod(5,1,2,1)**2*Prod(6,1,3,1)*
     -   ProdC(2,1,2,2)*ProdC(4,1,1,1)*ProdC(6,2,5,1)
	WProd1516 = ProdC(6,2,6,1)
	WProd1449 = Prod(2,2,2,1)*ProdC(6,2,6,1)
	WProd1528 = ProdC(2,2,2,1)*ProdC(6,2,6,1)
      
      
	Prp169 = 1.D0/WProd1437
	Prp170 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1443)
	Prp174 = 1.D0/(-WProd1439 - WProd1441 + WProd1443)
        Prp166 = 1.D0/
     -   ((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1437 - WProd1444 - 
     -     WProd1445)
	Prp167 = -(1.D0/WProd1445)
	Prp168 = 1.D0/WProd1450
	Prp173 = -(1.D0/WProd1478)
	Prp171 = 1.D0/(WProd1437 - WProd1478 - WProd1479)
	Prp172 = 1.D0/WProd1480
      
      
      
      endif

      end


