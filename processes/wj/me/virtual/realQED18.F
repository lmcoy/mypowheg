      subroutine realQED18(amp,resnr)
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

        tmp1148 = (32*Alfas*Cpl122*Prp176*Prp177*
     -      (WProd1572*WProd1573 - WProd1572*WProd1574 - 
     -        WProd1575*WProd1576 - WProd1576*WProd1577 + 
     -        WProd1576*WProd1578 - 2*WProd1571*WProd1580 - 
     -        2*WProd1579*WProd1580 + 2*WProd1580*WProd1581 + 
     -        2*WProd1580*WProd1582 + 2*WProd1580*WProd1583 + 
     -        2*WProd1580*WProd1584 + 2*WProd1580*WProd1585 + 
     -        3*WProd1580*WProd1586 - 2*WProd1580*WProd1587 - 
     -        2*WProd1580*WProd1588 + WProd1589 + WProd1590 - 
     -        WProd1591 - WProd1572*WProd1592*WProd1593 - 
     -        WProd1572*WProd1592*WProd1594 + 
     -        2*WProd1571*WProd1595 + 2*WProd1579*WProd1595 - 
     -        2*WProd1581*WProd1595 - 2*WProd1582*WProd1595 - 
     -        2*WProd1583*WProd1595 - 2*WProd1584*WProd1595 - 
     -        2*WProd1585*WProd1595 - 3*WProd1586*WProd1595 + 
     -        2*WProd1587*WProd1595 + 2*WProd1588*WProd1595 + 
     -        WProd1572*WProd1596 + WProd1572*WProd1597 + 
     -        WProd1572*WProd1598))/WProd1599 + 
     -   (32*Alfas*Cpl122*Prp177*Prp178*
     -      (-(WProd1572*WProd1574) - WProd1576*WProd1577 + 
     -        WProd1576*WProd1578 - WProd1571*WProd1580 - 
     -        2*WProd1579*WProd1580 + 2*WProd1580*WProd1581 + 
     -        WProd1580*WProd1582 + 2*WProd1580*WProd1583 + 
     -        WProd1580*WProd1584 + 2*WProd1580*WProd1585 + 
     -        2*WProd1580*WProd1586 - WProd1580*WProd1587 - 
     -        WProd1580*WProd1588 + WProd1590 - WProd1591 - 
     -        WProd1575*WProd1593 - 
     -        WProd1572*WProd1592*WProd1593 - 
     -        WProd1575*WProd1594 - 
     -        WProd1572*WProd1592*WProd1594 + 
     -        WProd1571*WProd1595 + 2*WProd1579*WProd1595 - 
     -        2*WProd1581*WProd1595 - WProd1582*WProd1595 - 
     -        2*WProd1583*WProd1595 - WProd1584*WProd1595 - 
     -        2*WProd1585*WProd1595 - 2*WProd1586*WProd1595 + 
     -        WProd1587*WProd1595 + WProd1588*WProd1595 + 
     -        WProd1572*WProd1596 + WProd1572*WProd1597 + 
     -        WProd1572*WProd1598 + WProd1573*WProd1601 - 
     -        WProd1574*WProd1601 + 
     -        WProd1576*WProd1592*WProd1601 - 
     -        WProd1592*WProd1593*WProd1601 - 
     -        WProd1592*WProd1594*WProd1601 + 
     -        WProd1596*WProd1601 + WProd1597*WProd1601 + 
     -        WProd1598*WProd1601 + WProd1573*WProd1602 - 
     -        WProd1574*WProd1602 + 
     -        WProd1576*WProd1592*WProd1602 - 
     -        WProd1592*WProd1593*WProd1602 - 
     -        WProd1592*WProd1594*WProd1602 + 
     -        WProd1596*WProd1602 + WProd1597*WProd1602 + 
     -        WProd1598*WProd1602 - WProd1571*WProd1603 - 
     -        2*WProd1579*WProd1603 + 2*WProd1581*WProd1603 + 
     -        WProd1582*WProd1603 + 2*WProd1583*WProd1603 + 
     -        WProd1584*WProd1603 + 2*WProd1585*WProd1603 + 
     -        2*WProd1586*WProd1603 - WProd1587*WProd1603 - 
     -        WProd1588*WProd1603 - WProd1576*WProd1604 - 
     -        WProd1576*WProd1605 - WProd1601*WProd1606 - 
     -        WProd1602*WProd1606 + WProd1607 + WProd1608 + 
     -        WProd1609 + WProd1610 + WProd1593*WProd1611 + 
     -        WProd1594*WProd1611 - WProd1601*WProd1612 - 
     -        WProd1602*WProd1612))/WProd1599 + 
     -   (16*Alfas*Cpl122*Prp178*Prp179*
     -      (-2*WProd1572*WProd1574 - 2*WProd1576*WProd1577 + 
     -        2*WProd1576*WProd1578 - 
     -        2*WProd1572*WProd1592*WProd1593 - 
     -        2*WProd1572*WProd1592*WProd1594 + 
     -        2*WProd1572*WProd1596 + 2*WProd1572*WProd1597 + 
     -        2*WProd1572*WProd1598 - 2*WProd1574*WProd1601 + 
     -        2*WProd1576*WProd1592*WProd1601 - 
     -        2*WProd1592*WProd1593*WProd1601 - 
     -        2*WProd1592*WProd1594*WProd1601 + 
     -        2*WProd1596*WProd1601 + 2*WProd1597*WProd1601 + 
     -        2*WProd1598*WProd1601 - 2*WProd1574*WProd1602 + 
     -        2*WProd1576*WProd1592*WProd1602 - 
     -        2*WProd1592*WProd1593*WProd1602 - 
     -        2*WProd1592*WProd1594*WProd1602 + 
     -        2*WProd1596*WProd1602 + 2*WProd1597*WProd1602 + 
     -        2*WProd1598*WProd1602 - 2*WProd1576*WProd1604 - 
     -        2*WProd1576*WProd1605 + 2*WProd1574*WProd1614 + 
     -        2*WProd1592*WProd1593*WProd1614 + 
     -        2*WProd1592*WProd1594*WProd1614 - 
     -        2*WProd1596*WProd1614 - 2*WProd1597*WProd1614 - 
     -        2*WProd1598*WProd1614 - 
     -        3*WProd1571*WProd1592*WProd1615 - 
     -        4*WProd1579*WProd1592*WProd1615 + 
     -        4*WProd1581*WProd1592*WProd1615 + 
     -        2*WProd1582*WProd1592*WProd1615 + 
     -        4*WProd1583*WProd1592*WProd1615 + 
     -        2*WProd1584*WProd1592*WProd1615 + 
     -        4*WProd1585*WProd1592*WProd1615 + 
     -        4*WProd1586*WProd1592*WProd1615 - 
     -        3*WProd1587*WProd1592*WProd1615 - 
     -        3*WProd1588*WProd1592*WProd1615 + 
     -        2*WProd1615*WProd1616 + 2*WProd1615*WProd1617 + 
     -        2*WProd1615*WProd1618 - 2*WProd1615*WProd1619))/
     -    WProd1599
        tmp1148 = tmp1148 + 
     -   Prp180*((-32*Alfas*Cpl122*Prp176*
     -         (WProd1576*WProd1577 + WProd1571*WProd1580 + 
     -           WProd1579*WProd1580 - WProd1580*WProd1581 - 
     -           WProd1580*WProd1582 - WProd1580*WProd1583 - 
     -           WProd1580*WProd1584 - WProd1580*WProd1585 - 
     -           WProd1580*WProd1586 + WProd1580*WProd1587 + 
     -           WProd1580*WProd1588 - WProd1590 - 
     -           WProd1572*WProd1598 + WProd1572*WProd1620 - 
     -           WProd1576*WProd1621 + WProd1622))/WProd1599 + 
     -      (16*Alfas*Cpl122*Prp179*
     -         (-2*WProd1576*WProd1577 + 2*WProd1572*WProd1598 + 
     -           2*WProd1598*WProd1601 + 2*WProd1598*WProd1602 - 
     -           2*WProd1598*WProd1614 - 
     -           2*WProd1571*WProd1592*WProd1615 - 
     -           WProd1579*WProd1592*WProd1615 + 
     -           2*WProd1581*WProd1592*WProd1615 + 
     -           WProd1582*WProd1592*WProd1615 + 
     -           2*WProd1583*WProd1592*WProd1615 + 
     -           WProd1584*WProd1592*WProd1615 + 
     -           WProd1585*WProd1592*WProd1615 + 
     -           WProd1586*WProd1592*WProd1615 - 
     -           WProd1587*WProd1592*WProd1615 - 
     -           WProd1588*WProd1592*WProd1615 + 
     -           2*WProd1615*WProd1617 - 2*WProd1572*WProd1620 - 
     -           2*WProd1601*WProd1620 - 2*WProd1602*WProd1620 + 
     -           2*WProd1614*WProd1620 + 2*WProd1576*WProd1621 - 
     -           2*WProd1615*WProd1623))/WProd1599)
	amp(18,5,1) = Prp175*tmp1148

#ifdef DEBUG
	print *, 'amp(18,5,1) =', amp(18,5,1)
#endif

        tmp1149 = (32*Alfas*Cpl122*Prp176*Prp177*
     -      (WProd1572*WProd1625 - WProd1572*WProd1626 - 
     -        WProd1575*WProd1627 - WProd1577*WProd1627 + 
     -        WProd1578*WProd1627 - 2*WProd1571*WProd1628 - 
     -        2*WProd1579*WProd1628 + 2*WProd1581*WProd1628 + 
     -        2*WProd1582*WProd1628 + 2*WProd1583*WProd1628 + 
     -        2*WProd1584*WProd1628 + 2*WProd1585*WProd1628 + 
     -        3*WProd1586*WProd1628 - 2*WProd1587*WProd1628 - 
     -        2*WProd1588*WProd1628 + WProd1629 + WProd1630 - 
     -        WProd1631 - WProd1572*WProd1592*WProd1632 - 
     -        WProd1572*WProd1592*WProd1633 + 
     -        2*WProd1571*WProd1634 + 2*WProd1579*WProd1634 - 
     -        2*WProd1581*WProd1634 - 2*WProd1582*WProd1634 - 
     -        2*WProd1583*WProd1634 - 2*WProd1584*WProd1634 - 
     -        2*WProd1585*WProd1634 - 3*WProd1586*WProd1634 + 
     -        2*WProd1587*WProd1634 + 2*WProd1588*WProd1634 + 
     -        WProd1572*WProd1635 + WProd1572*WProd1636 + 
     -        WProd1572*WProd1637))/WProd1624 + 
     -   (32*Alfas*Cpl122*Prp177*Prp178*
     -      (WProd1601*WProd1625 + WProd1602*WProd1625 - 
     -        WProd1572*WProd1626 - WProd1601*WProd1626 - 
     -        WProd1602*WProd1626 - WProd1577*WProd1627 + 
     -        WProd1578*WProd1627 + 
     -        WProd1592*WProd1601*WProd1627 + 
     -        WProd1592*WProd1602*WProd1627 - 
     -        WProd1604*WProd1627 - WProd1605*WProd1627 - 
     -        WProd1571*WProd1628 - 2*WProd1579*WProd1628 + 
     -        2*WProd1581*WProd1628 + WProd1582*WProd1628 + 
     -        2*WProd1583*WProd1628 + WProd1584*WProd1628 + 
     -        2*WProd1585*WProd1628 + 2*WProd1586*WProd1628 - 
     -        WProd1587*WProd1628 - WProd1588*WProd1628 + 
     -        WProd1630 - WProd1631 - WProd1575*WProd1632 - 
     -        WProd1572*WProd1592*WProd1632 - 
     -        WProd1592*WProd1601*WProd1632 - 
     -        WProd1592*WProd1602*WProd1632 + 
     -        WProd1611*WProd1632 - WProd1575*WProd1633 - 
     -        WProd1572*WProd1592*WProd1633 - 
     -        WProd1592*WProd1601*WProd1633 - 
     -        WProd1592*WProd1602*WProd1633 + 
     -        WProd1611*WProd1633 + WProd1571*WProd1634 + 
     -        2*WProd1579*WProd1634 - 2*WProd1581*WProd1634 - 
     -        WProd1582*WProd1634 - 2*WProd1583*WProd1634 - 
     -        WProd1584*WProd1634 - 2*WProd1585*WProd1634 - 
     -        2*WProd1586*WProd1634 + WProd1587*WProd1634 + 
     -        WProd1588*WProd1634 + WProd1572*WProd1635 + 
     -        WProd1601*WProd1635 + WProd1602*WProd1635 + 
     -        WProd1572*WProd1636 + WProd1601*WProd1636 + 
     -        WProd1602*WProd1636 + WProd1572*WProd1637 + 
     -        WProd1601*WProd1637 + WProd1602*WProd1637 - 
     -        WProd1571*WProd1638 - 2*WProd1579*WProd1638 + 
     -        2*WProd1581*WProd1638 + WProd1582*WProd1638 + 
     -        2*WProd1583*WProd1638 + WProd1584*WProd1638 + 
     -        2*WProd1585*WProd1638 + 2*WProd1586*WProd1638 - 
     -        WProd1587*WProd1638 - WProd1588*WProd1638 - 
     -        WProd1601*WProd1639 - WProd1602*WProd1639 + 
     -        WProd1640 + WProd1641 + WProd1642 + WProd1643 - 
     -        WProd1601*WProd1644 - WProd1602*WProd1644))/WProd1624
     -     + (16*Alfas*Cpl122*Prp178*Prp179*
     -      (-2*WProd1572*WProd1626 - 2*WProd1601*WProd1626 - 
     -        2*WProd1602*WProd1626 + 2*WProd1614*WProd1626 - 
     -        2*WProd1577*WProd1627 + 2*WProd1578*WProd1627 + 
     -        2*WProd1592*WProd1601*WProd1627 + 
     -        2*WProd1592*WProd1602*WProd1627 - 
     -        2*WProd1604*WProd1627 - 2*WProd1605*WProd1627 - 
     -        2*WProd1572*WProd1592*WProd1632 - 
     -        2*WProd1592*WProd1601*WProd1632 - 
     -        2*WProd1592*WProd1602*WProd1632 + 
     -        2*WProd1592*WProd1614*WProd1632 - 
     -        2*WProd1572*WProd1592*WProd1633 - 
     -        2*WProd1592*WProd1601*WProd1633 - 
     -        2*WProd1592*WProd1602*WProd1633 + 
     -        2*WProd1592*WProd1614*WProd1633 + 
     -        2*WProd1572*WProd1635 + 2*WProd1601*WProd1635 + 
     -        2*WProd1602*WProd1635 - 2*WProd1614*WProd1635 + 
     -        2*WProd1572*WProd1636 + 2*WProd1601*WProd1636 + 
     -        2*WProd1602*WProd1636 - 2*WProd1614*WProd1636 + 
     -        2*WProd1572*WProd1637 + 2*WProd1601*WProd1637 + 
     -        2*WProd1602*WProd1637 - 2*WProd1614*WProd1637 - 
     -        3*WProd1571*WProd1592*WProd1645 - 
     -        4*WProd1579*WProd1592*WProd1645 + 
     -        4*WProd1581*WProd1592*WProd1645 + 
     -        2*WProd1582*WProd1592*WProd1645 + 
     -        4*WProd1583*WProd1592*WProd1645 + 
     -        2*WProd1584*WProd1592*WProd1645 + 
     -        4*WProd1585*WProd1592*WProd1645 + 
     -        4*WProd1586*WProd1592*WProd1645 - 
     -        3*WProd1587*WProd1592*WProd1645 - 
     -        3*WProd1588*WProd1592*WProd1645 + 
     -        2*WProd1616*WProd1645 + 2*WProd1617*WProd1645 + 
     -        2*WProd1618*WProd1645 - 2*WProd1619*WProd1645))/
     -    WProd1624
        tmp1149 = tmp1149 + 
     -   Prp180*((16*Alfas*Cpl122*Prp179*
     -         (-2*WProd1577*WProd1627 + 2*WProd1621*WProd1627 + 
     -           2*WProd1572*WProd1637 + 2*WProd1601*WProd1637 + 
     -           2*WProd1602*WProd1637 - 2*WProd1614*WProd1637 - 
     -           2*WProd1571*WProd1592*WProd1645 - 
     -           WProd1579*WProd1592*WProd1645 + 
     -           2*WProd1581*WProd1592*WProd1645 + 
     -           WProd1582*WProd1592*WProd1645 + 
     -           2*WProd1583*WProd1592*WProd1645 + 
     -           WProd1584*WProd1592*WProd1645 + 
     -           WProd1585*WProd1592*WProd1645 + 
     -           WProd1586*WProd1592*WProd1645 - 
     -           WProd1587*WProd1592*WProd1645 - 
     -           WProd1588*WProd1592*WProd1645 + 
     -           2*WProd1617*WProd1645 - 2*WProd1623*WProd1645 - 
     -           2*WProd1572*WProd1646 - 2*WProd1601*WProd1646 - 
     -           2*WProd1602*WProd1646 + 2*WProd1614*WProd1646))/
     -       WProd1624 - (32*Alfas*Cpl122*Prp176*
     -         (WProd1577*WProd1627 - WProd1621*WProd1627 + 
     -           WProd1571*WProd1628 + WProd1579*WProd1628 - 
     -           WProd1581*WProd1628 - WProd1582*WProd1628 - 
     -           WProd1583*WProd1628 - WProd1584*WProd1628 - 
     -           WProd1585*WProd1628 - WProd1586*WProd1628 + 
     -           WProd1587*WProd1628 + WProd1588*WProd1628 - 
     -           WProd1630 - WProd1572*WProd1637 + 
     -           WProd1572*WProd1646 + WProd1647))/WProd1624)
	amp(18,5,8) = Prp175*tmp1149

#ifdef DEBUG
	print *, 'amp(18,5,8) =', amp(18,5,8)
#endif

        tmp1150 = Prp176*Prp177*
     -    ((-32*Alfas*Cpl122*(WProd1648 + WProd1649))/WProd1650 + 
     -      (32*Alfas*Cpl122*
     -         (WProd1573*WProd1651 - WProd1574*WProd1651 - 
     -           WProd1592*WProd1593*WProd1651 - 
     -           WProd1592*WProd1594*WProd1651 + 
     -           WProd1596*WProd1651 + WProd1597*WProd1651 + 
     -           WProd1598*WProd1651 - WProd1576*WProd1652 - 
     -           WProd1576*WProd1653 + WProd1576*WProd1654 - 
     -           2*WProd1571*WProd1655 - 2*WProd1579*WProd1655 + 
     -           2*WProd1581*WProd1655 + 2*WProd1582*WProd1655 + 
     -           2*WProd1583*WProd1655 + 2*WProd1584*WProd1655 + 
     -           2*WProd1585*WProd1655 + 3*WProd1586*WProd1655 - 
     -           2*WProd1587*WProd1655 - 2*WProd1588*WProd1655 + 
     -           WProd1656 + WProd1657 - WProd1658 + 
     -           WProd1593*WProd1659 + WProd1594*WProd1659 + 
     -           2*WProd1571*WProd1660 + 2*WProd1579*WProd1660 - 
     -           2*WProd1581*WProd1660 - 2*WProd1582*WProd1660 - 
     -           2*WProd1583*WProd1660 - 2*WProd1584*WProd1660 - 
     -           2*WProd1585*WProd1660 - 3*WProd1586*WProd1660 + 
     -           2*WProd1587*WProd1660 + 2*WProd1588*WProd1660))/
     -       WProd1661) + Prp177*Prp178*
     -    ((-32*Alfas*Cpl122*(WProd1648 + WProd1649))/WProd1650 + 
     -      (32*Alfas*Cpl122*
     -         (-(WProd1574*WProd1651) - 
     -           WProd1592*WProd1593*WProd1651 - 
     -           WProd1592*WProd1594*WProd1651 + 
     -           WProd1596*WProd1651 + WProd1597*WProd1651 + 
     -           WProd1598*WProd1651 - WProd1593*WProd1652 - 
     -           WProd1594*WProd1652 - WProd1576*WProd1653 + 
     -           WProd1576*WProd1654 - WProd1571*WProd1655 - 
     -           2*WProd1579*WProd1655 + 2*WProd1581*WProd1655 + 
     -           WProd1582*WProd1655 + 2*WProd1583*WProd1655 + 
     -           WProd1584*WProd1655 + 2*WProd1585*WProd1655 + 
     -           2*WProd1586*WProd1655 - WProd1587*WProd1655 - 
     -           WProd1588*WProd1655 + WProd1657 - WProd1658 + 
     -           WProd1593*WProd1659 + WProd1594*WProd1659 + 
     -           WProd1571*WProd1660 + 2*WProd1579*WProd1660 - 
     -           2*WProd1581*WProd1660 - WProd1582*WProd1660 - 
     -           2*WProd1583*WProd1660 - WProd1584*WProd1660 - 
     -           2*WProd1585*WProd1660 - 2*WProd1586*WProd1660 + 
     -           WProd1587*WProd1660 + WProd1588*WProd1660 + 
     -           WProd1573*WProd1662 - WProd1574*WProd1662 + 
     -           WProd1576*WProd1592*WProd1662 - 
     -           WProd1592*WProd1593*WProd1662 - 
     -           WProd1592*WProd1594*WProd1662 + 
     -           WProd1596*WProd1662 + WProd1597*WProd1662 + 
     -           WProd1598*WProd1662 - WProd1606*WProd1662 - 
     -           WProd1612*WProd1662 + WProd1573*WProd1663 - 
     -           WProd1574*WProd1663 + 
     -           WProd1576*WProd1592*WProd1663 - 
     -           WProd1592*WProd1593*WProd1663 - 
     -           WProd1592*WProd1594*WProd1663 + 
     -           WProd1596*WProd1663 + WProd1597*WProd1663 + 
     -           WProd1598*WProd1663 - WProd1606*WProd1663 - 
     -           WProd1612*WProd1663 - WProd1571*WProd1664 - 
     -           2*WProd1579*WProd1664 + 2*WProd1581*WProd1664 + 
     -           WProd1582*WProd1664 + 2*WProd1583*WProd1664 + 
     -           WProd1584*WProd1664 + 2*WProd1585*WProd1664 + 
     -           2*WProd1586*WProd1664 - WProd1587*WProd1664 - 
     -           WProd1588*WProd1664 - WProd1576*WProd1665 - 
     -           WProd1576*WProd1666 + WProd1667 + WProd1668 + 
     -           WProd1669 + WProd1670 + WProd1593*WProd1671 + 
     -           WProd1594*WProd1671))/WProd1661) + 
     -   (16*Alfas*Cpl122*Prp178*Prp179*
     -      (-2*WProd1574*WProd1651 - 
     -        2*WProd1592*WProd1593*WProd1651 - 
     -        2*WProd1592*WProd1594*WProd1651 + 
     -        2*WProd1596*WProd1651 + 2*WProd1597*WProd1651 + 
     -        2*WProd1598*WProd1651 - 2*WProd1576*WProd1653 + 
     -        2*WProd1576*WProd1654 - 2*WProd1574*WProd1662 + 
     -        2*WProd1576*WProd1592*WProd1662 - 
     -        2*WProd1592*WProd1593*WProd1662 - 
     -        2*WProd1592*WProd1594*WProd1662 + 
     -        2*WProd1596*WProd1662 + 2*WProd1597*WProd1662 + 
     -        2*WProd1598*WProd1662 - 2*WProd1574*WProd1663 + 
     -        2*WProd1576*WProd1592*WProd1663 - 
     -        2*WProd1592*WProd1593*WProd1663 - 
     -        2*WProd1592*WProd1594*WProd1663 + 
     -        2*WProd1596*WProd1663 + 2*WProd1597*WProd1663 + 
     -        2*WProd1598*WProd1663 - 2*WProd1576*WProd1665 - 
     -        2*WProd1576*WProd1666 + 2*WProd1574*WProd1672 + 
     -        2*WProd1592*WProd1593*WProd1672 + 
     -        2*WProd1592*WProd1594*WProd1672 - 
     -        2*WProd1596*WProd1672 - 2*WProd1597*WProd1672 - 
     -        2*WProd1598*WProd1672 - 
     -        3*WProd1571*WProd1592*WProd1673 - 
     -        4*WProd1579*WProd1592*WProd1673 + 
     -        4*WProd1581*WProd1592*WProd1673 + 
     -        2*WProd1582*WProd1592*WProd1673 + 
     -        4*WProd1583*WProd1592*WProd1673 + 
     -        2*WProd1584*WProd1592*WProd1673 + 
     -        4*WProd1585*WProd1592*WProd1673 + 
     -        4*WProd1586*WProd1592*WProd1673 - 
     -        3*WProd1587*WProd1592*WProd1673 - 
     -        3*WProd1588*WProd1592*WProd1673 + 
     -        2*WProd1616*WProd1673 + 2*WProd1617*WProd1673 + 
     -        2*WProd1618*WProd1673 - 2*WProd1619*WProd1673))/
     -    WProd1661
        tmp1150 = tmp1150 + 
     -   Prp180*((16*Alfas*Cpl122*Prp179*
     -         (2*WProd1598*WProd1651 - 2*WProd1620*WProd1651 - 
     -           2*WProd1576*WProd1653 + 2*WProd1598*WProd1662 - 
     -           2*WProd1620*WProd1662 + 2*WProd1598*WProd1663 - 
     -           2*WProd1620*WProd1663 - 2*WProd1598*WProd1672 + 
     -           2*WProd1620*WProd1672 - 
     -           2*WProd1571*WProd1592*WProd1673 - 
     -           WProd1579*WProd1592*WProd1673 + 
     -           2*WProd1581*WProd1592*WProd1673 + 
     -           WProd1582*WProd1592*WProd1673 + 
     -           2*WProd1583*WProd1592*WProd1673 + 
     -           WProd1584*WProd1592*WProd1673 + 
     -           WProd1585*WProd1592*WProd1673 + 
     -           WProd1586*WProd1592*WProd1673 - 
     -           WProd1587*WProd1592*WProd1673 - 
     -           WProd1588*WProd1592*WProd1673 + 
     -           2*WProd1617*WProd1673 - 2*WProd1623*WProd1673 + 
     -           2*WProd1576*WProd1674))/WProd1661 - 
     -      (32*Alfas*Cpl122*Prp176*
     -         (-(WProd1598*WProd1651) + WProd1620*WProd1651 + 
     -           WProd1576*WProd1653 + WProd1571*WProd1655 + 
     -           WProd1579*WProd1655 - WProd1581*WProd1655 - 
     -           WProd1582*WProd1655 - WProd1583*WProd1655 - 
     -           WProd1584*WProd1655 - WProd1585*WProd1655 - 
     -           WProd1586*WProd1655 + WProd1587*WProd1655 + 
     -           WProd1588*WProd1655 - WProd1657 - 
     -           WProd1576*WProd1674 + WProd1675))/WProd1661)
	amp(18,5,6) = Prp175*tmp1150

#ifdef DEBUG
	print *, 'amp(18,5,6) =', amp(18,5,6)
#endif

        tmp1151 = Prp176*Prp177*
     -    ((-32*Alfas*Cpl122*(WProd1677 + WProd1678))/WProd1676 + 
     -      (32*Alfas*Cpl122*
     -         (WProd1625*WProd1651 - WProd1626*WProd1651 - 
     -           WProd1592*WProd1632*WProd1651 - 
     -           WProd1592*WProd1633*WProd1651 + 
     -           WProd1635*WProd1651 + WProd1636*WProd1651 + 
     -           WProd1637*WProd1651 - WProd1627*WProd1652 - 
     -           WProd1627*WProd1653 + WProd1627*WProd1654 + 
     -           WProd1632*WProd1659 + WProd1633*WProd1659 - 
     -           2*WProd1571*WProd1680 - 2*WProd1579*WProd1680 + 
     -           2*WProd1581*WProd1680 + 2*WProd1582*WProd1680 + 
     -           2*WProd1583*WProd1680 + 2*WProd1584*WProd1680 + 
     -           2*WProd1585*WProd1680 + 3*WProd1586*WProd1680 - 
     -           2*WProd1587*WProd1680 - 2*WProd1588*WProd1680 + 
     -           WProd1681 + WProd1682 - WProd1683 + 
     -           2*WProd1571*WProd1684 + 2*WProd1579*WProd1684 - 
     -           2*WProd1581*WProd1684 - 2*WProd1582*WProd1684 - 
     -           2*WProd1583*WProd1684 - 2*WProd1584*WProd1684 - 
     -           2*WProd1585*WProd1684 - 3*WProd1586*WProd1684 + 
     -           2*WProd1587*WProd1684 + 2*WProd1588*WProd1684))/
     -       WProd1679) + Prp177*Prp178*
     -    ((-32*Alfas*Cpl122*(WProd1677 + WProd1678))/WProd1676 + 
     -      (32*Alfas*Cpl122*
     -         (-(WProd1626*WProd1651) - 
     -           WProd1592*WProd1632*WProd1651 - 
     -           WProd1592*WProd1633*WProd1651 + 
     -           WProd1635*WProd1651 + WProd1636*WProd1651 + 
     -           WProd1637*WProd1651 - WProd1632*WProd1652 - 
     -           WProd1633*WProd1652 - WProd1627*WProd1653 + 
     -           WProd1627*WProd1654 + WProd1632*WProd1659 + 
     -           WProd1633*WProd1659 + WProd1625*WProd1662 - 
     -           WProd1626*WProd1662 + 
     -           WProd1592*WProd1627*WProd1662 - 
     -           WProd1592*WProd1632*WProd1662 - 
     -           WProd1592*WProd1633*WProd1662 + 
     -           WProd1635*WProd1662 + WProd1636*WProd1662 + 
     -           WProd1637*WProd1662 - WProd1639*WProd1662 - 
     -           WProd1644*WProd1662 + WProd1625*WProd1663 - 
     -           WProd1626*WProd1663 + 
     -           WProd1592*WProd1627*WProd1663 - 
     -           WProd1592*WProd1632*WProd1663 - 
     -           WProd1592*WProd1633*WProd1663 + 
     -           WProd1635*WProd1663 + WProd1636*WProd1663 + 
     -           WProd1637*WProd1663 - WProd1639*WProd1663 - 
     -           WProd1644*WProd1663 - WProd1627*WProd1665 - 
     -           WProd1627*WProd1666 + WProd1632*WProd1671 + 
     -           WProd1633*WProd1671 - WProd1571*WProd1680 - 
     -           2*WProd1579*WProd1680 + 2*WProd1581*WProd1680 + 
     -           WProd1582*WProd1680 + 2*WProd1583*WProd1680 + 
     -           WProd1584*WProd1680 + 2*WProd1585*WProd1680 + 
     -           2*WProd1586*WProd1680 - WProd1587*WProd1680 - 
     -           WProd1588*WProd1680 + WProd1682 - WProd1683 + 
     -           WProd1571*WProd1684 + 2*WProd1579*WProd1684 - 
     -           2*WProd1581*WProd1684 - WProd1582*WProd1684 - 
     -           2*WProd1583*WProd1684 - WProd1584*WProd1684 - 
     -           2*WProd1585*WProd1684 - 2*WProd1586*WProd1684 + 
     -           WProd1587*WProd1684 + WProd1588*WProd1684 - 
     -           WProd1571*WProd1685 - 2*WProd1579*WProd1685 + 
     -           2*WProd1581*WProd1685 + WProd1582*WProd1685 + 
     -           2*WProd1583*WProd1685 + WProd1584*WProd1685 + 
     -           2*WProd1585*WProd1685 + 2*WProd1586*WProd1685 - 
     -           WProd1587*WProd1685 - WProd1588*WProd1685 + 
     -           WProd1686 + WProd1687 + WProd1688 + WProd1689))/
     -       WProd1679) + (16*Alfas*Cpl122*Prp178*Prp179*
     -      (-2*WProd1626*WProd1651 - 
     -        2*WProd1592*WProd1632*WProd1651 - 
     -        2*WProd1592*WProd1633*WProd1651 + 
     -        2*WProd1635*WProd1651 + 2*WProd1636*WProd1651 + 
     -        2*WProd1637*WProd1651 - 2*WProd1627*WProd1653 + 
     -        2*WProd1627*WProd1654 - 2*WProd1626*WProd1662 + 
     -        2*WProd1592*WProd1627*WProd1662 - 
     -        2*WProd1592*WProd1632*WProd1662 - 
     -        2*WProd1592*WProd1633*WProd1662 + 
     -        2*WProd1635*WProd1662 + 2*WProd1636*WProd1662 + 
     -        2*WProd1637*WProd1662 - 2*WProd1626*WProd1663 + 
     -        2*WProd1592*WProd1627*WProd1663 - 
     -        2*WProd1592*WProd1632*WProd1663 - 
     -        2*WProd1592*WProd1633*WProd1663 + 
     -        2*WProd1635*WProd1663 + 2*WProd1636*WProd1663 + 
     -        2*WProd1637*WProd1663 - 2*WProd1627*WProd1665 - 
     -        2*WProd1627*WProd1666 + 2*WProd1626*WProd1672 + 
     -        2*WProd1592*WProd1632*WProd1672 + 
     -        2*WProd1592*WProd1633*WProd1672 - 
     -        2*WProd1635*WProd1672 - 2*WProd1636*WProd1672 - 
     -        2*WProd1637*WProd1672 - 
     -        3*WProd1571*WProd1592*WProd1690 - 
     -        4*WProd1579*WProd1592*WProd1690 + 
     -        4*WProd1581*WProd1592*WProd1690 + 
     -        2*WProd1582*WProd1592*WProd1690 + 
     -        4*WProd1583*WProd1592*WProd1690 + 
     -        2*WProd1584*WProd1592*WProd1690 + 
     -        4*WProd1585*WProd1592*WProd1690 + 
     -        4*WProd1586*WProd1592*WProd1690 - 
     -        3*WProd1587*WProd1592*WProd1690 - 
     -        3*WProd1588*WProd1592*WProd1690 + 
     -        2*WProd1616*WProd1690 + 2*WProd1617*WProd1690 + 
     -        2*WProd1618*WProd1690 - 2*WProd1619*WProd1690))/
     -    WProd1679
        tmp1151 = tmp1151 + 
     -   Prp180*((16*Alfas*Cpl122*Prp179*
     -         (2*WProd1637*WProd1651 - 2*WProd1646*WProd1651 - 
     -           2*WProd1627*WProd1653 + 2*WProd1637*WProd1662 - 
     -           2*WProd1646*WProd1662 + 2*WProd1637*WProd1663 - 
     -           2*WProd1646*WProd1663 - 2*WProd1637*WProd1672 + 
     -           2*WProd1646*WProd1672 + 2*WProd1627*WProd1674 - 
     -           2*WProd1571*WProd1592*WProd1690 - 
     -           WProd1579*WProd1592*WProd1690 + 
     -           2*WProd1581*WProd1592*WProd1690 + 
     -           WProd1582*WProd1592*WProd1690 + 
     -           2*WProd1583*WProd1592*WProd1690 + 
     -           WProd1584*WProd1592*WProd1690 + 
     -           WProd1585*WProd1592*WProd1690 + 
     -           WProd1586*WProd1592*WProd1690 - 
     -           WProd1587*WProd1592*WProd1690 - 
     -           WProd1588*WProd1592*WProd1690 + 
     -           2*WProd1617*WProd1690 - 2*WProd1623*WProd1690))/
     -       WProd1679 - (32*Alfas*Cpl122*Prp176*
     -         (-(WProd1637*WProd1651) + WProd1646*WProd1651 + 
     -           WProd1627*WProd1653 - WProd1627*WProd1674 + 
     -           WProd1571*WProd1680 + WProd1579*WProd1680 - 
     -           WProd1581*WProd1680 - WProd1582*WProd1680 - 
     -           WProd1583*WProd1680 - WProd1584*WProd1680 - 
     -           WProd1585*WProd1680 - WProd1586*WProd1680 + 
     -           WProd1587*WProd1680 + WProd1588*WProd1680 - 
     -           WProd1682 + WProd1691))/WProd1679)
	amp(18,5,4) = Prp175*tmp1151

#ifdef DEBUG
	print *, 'amp(18,5,4) =', amp(18,5,4)
#endif

        amp(18,6,1) = Prp175*
     -   ((-16*Alfas*Cpl122*Prp178*Prp179*
     -        (-2*WProd1572*WProd1574 - 2*WProd1576*WProd1577 + 
     -          2*WProd1576*WProd1578 - 
     -          2*WProd1572*WProd1592*WProd1593 - 
     -          2*WProd1572*WProd1592*WProd1594 + 
     -          2*WProd1572*WProd1596 + 2*WProd1572*WProd1597 + 
     -          2*WProd1572*WProd1598 - 2*WProd1574*WProd1601 + 
     -          2*WProd1576*WProd1592*WProd1601 - 
     -          2*WProd1592*WProd1593*WProd1601 - 
     -          2*WProd1592*WProd1594*WProd1601 + 
     -          2*WProd1596*WProd1601 + 2*WProd1597*WProd1601 + 
     -          2*WProd1598*WProd1601 - 2*WProd1574*WProd1602 + 
     -          2*WProd1576*WProd1592*WProd1602 - 
     -          2*WProd1592*WProd1593*WProd1602 - 
     -          2*WProd1592*WProd1594*WProd1602 + 
     -          2*WProd1596*WProd1602 + 2*WProd1597*WProd1602 + 
     -          2*WProd1598*WProd1602 - 2*WProd1576*WProd1604 - 
     -          2*WProd1576*WProd1605 + 2*WProd1574*WProd1614 + 
     -          2*WProd1592*WProd1593*WProd1614 + 
     -          2*WProd1592*WProd1594*WProd1614 - 
     -          2*WProd1596*WProd1614 - 2*WProd1597*WProd1614 - 
     -          2*WProd1598*WProd1614 - 
     -          3*WProd1571*WProd1592*WProd1615 - 
     -          4*WProd1579*WProd1592*WProd1615 + 
     -          4*WProd1581*WProd1592*WProd1615 + 
     -          2*WProd1582*WProd1592*WProd1615 + 
     -          4*WProd1583*WProd1592*WProd1615 + 
     -          2*WProd1584*WProd1592*WProd1615 + 
     -          4*WProd1585*WProd1592*WProd1615 + 
     -          4*WProd1586*WProd1592*WProd1615 - 
     -          3*WProd1587*WProd1592*WProd1615 - 
     -          3*WProd1588*WProd1592*WProd1615 + 
     -          2*WProd1615*WProd1616 + 2*WProd1615*WProd1617 + 
     -          2*WProd1615*WProd1618 - 2*WProd1615*WProd1619))/
     -      WProd1599 + Prp180*
     -      ((-16*Alfas*Cpl122*Prp179*
     -           (-2*WProd1576*WProd1577 + 2*WProd1572*WProd1598 + 
     -             2*WProd1598*WProd1601 + 2*WProd1598*WProd1602 - 
     -             2*WProd1598*WProd1614 - 
     -             2*WProd1571*WProd1592*WProd1615 - 
     -             WProd1579*WProd1592*WProd1615 + 
     -             2*WProd1581*WProd1592*WProd1615 + 
     -             WProd1582*WProd1592*WProd1615 + 
     -             2*WProd1583*WProd1592*WProd1615 + 
     -             WProd1584*WProd1592*WProd1615 + 
     -             WProd1585*WProd1592*WProd1615 + 
     -             WProd1586*WProd1592*WProd1615 - 
     -             WProd1587*WProd1592*WProd1615 - 
     -             WProd1588*WProd1592*WProd1615 + 
     -             2*WProd1615*WProd1617 - 2*WProd1572*WProd1620 - 
     -             2*WProd1601*WProd1620 - 2*WProd1602*WProd1620 + 
     -             2*WProd1614*WProd1620 + 2*WProd1576*WProd1621 - 
     -             2*WProd1615*WProd1623))/WProd1599 - 
     -        (32*Alfas*Cpl122*Prp182*
     -           (WProd1571*WProd1580 - WProd1580*WProd1581 - 
     -             WProd1580*WProd1583 - WProd1590 + 
     -             WProd1598*WProd1601 + WProd1598*WProd1602 - 
     -             WProd1598*WProd1614 - 
     -             WProd1571*WProd1592*WProd1615 + 
     -             WProd1581*WProd1592*WProd1615 + 
     -             WProd1583*WProd1592*WProd1615 + 
     -             WProd1615*WProd1617 - WProd1601*WProd1620 - 
     -             WProd1602*WProd1620 + WProd1614*WProd1620 + 
     -             WProd1622 - WProd1615*WProd1623))/WProd1599))
        amp(18,6,1) = Prp175*Prp181*
     -    ((32*Alfas*Cpl122*Prp178*
     -         (-2*WProd1571*WProd1580 - 2*WProd1579*WProd1580 + 
     -           2*WProd1580*WProd1581 + WProd1580*WProd1582 + 
     -           2*WProd1580*WProd1583 + WProd1580*WProd1584 + 
     -           2*WProd1580*WProd1585 + 2*WProd1580*WProd1586 - 
     -           2*WProd1580*WProd1587 - 2*WProd1580*WProd1588 + 
     -           WProd1590 - WProd1591 - WProd1575*WProd1593 - 
     -           WProd1575*WProd1594 + 2*WProd1571*WProd1595 + 
     -           2*WProd1579*WProd1595 - 2*WProd1581*WProd1595 - 
     -           WProd1582*WProd1595 - 2*WProd1583*WProd1595 - 
     -           WProd1584*WProd1595 - 2*WProd1585*WProd1595 - 
     -           2*WProd1586*WProd1595 + 2*WProd1587*WProd1595 + 
     -           2*WProd1588*WProd1595 + WProd1573*WProd1601 + 
     -           WProd1573*WProd1602 - 2*WProd1571*WProd1603 - 
     -           2*WProd1579*WProd1603 + 2*WProd1581*WProd1603 + 
     -           WProd1582*WProd1603 + 2*WProd1583*WProd1603 + 
     -           WProd1584*WProd1603 + 2*WProd1585*WProd1603 + 
     -           2*WProd1586*WProd1603 - 2*WProd1587*WProd1603 - 
     -           2*WProd1588*WProd1603 - WProd1601*WProd1606 - 
     -           WProd1602*WProd1606 + WProd1607 + WProd1608 + 
     -           WProd1609 + WProd1610 + WProd1593*WProd1611 + 
     -           WProd1594*WProd1611 - WProd1601*WProd1612 - 
     -           WProd1602*WProd1612 - WProd1574*WProd1614 - 
     -           WProd1592*WProd1593*WProd1614 - 
     -           WProd1592*WProd1594*WProd1614 + 
     -           WProd1596*WProd1614 + WProd1597*WProd1614 + 
     -           WProd1598*WProd1614 + 
     -           2*WProd1571*WProd1592*WProd1615 + 
     -           2*WProd1579*WProd1592*WProd1615 - 
     -           2*WProd1581*WProd1592*WProd1615 - 
     -           WProd1582*WProd1592*WProd1615 - 
     -           2*WProd1583*WProd1592*WProd1615 - 
     -           WProd1584*WProd1592*WProd1615 - 
     -           2*WProd1585*WProd1592*WProd1615 - 
     -           2*WProd1586*WProd1592*WProd1615 + 
     -           2*WProd1587*WProd1592*WProd1615 + 
     -           2*WProd1588*WProd1592*WProd1615 - 
     -           WProd1615*WProd1616 - WProd1615*WProd1617 - 
     -           WProd1615*WProd1618 + WProd1615*WProd1619))/
     -       WProd1599 - (32*Alfas*Cpl122*Prp182*
     -         (3*WProd1571*WProd1580 + 2*WProd1579*WProd1580 - 
     -           3*WProd1580*WProd1581 - 2*WProd1580*WProd1582 - 
     -           3*WProd1580*WProd1583 - 2*WProd1580*WProd1584 - 
     -           3*WProd1580*WProd1585 - 3*WProd1580*WProd1586 + 
     -           3*WProd1580*WProd1587 + 3*WProd1580*WProd1588 + 
     -           WProd1589 - WProd1590 - WProd1573*WProd1601 - 
     -           WProd1573*WProd1602 + 3*WProd1571*WProd1603 + 
     -           2*WProd1579*WProd1603 - 3*WProd1581*WProd1603 - 
     -           2*WProd1582*WProd1603 - 3*WProd1583*WProd1603 - 
     -           2*WProd1584*WProd1603 - 3*WProd1585*WProd1603 - 
     -           3*WProd1586*WProd1603 + 3*WProd1587*WProd1603 + 
     -           3*WProd1588*WProd1603 + WProd1601*WProd1606 + 
     -           WProd1602*WProd1606 - WProd1607 - WProd1608 - 
     -           WProd1609 - WProd1610 + WProd1573*WProd1614 - 
     -           WProd1598*WProd1614 - 
     -           3*WProd1571*WProd1592*WProd1615 - 
     -           2*WProd1579*WProd1592*WProd1615 + 
     -           3*WProd1581*WProd1592*WProd1615 + 
     -           2*WProd1582*WProd1592*WProd1615 + 
     -           3*WProd1583*WProd1592*WProd1615 + 
     -           2*WProd1584*WProd1592*WProd1615 + 
     -           3*WProd1585*WProd1592*WProd1615 + 
     -           3*WProd1586*WProd1592*WProd1615 - 
     -           3*WProd1587*WProd1592*WProd1615 - 
     -           3*WProd1588*WProd1592*WProd1615 + 
     -           WProd1615*WProd1616 + WProd1615*WProd1617 + 
     -           WProd1615*WProd1618 + WProd1614*WProd1620 + 
     -           WProd1622 - WProd1615*WProd1623 - 
     -           WProd1615*WProd1692))/WProd1599) + amp(18,6,1)

#ifdef DEBUG
	print *, 'amp(18,6,1) =', amp(18,6,1)
#endif

        amp(18,6,8) = Prp175*
     -   ((-16*Alfas*Cpl122*Prp178*Prp179*
     -        (-2*WProd1572*WProd1626 - 2*WProd1601*WProd1626 - 
     -          2*WProd1602*WProd1626 + 2*WProd1614*WProd1626 - 
     -          2*WProd1577*WProd1627 + 2*WProd1578*WProd1627 + 
     -          2*WProd1592*WProd1601*WProd1627 + 
     -          2*WProd1592*WProd1602*WProd1627 - 
     -          2*WProd1604*WProd1627 - 2*WProd1605*WProd1627 - 
     -          2*WProd1572*WProd1592*WProd1632 - 
     -          2*WProd1592*WProd1601*WProd1632 - 
     -          2*WProd1592*WProd1602*WProd1632 + 
     -          2*WProd1592*WProd1614*WProd1632 - 
     -          2*WProd1572*WProd1592*WProd1633 - 
     -          2*WProd1592*WProd1601*WProd1633 - 
     -          2*WProd1592*WProd1602*WProd1633 + 
     -          2*WProd1592*WProd1614*WProd1633 + 
     -          2*WProd1572*WProd1635 + 2*WProd1601*WProd1635 + 
     -          2*WProd1602*WProd1635 - 2*WProd1614*WProd1635 + 
     -          2*WProd1572*WProd1636 + 2*WProd1601*WProd1636 + 
     -          2*WProd1602*WProd1636 - 2*WProd1614*WProd1636 + 
     -          2*WProd1572*WProd1637 + 2*WProd1601*WProd1637 + 
     -          2*WProd1602*WProd1637 - 2*WProd1614*WProd1637 - 
     -          3*WProd1571*WProd1592*WProd1645 - 
     -          4*WProd1579*WProd1592*WProd1645 + 
     -          4*WProd1581*WProd1592*WProd1645 + 
     -          2*WProd1582*WProd1592*WProd1645 + 
     -          4*WProd1583*WProd1592*WProd1645 + 
     -          2*WProd1584*WProd1592*WProd1645 + 
     -          4*WProd1585*WProd1592*WProd1645 + 
     -          4*WProd1586*WProd1592*WProd1645 - 
     -          3*WProd1587*WProd1592*WProd1645 - 
     -          3*WProd1588*WProd1592*WProd1645 + 
     -          2*WProd1616*WProd1645 + 2*WProd1617*WProd1645 + 
     -          2*WProd1618*WProd1645 - 2*WProd1619*WProd1645))/
     -      WProd1624 + Prp180*
     -      ((-16*Alfas*Cpl122*Prp179*
     -           (-2*WProd1577*WProd1627 + 2*WProd1621*WProd1627 + 
     -             2*WProd1572*WProd1637 + 2*WProd1601*WProd1637 + 
     -             2*WProd1602*WProd1637 - 2*WProd1614*WProd1637 - 
     -             2*WProd1571*WProd1592*WProd1645 - 
     -             WProd1579*WProd1592*WProd1645 + 
     -             2*WProd1581*WProd1592*WProd1645 + 
     -             WProd1582*WProd1592*WProd1645 + 
     -             2*WProd1583*WProd1592*WProd1645 + 
     -             WProd1584*WProd1592*WProd1645 + 
     -             WProd1585*WProd1592*WProd1645 + 
     -             WProd1586*WProd1592*WProd1645 - 
     -             WProd1587*WProd1592*WProd1645 - 
     -             WProd1588*WProd1592*WProd1645 + 
     -             2*WProd1617*WProd1645 - 2*WProd1623*WProd1645 - 
     -             2*WProd1572*WProd1646 - 2*WProd1601*WProd1646 - 
     -             2*WProd1602*WProd1646 + 2*WProd1614*WProd1646))/
     -         WProd1624 - 
     -        (32*Alfas*Cpl122*Prp182*
     -           (WProd1571*WProd1628 - WProd1581*WProd1628 - 
     -             WProd1583*WProd1628 - WProd1630 + 
     -             WProd1601*WProd1637 + WProd1602*WProd1637 - 
     -             WProd1614*WProd1637 - 
     -             WProd1571*WProd1592*WProd1645 + 
     -             WProd1581*WProd1592*WProd1645 + 
     -             WProd1583*WProd1592*WProd1645 + 
     -             WProd1617*WProd1645 - WProd1623*WProd1645 - 
     -             WProd1601*WProd1646 - WProd1602*WProd1646 + 
     -             WProd1614*WProd1646 + WProd1647))/WProd1624))
        amp(18,6,8) = Prp175*Prp181*
     -    ((32*Alfas*Cpl122*Prp178*
     -         (WProd1601*WProd1625 + WProd1602*WProd1625 - 
     -           WProd1614*WProd1626 - 2*WProd1571*WProd1628 - 
     -           2*WProd1579*WProd1628 + 2*WProd1581*WProd1628 + 
     -           WProd1582*WProd1628 + 2*WProd1583*WProd1628 + 
     -           WProd1584*WProd1628 + 2*WProd1585*WProd1628 + 
     -           2*WProd1586*WProd1628 - 2*WProd1587*WProd1628 - 
     -           2*WProd1588*WProd1628 + WProd1630 - WProd1631 - 
     -           WProd1575*WProd1632 + WProd1611*WProd1632 - 
     -           WProd1592*WProd1614*WProd1632 - 
     -           WProd1575*WProd1633 + WProd1611*WProd1633 - 
     -           WProd1592*WProd1614*WProd1633 + 
     -           2*WProd1571*WProd1634 + 2*WProd1579*WProd1634 - 
     -           2*WProd1581*WProd1634 - WProd1582*WProd1634 - 
     -           2*WProd1583*WProd1634 - WProd1584*WProd1634 - 
     -           2*WProd1585*WProd1634 - 2*WProd1586*WProd1634 + 
     -           2*WProd1587*WProd1634 + 2*WProd1588*WProd1634 + 
     -           WProd1614*WProd1635 + WProd1614*WProd1636 + 
     -           WProd1614*WProd1637 - 2*WProd1571*WProd1638 - 
     -           2*WProd1579*WProd1638 + 2*WProd1581*WProd1638 + 
     -           WProd1582*WProd1638 + 2*WProd1583*WProd1638 + 
     -           WProd1584*WProd1638 + 2*WProd1585*WProd1638 + 
     -           2*WProd1586*WProd1638 - 2*WProd1587*WProd1638 - 
     -           2*WProd1588*WProd1638 - WProd1601*WProd1639 - 
     -           WProd1602*WProd1639 + WProd1640 + WProd1641 + 
     -           WProd1642 + WProd1643 - WProd1601*WProd1644 - 
     -           WProd1602*WProd1644 + 
     -           2*WProd1571*WProd1592*WProd1645 + 
     -           2*WProd1579*WProd1592*WProd1645 - 
     -           2*WProd1581*WProd1592*WProd1645 - 
     -           WProd1582*WProd1592*WProd1645 - 
     -           2*WProd1583*WProd1592*WProd1645 - 
     -           WProd1584*WProd1592*WProd1645 - 
     -           2*WProd1585*WProd1592*WProd1645 - 
     -           2*WProd1586*WProd1592*WProd1645 + 
     -           2*WProd1587*WProd1592*WProd1645 + 
     -           2*WProd1588*WProd1592*WProd1645 - 
     -           WProd1616*WProd1645 - WProd1617*WProd1645 - 
     -           WProd1618*WProd1645 + WProd1619*WProd1645))/
     -       WProd1624 - (32*Alfas*Cpl122*Prp182*
     -         (-(WProd1601*WProd1625) - WProd1602*WProd1625 + 
     -           WProd1614*WProd1625 + 3*WProd1571*WProd1628 + 
     -           2*WProd1579*WProd1628 - 3*WProd1581*WProd1628 - 
     -           2*WProd1582*WProd1628 - 3*WProd1583*WProd1628 - 
     -           2*WProd1584*WProd1628 - 3*WProd1585*WProd1628 - 
     -           3*WProd1586*WProd1628 + 3*WProd1587*WProd1628 + 
     -           3*WProd1588*WProd1628 + WProd1629 - WProd1630 - 
     -           WProd1614*WProd1637 + 3*WProd1571*WProd1638 + 
     -           2*WProd1579*WProd1638 - 3*WProd1581*WProd1638 - 
     -           2*WProd1582*WProd1638 - 3*WProd1583*WProd1638 - 
     -           2*WProd1584*WProd1638 - 3*WProd1585*WProd1638 - 
     -           3*WProd1586*WProd1638 + 3*WProd1587*WProd1638 + 
     -           3*WProd1588*WProd1638 + WProd1601*WProd1639 + 
     -           WProd1602*WProd1639 - WProd1640 - WProd1641 - 
     -           WProd1642 - WProd1643 - 
     -           3*WProd1571*WProd1592*WProd1645 - 
     -           2*WProd1579*WProd1592*WProd1645 + 
     -           3*WProd1581*WProd1592*WProd1645 + 
     -           2*WProd1582*WProd1592*WProd1645 + 
     -           3*WProd1583*WProd1592*WProd1645 + 
     -           2*WProd1584*WProd1592*WProd1645 + 
     -           3*WProd1585*WProd1592*WProd1645 + 
     -           3*WProd1586*WProd1592*WProd1645 - 
     -           3*WProd1587*WProd1592*WProd1645 - 
     -           3*WProd1588*WProd1592*WProd1645 + 
     -           WProd1616*WProd1645 + WProd1617*WProd1645 + 
     -           WProd1618*WProd1645 - WProd1623*WProd1645 + 
     -           WProd1614*WProd1646 + WProd1647 - 
     -           WProd1645*WProd1692))/WProd1624) + amp(18,6,8)

#ifdef DEBUG
	print *, 'amp(18,6,8) =', amp(18,6,8)
#endif

        amp(18,6,6) = Prp175*
     -   ((-16*Alfas*Cpl122*Prp178*Prp179*
     -        (-2*WProd1574*WProd1651 - 
     -          2*WProd1592*WProd1593*WProd1651 - 
     -          2*WProd1592*WProd1594*WProd1651 + 
     -          2*WProd1596*WProd1651 + 2*WProd1597*WProd1651 + 
     -          2*WProd1598*WProd1651 - 2*WProd1576*WProd1653 + 
     -          2*WProd1576*WProd1654 - 2*WProd1574*WProd1662 + 
     -          2*WProd1576*WProd1592*WProd1662 - 
     -          2*WProd1592*WProd1593*WProd1662 - 
     -          2*WProd1592*WProd1594*WProd1662 + 
     -          2*WProd1596*WProd1662 + 2*WProd1597*WProd1662 + 
     -          2*WProd1598*WProd1662 - 2*WProd1574*WProd1663 + 
     -          2*WProd1576*WProd1592*WProd1663 - 
     -          2*WProd1592*WProd1593*WProd1663 - 
     -          2*WProd1592*WProd1594*WProd1663 + 
     -          2*WProd1596*WProd1663 + 2*WProd1597*WProd1663 + 
     -          2*WProd1598*WProd1663 - 2*WProd1576*WProd1665 - 
     -          2*WProd1576*WProd1666 + 2*WProd1574*WProd1672 + 
     -          2*WProd1592*WProd1593*WProd1672 + 
     -          2*WProd1592*WProd1594*WProd1672 - 
     -          2*WProd1596*WProd1672 - 2*WProd1597*WProd1672 - 
     -          2*WProd1598*WProd1672 - 
     -          3*WProd1571*WProd1592*WProd1673 - 
     -          4*WProd1579*WProd1592*WProd1673 + 
     -          4*WProd1581*WProd1592*WProd1673 + 
     -          2*WProd1582*WProd1592*WProd1673 + 
     -          4*WProd1583*WProd1592*WProd1673 + 
     -          2*WProd1584*WProd1592*WProd1673 + 
     -          4*WProd1585*WProd1592*WProd1673 + 
     -          4*WProd1586*WProd1592*WProd1673 - 
     -          3*WProd1587*WProd1592*WProd1673 - 
     -          3*WProd1588*WProd1592*WProd1673 + 
     -          2*WProd1616*WProd1673 + 2*WProd1617*WProd1673 + 
     -          2*WProd1618*WProd1673 - 2*WProd1619*WProd1673))/
     -      WProd1661 + Prp180*
     -      ((-16*Alfas*Cpl122*Prp179*
     -           (2*WProd1598*WProd1651 - 2*WProd1620*WProd1651 - 
     -             2*WProd1576*WProd1653 + 2*WProd1598*WProd1662 - 
     -             2*WProd1620*WProd1662 + 2*WProd1598*WProd1663 - 
     -             2*WProd1620*WProd1663 - 2*WProd1598*WProd1672 + 
     -             2*WProd1620*WProd1672 - 
     -             2*WProd1571*WProd1592*WProd1673 - 
     -             WProd1579*WProd1592*WProd1673 + 
     -             2*WProd1581*WProd1592*WProd1673 + 
     -             WProd1582*WProd1592*WProd1673 + 
     -             2*WProd1583*WProd1592*WProd1673 + 
     -             WProd1584*WProd1592*WProd1673 + 
     -             WProd1585*WProd1592*WProd1673 + 
     -             WProd1586*WProd1592*WProd1673 - 
     -             WProd1587*WProd1592*WProd1673 - 
     -             WProd1588*WProd1592*WProd1673 + 
     -             2*WProd1617*WProd1673 - 2*WProd1623*WProd1673 + 
     -             2*WProd1576*WProd1674))/WProd1661 - 
     -        (32*Alfas*Cpl122*Prp182*
     -           (WProd1571*WProd1655 - WProd1581*WProd1655 - 
     -             WProd1583*WProd1655 - WProd1657 + 
     -             WProd1598*WProd1662 - WProd1620*WProd1662 + 
     -             WProd1598*WProd1663 - WProd1620*WProd1663 - 
     -             WProd1598*WProd1672 + WProd1620*WProd1672 - 
     -             WProd1571*WProd1592*WProd1673 + 
     -             WProd1581*WProd1592*WProd1673 + 
     -             WProd1583*WProd1592*WProd1673 + 
     -             WProd1617*WProd1673 - WProd1623*WProd1673 + 
     -             WProd1675))/WProd1661))
        amp(18,6,6) = Prp175*Prp181*
     -    (Prp178*((-32*Alfas*Cpl122*(WProd1648 + WProd1649))/
     -          WProd1650 + 
     -         (32*Alfas*Cpl122*
     -            (-(WProd1593*WProd1652) - WProd1594*WProd1652 - 
     -              2*WProd1571*WProd1655 - 
     -              2*WProd1579*WProd1655 + 
     -              2*WProd1581*WProd1655 + WProd1582*WProd1655 + 
     -              2*WProd1583*WProd1655 + WProd1584*WProd1655 + 
     -              2*WProd1585*WProd1655 + 
     -              2*WProd1586*WProd1655 - 
     -              2*WProd1587*WProd1655 - 
     -              2*WProd1588*WProd1655 + WProd1657 - 
     -              WProd1658 + WProd1593*WProd1659 + 
     -              WProd1594*WProd1659 + 2*WProd1571*WProd1660 + 
     -              2*WProd1579*WProd1660 - 
     -              2*WProd1581*WProd1660 - WProd1582*WProd1660 - 
     -              2*WProd1583*WProd1660 - WProd1584*WProd1660 - 
     -              2*WProd1585*WProd1660 - 
     -              2*WProd1586*WProd1660 + 
     -              2*WProd1587*WProd1660 + 
     -              2*WProd1588*WProd1660 + WProd1573*WProd1662 - 
     -              WProd1606*WProd1662 - WProd1612*WProd1662 + 
     -              WProd1573*WProd1663 - WProd1606*WProd1663 - 
     -              WProd1612*WProd1663 - 2*WProd1571*WProd1664 - 
     -              2*WProd1579*WProd1664 + 
     -              2*WProd1581*WProd1664 + WProd1582*WProd1664 + 
     -              2*WProd1583*WProd1664 + WProd1584*WProd1664 + 
     -              2*WProd1585*WProd1664 + 
     -              2*WProd1586*WProd1664 - 
     -              2*WProd1587*WProd1664 - 
     -              2*WProd1588*WProd1664 + WProd1667 + 
     -              WProd1668 + WProd1669 + WProd1670 + 
     -              WProd1593*WProd1671 + WProd1594*WProd1671 - 
     -              WProd1574*WProd1672 - 
     -              WProd1592*WProd1593*WProd1672 - 
     -              WProd1592*WProd1594*WProd1672 + 
     -              WProd1596*WProd1672 + WProd1597*WProd1672 + 
     -              WProd1598*WProd1672 + 
     -              2*WProd1571*WProd1592*WProd1673 + 
     -              2*WProd1579*WProd1592*WProd1673 - 
     -              2*WProd1581*WProd1592*WProd1673 - 
     -              WProd1582*WProd1592*WProd1673 - 
     -              2*WProd1583*WProd1592*WProd1673 - 
     -              WProd1584*WProd1592*WProd1673 - 
     -              2*WProd1585*WProd1592*WProd1673 - 
     -              2*WProd1586*WProd1592*WProd1673 + 
     -              2*WProd1587*WProd1592*WProd1673 + 
     -              2*WProd1588*WProd1592*WProd1673 - 
     -              WProd1616*WProd1673 - WProd1617*WProd1673 - 
     -              WProd1618*WProd1673 + WProd1619*WProd1673))/
     -          WProd1661) - 
     -      (32*Alfas*Cpl122*Prp182*
     -         (3*WProd1571*WProd1655 + 2*WProd1579*WProd1655 - 
     -           3*WProd1581*WProd1655 - 2*WProd1582*WProd1655 - 
     -           3*WProd1583*WProd1655 - 2*WProd1584*WProd1655 - 
     -           3*WProd1585*WProd1655 - 3*WProd1586*WProd1655 + 
     -           3*WProd1587*WProd1655 + 3*WProd1588*WProd1655 + 
     -           WProd1656 - WProd1657 - WProd1573*WProd1662 + 
     -           WProd1606*WProd1662 - WProd1573*WProd1663 + 
     -           WProd1606*WProd1663 + 3*WProd1571*WProd1664 + 
     -           2*WProd1579*WProd1664 - 3*WProd1581*WProd1664 - 
     -           2*WProd1582*WProd1664 - 3*WProd1583*WProd1664 - 
     -           2*WProd1584*WProd1664 - 3*WProd1585*WProd1664 - 
     -           3*WProd1586*WProd1664 + 3*WProd1587*WProd1664 + 
     -           3*WProd1588*WProd1664 - WProd1667 - WProd1668 - 
     -           WProd1669 - WProd1670 + WProd1573*WProd1672 - 
     -           WProd1598*WProd1672 + WProd1620*WProd1672 - 
     -           3*WProd1571*WProd1592*WProd1673 - 
     -           2*WProd1579*WProd1592*WProd1673 + 
     -           3*WProd1581*WProd1592*WProd1673 + 
     -           2*WProd1582*WProd1592*WProd1673 + 
     -           3*WProd1583*WProd1592*WProd1673 + 
     -           2*WProd1584*WProd1592*WProd1673 + 
     -           3*WProd1585*WProd1592*WProd1673 + 
     -           3*WProd1586*WProd1592*WProd1673 - 
     -           3*WProd1587*WProd1592*WProd1673 - 
     -           3*WProd1588*WProd1592*WProd1673 + 
     -           WProd1616*WProd1673 + WProd1617*WProd1673 + 
     -           WProd1618*WProd1673 - WProd1623*WProd1673 + 
     -           WProd1675 - WProd1673*WProd1692))/WProd1661) + 
     -   amp(18,6,6)

#ifdef DEBUG
	print *, 'amp(18,6,6) =', amp(18,6,6)
#endif

        amp(18,6,4) = Prp175*
     -   ((-16*Alfas*Cpl122*Prp178*Prp179*
     -        (-2*WProd1626*WProd1651 - 
     -          2*WProd1592*WProd1632*WProd1651 - 
     -          2*WProd1592*WProd1633*WProd1651 + 
     -          2*WProd1635*WProd1651 + 2*WProd1636*WProd1651 + 
     -          2*WProd1637*WProd1651 - 2*WProd1627*WProd1653 + 
     -          2*WProd1627*WProd1654 - 2*WProd1626*WProd1662 + 
     -          2*WProd1592*WProd1627*WProd1662 - 
     -          2*WProd1592*WProd1632*WProd1662 - 
     -          2*WProd1592*WProd1633*WProd1662 + 
     -          2*WProd1635*WProd1662 + 2*WProd1636*WProd1662 + 
     -          2*WProd1637*WProd1662 - 2*WProd1626*WProd1663 + 
     -          2*WProd1592*WProd1627*WProd1663 - 
     -          2*WProd1592*WProd1632*WProd1663 - 
     -          2*WProd1592*WProd1633*WProd1663 + 
     -          2*WProd1635*WProd1663 + 2*WProd1636*WProd1663 + 
     -          2*WProd1637*WProd1663 - 2*WProd1627*WProd1665 - 
     -          2*WProd1627*WProd1666 + 2*WProd1626*WProd1672 + 
     -          2*WProd1592*WProd1632*WProd1672 + 
     -          2*WProd1592*WProd1633*WProd1672 - 
     -          2*WProd1635*WProd1672 - 2*WProd1636*WProd1672 - 
     -          2*WProd1637*WProd1672 - 
     -          3*WProd1571*WProd1592*WProd1690 - 
     -          4*WProd1579*WProd1592*WProd1690 + 
     -          4*WProd1581*WProd1592*WProd1690 + 
     -          2*WProd1582*WProd1592*WProd1690 + 
     -          4*WProd1583*WProd1592*WProd1690 + 
     -          2*WProd1584*WProd1592*WProd1690 + 
     -          4*WProd1585*WProd1592*WProd1690 + 
     -          4*WProd1586*WProd1592*WProd1690 - 
     -          3*WProd1587*WProd1592*WProd1690 - 
     -          3*WProd1588*WProd1592*WProd1690 + 
     -          2*WProd1616*WProd1690 + 2*WProd1617*WProd1690 + 
     -          2*WProd1618*WProd1690 - 2*WProd1619*WProd1690))/
     -      WProd1679 + Prp180*
     -      ((-16*Alfas*Cpl122*Prp179*
     -           (2*WProd1637*WProd1651 - 2*WProd1646*WProd1651 - 
     -             2*WProd1627*WProd1653 + 2*WProd1637*WProd1662 - 
     -             2*WProd1646*WProd1662 + 2*WProd1637*WProd1663 - 
     -             2*WProd1646*WProd1663 - 2*WProd1637*WProd1672 + 
     -             2*WProd1646*WProd1672 + 2*WProd1627*WProd1674 - 
     -             2*WProd1571*WProd1592*WProd1690 - 
     -             WProd1579*WProd1592*WProd1690 + 
     -             2*WProd1581*WProd1592*WProd1690 + 
     -             WProd1582*WProd1592*WProd1690 + 
     -             2*WProd1583*WProd1592*WProd1690 + 
     -             WProd1584*WProd1592*WProd1690 + 
     -             WProd1585*WProd1592*WProd1690 + 
     -             WProd1586*WProd1592*WProd1690 - 
     -             WProd1587*WProd1592*WProd1690 - 
     -             WProd1588*WProd1592*WProd1690 + 
     -             2*WProd1617*WProd1690 - 2*WProd1623*WProd1690))/
     -         WProd1679 - 
     -        (32*Alfas*Cpl122*Prp182*
     -           (WProd1637*WProd1662 - WProd1646*WProd1662 + 
     -             WProd1637*WProd1663 - WProd1646*WProd1663 - 
     -             WProd1637*WProd1672 + WProd1646*WProd1672 + 
     -             WProd1571*WProd1680 - WProd1581*WProd1680 - 
     -             WProd1583*WProd1680 - WProd1682 - 
     -             WProd1571*WProd1592*WProd1690 + 
     -             WProd1581*WProd1592*WProd1690 + 
     -             WProd1583*WProd1592*WProd1690 + 
     -             WProd1617*WProd1690 - WProd1623*WProd1690 + 
     -             WProd1691))/WProd1679))
        amp(18,6,4) = Prp175*Prp181*
     -    (Prp178*((-32*Alfas*Cpl122*(WProd1677 + WProd1678))/
     -          WProd1676 + 
     -         (32*Alfas*Cpl122*
     -            (-(WProd1632*WProd1652) - WProd1633*WProd1652 + 
     -              WProd1632*WProd1659 + WProd1633*WProd1659 + 
     -              WProd1625*WProd1662 - WProd1639*WProd1662 - 
     -              WProd1644*WProd1662 + WProd1625*WProd1663 - 
     -              WProd1639*WProd1663 - WProd1644*WProd1663 + 
     -              WProd1632*WProd1671 + WProd1633*WProd1671 - 
     -              WProd1626*WProd1672 - 
     -              WProd1592*WProd1632*WProd1672 - 
     -              WProd1592*WProd1633*WProd1672 + 
     -              WProd1635*WProd1672 + WProd1636*WProd1672 + 
     -              WProd1637*WProd1672 - 2*WProd1571*WProd1680 - 
     -              2*WProd1579*WProd1680 + 
     -              2*WProd1581*WProd1680 + WProd1582*WProd1680 + 
     -              2*WProd1583*WProd1680 + WProd1584*WProd1680 + 
     -              2*WProd1585*WProd1680 + 
     -              2*WProd1586*WProd1680 - 
     -              2*WProd1587*WProd1680 - 
     -              2*WProd1588*WProd1680 + WProd1682 - 
     -              WProd1683 + 2*WProd1571*WProd1684 + 
     -              2*WProd1579*WProd1684 - 
     -              2*WProd1581*WProd1684 - WProd1582*WProd1684 - 
     -              2*WProd1583*WProd1684 - WProd1584*WProd1684 - 
     -              2*WProd1585*WProd1684 - 
     -              2*WProd1586*WProd1684 + 
     -              2*WProd1587*WProd1684 + 
     -              2*WProd1588*WProd1684 - 
     -              2*WProd1571*WProd1685 - 
     -              2*WProd1579*WProd1685 + 
     -              2*WProd1581*WProd1685 + WProd1582*WProd1685 + 
     -              2*WProd1583*WProd1685 + WProd1584*WProd1685 + 
     -              2*WProd1585*WProd1685 + 
     -              2*WProd1586*WProd1685 - 
     -              2*WProd1587*WProd1685 - 
     -              2*WProd1588*WProd1685 + WProd1686 + 
     -              WProd1687 + WProd1688 + WProd1689 + 
     -              2*WProd1571*WProd1592*WProd1690 + 
     -              2*WProd1579*WProd1592*WProd1690 - 
     -              2*WProd1581*WProd1592*WProd1690 - 
     -              WProd1582*WProd1592*WProd1690 - 
     -              2*WProd1583*WProd1592*WProd1690 - 
     -              WProd1584*WProd1592*WProd1690 - 
     -              2*WProd1585*WProd1592*WProd1690 - 
     -              2*WProd1586*WProd1592*WProd1690 + 
     -              2*WProd1587*WProd1592*WProd1690 + 
     -              2*WProd1588*WProd1592*WProd1690 - 
     -              WProd1616*WProd1690 - WProd1617*WProd1690 - 
     -              WProd1618*WProd1690 + WProd1619*WProd1690))/
     -          WProd1679) - 
     -      (32*Alfas*Cpl122*Prp182*
     -         (-(WProd1625*WProd1662) + WProd1639*WProd1662 - 
     -           WProd1625*WProd1663 + WProd1639*WProd1663 + 
     -           WProd1625*WProd1672 - WProd1637*WProd1672 + 
     -           WProd1646*WProd1672 + 3*WProd1571*WProd1680 + 
     -           2*WProd1579*WProd1680 - 3*WProd1581*WProd1680 - 
     -           2*WProd1582*WProd1680 - 3*WProd1583*WProd1680 - 
     -           2*WProd1584*WProd1680 - 3*WProd1585*WProd1680 - 
     -           3*WProd1586*WProd1680 + 3*WProd1587*WProd1680 + 
     -           3*WProd1588*WProd1680 + WProd1681 - WProd1682 + 
     -           3*WProd1571*WProd1685 + 2*WProd1579*WProd1685 - 
     -           3*WProd1581*WProd1685 - 2*WProd1582*WProd1685 - 
     -           3*WProd1583*WProd1685 - 2*WProd1584*WProd1685 - 
     -           3*WProd1585*WProd1685 - 3*WProd1586*WProd1685 + 
     -           3*WProd1587*WProd1685 + 3*WProd1588*WProd1685 - 
     -           WProd1686 - WProd1687 - WProd1688 - WProd1689 - 
     -           3*WProd1571*WProd1592*WProd1690 - 
     -           2*WProd1579*WProd1592*WProd1690 + 
     -           3*WProd1581*WProd1592*WProd1690 + 
     -           2*WProd1582*WProd1592*WProd1690 + 
     -           3*WProd1583*WProd1592*WProd1690 + 
     -           2*WProd1584*WProd1592*WProd1690 + 
     -           3*WProd1585*WProd1592*WProd1690 + 
     -           3*WProd1586*WProd1592*WProd1690 - 
     -           3*WProd1587*WProd1592*WProd1690 - 
     -           3*WProd1588*WProd1592*WProd1690 + 
     -           WProd1616*WProd1690 + WProd1617*WProd1690 + 
     -           WProd1618*WProd1690 - WProd1623*WProd1690 + 
     -           WProd1691 - WProd1690*WProd1692))/WProd1679) + 
     -   amp(18,6,4)

#ifdef DEBUG
	print *, 'amp(18,6,4) =', amp(18,6,4)
#endif

      endif

      end


