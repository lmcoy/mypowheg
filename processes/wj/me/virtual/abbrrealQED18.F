      subroutine abbrrealQED18(resnr)
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

      
	WProd1676 = Prod(6,2,6,1)
	WProd1624 = Prod(2,2,2,1)*Prod(6,2,6,1)
	WProd1579 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1614 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd1672 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd1679 = Prod(6,2,6,1)*ProdC(2,2,2,1)
	WProd1581 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1582 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd1601 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd1662 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd1583 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1592 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd1623 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1621 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1659 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
        WProd1674 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,1,1))
	WProd1584 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd1602 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd1692 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1611 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1575 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1652 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
        WProd1616 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1604 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1617 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1577 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1618 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1605 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
	WProd1663 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd1671 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1665 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1653 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1666 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd1571 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1585 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd1586 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd1572 = Prod(5,1,2,2)*ProdC(5,1,2,1)
	WProd1651 = Prod(5,1,2,1)*ProdC(5,1,2,2)
	WProd1587 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd1588 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd1619 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1578 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1654 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(2,2,1,1)*
     -     ProdC(5,1,4,1))
	WProd1600 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1646 = Prod(3,1,1,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd1625 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1638 = Prod(3,1,2,2)*Prod(6,2,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1685 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,2,2)*
     -   ProdC(6,1,1,1)
        WProd1626 = Prod(5,1,3,1)*Prod(6,2,5,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
	WProd1613 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd1627 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd1645 = Prod(6,2,2,2)*ProdC(6,1,2,1)
        WProd1628 = -(Prod(5,1,2,2)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,1))
        WProd1639 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd1647 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,1))
        WProd1629 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1640 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1641 = Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,2,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1)
        WProd1630 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1642 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1643 = Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1)
        WProd1686 = Prod(3,1,2,1)**2*Prod(6,2,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,2)*ProdC(6,1,2,1)
        WProd1687 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,2)*ProdC(6,1,2,1)
        WProd1631 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
	WProd1690 = Prod(6,2,2,1)*ProdC(6,1,2,2)
        WProd1680 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,2))
        WProd1691 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,2))
        WProd1681 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1688 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1682 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1689 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1683 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
	WProd1632 = Prod(6,2,3,1)*ProdC(6,1,3,1)
	WProd1633 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd1644 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1634 = -(Prod(5,1,2,2)*Prod(6,2,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1684 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(2,2,1,1)*
     -     ProdC(6,1,4,1))
        WProd1677 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(6,1,4,1))
        WProd1635 = -(Prod(5,1,3,1)*Prod(6,2,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1678 = -(Prod(4,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,1,4,1))
        WProd1636 = -(Prod(5,1,4,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1637 = -(Prod(4,1,3,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
	WProd1693 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1620 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,1,1)
        WProd1573 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1603 = Prod(3,1,2,2)*Prod(6,1,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1664 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,2,2)*
     -   ProdC(6,2,1,1)
        WProd1574 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
	WProd1576 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd1615 = Prod(6,1,2,2)*ProdC(6,2,2,1)
        WProd1580 = -(Prod(5,1,2,2)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,1))
        WProd1606 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,2,1)
        WProd1622 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,1))
        WProd1589 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1607 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1608 = Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,1,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1)
        WProd1590 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1609 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1610 = Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1)
        WProd1667 = Prod(3,1,2,1)**2*Prod(6,1,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,2)*ProdC(6,2,2,1)
        WProd1668 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,2)*ProdC(6,2,2,1)
        WProd1591 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
	WProd1673 = Prod(6,1,2,1)*ProdC(6,2,2,2)
        WProd1655 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,2))
        WProd1675 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,2))
        WProd1656 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1669 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1657 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1670 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1658 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
	WProd1593 = Prod(6,1,3,1)*ProdC(6,2,3,1)
	WProd1594 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd1612 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1595 = -(Prod(5,1,2,2)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1660 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(2,2,1,1)*
     -     ProdC(6,2,4,1))
        WProd1648 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(6,2,4,1))
        WProd1596 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1649 = -(Prod(4,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,2,4,1))
        WProd1597 = -(Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1598 = -(Prod(4,1,3,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
	WProd1650 = ProdC(6,2,6,1)
	WProd1599 = Prod(2,2,2,1)*ProdC(6,2,6,1)
	WProd1661 = ProdC(2,2,2,1)*ProdC(6,2,6,1)
      
      
	Prp175 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1571)
	Prp181 = 1.D0/WProd1579
	Prp180 = 1.D0/(WProd1571 - WProd1581 - WProd1583)
	Prp176 = -(1.D0/WProd1586)
	Prp177 = -(1.D0/WProd1600)
	Prp178 = 1.D0/(WProd1579 - WProd1600 - WProd1613)
	Prp179 = -(1.D0/WProd1613)
	Prp182 = 1.D0/WProd1693
      
      
      
      endif

      end


