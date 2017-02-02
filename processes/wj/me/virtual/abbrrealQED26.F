      subroutine abbrrealQED26(resnr)
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

      
	WProd2010 = Prod(6,2,6,1)
	WProd1947 = Prod(2,2,2,1)*Prod(6,2,6,1)
	WProd1926 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1925 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd1995 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd2012 = Prod(6,2,6,1)*ProdC(2,2,2,1)
	WProd1927 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1903 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd1900 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd1986 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd1928 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1905 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd1895 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd1978 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd1929 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1910 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd1906 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd1902 = Prod(5,1,2,2)*ProdC(5,1,2,1)
        WProd1944 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd1941 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd1888 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd1912 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd1914 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd1990 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
	WProd1987 = Prod(5,1,2,1)*ProdC(5,1,2,2)
        WProd2004 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,2))
        WProd1991 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd1930 = Prod(5,1,3,1)*ProdC(5,1,3,1)
        WProd1890 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd1945 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd2008 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd1915 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
	WProd1931 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd1992 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2006 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd1891 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1946 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1933 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1938 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2009 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2000 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd1934 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1939 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2001 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,2)*
     -     ProdC(5,1,4,1))
	WProd1909 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd1940 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd1959 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd1948 = Prod(6,2,2,2)*ProdC(6,1,2,1)
        WProd1973 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd1949 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd1960 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd1950 = -(Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd1961 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd1951 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2013 = -(Prod(3,1,2,1)**2*Prod(6,2,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd1952 = -(Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd1970 = -(Prod(2,2,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd1953 = -(Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd1962 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd1971 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd1954 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2014 = -(Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd1972 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
	WProd2015 = Prod(6,2,2,1)*ProdC(6,1,2,2)
        WProd2030 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2016 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2022 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,2))
        WProd2017 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,2))
        WProd2025 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,2))
        WProd2023 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2026 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2018 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2027 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
	WProd1963 = Prod(6,2,3,1)*ProdC(6,1,3,1)
        WProd2019 = Prod(3,1,2,1)**2*Prod(6,2,1,1)*
     -   ProdC(2,1,2,2)*ProdC(5,1,4,1)*ProdC(6,1,3,1)
	WProd1964 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd1965 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,1,4,1))
        WProd1974 = -(Prod(2,2,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,1,4,1))
        WProd2031 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,2)*
     -     ProdC(6,1,4,1))
        WProd1966 = -(Prod(3,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,3,1)*
     -     ProdC(6,1,4,1))
        WProd2028 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,3,1)*ProdC(6,1,4,1))
        WProd1955 = Prod(4,1,3,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1967 = -(Prod(4,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd2029 = -(Prod(2,1,1,1)*Prod(4,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,1,4,1))
        WProd2020 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -   ProdC(2,1,2,2)*ProdC(5,1,4,1)*ProdC(6,1,4,1)
	WProd1911 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1968 = Prod(3,1,1,1)*Prod(6,2,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1956 = Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd1957 = Prod(3,1,2,2)*Prod(6,2,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd1969 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,2,1,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,5,1)
        WProd2024 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,1,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(6,1,5,1)
        WProd2021 = Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(4,1,2,2)*
     -   ProdC(6,1,5,1)
        WProd2011 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(5,1,4,1)*ProdC(6,1,5,1))
        WProd1958 = Prod(5,1,3,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
	WProd1913 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd1889 = Prod(6,1,2,2)*ProdC(6,2,2,1)
        WProd1942 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd1892 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd1916 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd1893 = -(Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd1917 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd1894 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd1977 = -(Prod(3,1,2,1)**2*Prod(6,1,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd1896 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd1935 = -(Prod(2,2,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd1897 = -(Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1918 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1936 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1898 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1979 = -(Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd1937 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
	WProd1980 = Prod(6,1,2,1)*ProdC(6,2,2,2)
        WProd2005 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd1981 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd1993 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,2))
        WProd1982 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,2))
        WProd1997 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,2))
        WProd1994 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd1998 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd1983 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd1999 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
	WProd1919 = Prod(6,1,3,1)*ProdC(6,2,3,1)
        WProd1984 = Prod(3,1,2,1)**2*Prod(6,1,1,1)*
     -   ProdC(2,1,2,2)*ProdC(5,1,4,1)*ProdC(6,2,3,1)
	WProd1920 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd1921 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,2,4,1))
        WProd1943 = -(Prod(2,2,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,2,4,1))
        WProd2007 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,2)*
     -     ProdC(6,2,4,1))
        WProd1922 = -(Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,3,1)*
     -     ProdC(6,2,4,1))
        WProd2002 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,3,1)*ProdC(6,2,4,1))
        WProd1899 = Prod(4,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,4,1)
        WProd1923 = -(Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,4,1))
        WProd2003 = -(Prod(2,1,1,1)*Prod(4,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,2,4,1))
        WProd1985 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(2,1,2,2)*ProdC(5,1,4,1)*ProdC(6,2,4,1)
        WProd1924 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,1)
        WProd1901 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,5,1)
        WProd1904 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,5,1)
        WProd1932 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,5,1)
        WProd1996 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(6,2,5,1)
        WProd1988 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(6,2,5,1)
        WProd1975 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(5,1,4,1)*ProdC(6,2,5,1))
        WProd1907 = Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,5,1)
	WProd1976 = ProdC(6,2,6,1)
	WProd1908 = Prod(2,2,2,1)*ProdC(6,2,6,1)
	WProd1989 = ProdC(2,2,2,1)*ProdC(6,2,6,1)
      
      
	Prp215 = -(1.D0/WProd1905)
	Prp220 = -(1.D0/WProd1906)
	Prp216 = -(1.D0/WProd1909)
	Prp218 = 1.D0/WProd1911
        Prp217 = 1.D0/
     -   ((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd1909 - WProd1910 + 
     -     WProd1911)
	Prp223 = 1.D0/WProd1926
	Prp219 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1929)
	Prp221 = 1.D0/(-WProd1927 - WProd1928 + WProd1929)
	Prp222 = 1.D0/(-WProd1909 + WProd1926 - WProd1940)
      
      
      
      endif

      end


