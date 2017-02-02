      subroutine abbrrealQED27(resnr)
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

      
	WProd2168 = Prod(6,2,6,1)
	WProd2095 = Prod(2,2,2,1)*Prod(6,2,6,1)
	WProd2050 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd2047 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd2140 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd2160 = Prod(6,2,6,1)*ProdC(2,2,2,1)
	WProd2052 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd2070 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd2049 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd2141 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd2053 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd2054 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd2083 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd2152 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd2055 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd2056 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd2057 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd2037 = Prod(5,1,2,2)*ProdC(5,1,2,1)
        WProd2079 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd2067 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd2081 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd2042 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd2044 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd2137 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
	WProd2131 = Prod(5,1,2,1)*ProdC(5,1,2,2)
        WProd2144 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,2))
        WProd2138 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd2058 = Prod(5,1,3,1)*ProdC(5,1,3,1)
        WProd2084 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd2093 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd2158 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd2045 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
	WProd2059 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd2139 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd2085 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2094 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2086 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2072 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2159 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2148 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2082 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2069 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd2145 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,2,2)*
     -     ProdC(5,1,4,1))
	WProd2076 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd2096 = Prod(2,1,1,1)*Prod(2,2,1,1)*Prod(6,2,3,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,1,1)
        WProd2161 = Prod(2,1,1,1)**2*Prod(6,2,3,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,2)*ProdC(6,1,1,1)
        WProd2169 = Prod(2,1,1,1)**2*Prod(6,2,3,1)*
     -   ProdC(5,1,4,1)*ProdC(6,1,1,1)
	WProd2077 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd2114 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd2119 = Prod(6,2,2,2)*ProdC(6,1,2,1)
        WProd2097 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,2,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd2098 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd2099 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd2100 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,1))
        WProd2120 = -(Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2121 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2177 = -(Prod(3,1,2,1)**2*Prod(6,2,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,3,1)*ProdC(6,1,2,1))
        WProd2122 = -(Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd2101 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd2116 = -(Prod(2,2,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,1))
        WProd2123 = -(Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2117 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2124 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2178 = -(Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2102 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2103 = -(Prod(2,2,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
        WProd2104 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,2,5,1)*
     -   ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1)
        WProd2162 = -(Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,2,2)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
	WProd2176 = Prod(6,2,2,1)*ProdC(6,1,2,2)
        WProd2163 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2164 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2172 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,1,2,2))
        WProd2179 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,1,2,2))
        WProd2174 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,2,2))
        WProd2175 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2180 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2165 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
        WProd2166 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2)
	WProd2105 = Prod(6,2,3,1)*ProdC(6,1,3,1)
	WProd2066 = Prod(6,1,4,1)*ProdC(6,1,4,1)
	WProd2106 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd2107 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,1,4,1))
        WProd2108 = -(Prod(2,2,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,1,4,1))
        WProd2167 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,2,2)*
     -     ProdC(6,1,4,1))
        WProd2109 = -(Prod(3,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,3,1)*
     -     ProdC(6,1,4,1))
        WProd2170 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,3,1)*ProdC(6,1,4,1))
        WProd2118 = Prod(4,1,3,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd2115 = -(Prod(4,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd2173 = -(Prod(2,1,1,1)*Prod(4,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,1,4,1))
	WProd2078 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd2110 = Prod(3,1,1,1)*Prod(6,2,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd2111 = Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd2125 = Prod(3,1,2,2)*Prod(6,2,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd2181 = Prod(3,1,2,1)*Prod(6,2,1,1)*ProdC(4,1,2,2)*
     -   ProdC(6,1,5,1)
        WProd2112 = Prod(5,1,3,1)*Prod(6,2,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd2113 = -(Prod(5,1,1,1)*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,5,1))
        WProd2171 = -(Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,1,5,1))
        WProd2032 = Prod(2,1,1,1)*Prod(2,2,1,1)*Prod(6,1,3,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,1,1)
        WProd2130 = Prod(2,1,1,1)**2*Prod(6,1,3,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,2)*ProdC(6,2,1,1)
        WProd2126 = Prod(2,1,1,1)**2*Prod(6,1,3,1)*
     -   ProdC(5,1,4,1)*ProdC(6,2,1,1)
	WProd2068 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd2080 = Prod(6,1,2,2)*ProdC(6,2,2,1)
        WProd2033 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd2034 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd2035 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd2036 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,1))
        WProd2087 = -(Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd2088 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd2153 = -(Prod(3,1,2,1)**2*Prod(6,1,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,3,1)*ProdC(6,2,2,1))
        WProd2089 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd2038 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd2073 = -(Prod(2,2,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,1))
        WProd2090 = -(Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2074 = -(Prod(2,2,1,1)*Prod(4,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2091 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2154 = -(Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,2,2)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2039 = -(Prod(2,2,1,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2040 = -(Prod(2,2,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
        WProd2041 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,5,1)*
     -   ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1)
        WProd2132 = -(Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,2,2)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
	WProd2151 = Prod(6,1,2,1)*ProdC(6,2,2,2)
        WProd2133 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd2134 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd2146 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,2,1)*ProdC(6,2,2,2))
        WProd2155 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,3,1)*ProdC(6,2,2,2))
        WProd2149 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,2,2))
        WProd2150 = -(Prod(2,1,1,1)*Prod(4,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd2156 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd2135 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
        WProd2136 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(5,1,2,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2)
	WProd2043 = Prod(6,1,3,1)*ProdC(6,2,3,1)
	WProd2046 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd2048 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,2,4,1))
        WProd2051 = -(Prod(2,2,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -     ProdC(6,2,4,1))
        WProd2142 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,2)*
     -     ProdC(6,2,4,1))
        WProd2060 = -(Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,3,1)*
     -     ProdC(6,2,4,1))
        WProd2127 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,3,1)*ProdC(6,2,4,1))
        WProd2075 = Prod(4,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,4,1)
        WProd2071 = -(Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,4,1))
        WProd2147 = -(Prod(2,1,1,1)*Prod(4,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,2,4,1))
        WProd2061 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,1)
        WProd2062 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,5,1)
        WProd2092 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,5,1)
        WProd2157 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(6,2,5,1)
        WProd2063 = Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,5,1)
        WProd2064 = -(Prod(5,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,2,5,1))
        WProd2128 = -(Prod(2,1,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(5,1,4,1)*ProdC(6,2,5,1))
	WProd2129 = ProdC(6,2,6,1)
	WProd2065 = Prod(2,2,2,1)*ProdC(6,2,6,1)
	WProd2143 = ProdC(2,2,2,1)*ProdC(6,2,6,1)
      
      
	Prp227 = 1.D0/WProd2050
	Prp228 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd2055)
	Prp232 = 1.D0/(-WProd2052 - WProd2053 + WProd2055)
        Prp224 = 1.D0/
     -   ((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd2050 - WProd2056 - 
     -     WProd2057)
	Prp225 = -(1.D0/WProd2057)
	Prp226 = 1.D0/WProd2066
	Prp231 = -(1.D0/WProd2076)
	Prp229 = 1.D0/(WProd2050 - WProd2076 - WProd2077)
	Prp230 = 1.D0/WProd2078
      
      
      
      endif

      end


