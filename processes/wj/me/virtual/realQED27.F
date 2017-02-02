      subroutine realQED27(amp,resnr)
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

        tmp1156 = (4*Cpl2*GS*Prp224*Prp225*
     -      (-2*WProd2034 - 2*WProd2035 + 2*WProd2039 + 
     -        2*WProd2042*WProd2043 - 2*WProd2043*WProd2044 - 
     -        2*WProd2037*WProd2043*WProd2045 + 
     -        2*WProd2042*WProd2046 - 2*WProd2044*WProd2046 - 
     -        2*WProd2037*WProd2045*WProd2046 + 
     -        2*WProd2048*WProd2049 - 3*WProd2050*WProd2051 + 
     -        3*WProd2051*WProd2052 + 3*WProd2051*WProd2053 + 
     -        WProd2051*WProd2054 - 3*WProd2051*WProd2055 + 
     -        3*WProd2051*WProd2056 + 2*WProd2051*WProd2057 - 
     -        3*WProd2051*WProd2058 - 3*WProd2051*WProd2059 + 
     -        2*WProd2037*WProd2060 + 2*WProd2037*WProd2061 + 
     -        2*WProd2037*WProd2062 - 2*WProd2037*WProd2063 + 
     -        WProd2051*WProd2070 + 2*WProd2037*WProd2071 + 
     -        2*WProd2079*WProd2080 + 2*WProd2080*WProd2081 - 
     -        2*WProd2080*WProd2082 + 2*WProd2048*WProd2083))/
     -    WProd2065 + 2.6666666666666665D0*
     -    (Cpl2*GS*Prp229*Prp231*
     -       (WProd2039 - WProd2037*WProd2063 - 
     -         WProd2049*WProd2063 - 2*WProd2050*WProd2073 + 
     -         2*WProd2052*WProd2073 + 2*WProd2053*WProd2073 + 
     -         2*WProd2054*WProd2073 - 3*WProd2055*WProd2073 + 
     -         2*WProd2056*WProd2073 + 2*WProd2057*WProd2073 - 
     -         3*WProd2058*WProd2073 - 3*WProd2059*WProd2073 + 
     -         2*WProd2070*WProd2073 + WProd2074 - 
     -         WProd2037*WProd2075 - WProd2049*WProd2075 + 
     -         2*WProd2045*WProd2050*WProd2080 - 
     -         2*WProd2045*WProd2052*WProd2080 - 
     -         2*WProd2045*WProd2053*WProd2080 - 
     -         2*WProd2045*WProd2054*WProd2080 + 
     -         3*WProd2045*WProd2055*WProd2080 - 
     -         2*WProd2045*WProd2056*WProd2080 - 
     -         2*WProd2045*WProd2057*WProd2080 + 
     -         3*WProd2045*WProd2058*WProd2080 + 
     -         3*WProd2045*WProd2059*WProd2080 - 
     -         2*WProd2045*WProd2070*WProd2080 - 
     -         WProd2080*WProd2082 - WProd2063*WProd2083 - 
     -         WProd2075*WProd2083 - WProd2080*WProd2086))/
     -     WProd2065 + 2.6666666666666665D0*
     -    (Cpl2*GS*Prp225*Prp231*
     -       (-WProd2035 + WProd2039 + WProd2037*WProd2062 + 
     -         WProd2049*WProd2062 - WProd2037*WProd2063 - 
     -         2*WProd2050*WProd2073 + 2*WProd2052*WProd2073 + 
     -         2*WProd2053*WProd2073 + WProd2054*WProd2073 - 
     -         2*WProd2055*WProd2073 + 2*WProd2056*WProd2073 + 
     -         WProd2057*WProd2073 - 2*WProd2058*WProd2073 - 
     -         2*WProd2059*WProd2073 + WProd2070*WProd2073 + 
     -         WProd2074 - WProd2037*WProd2075 + 
     -         2*WProd2045*WProd2050*WProd2080 - 
     -         2*WProd2045*WProd2052*WProd2080 - 
     -         2*WProd2045*WProd2053*WProd2080 - 
     -         WProd2045*WProd2054*WProd2080 + 
     -         2*WProd2045*WProd2055*WProd2080 - 
     -         2*WProd2045*WProd2056*WProd2080 - 
     -         WProd2045*WProd2057*WProd2080 + 
     -         2*WProd2045*WProd2058*WProd2080 + 
     -         2*WProd2045*WProd2059*WProd2080 - 
     -         WProd2045*WProd2070*WProd2080 + 
     -         WProd2080*WProd2081 - WProd2080*WProd2082 + 
     -         WProd2062*WProd2083 - WProd2080*WProd2084 - 
     -         WProd2080*WProd2085 - WProd2080*WProd2086 + 
     -         WProd2087 + WProd2088 - WProd2049*WProd2089 - 
     -         WProd2083*WProd2089 + WProd2090 + WProd2091 - 
     -         2*WProd2050*WProd2092 + 2*WProd2052*WProd2092 + 
     -         2*WProd2053*WProd2092 + WProd2054*WProd2092 - 
     -         2*WProd2055*WProd2092 + 2*WProd2056*WProd2092 + 
     -         WProd2057*WProd2092 - 2*WProd2058*WProd2092 - 
     -         2*WProd2059*WProd2092 + WProd2070*WProd2092))/
     -     WProd2065
        tmp1156 = tmp1156 + 
     -   Prp232*(5.333333333333333D0*
     -       (Cpl2*GS*Prp225*
     -          (-WProd2034 + WProd2042*WProd2043 - 
     -            WProd2043*WProd2044 - 
     -            WProd2037*WProd2043*WProd2045 + 
     -            WProd2042*WProd2046 - WProd2044*WProd2046 - 
     -            WProd2037*WProd2045*WProd2046 + 
     -            WProd2048*WProd2049 - WProd2050*WProd2051 + 
     -            WProd2051*WProd2052 + WProd2051*WProd2053 - 
     -            WProd2051*WProd2055 + WProd2051*WProd2056 + 
     -            WProd2051*WProd2057 - WProd2051*WProd2058 - 
     -            WProd2051*WProd2059 + WProd2037*WProd2060 + 
     -            WProd2037*WProd2061 - WProd2049*WProd2062 + 
     -            WProd2037*WProd2071 + WProd2050*WProd2073 - 
     -            WProd2052*WProd2073 - WProd2053*WProd2073 + 
     -            WProd2055*WProd2073 - WProd2056*WProd2073 - 
     -            WProd2057*WProd2073 + WProd2058*WProd2073 + 
     -            WProd2059*WProd2073 - WProd2074 + 
     -            WProd2037*WProd2075 - 
     -            WProd2045*WProd2050*WProd2080 + 
     -            WProd2045*WProd2052*WProd2080 + 
     -            WProd2045*WProd2053*WProd2080 - 
     -            WProd2045*WProd2055*WProd2080 + 
     -            WProd2045*WProd2056*WProd2080 + 
     -            WProd2045*WProd2057*WProd2080 - 
     -            WProd2045*WProd2058*WProd2080 - 
     -            WProd2045*WProd2059*WProd2080 + 
     -            WProd2079*WProd2080 + WProd2048*WProd2083 - 
     -            WProd2062*WProd2083 + WProd2080*WProd2084 + 
     -            WProd2080*WProd2085 + WProd2080*WProd2086 - 
     -            WProd2087 - WProd2088 + WProd2049*WProd2089 + 
     -            WProd2083*WProd2089 - WProd2090 - WProd2091 + 
     -            WProd2050*WProd2092 - WProd2052*WProd2092 - 
     -            WProd2053*WProd2092 + WProd2055*WProd2092 - 
     -            WProd2056*WProd2092 - WProd2057*WProd2092 + 
     -            WProd2058*WProd2092 + WProd2059*WProd2092))/
     -        WProd2065 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp230*
     -          (-WProd2034 + WProd2042*WProd2043 - 
     -            WProd2043*WProd2044 + WProd2042*WProd2046 - 
     -            WProd2044*WProd2046 - 
     -            WProd2043*WProd2045*WProd2047 - 
     -            WProd2045*WProd2046*WProd2047 + 
     -            WProd2043*WProd2045*WProd2049 + 
     -            WProd2045*WProd2046*WProd2049 + 
     -            WProd2048*WProd2049 - 2*WProd2050*WProd2051 + 
     -            WProd2051*WProd2052 + WProd2051*WProd2053 + 
     -            WProd2051*WProd2054 - WProd2051*WProd2055 + 
     -            2*WProd2051*WProd2056 + 2*WProd2051*WProd2057 - 
     -            2*WProd2051*WProd2058 - 2*WProd2051*WProd2059 + 
     -            WProd2047*WProd2060 - WProd2049*WProd2060 + 
     -            WProd2047*WProd2061 - WProd2049*WProd2061 - 
     -            WProd2049*WProd2062 - 
     -            WProd2045*WProd2049*WProd2068 + 
     -            WProd2067*WProd2068 + WProd2051*WProd2070 + 
     -            WProd2047*WProd2071 - WProd2049*WProd2071 + 
     -            WProd2068*WProd2072 + 2*WProd2050*WProd2073 - 
     -            WProd2052*WProd2073 - WProd2053*WProd2073 - 
     -            WProd2054*WProd2073 + WProd2055*WProd2073 - 
     -            2*WProd2056*WProd2073 - 2*WProd2057*WProd2073 + 
     -            2*WProd2058*WProd2073 + 2*WProd2059*WProd2073 - 
     -            WProd2070*WProd2073 - WProd2074 + 
     -            WProd2047*WProd2075 - WProd2049*WProd2075 + 
     -            WProd2043*WProd2045*WProd2083 + 
     -            WProd2045*WProd2046*WProd2083 + 
     -            WProd2048*WProd2083 - WProd2060*WProd2083 - 
     -            WProd2061*WProd2083 - WProd2062*WProd2083 - 
     -            WProd2045*WProd2068*WProd2083 - 
     -            WProd2071*WProd2083 - WProd2075*WProd2083 - 
     -            WProd2087 - WProd2088 + WProd2049*WProd2089 + 
     -            WProd2083*WProd2089 - WProd2090 - WProd2091 + 
     -            2*WProd2050*WProd2092 - WProd2052*WProd2092 - 
     -            WProd2053*WProd2092 - WProd2054*WProd2092 + 
     -            WProd2055*WProd2092 - 2*WProd2056*WProd2092 - 
     -            2*WProd2057*WProd2092 + 2*WProd2058*WProd2092 + 
     -            2*WProd2059*WProd2092 - WProd2070*WProd2092 + 
     -            WProd2068*WProd2093 + WProd2068*WProd2094))/
     -        WProd2065)
        amp(27,7,3) = Prp228*tmp1156 - 
     -   (8*Cpl2*GS*Prp224*Prp225*Prp226*
     -      (WProd2032 - WProd2033 + WProd2034 + WProd2035 - 
     -        WProd2036 - WProd2037*WProd2038 - WProd2039 + 
     -        WProd2040 + WProd2041 - WProd2042*WProd2043 + 
     -        WProd2043*WProd2044 + 
     -        WProd2037*WProd2043*WProd2045 - 
     -        WProd2042*WProd2046 + WProd2044*WProd2046 + 
     -        WProd2037*WProd2045*WProd2046 + 
     -        WProd2037*WProd2048 - WProd2047*WProd2048 - 
     -        WProd2048*WProd2049 + WProd2050*WProd2051 - 
     -        WProd2051*WProd2052 - WProd2051*WProd2053 - 
     -        WProd2051*WProd2054 + WProd2051*WProd2055 - 
     -        WProd2051*WProd2056 - WProd2051*WProd2057 + 
     -        WProd2051*WProd2058 + WProd2051*WProd2059 - 
     -        WProd2037*WProd2060 - WProd2037*WProd2061 - 
     -        WProd2037*WProd2062 + WProd2037*WProd2063 + 
     -        WProd2037*WProd2064))/WProd2065 + 
     -   Prp227*((8*Cpl2*GS*Prp224*Prp226*
     -         (-WProd2034 - WProd2035 + WProd2039 + 
     -           WProd2038*WProd2047 - 
     -           WProd2043*WProd2045*WProd2047 - 
     -           WProd2045*WProd2046*WProd2047 - 
     -           2*WProd2050*WProd2051 + WProd2051*WProd2052 + 
     -           WProd2051*WProd2053 + WProd2051*WProd2054 - 
     -           WProd2051*WProd2055 + WProd2051*WProd2056 + 
     -           WProd2051*WProd2057 - WProd2051*WProd2058 - 
     -           WProd2051*WProd2059 + WProd2047*WProd2060 + 
     -           WProd2047*WProd2061 + WProd2047*WProd2062 - 
     -           WProd2047*WProd2063 - WProd2047*WProd2064 + 
     -           WProd2042*WProd2068 + WProd2067*WProd2068 - 
     -           WProd2068*WProd2069 + WProd2051*WProd2070))/
     -       WProd2065 + Prp228*
     -       ((4*Cpl2*GS*Prp224*
     -            (-2*WProd2034 - 2*WProd2035 + 2*WProd2039 - 
     -              2*WProd2043*WProd2045*WProd2047 - 
     -              2*WProd2045*WProd2046*WProd2047 - 
     -              4*WProd2050*WProd2051 + 
     -              3*WProd2051*WProd2052 + 
     -              3*WProd2051*WProd2053 + 
     -              3*WProd2051*WProd2054 - 
     -              3*WProd2051*WProd2055 + 
     -              3*WProd2051*WProd2056 + 
     -              3*WProd2051*WProd2057 - 
     -              3*WProd2051*WProd2058 - 
     -              3*WProd2051*WProd2059 + 
     -              2*WProd2047*WProd2060 + 
     -              2*WProd2047*WProd2061 + 
     -              2*WProd2047*WProd2062 - 
     -              2*WProd2047*WProd2063 + 
     -              2*WProd2042*WProd2068 + 
     -              2*WProd2067*WProd2068 - 
     -              2*WProd2068*WProd2069 + 
     -              3*WProd2051*WProd2070 + 2*WProd2047*WProd2071))
     -           /WProd2065 - 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp229*
     -             (-WProd2039 + WProd2047*WProd2063 + 
     -               WProd2068*WProd2069 + WProd2068*WProd2072 + 
     -               2*WProd2050*WProd2073 - 
     -               2*WProd2052*WProd2073 - 
     -               2*WProd2053*WProd2073 - 
     -               2*WProd2054*WProd2073 + 
     -               2*WProd2055*WProd2073 - 
     -               2*WProd2056*WProd2073 - 
     -               2*WProd2057*WProd2073 + 
     -               2*WProd2058*WProd2073 + 
     -               2*WProd2059*WProd2073 - 
     -               2*WProd2070*WProd2073 - WProd2074 + 
     -               WProd2047*WProd2075))/WProd2065 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*Prp230*
     -             (-WProd2034 - WProd2035 - 
     -               WProd2043*WProd2045*WProd2047 - 
     -               WProd2045*WProd2046*WProd2047 - 
     -               2*WProd2050*WProd2051 + WProd2051*WProd2052 + 
     -               WProd2051*WProd2053 + WProd2051*WProd2054 - 
     -               WProd2051*WProd2055 + WProd2051*WProd2056 + 
     -               WProd2051*WProd2057 - WProd2051*WProd2058 - 
     -               WProd2051*WProd2059 + WProd2047*WProd2060 + 
     -               WProd2047*WProd2061 + WProd2047*WProd2062 + 
     -               WProd2042*WProd2068 + WProd2067*WProd2068 + 
     -               WProd2051*WProd2070 + WProd2047*WProd2071 + 
     -               WProd2068*WProd2072 + 2*WProd2050*WProd2073 - 
     -               WProd2052*WProd2073 - WProd2053*WProd2073 - 
     -               WProd2054*WProd2073 + WProd2055*WProd2073 - 
     -               WProd2056*WProd2073 - WProd2057*WProd2073 + 
     -               WProd2058*WProd2073 + WProd2059*WProd2073 - 
     -               WProd2070*WProd2073 - WProd2074 + 
     -               WProd2047*WProd2075))/WProd2065))

#ifdef DEBUG
	print *, 'amp(27,7,3) =', amp(27,7,3)
#endif

        tmp1157 = (4*Cpl2*GS*Prp224*Prp225*
     -      (-2*WProd2098 - 2*WProd2099 + 2*WProd2102 + 
     -        2*WProd2042*WProd2105 - 2*WProd2044*WProd2105 - 
     -        2*WProd2037*WProd2045*WProd2105 + 
     -        2*WProd2042*WProd2106 - 2*WProd2044*WProd2106 - 
     -        2*WProd2037*WProd2045*WProd2106 + 
     -        2*WProd2049*WProd2107 + 2*WProd2083*WProd2107 - 
     -        3*WProd2050*WProd2108 + 3*WProd2052*WProd2108 + 
     -        3*WProd2053*WProd2108 + WProd2054*WProd2108 - 
     -        3*WProd2055*WProd2108 + 3*WProd2056*WProd2108 + 
     -        2*WProd2057*WProd2108 - 3*WProd2058*WProd2108 - 
     -        3*WProd2059*WProd2108 + WProd2070*WProd2108 + 
     -        2*WProd2037*WProd2109 + 2*WProd2037*WProd2110 + 
     -        2*WProd2037*WProd2111 - 2*WProd2037*WProd2112 + 
     -        2*WProd2037*WProd2115 + 2*WProd2079*WProd2119 + 
     -        2*WProd2081*WProd2119 - 2*WProd2082*WProd2119))/
     -    WProd2095 + 2.6666666666666665D0*
     -    (Cpl2*GS*Prp229*Prp231*
     -       (WProd2102 - WProd2037*WProd2112 - 
     -         WProd2049*WProd2112 - WProd2083*WProd2112 - 
     -         2*WProd2050*WProd2116 + 2*WProd2052*WProd2116 + 
     -         2*WProd2053*WProd2116 + 2*WProd2054*WProd2116 - 
     -         3*WProd2055*WProd2116 + 2*WProd2056*WProd2116 + 
     -         2*WProd2057*WProd2116 - 3*WProd2058*WProd2116 - 
     -         3*WProd2059*WProd2116 + 2*WProd2070*WProd2116 + 
     -         WProd2117 - WProd2037*WProd2118 - 
     -         WProd2049*WProd2118 - WProd2083*WProd2118 + 
     -         2*WProd2045*WProd2050*WProd2119 - 
     -         2*WProd2045*WProd2052*WProd2119 - 
     -         2*WProd2045*WProd2053*WProd2119 - 
     -         2*WProd2045*WProd2054*WProd2119 + 
     -         3*WProd2045*WProd2055*WProd2119 - 
     -         2*WProd2045*WProd2056*WProd2119 - 
     -         2*WProd2045*WProd2057*WProd2119 + 
     -         3*WProd2045*WProd2058*WProd2119 + 
     -         3*WProd2045*WProd2059*WProd2119 - 
     -         2*WProd2045*WProd2070*WProd2119 - 
     -         WProd2082*WProd2119 - WProd2086*WProd2119))/
     -     WProd2095 + 2.6666666666666665D0*
     -    (Cpl2*GS*Prp225*Prp231*
     -       (-WProd2099 + WProd2102 + WProd2037*WProd2111 + 
     -         WProd2049*WProd2111 + WProd2083*WProd2111 - 
     -         WProd2037*WProd2112 - 2*WProd2050*WProd2116 + 
     -         2*WProd2052*WProd2116 + 2*WProd2053*WProd2116 + 
     -         WProd2054*WProd2116 - 2*WProd2055*WProd2116 + 
     -         2*WProd2056*WProd2116 + WProd2057*WProd2116 - 
     -         2*WProd2058*WProd2116 - 2*WProd2059*WProd2116 + 
     -         WProd2070*WProd2116 + WProd2117 - 
     -         WProd2037*WProd2118 + 
     -         2*WProd2045*WProd2050*WProd2119 - 
     -         2*WProd2045*WProd2052*WProd2119 - 
     -         2*WProd2045*WProd2053*WProd2119 - 
     -         WProd2045*WProd2054*WProd2119 + 
     -         2*WProd2045*WProd2055*WProd2119 - 
     -         2*WProd2045*WProd2056*WProd2119 - 
     -         WProd2045*WProd2057*WProd2119 + 
     -         2*WProd2045*WProd2058*WProd2119 + 
     -         2*WProd2045*WProd2059*WProd2119 - 
     -         WProd2045*WProd2070*WProd2119 + 
     -         WProd2081*WProd2119 - WProd2082*WProd2119 - 
     -         WProd2084*WProd2119 - WProd2085*WProd2119 - 
     -         WProd2086*WProd2119 + WProd2120 + WProd2121 - 
     -         WProd2049*WProd2122 - WProd2083*WProd2122 + 
     -         WProd2123 + WProd2124 - 2*WProd2050*WProd2125 + 
     -         2*WProd2052*WProd2125 + 2*WProd2053*WProd2125 + 
     -         WProd2054*WProd2125 - 2*WProd2055*WProd2125 + 
     -         2*WProd2056*WProd2125 + WProd2057*WProd2125 - 
     -         2*WProd2058*WProd2125 - 2*WProd2059*WProd2125 + 
     -         WProd2070*WProd2125))/WProd2095
        tmp1157 = tmp1157 + 
     -   Prp232*(5.333333333333333D0*
     -       (Cpl2*GS*Prp225*
     -          (-WProd2098 + WProd2042*WProd2105 - 
     -            WProd2044*WProd2105 - 
     -            WProd2037*WProd2045*WProd2105 + 
     -            WProd2042*WProd2106 - WProd2044*WProd2106 - 
     -            WProd2037*WProd2045*WProd2106 + 
     -            WProd2049*WProd2107 + WProd2083*WProd2107 - 
     -            WProd2050*WProd2108 + WProd2052*WProd2108 + 
     -            WProd2053*WProd2108 - WProd2055*WProd2108 + 
     -            WProd2056*WProd2108 + WProd2057*WProd2108 - 
     -            WProd2058*WProd2108 - WProd2059*WProd2108 + 
     -            WProd2037*WProd2109 + WProd2037*WProd2110 - 
     -            WProd2049*WProd2111 - WProd2083*WProd2111 + 
     -            WProd2037*WProd2115 + WProd2050*WProd2116 - 
     -            WProd2052*WProd2116 - WProd2053*WProd2116 + 
     -            WProd2055*WProd2116 - WProd2056*WProd2116 - 
     -            WProd2057*WProd2116 + WProd2058*WProd2116 + 
     -            WProd2059*WProd2116 - WProd2117 + 
     -            WProd2037*WProd2118 - 
     -            WProd2045*WProd2050*WProd2119 + 
     -            WProd2045*WProd2052*WProd2119 + 
     -            WProd2045*WProd2053*WProd2119 - 
     -            WProd2045*WProd2055*WProd2119 + 
     -            WProd2045*WProd2056*WProd2119 + 
     -            WProd2045*WProd2057*WProd2119 - 
     -            WProd2045*WProd2058*WProd2119 - 
     -            WProd2045*WProd2059*WProd2119 + 
     -            WProd2079*WProd2119 + WProd2084*WProd2119 + 
     -            WProd2085*WProd2119 + WProd2086*WProd2119 - 
     -            WProd2120 - WProd2121 + WProd2049*WProd2122 + 
     -            WProd2083*WProd2122 - WProd2123 - WProd2124 + 
     -            WProd2050*WProd2125 - WProd2052*WProd2125 - 
     -            WProd2053*WProd2125 + WProd2055*WProd2125 - 
     -            WProd2056*WProd2125 - WProd2057*WProd2125 + 
     -            WProd2058*WProd2125 + WProd2059*WProd2125))/
     -        WProd2095 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp230*
     -          (-WProd2098 + WProd2042*WProd2105 - 
     -            WProd2044*WProd2105 - 
     -            WProd2045*WProd2047*WProd2105 + 
     -            WProd2045*WProd2049*WProd2105 + 
     -            WProd2045*WProd2083*WProd2105 + 
     -            WProd2042*WProd2106 - WProd2044*WProd2106 - 
     -            WProd2045*WProd2047*WProd2106 + 
     -            WProd2045*WProd2049*WProd2106 + 
     -            WProd2045*WProd2083*WProd2106 + 
     -            WProd2049*WProd2107 + WProd2083*WProd2107 - 
     -            2*WProd2050*WProd2108 + WProd2052*WProd2108 + 
     -            WProd2053*WProd2108 + WProd2054*WProd2108 - 
     -            WProd2055*WProd2108 + 2*WProd2056*WProd2108 + 
     -            2*WProd2057*WProd2108 - 2*WProd2058*WProd2108 - 
     -            2*WProd2059*WProd2108 + WProd2070*WProd2108 + 
     -            WProd2047*WProd2109 - WProd2049*WProd2109 - 
     -            WProd2083*WProd2109 + WProd2047*WProd2110 - 
     -            WProd2049*WProd2110 - WProd2083*WProd2110 - 
     -            WProd2049*WProd2111 - WProd2083*WProd2111 - 
     -            WProd2045*WProd2049*WProd2114 + 
     -            WProd2067*WProd2114 + WProd2072*WProd2114 - 
     -            WProd2045*WProd2083*WProd2114 + 
     -            WProd2093*WProd2114 + WProd2094*WProd2114 + 
     -            WProd2047*WProd2115 - WProd2049*WProd2115 - 
     -            WProd2083*WProd2115 + 2*WProd2050*WProd2116 - 
     -            WProd2052*WProd2116 - WProd2053*WProd2116 - 
     -            WProd2054*WProd2116 + WProd2055*WProd2116 - 
     -            2*WProd2056*WProd2116 - 2*WProd2057*WProd2116 + 
     -            2*WProd2058*WProd2116 + 2*WProd2059*WProd2116 - 
     -            WProd2070*WProd2116 - WProd2117 + 
     -            WProd2047*WProd2118 - WProd2049*WProd2118 - 
     -            WProd2083*WProd2118 - WProd2120 - WProd2121 + 
     -            WProd2049*WProd2122 + WProd2083*WProd2122 - 
     -            WProd2123 - WProd2124 + 2*WProd2050*WProd2125 - 
     -            WProd2052*WProd2125 - WProd2053*WProd2125 - 
     -            WProd2054*WProd2125 + WProd2055*WProd2125 - 
     -            2*WProd2056*WProd2125 - 2*WProd2057*WProd2125 + 
     -            2*WProd2058*WProd2125 + 2*WProd2059*WProd2125 - 
     -            WProd2070*WProd2125))/WProd2095)
        amp(27,7,9) = Prp228*tmp1157 - 
     -   (8*Cpl2*GS*Prp224*Prp225*Prp226*
     -      (WProd2096 - WProd2097 + WProd2098 + WProd2099 - 
     -        WProd2100 - WProd2037*WProd2101 - WProd2102 + 
     -        WProd2103 + WProd2104 - WProd2042*WProd2105 + 
     -        WProd2044*WProd2105 + 
     -        WProd2037*WProd2045*WProd2105 - 
     -        WProd2042*WProd2106 + WProd2044*WProd2106 + 
     -        WProd2037*WProd2045*WProd2106 + 
     -        WProd2037*WProd2107 - WProd2047*WProd2107 - 
     -        WProd2049*WProd2107 + WProd2050*WProd2108 - 
     -        WProd2052*WProd2108 - WProd2053*WProd2108 - 
     -        WProd2054*WProd2108 + WProd2055*WProd2108 - 
     -        WProd2056*WProd2108 - WProd2057*WProd2108 + 
     -        WProd2058*WProd2108 + WProd2059*WProd2108 - 
     -        WProd2037*WProd2109 - WProd2037*WProd2110 - 
     -        WProd2037*WProd2111 + WProd2037*WProd2112 + 
     -        WProd2037*WProd2113))/WProd2095 + 
     -   Prp227*((8*Cpl2*GS*Prp224*Prp226*
     -         (-WProd2098 - WProd2099 + WProd2047*WProd2101 + 
     -           WProd2102 - WProd2045*WProd2047*WProd2105 - 
     -           WProd2045*WProd2047*WProd2106 - 
     -           2*WProd2050*WProd2108 + WProd2052*WProd2108 + 
     -           WProd2053*WProd2108 + WProd2054*WProd2108 - 
     -           WProd2055*WProd2108 + WProd2056*WProd2108 + 
     -           WProd2057*WProd2108 - WProd2058*WProd2108 - 
     -           WProd2059*WProd2108 + WProd2070*WProd2108 + 
     -           WProd2047*WProd2109 + WProd2047*WProd2110 + 
     -           WProd2047*WProd2111 - WProd2047*WProd2112 - 
     -           WProd2047*WProd2113 + WProd2042*WProd2114 + 
     -           WProd2067*WProd2114 - WProd2069*WProd2114))/
     -       WProd2095 + Prp228*
     -       ((4*Cpl2*GS*Prp224*
     -            (-2*WProd2098 - 2*WProd2099 + 2*WProd2102 - 
     -              2*WProd2045*WProd2047*WProd2105 - 
     -              2*WProd2045*WProd2047*WProd2106 - 
     -              4*WProd2050*WProd2108 + 
     -              3*WProd2052*WProd2108 + 
     -              3*WProd2053*WProd2108 + 
     -              3*WProd2054*WProd2108 - 
     -              3*WProd2055*WProd2108 + 
     -              3*WProd2056*WProd2108 + 
     -              3*WProd2057*WProd2108 - 
     -              3*WProd2058*WProd2108 - 
     -              3*WProd2059*WProd2108 + 
     -              3*WProd2070*WProd2108 + 
     -              2*WProd2047*WProd2109 + 
     -              2*WProd2047*WProd2110 + 
     -              2*WProd2047*WProd2111 - 
     -              2*WProd2047*WProd2112 + 
     -              2*WProd2042*WProd2114 + 
     -              2*WProd2067*WProd2114 - 
     -              2*WProd2069*WProd2114 + 2*WProd2047*WProd2115))
     -           /WProd2095 - 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp229*
     -             (-WProd2102 + WProd2047*WProd2112 + 
     -               WProd2069*WProd2114 + WProd2072*WProd2114 + 
     -               2*WProd2050*WProd2116 - 
     -               2*WProd2052*WProd2116 - 
     -               2*WProd2053*WProd2116 - 
     -               2*WProd2054*WProd2116 + 
     -               2*WProd2055*WProd2116 - 
     -               2*WProd2056*WProd2116 - 
     -               2*WProd2057*WProd2116 + 
     -               2*WProd2058*WProd2116 + 
     -               2*WProd2059*WProd2116 - 
     -               2*WProd2070*WProd2116 - WProd2117 + 
     -               WProd2047*WProd2118))/WProd2095 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*Prp230*
     -             (-WProd2098 - WProd2099 - 
     -               WProd2045*WProd2047*WProd2105 - 
     -               WProd2045*WProd2047*WProd2106 - 
     -               2*WProd2050*WProd2108 + WProd2052*WProd2108 + 
     -               WProd2053*WProd2108 + WProd2054*WProd2108 - 
     -               WProd2055*WProd2108 + WProd2056*WProd2108 + 
     -               WProd2057*WProd2108 - WProd2058*WProd2108 - 
     -               WProd2059*WProd2108 + WProd2070*WProd2108 + 
     -               WProd2047*WProd2109 + WProd2047*WProd2110 + 
     -               WProd2047*WProd2111 + WProd2042*WProd2114 + 
     -               WProd2067*WProd2114 + WProd2072*WProd2114 + 
     -               WProd2047*WProd2115 + 2*WProd2050*WProd2116 - 
     -               WProd2052*WProd2116 - WProd2053*WProd2116 - 
     -               WProd2054*WProd2116 + WProd2055*WProd2116 - 
     -               WProd2056*WProd2116 - WProd2057*WProd2116 + 
     -               WProd2058*WProd2116 + WProd2059*WProd2116 - 
     -               WProd2070*WProd2116 - WProd2117 + 
     -               WProd2047*WProd2118))/WProd2095))

#ifdef DEBUG
	print *, 'amp(27,7,9) =', amp(27,7,9)
#endif

        tmp1158 = Prp224*Prp225*
     -    ((8*Cpl2*GS*(WProd2127 + WProd2147))/WProd2129 + 
     -      (4*Cpl2*GS*(-2*WProd2043*WProd2045*WProd2131 - 
     -           2*WProd2045*WProd2046*WProd2131 + 
     -           2*WProd2060*WProd2131 + 2*WProd2061*WProd2131 + 
     -           2*WProd2062*WProd2131 - 2*WProd2063*WProd2131 + 
     -           2*WProd2071*WProd2131 - 2*WProd2134 + 
     -           2*WProd2135 - 2*WProd2043*WProd2137 - 
     -           2*WProd2046*WProd2137 + 2*WProd2043*WProd2138 + 
     -           2*WProd2046*WProd2138 - 2*WProd2043*WProd2139 - 
     -           2*WProd2046*WProd2139 + 2*WProd2048*WProd2141 - 
     -           3*WProd2050*WProd2142 + 3*WProd2052*WProd2142 + 
     -           3*WProd2053*WProd2142 + WProd2054*WProd2142 - 
     -           3*WProd2055*WProd2142 + 3*WProd2056*WProd2142 + 
     -           2*WProd2057*WProd2142 - 3*WProd2058*WProd2142 - 
     -           3*WProd2059*WProd2142 + WProd2070*WProd2142 - 
     -           2*WProd2146 + 2*WProd2079*WProd2151 + 
     -           2*WProd2081*WProd2151 - 2*WProd2082*WProd2151 + 
     -           2*WProd2048*WProd2152))/WProd2143) + 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp229*Prp231*
     -       (-(WProd2063*WProd2131) - WProd2075*WProd2131 + 
     -         WProd2135 - WProd2063*WProd2141 - 
     -         WProd2075*WProd2141 - 2*WProd2050*WProd2149 + 
     -         2*WProd2052*WProd2149 + 2*WProd2053*WProd2149 + 
     -         2*WProd2054*WProd2149 - 3*WProd2055*WProd2149 + 
     -         2*WProd2056*WProd2149 + 2*WProd2057*WProd2149 - 
     -         3*WProd2058*WProd2149 - 3*WProd2059*WProd2149 + 
     -         2*WProd2070*WProd2149 + WProd2150 + 
     -         2*WProd2045*WProd2050*WProd2151 - 
     -         2*WProd2045*WProd2052*WProd2151 - 
     -         2*WProd2045*WProd2053*WProd2151 - 
     -         2*WProd2045*WProd2054*WProd2151 + 
     -         3*WProd2045*WProd2055*WProd2151 - 
     -         2*WProd2045*WProd2056*WProd2151 - 
     -         2*WProd2045*WProd2057*WProd2151 + 
     -         3*WProd2045*WProd2058*WProd2151 + 
     -         3*WProd2045*WProd2059*WProd2151 - 
     -         2*WProd2045*WProd2070*WProd2151 - 
     -         WProd2082*WProd2151 - WProd2086*WProd2151 - 
     -         WProd2063*WProd2152 - WProd2075*WProd2152))/
     -     WProd2143 + 2.6666666666666665D0*
     -    (Cpl2*GS*Prp225*Prp231*
     -       (WProd2062*WProd2131 - WProd2063*WProd2131 - 
     -         WProd2075*WProd2131 + WProd2135 + 
     -         WProd2062*WProd2141 - WProd2089*WProd2141 - 
     -         WProd2146 - 2*WProd2050*WProd2149 + 
     -         2*WProd2052*WProd2149 + 2*WProd2053*WProd2149 + 
     -         WProd2054*WProd2149 - 2*WProd2055*WProd2149 + 
     -         2*WProd2056*WProd2149 + WProd2057*WProd2149 - 
     -         2*WProd2058*WProd2149 - 2*WProd2059*WProd2149 + 
     -         WProd2070*WProd2149 + WProd2150 + 
     -         2*WProd2045*WProd2050*WProd2151 - 
     -         2*WProd2045*WProd2052*WProd2151 - 
     -         2*WProd2045*WProd2053*WProd2151 - 
     -         WProd2045*WProd2054*WProd2151 + 
     -         2*WProd2045*WProd2055*WProd2151 - 
     -         2*WProd2045*WProd2056*WProd2151 - 
     -         WProd2045*WProd2057*WProd2151 + 
     -         2*WProd2045*WProd2058*WProd2151 + 
     -         2*WProd2045*WProd2059*WProd2151 - 
     -         WProd2045*WProd2070*WProd2151 + 
     -         WProd2081*WProd2151 - WProd2082*WProd2151 - 
     -         WProd2084*WProd2151 - WProd2085*WProd2151 - 
     -         WProd2086*WProd2151 + WProd2062*WProd2152 - 
     -         WProd2089*WProd2152 + WProd2153 + WProd2154 + 
     -         WProd2155 + WProd2156 - 2*WProd2050*WProd2157 + 
     -         2*WProd2052*WProd2157 + 2*WProd2053*WProd2157 + 
     -         WProd2054*WProd2157 - 2*WProd2055*WProd2157 + 
     -         2*WProd2056*WProd2157 + WProd2057*WProd2157 - 
     -         2*WProd2058*WProd2157 - 2*WProd2059*WProd2157 + 
     -         WProd2070*WProd2157))/WProd2143
        tmp1158 = tmp1158 + 
     -   Prp232*(Prp225*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2127 + WProd2147))/WProd2129 - 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd2043*WProd2045*WProd2131 + 
     -               WProd2045*WProd2046*WProd2131 - 
     -               WProd2060*WProd2131 - WProd2061*WProd2131 - 
     -               WProd2071*WProd2131 - WProd2075*WProd2131 + 
     -               WProd2134 + WProd2043*WProd2137 + 
     -               WProd2046*WProd2137 - WProd2043*WProd2138 - 
     -               WProd2046*WProd2138 + WProd2043*WProd2139 + 
     -               WProd2046*WProd2139 - WProd2048*WProd2141 + 
     -               WProd2062*WProd2141 - WProd2089*WProd2141 + 
     -               WProd2050*WProd2142 - WProd2052*WProd2142 - 
     -               WProd2053*WProd2142 + WProd2055*WProd2142 - 
     -               WProd2056*WProd2142 - WProd2057*WProd2142 + 
     -               WProd2058*WProd2142 + WProd2059*WProd2142 - 
     -               WProd2050*WProd2149 + WProd2052*WProd2149 + 
     -               WProd2053*WProd2149 - WProd2055*WProd2149 + 
     -               WProd2056*WProd2149 + WProd2057*WProd2149 - 
     -               WProd2058*WProd2149 - WProd2059*WProd2149 + 
     -               WProd2150 + WProd2045*WProd2050*WProd2151 - 
     -               WProd2045*WProd2052*WProd2151 - 
     -               WProd2045*WProd2053*WProd2151 + 
     -               WProd2045*WProd2055*WProd2151 - 
     -               WProd2045*WProd2056*WProd2151 - 
     -               WProd2045*WProd2057*WProd2151 + 
     -               WProd2045*WProd2058*WProd2151 + 
     -               WProd2045*WProd2059*WProd2151 - 
     -               WProd2079*WProd2151 - WProd2084*WProd2151 - 
     -               WProd2085*WProd2151 - WProd2086*WProd2151 - 
     -               WProd2048*WProd2152 + WProd2062*WProd2152 - 
     -               WProd2089*WProd2152 + WProd2153 + WProd2154 + 
     -               WProd2155 + WProd2156 - WProd2050*WProd2157 + 
     -               WProd2052*WProd2157 + WProd2053*WProd2157 - 
     -               WProd2055*WProd2157 + WProd2056*WProd2157 + 
     -               WProd2057*WProd2157 - WProd2058*WProd2157 - 
     -               WProd2059*WProd2157))/WProd2143) + 
     -      Prp230*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2127 + WProd2147))/WProd2129 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-WProd2134 - WProd2043*WProd2137 - 
     -               WProd2046*WProd2137 + WProd2043*WProd2138 + 
     -               WProd2046*WProd2138 - WProd2043*WProd2139 - 
     -               WProd2046*WProd2139 - 
     -               WProd2043*WProd2045*WProd2140 - 
     -               WProd2045*WProd2046*WProd2140 + 
     -               WProd2060*WProd2140 + WProd2061*WProd2140 + 
     -               WProd2071*WProd2140 + WProd2075*WProd2140 + 
     -               WProd2043*WProd2045*WProd2141 + 
     -               WProd2045*WProd2046*WProd2141 + 
     -               WProd2048*WProd2141 - WProd2060*WProd2141 - 
     -               WProd2061*WProd2141 - WProd2062*WProd2141 - 
     -               WProd2045*WProd2068*WProd2141 - 
     -               WProd2071*WProd2141 - WProd2075*WProd2141 + 
     -               WProd2089*WProd2141 - 2*WProd2050*WProd2142 + 
     -               WProd2052*WProd2142 + WProd2053*WProd2142 + 
     -               WProd2054*WProd2142 - WProd2055*WProd2142 + 
     -               2*WProd2056*WProd2142 + 
     -               2*WProd2057*WProd2142 - 
     -               2*WProd2058*WProd2142 - 
     -               2*WProd2059*WProd2142 + WProd2070*WProd2142 + 
     -               WProd2068*WProd2144 + WProd2068*WProd2148 + 
     -               2*WProd2050*WProd2149 - WProd2052*WProd2149 - 
     -               WProd2053*WProd2149 - WProd2054*WProd2149 + 
     -               WProd2055*WProd2149 - 2*WProd2056*WProd2149 - 
     -               2*WProd2057*WProd2149 + 
     -               2*WProd2058*WProd2149 + 
     -               2*WProd2059*WProd2149 - WProd2070*WProd2149 - 
     -               WProd2150 + WProd2043*WProd2045*WProd2152 + 
     -               WProd2045*WProd2046*WProd2152 + 
     -               WProd2048*WProd2152 - WProd2060*WProd2152 - 
     -               WProd2061*WProd2152 - WProd2062*WProd2152 - 
     -               WProd2045*WProd2068*WProd2152 - 
     -               WProd2071*WProd2152 - WProd2075*WProd2152 + 
     -               WProd2089*WProd2152 - WProd2153 - WProd2154 - 
     -               WProd2155 - WProd2156 + 
     -               2*WProd2050*WProd2157 - WProd2052*WProd2157 - 
     -               WProd2053*WProd2157 - WProd2054*WProd2157 + 
     -               WProd2055*WProd2157 - 2*WProd2056*WProd2157 - 
     -               2*WProd2057*WProd2157 + 
     -               2*WProd2058*WProd2157 + 
     -               2*WProd2059*WProd2157 - WProd2070*WProd2157 + 
     -               WProd2068*WProd2158 + WProd2068*WProd2159))/
     -           WProd2143))
        amp(27,7,10) = Prp228*tmp1158 + 
     -   Prp224*Prp225*Prp226*
     -    ((8*Cpl2*GS*(WProd2126 + WProd2127 - WProd2128))/
     -       WProd2129 - (8*Cpl2*GS*
     -         (WProd2130 - WProd2038*WProd2131 + 
     -           WProd2043*WProd2045*WProd2131 + 
     -           WProd2045*WProd2046*WProd2131 + 
     -           WProd2048*WProd2131 - WProd2060*WProd2131 - 
     -           WProd2061*WProd2131 - WProd2062*WProd2131 + 
     -           WProd2063*WProd2131 + WProd2064*WProd2131 + 
     -           WProd2132 - WProd2133 + WProd2134 - WProd2135 + 
     -           WProd2136 + WProd2043*WProd2137 + 
     -           WProd2046*WProd2137 - WProd2043*WProd2138 - 
     -           WProd2046*WProd2138 + WProd2043*WProd2139 + 
     -           WProd2046*WProd2139 - WProd2048*WProd2140 - 
     -           WProd2048*WProd2141 + WProd2050*WProd2142 - 
     -           WProd2052*WProd2142 - WProd2053*WProd2142 - 
     -           WProd2054*WProd2142 + WProd2055*WProd2142 - 
     -           WProd2056*WProd2142 - WProd2057*WProd2142 + 
     -           WProd2058*WProd2142 + WProd2059*WProd2142))/
     -       WProd2143) + Prp227*
     -    (Prp224*Prp226*((8*Cpl2*GS*
     -            (WProd2126 + WProd2127 - WProd2128))/WProd2129 + 
     -         (8*Cpl2*GS*(-WProd2134 + WProd2135 + 
     -              WProd2068*WProd2138 - WProd2043*WProd2139 - 
     -              WProd2046*WProd2139 + WProd2038*WProd2140 - 
     -              WProd2043*WProd2045*WProd2140 - 
     -              WProd2045*WProd2046*WProd2140 + 
     -              WProd2060*WProd2140 + WProd2061*WProd2140 + 
     -              WProd2062*WProd2140 - WProd2063*WProd2140 - 
     -              WProd2064*WProd2140 - 2*WProd2050*WProd2142 + 
     -              WProd2052*WProd2142 + WProd2053*WProd2142 + 
     -              WProd2054*WProd2142 - WProd2055*WProd2142 + 
     -              WProd2056*WProd2142 + WProd2057*WProd2142 - 
     -              WProd2058*WProd2142 - WProd2059*WProd2142 + 
     -              WProd2070*WProd2142 + WProd2068*WProd2144 - 
     -              WProd2068*WProd2145 - WProd2146))/WProd2143) + 
     -      Prp228*(Prp224*
     -          ((4*Cpl2*GS*
     -               (-2*WProd2134 + 2*WProd2135 + 
     -                 2*WProd2068*WProd2138 - 
     -                 2*WProd2043*WProd2139 - 
     -                 2*WProd2046*WProd2139 - 
     -                 2*WProd2043*WProd2045*WProd2140 - 
     -                 2*WProd2045*WProd2046*WProd2140 + 
     -                 2*WProd2060*WProd2140 + 
     -                 2*WProd2061*WProd2140 + 
     -                 2*WProd2062*WProd2140 - 
     -                 2*WProd2063*WProd2140 + 
     -                 2*WProd2071*WProd2140 - 
     -                 4*WProd2050*WProd2142 + 
     -                 3*WProd2052*WProd2142 + 
     -                 3*WProd2053*WProd2142 + 
     -                 3*WProd2054*WProd2142 - 
     -                 3*WProd2055*WProd2142 + 
     -                 3*WProd2056*WProd2142 + 
     -                 3*WProd2057*WProd2142 - 
     -                 3*WProd2058*WProd2142 - 
     -                 3*WProd2059*WProd2142 + 
     -                 3*WProd2070*WProd2142 + 
     -                 2*WProd2068*WProd2144 - 
     -                 2*WProd2068*WProd2145 - 2*WProd2146))/
     -             WProd2143 + 
     -            (8*Cpl2*GS*(WProd2127 + WProd2147))/WProd2129) - 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp229*
     -             (-WProd2135 + WProd2063*WProd2140 + 
     -               WProd2075*WProd2140 + WProd2068*WProd2145 + 
     -               WProd2068*WProd2148 + 2*WProd2050*WProd2149 - 
     -               2*WProd2052*WProd2149 - 
     -               2*WProd2053*WProd2149 - 
     -               2*WProd2054*WProd2149 + 
     -               2*WProd2055*WProd2149 - 
     -               2*WProd2056*WProd2149 - 
     -               2*WProd2057*WProd2149 + 
     -               2*WProd2058*WProd2149 + 
     -               2*WProd2059*WProd2149 - 
     -               2*WProd2070*WProd2149 - WProd2150))/WProd2143+
     -           Prp230*(5.333333333333333D0*
     -             (Cpl2*GS*(WProd2127 + WProd2147))/WProd2129 + 
     -            5.333333333333333D0*
     -             (Cpl2*GS*
     -                (-WProd2134 + WProd2068*WProd2138 - 
     -                  WProd2043*WProd2139 - 
     -                  WProd2046*WProd2139 - 
     -                  WProd2043*WProd2045*WProd2140 - 
     -                  WProd2045*WProd2046*WProd2140 + 
     -                  WProd2060*WProd2140 + 
     -                  WProd2061*WProd2140 + 
     -                  WProd2062*WProd2140 + 
     -                  WProd2071*WProd2140 + 
     -                  WProd2075*WProd2140 - 
     -                  2*WProd2050*WProd2142 + 
     -                  WProd2052*WProd2142 + 
     -                  WProd2053*WProd2142 + 
     -                  WProd2054*WProd2142 - 
     -                  WProd2055*WProd2142 + 
     -                  WProd2056*WProd2142 + 
     -                  WProd2057*WProd2142 - 
     -                  WProd2058*WProd2142 - 
     -                  WProd2059*WProd2142 + 
     -                  WProd2070*WProd2142 + 
     -                  WProd2068*WProd2144 - WProd2146 + 
     -                  WProd2068*WProd2148 + 
     -                  2*WProd2050*WProd2149 - 
     -                  WProd2052*WProd2149 - 
     -                  WProd2053*WProd2149 - 
     -                  WProd2054*WProd2149 + 
     -                  WProd2055*WProd2149 - 
     -                  WProd2056*WProd2149 - 
     -                  WProd2057*WProd2149 + 
     -                  WProd2058*WProd2149 + 
     -                  WProd2059*WProd2149 - 
     -                  WProd2070*WProd2149 - WProd2150))/WProd2143
     -            )))

#ifdef DEBUG
	print *, 'amp(27,7,10) =', amp(27,7,10)
#endif

        tmp1159 = Prp224*Prp225*
     -    ((8*Cpl2*GS*(WProd2170 + WProd2173))/WProd2168 + 
     -      (4*Cpl2*GS*(-2*WProd2045*WProd2105*WProd2131 - 
     -           2*WProd2045*WProd2106*WProd2131 + 
     -           2*WProd2109*WProd2131 + 2*WProd2110*WProd2131 + 
     -           2*WProd2111*WProd2131 - 2*WProd2112*WProd2131 + 
     -           2*WProd2115*WProd2131 - 2*WProd2105*WProd2137 - 
     -           2*WProd2106*WProd2137 + 2*WProd2105*WProd2138 + 
     -           2*WProd2106*WProd2138 - 2*WProd2105*WProd2139 - 
     -           2*WProd2106*WProd2139 + 2*WProd2107*WProd2141 + 
     -           2*WProd2107*WProd2152 - 2*WProd2164 + 
     -           2*WProd2165 - 3*WProd2050*WProd2167 + 
     -           3*WProd2052*WProd2167 + 3*WProd2053*WProd2167 + 
     -           WProd2054*WProd2167 - 3*WProd2055*WProd2167 + 
     -           3*WProd2056*WProd2167 + 2*WProd2057*WProd2167 - 
     -           3*WProd2058*WProd2167 - 3*WProd2059*WProd2167 + 
     -           WProd2070*WProd2167 - 2*WProd2172 + 
     -           2*WProd2079*WProd2176 + 2*WProd2081*WProd2176 - 
     -           2*WProd2082*WProd2176))/WProd2160) + 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp229*Prp231*
     -       (-(WProd2112*WProd2131) - WProd2118*WProd2131 - 
     -         WProd2112*WProd2141 - WProd2118*WProd2141 - 
     -         WProd2112*WProd2152 - WProd2118*WProd2152 + 
     -         WProd2165 - 2*WProd2050*WProd2174 + 
     -         2*WProd2052*WProd2174 + 2*WProd2053*WProd2174 + 
     -         2*WProd2054*WProd2174 - 3*WProd2055*WProd2174 + 
     -         2*WProd2056*WProd2174 + 2*WProd2057*WProd2174 - 
     -         3*WProd2058*WProd2174 - 3*WProd2059*WProd2174 + 
     -         2*WProd2070*WProd2174 + WProd2175 + 
     -         2*WProd2045*WProd2050*WProd2176 - 
     -         2*WProd2045*WProd2052*WProd2176 - 
     -         2*WProd2045*WProd2053*WProd2176 - 
     -         2*WProd2045*WProd2054*WProd2176 + 
     -         3*WProd2045*WProd2055*WProd2176 - 
     -         2*WProd2045*WProd2056*WProd2176 - 
     -         2*WProd2045*WProd2057*WProd2176 + 
     -         3*WProd2045*WProd2058*WProd2176 + 
     -         3*WProd2045*WProd2059*WProd2176 - 
     -         2*WProd2045*WProd2070*WProd2176 - 
     -         WProd2082*WProd2176 - WProd2086*WProd2176))/
     -     WProd2160 + 2.6666666666666665D0*
     -    (Cpl2*GS*Prp225*Prp231*
     -       (WProd2111*WProd2131 - WProd2112*WProd2131 - 
     -         WProd2118*WProd2131 + WProd2111*WProd2141 - 
     -         WProd2122*WProd2141 + WProd2111*WProd2152 - 
     -         WProd2122*WProd2152 + WProd2165 - WProd2172 - 
     -         2*WProd2050*WProd2174 + 2*WProd2052*WProd2174 + 
     -         2*WProd2053*WProd2174 + WProd2054*WProd2174 - 
     -         2*WProd2055*WProd2174 + 2*WProd2056*WProd2174 + 
     -         WProd2057*WProd2174 - 2*WProd2058*WProd2174 - 
     -         2*WProd2059*WProd2174 + WProd2070*WProd2174 + 
     -         WProd2175 + 2*WProd2045*WProd2050*WProd2176 - 
     -         2*WProd2045*WProd2052*WProd2176 - 
     -         2*WProd2045*WProd2053*WProd2176 - 
     -         WProd2045*WProd2054*WProd2176 + 
     -         2*WProd2045*WProd2055*WProd2176 - 
     -         2*WProd2045*WProd2056*WProd2176 - 
     -         WProd2045*WProd2057*WProd2176 + 
     -         2*WProd2045*WProd2058*WProd2176 + 
     -         2*WProd2045*WProd2059*WProd2176 - 
     -         WProd2045*WProd2070*WProd2176 + 
     -         WProd2081*WProd2176 - WProd2082*WProd2176 - 
     -         WProd2084*WProd2176 - WProd2085*WProd2176 - 
     -         WProd2086*WProd2176 + WProd2177 + WProd2178 + 
     -         WProd2179 + WProd2180 - 2*WProd2050*WProd2181 + 
     -         2*WProd2052*WProd2181 + 2*WProd2053*WProd2181 + 
     -         WProd2054*WProd2181 - 2*WProd2055*WProd2181 + 
     -         2*WProd2056*WProd2181 + WProd2057*WProd2181 - 
     -         2*WProd2058*WProd2181 - 2*WProd2059*WProd2181 + 
     -         WProd2070*WProd2181))/WProd2160
        tmp1159 = tmp1159 + 
     -   Prp232*(Prp225*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2170 + WProd2173))/WProd2168 - 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd2045*WProd2105*WProd2131 + 
     -               WProd2045*WProd2106*WProd2131 - 
     -               WProd2109*WProd2131 - WProd2110*WProd2131 - 
     -               WProd2115*WProd2131 - WProd2118*WProd2131 + 
     -               WProd2105*WProd2137 + WProd2106*WProd2137 - 
     -               WProd2105*WProd2138 - WProd2106*WProd2138 + 
     -               WProd2105*WProd2139 + WProd2106*WProd2139 - 
     -               WProd2107*WProd2141 + WProd2111*WProd2141 - 
     -               WProd2122*WProd2141 - WProd2107*WProd2152 + 
     -               WProd2111*WProd2152 - WProd2122*WProd2152 + 
     -               WProd2164 + WProd2050*WProd2167 - 
     -               WProd2052*WProd2167 - WProd2053*WProd2167 + 
     -               WProd2055*WProd2167 - WProd2056*WProd2167 - 
     -               WProd2057*WProd2167 + WProd2058*WProd2167 + 
     -               WProd2059*WProd2167 - WProd2050*WProd2174 + 
     -               WProd2052*WProd2174 + WProd2053*WProd2174 - 
     -               WProd2055*WProd2174 + WProd2056*WProd2174 + 
     -               WProd2057*WProd2174 - WProd2058*WProd2174 - 
     -               WProd2059*WProd2174 + WProd2175 + 
     -               WProd2045*WProd2050*WProd2176 - 
     -               WProd2045*WProd2052*WProd2176 - 
     -               WProd2045*WProd2053*WProd2176 + 
     -               WProd2045*WProd2055*WProd2176 - 
     -               WProd2045*WProd2056*WProd2176 - 
     -               WProd2045*WProd2057*WProd2176 + 
     -               WProd2045*WProd2058*WProd2176 + 
     -               WProd2045*WProd2059*WProd2176 - 
     -               WProd2079*WProd2176 - WProd2084*WProd2176 - 
     -               WProd2085*WProd2176 - WProd2086*WProd2176 + 
     -               WProd2177 + WProd2178 + WProd2179 + 
     -               WProd2180 - WProd2050*WProd2181 + 
     -               WProd2052*WProd2181 + WProd2053*WProd2181 - 
     -               WProd2055*WProd2181 + WProd2056*WProd2181 + 
     -               WProd2057*WProd2181 - WProd2058*WProd2181 - 
     -               WProd2059*WProd2181))/WProd2160) + 
     -      Prp230*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2170 + WProd2173))/WProd2168 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd2105*WProd2137) - 
     -               WProd2106*WProd2137 + WProd2105*WProd2138 + 
     -               WProd2106*WProd2138 - WProd2105*WProd2139 - 
     -               WProd2106*WProd2139 - 
     -               WProd2045*WProd2105*WProd2140 - 
     -               WProd2045*WProd2106*WProd2140 + 
     -               WProd2109*WProd2140 + WProd2110*WProd2140 + 
     -               WProd2115*WProd2140 + WProd2118*WProd2140 + 
     -               WProd2045*WProd2105*WProd2141 + 
     -               WProd2045*WProd2106*WProd2141 + 
     -               WProd2107*WProd2141 - WProd2109*WProd2141 - 
     -               WProd2110*WProd2141 - WProd2111*WProd2141 - 
     -               WProd2045*WProd2114*WProd2141 - 
     -               WProd2115*WProd2141 - WProd2118*WProd2141 + 
     -               WProd2122*WProd2141 + WProd2114*WProd2144 + 
     -               WProd2114*WProd2148 + 
     -               WProd2045*WProd2105*WProd2152 + 
     -               WProd2045*WProd2106*WProd2152 + 
     -               WProd2107*WProd2152 - WProd2109*WProd2152 - 
     -               WProd2110*WProd2152 - WProd2111*WProd2152 - 
     -               WProd2045*WProd2114*WProd2152 - 
     -               WProd2115*WProd2152 - WProd2118*WProd2152 + 
     -               WProd2122*WProd2152 + WProd2114*WProd2158 + 
     -               WProd2114*WProd2159 - WProd2164 - 
     -               2*WProd2050*WProd2167 + WProd2052*WProd2167 + 
     -               WProd2053*WProd2167 + WProd2054*WProd2167 - 
     -               WProd2055*WProd2167 + 2*WProd2056*WProd2167 + 
     -               2*WProd2057*WProd2167 - 
     -               2*WProd2058*WProd2167 - 
     -               2*WProd2059*WProd2167 + WProd2070*WProd2167 + 
     -               2*WProd2050*WProd2174 - WProd2052*WProd2174 - 
     -               WProd2053*WProd2174 - WProd2054*WProd2174 + 
     -               WProd2055*WProd2174 - 2*WProd2056*WProd2174 - 
     -               2*WProd2057*WProd2174 + 
     -               2*WProd2058*WProd2174 + 
     -               2*WProd2059*WProd2174 - WProd2070*WProd2174 - 
     -               WProd2175 - WProd2177 - WProd2178 - 
     -               WProd2179 - WProd2180 + 
     -               2*WProd2050*WProd2181 - WProd2052*WProd2181 - 
     -               WProd2053*WProd2181 - WProd2054*WProd2181 + 
     -               WProd2055*WProd2181 - 2*WProd2056*WProd2181 - 
     -               2*WProd2057*WProd2181 + 
     -               2*WProd2058*WProd2181 + 
     -               2*WProd2059*WProd2181 - WProd2070*WProd2181))/
     -           WProd2160))
        amp(27,7,5) = Prp228*tmp1159 + 
     -   Prp224*Prp225*Prp226*
     -    ((-8*Cpl2*GS*(-(WProd2101*WProd2131) + 
     -           WProd2045*WProd2105*WProd2131 + 
     -           WProd2045*WProd2106*WProd2131 + 
     -           WProd2107*WProd2131 - WProd2109*WProd2131 - 
     -           WProd2110*WProd2131 - WProd2111*WProd2131 + 
     -           WProd2112*WProd2131 + WProd2113*WProd2131 + 
     -           WProd2105*WProd2137 + WProd2106*WProd2137 - 
     -           WProd2105*WProd2138 - WProd2106*WProd2138 + 
     -           WProd2105*WProd2139 + WProd2106*WProd2139 - 
     -           WProd2107*WProd2140 - WProd2107*WProd2141 + 
     -           WProd2161 + WProd2162 - WProd2163 + WProd2164 - 
     -           WProd2165 + WProd2166 + WProd2050*WProd2167 - 
     -           WProd2052*WProd2167 - WProd2053*WProd2167 - 
     -           WProd2054*WProd2167 + WProd2055*WProd2167 - 
     -           WProd2056*WProd2167 - WProd2057*WProd2167 + 
     -           WProd2058*WProd2167 + WProd2059*WProd2167))/
     -       WProd2160 + (8*Cpl2*GS*
     -         (WProd2169 + WProd2170 - WProd2171))/WProd2168) + 
     -   Prp227*(Prp224*Prp226*
     -       ((8*Cpl2*GS*(WProd2169 + WProd2170 - WProd2171))/
     -          WProd2168 + 
     -         (8*Cpl2*GS*(WProd2114*WProd2138 - 
     -              WProd2105*WProd2139 - WProd2106*WProd2139 + 
     -              WProd2101*WProd2140 - 
     -              WProd2045*WProd2105*WProd2140 - 
     -              WProd2045*WProd2106*WProd2140 + 
     -              WProd2109*WProd2140 + WProd2110*WProd2140 + 
     -              WProd2111*WProd2140 - WProd2112*WProd2140 - 
     -              WProd2113*WProd2140 + WProd2114*WProd2144 - 
     -              WProd2114*WProd2145 - WProd2164 + WProd2165 - 
     -              2*WProd2050*WProd2167 + WProd2052*WProd2167 + 
     -              WProd2053*WProd2167 + WProd2054*WProd2167 - 
     -              WProd2055*WProd2167 + WProd2056*WProd2167 + 
     -              WProd2057*WProd2167 - WProd2058*WProd2167 - 
     -              WProd2059*WProd2167 + WProd2070*WProd2167 - 
     -              WProd2172))/WProd2160) + 
     -      Prp228*(Prp224*
     -          ((4*Cpl2*GS*
     -               (2*WProd2114*WProd2138 - 
     -                 2*WProd2105*WProd2139 - 
     -                 2*WProd2106*WProd2139 - 
     -                 2*WProd2045*WProd2105*WProd2140 - 
     -                 2*WProd2045*WProd2106*WProd2140 + 
     -                 2*WProd2109*WProd2140 + 
     -                 2*WProd2110*WProd2140 + 
     -                 2*WProd2111*WProd2140 - 
     -                 2*WProd2112*WProd2140 + 
     -                 2*WProd2115*WProd2140 + 
     -                 2*WProd2114*WProd2144 - 
     -                 2*WProd2114*WProd2145 - 2*WProd2164 + 
     -                 2*WProd2165 - 4*WProd2050*WProd2167 + 
     -                 3*WProd2052*WProd2167 + 
     -                 3*WProd2053*WProd2167 + 
     -                 3*WProd2054*WProd2167 - 
     -                 3*WProd2055*WProd2167 + 
     -                 3*WProd2056*WProd2167 + 
     -                 3*WProd2057*WProd2167 - 
     -                 3*WProd2058*WProd2167 - 
     -                 3*WProd2059*WProd2167 + 
     -                 3*WProd2070*WProd2167 - 2*WProd2172))/
     -             WProd2160 + 
     -            (8*Cpl2*GS*(WProd2170 + WProd2173))/WProd2168) - 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*Prp229*
     -             (WProd2112*WProd2140 + WProd2118*WProd2140 + 
     -               WProd2114*WProd2145 + WProd2114*WProd2148 - 
     -               WProd2165 + 2*WProd2050*WProd2174 - 
     -               2*WProd2052*WProd2174 - 
     -               2*WProd2053*WProd2174 - 
     -               2*WProd2054*WProd2174 + 
     -               2*WProd2055*WProd2174 - 
     -               2*WProd2056*WProd2174 - 
     -               2*WProd2057*WProd2174 + 
     -               2*WProd2058*WProd2174 + 
     -               2*WProd2059*WProd2174 - 
     -               2*WProd2070*WProd2174 - WProd2175))/WProd2160+
     -           Prp230*(5.333333333333333D0*
     -             (Cpl2*GS*(WProd2170 + WProd2173))/WProd2168 + 
     -            5.333333333333333D0*
     -             (Cpl2*GS*
     -                (WProd2114*WProd2138 - WProd2105*WProd2139 - 
     -                  WProd2106*WProd2139 - 
     -                  WProd2045*WProd2105*WProd2140 - 
     -                  WProd2045*WProd2106*WProd2140 + 
     -                  WProd2109*WProd2140 + 
     -                  WProd2110*WProd2140 + 
     -                  WProd2111*WProd2140 + 
     -                  WProd2115*WProd2140 + 
     -                  WProd2118*WProd2140 + 
     -                  WProd2114*WProd2144 + 
     -                  WProd2114*WProd2148 - WProd2164 - 
     -                  2*WProd2050*WProd2167 + 
     -                  WProd2052*WProd2167 + 
     -                  WProd2053*WProd2167 + 
     -                  WProd2054*WProd2167 - 
     -                  WProd2055*WProd2167 + 
     -                  WProd2056*WProd2167 + 
     -                  WProd2057*WProd2167 - 
     -                  WProd2058*WProd2167 - 
     -                  WProd2059*WProd2167 + 
     -                  WProd2070*WProd2167 - WProd2172 + 
     -                  2*WProd2050*WProd2174 - 
     -                  WProd2052*WProd2174 - 
     -                  WProd2053*WProd2174 - 
     -                  WProd2054*WProd2174 + 
     -                  WProd2055*WProd2174 - 
     -                  WProd2056*WProd2174 - 
     -                  WProd2057*WProd2174 + 
     -                  WProd2058*WProd2174 + 
     -                  WProd2059*WProd2174 - 
     -                  WProd2070*WProd2174 - WProd2175))/WProd2160
     -            )))

#ifdef DEBUG
	print *, 'amp(27,7,5) =', amp(27,7,5)
#endif

      endif

      end


