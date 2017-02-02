      subroutine realQED26(amp,resnr)
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

        tmp1152 = Prp216*
     -    ((-4*Cpl2*GS*Prp217*
     -         (2*WProd1888*WProd1889 - 2*WProd1889*WProd1890 - 
     -           2*WProd1889*WProd1891 - 2*WProd1892 + 
     -           2*WProd1893 + 2*WProd1894 - 
     -           2*WProd1895*WProd1896 + 2*WProd1897 + 
     -           2*WProd1898 + 2*WProd1895*WProd1899 - 
     -           2*WProd1896*WProd1900 + 2*WProd1899*WProd1900 + 
     -           2*WProd1895*WProd1901 + 2*WProd1900*WProd1901 + 
     -           2*WProd1901*WProd1902 + 3*WProd1903*WProd1904 + 
     -           3*WProd1904*WProd1905 + 3*WProd1904*WProd1906 + 
     -           2*WProd1895*WProd1907 + 2*WProd1900*WProd1907 + 
     -           4*WProd1904*WProd1910 - 4*WProd1904*WProd1926 + 
     -           4*WProd1904*WProd1927 + 4*WProd1904*WProd1928 - 
     -           5*WProd1904*WProd1929 - 5*WProd1904*WProd1930 - 
     -           5*WProd1904*WProd1931))/WProd1908 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp222*
     -          (-(WProd1895*WProd1899) - WProd1899*WProd1900 - 
     -            WProd1899*WProd1902 - WProd1895*WProd1907 - 
     -            WProd1900*WProd1907 - WProd1902*WProd1907 - 
     -            2*WProd1889*WProd1903*WProd1915 - 
     -            2*WProd1889*WProd1905*WProd1915 - 
     -            2*WProd1889*WProd1906*WProd1915 - 
     -            2*WProd1889*WProd1910*WProd1915 + 
     -            2*WProd1889*WProd1915*WProd1926 - 
     -            2*WProd1889*WProd1915*WProd1927 - 
     -            2*WProd1889*WProd1915*WProd1928 + 
     -            3*WProd1889*WProd1915*WProd1929 + 
     -            3*WProd1889*WProd1915*WProd1930 + 
     -            3*WProd1889*WProd1915*WProd1931 - 
     -            WProd1889*WProd1933 - WProd1889*WProd1934 + 
     -            2*WProd1903*WProd1935 + 2*WProd1905*WProd1935 + 
     -            2*WProd1906*WProd1935 + 2*WProd1910*WProd1935 - 
     -            2*WProd1926*WProd1935 + 2*WProd1927*WProd1935 + 
     -            2*WProd1928*WProd1935 - 3*WProd1929*WProd1935 - 
     -            3*WProd1930*WProd1935 - 3*WProd1931*WProd1935 + 
     -            WProd1936 + WProd1937))/WProd1908) - 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp216*Prp220*
     -       (WProd1888*WProd1889 - WProd1889*WProd1890 - 
     -         WProd1889*WProd1891 - WProd1892 + WProd1893 + 
     -         WProd1894 - WProd1895*WProd1896 + WProd1897 + 
     -         WProd1898 - WProd1896*WProd1900 + 
     -         WProd1895*WProd1901 + WProd1900*WProd1901 - 
     -         WProd1899*WProd1902 + WProd1901*WProd1902 + 
     -         WProd1903*WProd1904 + WProd1904*WProd1905 + 
     -         WProd1904*WProd1906 - WProd1902*WProd1907 + 
     -         2*WProd1904*WProd1910 - 
     -         WProd1889*WProd1903*WProd1915 - 
     -         WProd1889*WProd1905*WProd1915 - 
     -         WProd1889*WProd1906*WProd1915 - 
     -         2*WProd1889*WProd1910*WProd1915 - 
     -         2*WProd1904*WProd1926 + 
     -         2*WProd1889*WProd1915*WProd1926 + 
     -         2*WProd1904*WProd1927 - 
     -         2*WProd1889*WProd1915*WProd1927 + 
     -         2*WProd1904*WProd1928 - 
     -         2*WProd1889*WProd1915*WProd1928 - 
     -         2*WProd1904*WProd1929 + 
     -         2*WProd1889*WProd1915*WProd1929 - 
     -         2*WProd1904*WProd1930 + 
     -         2*WProd1889*WProd1915*WProd1930 - 
     -         2*WProd1904*WProd1931 + 
     -         2*WProd1889*WProd1915*WProd1931 - 
     -         WProd1889*WProd1933 - WProd1889*WProd1934 + 
     -         WProd1903*WProd1935 + WProd1905*WProd1935 + 
     -         WProd1906*WProd1935 + 2*WProd1910*WProd1935 - 
     -         2*WProd1926*WProd1935 + 2*WProd1927*WProd1935 + 
     -         2*WProd1928*WProd1935 - 2*WProd1929*WProd1935 - 
     -         2*WProd1930*WProd1935 - 2*WProd1931*WProd1935 + 
     -         WProd1936 + WProd1937))/WProd1908 + 
     -   Prp223*(-(2.6666666666666665D0*
     -         (Cpl2*GS*Prp222*
     -            (WProd1899*WProd1925 + WProd1907*WProd1925 - 
     -              2*WProd1903*WProd1935 - 
     -              2*WProd1905*WProd1935 - 
     -              2*WProd1906*WProd1935 - 
     -              2*WProd1910*WProd1935 + 
     -              2*WProd1926*WProd1935 - 
     -              2*WProd1927*WProd1935 - 
     -              2*WProd1928*WProd1935 + 
     -              2*WProd1929*WProd1935 + 
     -              2*WProd1930*WProd1935 + 
     -              2*WProd1931*WProd1935 - WProd1936 - 
     -              WProd1937 + WProd1913*WProd1938 + 
     -              WProd1913*WProd1939))/WProd1908) - 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp218*
     -          (-WProd1892 + WProd1912*WProd1913 + 
     -            WProd1899*WProd1925 + WProd1901*WProd1925 - 
     -            WProd1915*WProd1919*WProd1925 - 
     -            WProd1915*WProd1920*WProd1925 + 
     -            WProd1922*WProd1925 + WProd1923*WProd1925 + 
     -            WProd1924*WProd1925 - WProd1903*WProd1935 - 
     -            WProd1905*WProd1935 - WProd1906*WProd1935 - 
     -            WProd1910*WProd1935 + 2*WProd1926*WProd1935 - 
     -            WProd1927*WProd1935 - WProd1928*WProd1935 + 
     -            WProd1929*WProd1935 + WProd1930*WProd1935 + 
     -            WProd1931*WProd1935 - WProd1936 + 
     -            WProd1913*WProd1938 + WProd1913*WProd1941 - 
     -            WProd1942 + WProd1903*WProd1943 + 
     -            WProd1905*WProd1943 + WProd1906*WProd1943 + 
     -            WProd1910*WProd1943 - 2*WProd1926*WProd1943 + 
     -            WProd1927*WProd1943 + WProd1928*WProd1943 - 
     -            WProd1929*WProd1943 - WProd1930*WProd1943 - 
     -            WProd1931*WProd1943))/WProd1908)
        tmp1152 = tmp1152 - 
     -   (4*Cpl2*GS*Prp217*Prp218*
     -      (-2*WProd1892 + 2*WProd1893 + 2*WProd1894 - 
     -        2*WProd1895*WProd1896 + 2*WProd1897 + 2*WProd1898 + 
     -        2*WProd1895*WProd1899 - 2*WProd1896*WProd1900 + 
     -        2*WProd1899*WProd1900 + 2*WProd1895*WProd1901 + 
     -        2*WProd1900*WProd1901 + 2*WProd1903*WProd1904 + 
     -        2*WProd1904*WProd1905 + 3*WProd1904*WProd1906 + 
     -        3*WProd1904*WProd1910 + 2*WProd1912*WProd1913 + 
     -        2*WProd1895*WProd1913*WProd1915 + 
     -        2*WProd1900*WProd1913*WProd1915 - 
     -        2*WProd1912*WProd1919 + 2*WProd1914*WProd1919 - 
     -        2*WProd1895*WProd1915*WProd1919 - 
     -        2*WProd1900*WProd1915*WProd1919 - 
     -        2*WProd1912*WProd1920 + 2*WProd1914*WProd1920 - 
     -        2*WProd1895*WProd1915*WProd1920 - 
     -        2*WProd1900*WProd1915*WProd1920 - 
     -        2*WProd1895*WProd1921 - 2*WProd1900*WProd1921 + 
     -        2*WProd1895*WProd1922 + 2*WProd1900*WProd1922 + 
     -        2*WProd1895*WProd1923 + 2*WProd1900*WProd1923 + 
     -        2*WProd1895*WProd1924 + 2*WProd1900*WProd1924 + 
     -        2*WProd1901*WProd1925 - 4*WProd1904*WProd1926 + 
     -        2*WProd1904*WProd1927 + 2*WProd1904*WProd1928 - 
     -        2*WProd1904*WProd1929 - 3*WProd1904*WProd1930 - 
     -        3*WProd1904*WProd1931 - 2*WProd1913*WProd1945 - 
     -        2*WProd1913*WProd1946))/WProd1908 + 
     -   Prp221*(5.333333333333333D0*
     -       (Cpl2*GS*Prp220*
     -          (WProd1889*WProd1890 + WProd1889*WProd1891 - 
     -            WProd1893 - WProd1894 + WProd1895*WProd1896 - 
     -            WProd1897 - WProd1898 + WProd1896*WProd1900 - 
     -            WProd1895*WProd1901 - WProd1900*WProd1901 + 
     -            WProd1899*WProd1902 - WProd1904*WProd1906 - 
     -            WProd1904*WProd1910 + 
     -            WProd1889*WProd1906*WProd1915 + 
     -            WProd1889*WProd1910*WProd1915 + 
     -            WProd1912*WProd1919 - WProd1914*WProd1919 - 
     -            WProd1902*WProd1915*WProd1919 + 
     -            WProd1912*WProd1920 - WProd1914*WProd1920 - 
     -            WProd1902*WProd1915*WProd1920 + 
     -            WProd1895*WProd1921 + WProd1900*WProd1921 + 
     -            WProd1902*WProd1922 + WProd1902*WProd1923 + 
     -            WProd1902*WProd1924 + WProd1904*WProd1926 - 
     -            WProd1889*WProd1915*WProd1926 - 
     -            WProd1904*WProd1927 + 
     -            WProd1889*WProd1915*WProd1927 - 
     -            WProd1904*WProd1928 + 
     -            WProd1889*WProd1915*WProd1928 + 
     -            WProd1904*WProd1929 - 
     -            WProd1889*WProd1915*WProd1929 + 
     -            WProd1904*WProd1930 - 
     -            WProd1889*WProd1915*WProd1930 + 
     -            WProd1904*WProd1931 - 
     -            WProd1889*WProd1915*WProd1931 + 
     -            WProd1889*WProd1933 - WProd1906*WProd1935 - 
     -            WProd1910*WProd1935 + WProd1926*WProd1935 - 
     -            WProd1927*WProd1935 - WProd1928*WProd1935 + 
     -            WProd1929*WProd1935 + WProd1930*WProd1935 + 
     -            WProd1931*WProd1935 - WProd1936 - WProd1942 + 
     -            WProd1906*WProd1943 + WProd1910*WProd1943 - 
     -            WProd1926*WProd1943 + WProd1927*WProd1943 + 
     -            WProd1928*WProd1943 - WProd1929*WProd1943 - 
     -            WProd1930*WProd1943 - WProd1931*WProd1943 + 
     -            WProd1889*WProd1944))/WProd1908 + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp218*
     -          (-WProd1893 - WProd1894 + WProd1895*WProd1896 - 
     -            WProd1897 - WProd1898 - WProd1895*WProd1899 + 
     -            WProd1896*WProd1900 - WProd1899*WProd1900 - 
     -            WProd1895*WProd1901 - WProd1900*WProd1901 - 
     -            WProd1903*WProd1904 - WProd1904*WProd1905 - 
     -            2*WProd1904*WProd1906 - 2*WProd1904*WProd1910 - 
     -            WProd1895*WProd1913*WProd1915 - 
     -            WProd1900*WProd1913*WProd1915 + 
     -            WProd1912*WProd1919 - WProd1914*WProd1919 + 
     -            WProd1895*WProd1915*WProd1919 + 
     -            WProd1900*WProd1915*WProd1919 + 
     -            WProd1912*WProd1920 - WProd1914*WProd1920 + 
     -            WProd1895*WProd1915*WProd1920 + 
     -            WProd1900*WProd1915*WProd1920 + 
     -            WProd1895*WProd1921 + WProd1900*WProd1921 - 
     -            WProd1895*WProd1922 - WProd1900*WProd1922 - 
     -            WProd1895*WProd1923 - WProd1900*WProd1923 - 
     -            WProd1895*WProd1924 - WProd1900*WProd1924 + 
     -            WProd1899*WProd1925 - 
     -            WProd1915*WProd1919*WProd1925 - 
     -            WProd1915*WProd1920*WProd1925 + 
     -            WProd1922*WProd1925 + WProd1923*WProd1925 + 
     -            WProd1924*WProd1925 + 2*WProd1904*WProd1926 - 
     -            WProd1904*WProd1927 - WProd1904*WProd1928 + 
     -            WProd1904*WProd1929 + 2*WProd1904*WProd1930 + 
     -            2*WProd1904*WProd1931 - WProd1903*WProd1935 - 
     -            WProd1905*WProd1935 - 2*WProd1906*WProd1935 - 
     -            2*WProd1910*WProd1935 + 2*WProd1926*WProd1935 - 
     -            WProd1927*WProd1935 - WProd1928*WProd1935 + 
     -            WProd1929*WProd1935 + 2*WProd1930*WProd1935 + 
     -            2*WProd1931*WProd1935 - WProd1936 + 
     -            WProd1913*WProd1938 + WProd1913*WProd1941 - 
     -            WProd1942 + WProd1903*WProd1943 + 
     -            WProd1905*WProd1943 + 2*WProd1906*WProd1943 + 
     -            2*WProd1910*WProd1943 - 2*WProd1926*WProd1943 + 
     -            WProd1927*WProd1943 + WProd1928*WProd1943 - 
     -            WProd1929*WProd1943 - 2*WProd1930*WProd1943 - 
     -            2*WProd1931*WProd1943 + WProd1913*WProd1945 + 
     -            WProd1913*WProd1946))/WProd1908)
        amp(26,7,3) = Prp219*tmp1152 - 
     -   (8*Cpl2*GS*Prp215*Prp216*Prp217*
     -      (WProd1888*WProd1889 - WProd1889*WProd1890 - 
     -        WProd1889*WProd1891 - WProd1892 + WProd1893 + 
     -        WProd1894 - WProd1895*WProd1896 + WProd1897 + 
     -        WProd1898 + WProd1895*WProd1899 + 
     -        WProd1895*WProd1901 + WProd1900*WProd1901 + 
     -        WProd1901*WProd1902 - WProd1903*WProd1904 - 
     -        WProd1904*WProd1905 - WProd1904*WProd1906 + 
     -        WProd1895*WProd1907))/WProd1908 - 
     -   (8*Cpl2*GS*Prp215*Prp217*Prp218*
     -      (-WProd1892 + WProd1893 + WProd1894 - 
     -        WProd1895*WProd1896 + WProd1897 + WProd1898 + 
     -        WProd1895*WProd1899 + WProd1895*WProd1901 + 
     -        WProd1900*WProd1901 - 2*WProd1903*WProd1904 - 
     -        WProd1904*WProd1905 - WProd1904*WProd1906 - 
     -        WProd1904*WProd1910 + WProd1912*WProd1913 - 
     -        WProd1913*WProd1914 + 
     -        WProd1895*WProd1913*WProd1915 + WProd1916 - 
     -        WProd1917 - WProd1918 - WProd1912*WProd1919 + 
     -        WProd1914*WProd1919 - 
     -        WProd1895*WProd1915*WProd1919 - 
     -        WProd1912*WProd1920 + WProd1914*WProd1920 - 
     -        WProd1895*WProd1915*WProd1920 - 
     -        WProd1895*WProd1921 + WProd1895*WProd1922 + 
     -        WProd1895*WProd1923 + WProd1895*WProd1924 + 
     -        WProd1901*WProd1925 + WProd1904*WProd1926 - 
     -        WProd1904*WProd1927 - WProd1904*WProd1928 + 
     -        WProd1904*WProd1929 + WProd1904*WProd1930 + 
     -        WProd1904*WProd1931 - WProd1932))/WProd1908

#ifdef DEBUG
	print *, 'amp(26,7,3) =', amp(26,7,3)
#endif

        tmp1153 = (-4*Cpl2*GS*Prp217*Prp218*
     -      (-2*WProd1949 + 2*WProd1950 + 2*WProd1951 - 
     -        2*WProd1895*WProd1952 - 2*WProd1900*WProd1952 + 
     -        2*WProd1953 + 2*WProd1954 + 2*WProd1895*WProd1955 + 
     -        2*WProd1900*WProd1955 + 2*WProd1895*WProd1956 + 
     -        2*WProd1900*WProd1956 + 2*WProd1925*WProd1956 + 
     -        2*WProd1903*WProd1957 + 2*WProd1905*WProd1957 + 
     -        3*WProd1906*WProd1957 + 3*WProd1910*WProd1957 - 
     -        4*WProd1926*WProd1957 + 2*WProd1927*WProd1957 + 
     -        2*WProd1928*WProd1957 - 2*WProd1929*WProd1957 - 
     -        3*WProd1930*WProd1957 - 3*WProd1931*WProd1957 + 
     -        2*WProd1912*WProd1959 + 
     -        2*WProd1895*WProd1915*WProd1959 + 
     -        2*WProd1900*WProd1915*WProd1959 - 
     -        2*WProd1945*WProd1959 - 2*WProd1946*WProd1959 - 
     -        2*WProd1912*WProd1963 + 2*WProd1914*WProd1963 - 
     -        2*WProd1895*WProd1915*WProd1963 - 
     -        2*WProd1900*WProd1915*WProd1963 - 
     -        2*WProd1912*WProd1964 + 2*WProd1914*WProd1964 - 
     -        2*WProd1895*WProd1915*WProd1964 - 
     -        2*WProd1900*WProd1915*WProd1964 - 
     -        2*WProd1895*WProd1965 - 2*WProd1900*WProd1965 + 
     -        2*WProd1895*WProd1966 + 2*WProd1900*WProd1966 + 
     -        2*WProd1895*WProd1967 + 2*WProd1900*WProd1967 + 
     -        2*WProd1895*WProd1968 + 2*WProd1900*WProd1968))/
     -    WProd1947 + Prp216*
     -    ((-4*Cpl2*GS*Prp217*
     -         (2*WProd1888*WProd1948 - 2*WProd1890*WProd1948 - 
     -           2*WProd1891*WProd1948 - 2*WProd1949 + 
     -           2*WProd1950 + 2*WProd1951 - 
     -           2*WProd1895*WProd1952 - 2*WProd1900*WProd1952 + 
     -           2*WProd1953 + 2*WProd1954 + 
     -           2*WProd1895*WProd1955 + 2*WProd1900*WProd1955 + 
     -           2*WProd1895*WProd1956 + 2*WProd1900*WProd1956 + 
     -           2*WProd1902*WProd1956 + 3*WProd1903*WProd1957 + 
     -           3*WProd1905*WProd1957 + 3*WProd1906*WProd1957 + 
     -           4*WProd1910*WProd1957 - 4*WProd1926*WProd1957 + 
     -           4*WProd1927*WProd1957 + 4*WProd1928*WProd1957 - 
     -           5*WProd1929*WProd1957 - 5*WProd1930*WProd1957 - 
     -           5*WProd1931*WProd1957 + 2*WProd1895*WProd1958 + 
     -           2*WProd1900*WProd1958))/WProd1947 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp222*
     -          (-2*WProd1903*WProd1915*WProd1948 - 
     -            2*WProd1905*WProd1915*WProd1948 - 
     -            2*WProd1906*WProd1915*WProd1948 - 
     -            2*WProd1910*WProd1915*WProd1948 + 
     -            2*WProd1915*WProd1926*WProd1948 - 
     -            2*WProd1915*WProd1927*WProd1948 - 
     -            2*WProd1915*WProd1928*WProd1948 + 
     -            3*WProd1915*WProd1929*WProd1948 + 
     -            3*WProd1915*WProd1930*WProd1948 + 
     -            3*WProd1915*WProd1931*WProd1948 - 
     -            WProd1933*WProd1948 - WProd1934*WProd1948 - 
     -            WProd1895*WProd1955 - WProd1900*WProd1955 - 
     -            WProd1902*WProd1955 - WProd1895*WProd1958 - 
     -            WProd1900*WProd1958 - WProd1902*WProd1958 + 
     -            2*WProd1903*WProd1970 + 2*WProd1905*WProd1970 + 
     -            2*WProd1906*WProd1970 + 2*WProd1910*WProd1970 - 
     -            2*WProd1926*WProd1970 + 2*WProd1927*WProd1970 + 
     -            2*WProd1928*WProd1970 - 3*WProd1929*WProd1970 - 
     -            3*WProd1930*WProd1970 - 3*WProd1931*WProd1970 + 
     -            WProd1971 + WProd1972))/WProd1947) - 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp216*Prp220*
     -       (WProd1888*WProd1948 - WProd1890*WProd1948 - 
     -         WProd1891*WProd1948 - 
     -         WProd1903*WProd1915*WProd1948 - 
     -         WProd1905*WProd1915*WProd1948 - 
     -         WProd1906*WProd1915*WProd1948 - 
     -         2*WProd1910*WProd1915*WProd1948 + 
     -         2*WProd1915*WProd1926*WProd1948 - 
     -         2*WProd1915*WProd1927*WProd1948 - 
     -         2*WProd1915*WProd1928*WProd1948 + 
     -         2*WProd1915*WProd1929*WProd1948 + 
     -         2*WProd1915*WProd1930*WProd1948 + 
     -         2*WProd1915*WProd1931*WProd1948 - 
     -         WProd1933*WProd1948 - WProd1934*WProd1948 - 
     -         WProd1949 + WProd1950 + WProd1951 - 
     -         WProd1895*WProd1952 - WProd1900*WProd1952 + 
     -         WProd1953 + WProd1954 - WProd1902*WProd1955 + 
     -         WProd1895*WProd1956 + WProd1900*WProd1956 + 
     -         WProd1902*WProd1956 + WProd1903*WProd1957 + 
     -         WProd1905*WProd1957 + WProd1906*WProd1957 + 
     -         2*WProd1910*WProd1957 - 2*WProd1926*WProd1957 + 
     -         2*WProd1927*WProd1957 + 2*WProd1928*WProd1957 - 
     -         2*WProd1929*WProd1957 - 2*WProd1930*WProd1957 - 
     -         2*WProd1931*WProd1957 - WProd1902*WProd1958 + 
     -         WProd1903*WProd1970 + WProd1905*WProd1970 + 
     -         WProd1906*WProd1970 + 2*WProd1910*WProd1970 - 
     -         2*WProd1926*WProd1970 + 2*WProd1927*WProd1970 + 
     -         2*WProd1928*WProd1970 - 2*WProd1929*WProd1970 - 
     -         2*WProd1930*WProd1970 - 2*WProd1931*WProd1970 + 
     -         WProd1971 + WProd1972))/WProd1947
        tmp1153 = tmp1153 + 
     -   Prp223*(-(2.6666666666666665D0*
     -         (Cpl2*GS*Prp222*
     -            (WProd1925*WProd1955 + WProd1925*WProd1958 + 
     -              WProd1938*WProd1959 + WProd1939*WProd1959 - 
     -              2*WProd1903*WProd1970 - 
     -              2*WProd1905*WProd1970 - 
     -              2*WProd1906*WProd1970 - 
     -              2*WProd1910*WProd1970 + 
     -              2*WProd1926*WProd1970 - 
     -              2*WProd1927*WProd1970 - 
     -              2*WProd1928*WProd1970 + 
     -              2*WProd1929*WProd1970 + 
     -              2*WProd1930*WProd1970 + 
     -              2*WProd1931*WProd1970 - WProd1971 - WProd1972))
     -           /WProd1947) - 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp218*
     -          (-WProd1949 + WProd1925*WProd1955 + 
     -            WProd1925*WProd1956 + WProd1912*WProd1959 + 
     -            WProd1938*WProd1959 + WProd1941*WProd1959 - 
     -            WProd1915*WProd1925*WProd1963 - 
     -            WProd1915*WProd1925*WProd1964 + 
     -            WProd1925*WProd1966 + WProd1925*WProd1967 + 
     -            WProd1925*WProd1968 - WProd1903*WProd1970 - 
     -            WProd1905*WProd1970 - WProd1906*WProd1970 - 
     -            WProd1910*WProd1970 + 2*WProd1926*WProd1970 - 
     -            WProd1927*WProd1970 - WProd1928*WProd1970 + 
     -            WProd1929*WProd1970 + WProd1930*WProd1970 + 
     -            WProd1931*WProd1970 - WProd1971 - WProd1973 + 
     -            WProd1903*WProd1974 + WProd1905*WProd1974 + 
     -            WProd1906*WProd1974 + WProd1910*WProd1974 - 
     -            2*WProd1926*WProd1974 + WProd1927*WProd1974 + 
     -            WProd1928*WProd1974 - WProd1929*WProd1974 - 
     -            WProd1930*WProd1974 - WProd1931*WProd1974))/
     -        WProd1947) + 
     -   Prp221*(5.333333333333333D0*
     -       (Cpl2*GS*Prp218*
     -          (-WProd1950 - WProd1951 + WProd1895*WProd1952 + 
     -            WProd1900*WProd1952 - WProd1953 - WProd1954 - 
     -            WProd1895*WProd1955 - WProd1900*WProd1955 + 
     -            WProd1925*WProd1955 - WProd1895*WProd1956 - 
     -            WProd1900*WProd1956 - WProd1903*WProd1957 - 
     -            WProd1905*WProd1957 - 2*WProd1906*WProd1957 - 
     -            2*WProd1910*WProd1957 + 2*WProd1926*WProd1957 - 
     -            WProd1927*WProd1957 - WProd1928*WProd1957 + 
     -            WProd1929*WProd1957 + 2*WProd1930*WProd1957 + 
     -            2*WProd1931*WProd1957 - 
     -            WProd1895*WProd1915*WProd1959 - 
     -            WProd1900*WProd1915*WProd1959 + 
     -            WProd1938*WProd1959 + WProd1941*WProd1959 + 
     -            WProd1945*WProd1959 + WProd1946*WProd1959 + 
     -            WProd1912*WProd1963 - WProd1914*WProd1963 + 
     -            WProd1895*WProd1915*WProd1963 + 
     -            WProd1900*WProd1915*WProd1963 - 
     -            WProd1915*WProd1925*WProd1963 + 
     -            WProd1912*WProd1964 - WProd1914*WProd1964 + 
     -            WProd1895*WProd1915*WProd1964 + 
     -            WProd1900*WProd1915*WProd1964 - 
     -            WProd1915*WProd1925*WProd1964 + 
     -            WProd1895*WProd1965 + WProd1900*WProd1965 - 
     -            WProd1895*WProd1966 - WProd1900*WProd1966 + 
     -            WProd1925*WProd1966 - WProd1895*WProd1967 - 
     -            WProd1900*WProd1967 + WProd1925*WProd1967 - 
     -            WProd1895*WProd1968 - WProd1900*WProd1968 + 
     -            WProd1925*WProd1968 - WProd1903*WProd1970 - 
     -            WProd1905*WProd1970 - 2*WProd1906*WProd1970 - 
     -            2*WProd1910*WProd1970 + 2*WProd1926*WProd1970 - 
     -            WProd1927*WProd1970 - WProd1928*WProd1970 + 
     -            WProd1929*WProd1970 + 2*WProd1930*WProd1970 + 
     -            2*WProd1931*WProd1970 - WProd1971 - WProd1973 + 
     -            WProd1903*WProd1974 + WProd1905*WProd1974 + 
     -            2*WProd1906*WProd1974 + 2*WProd1910*WProd1974 - 
     -            2*WProd1926*WProd1974 + WProd1927*WProd1974 + 
     -            WProd1928*WProd1974 - WProd1929*WProd1974 - 
     -            2*WProd1930*WProd1974 - 2*WProd1931*WProd1974))/
     -        WProd1947 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp220*
     -          (WProd1890*WProd1948 + WProd1891*WProd1948 + 
     -            WProd1906*WProd1915*WProd1948 + 
     -            WProd1910*WProd1915*WProd1948 - 
     -            WProd1915*WProd1926*WProd1948 + 
     -            WProd1915*WProd1927*WProd1948 + 
     -            WProd1915*WProd1928*WProd1948 - 
     -            WProd1915*WProd1929*WProd1948 - 
     -            WProd1915*WProd1930*WProd1948 - 
     -            WProd1915*WProd1931*WProd1948 + 
     -            WProd1933*WProd1948 + WProd1944*WProd1948 - 
     -            WProd1950 - WProd1951 + WProd1895*WProd1952 + 
     -            WProd1900*WProd1952 - WProd1953 - WProd1954 + 
     -            WProd1902*WProd1955 - WProd1895*WProd1956 - 
     -            WProd1900*WProd1956 - WProd1906*WProd1957 - 
     -            WProd1910*WProd1957 + WProd1926*WProd1957 - 
     -            WProd1927*WProd1957 - WProd1928*WProd1957 + 
     -            WProd1929*WProd1957 + WProd1930*WProd1957 + 
     -            WProd1931*WProd1957 + WProd1912*WProd1963 - 
     -            WProd1914*WProd1963 - 
     -            WProd1902*WProd1915*WProd1963 + 
     -            WProd1912*WProd1964 - WProd1914*WProd1964 - 
     -            WProd1902*WProd1915*WProd1964 + 
     -            WProd1895*WProd1965 + WProd1900*WProd1965 + 
     -            WProd1902*WProd1966 + WProd1902*WProd1967 + 
     -            WProd1902*WProd1968 - WProd1906*WProd1970 - 
     -            WProd1910*WProd1970 + WProd1926*WProd1970 - 
     -            WProd1927*WProd1970 - WProd1928*WProd1970 + 
     -            WProd1929*WProd1970 + WProd1930*WProd1970 + 
     -            WProd1931*WProd1970 - WProd1971 - WProd1973 + 
     -            WProd1906*WProd1974 + WProd1910*WProd1974 - 
     -            WProd1926*WProd1974 + WProd1927*WProd1974 + 
     -            WProd1928*WProd1974 - WProd1929*WProd1974 - 
     -            WProd1930*WProd1974 - WProd1931*WProd1974))/
     -        WProd1947)
        amp(26,7,9) = Prp219*tmp1153 - 
     -   (8*Cpl2*GS*Prp215*Prp216*Prp217*
     -      (WProd1888*WProd1948 - WProd1890*WProd1948 - 
     -        WProd1891*WProd1948 - WProd1949 + WProd1950 + 
     -        WProd1951 - WProd1895*WProd1952 + WProd1953 + 
     -        WProd1954 + WProd1895*WProd1955 + 
     -        WProd1895*WProd1956 + WProd1900*WProd1956 + 
     -        WProd1902*WProd1956 - WProd1903*WProd1957 - 
     -        WProd1905*WProd1957 - WProd1906*WProd1957 + 
     -        WProd1895*WProd1958))/WProd1947 - 
     -   (8*Cpl2*GS*Prp215*Prp217*Prp218*
     -      (-WProd1949 + WProd1950 + WProd1951 - 
     -        WProd1895*WProd1952 + WProd1953 + WProd1954 + 
     -        WProd1895*WProd1955 + WProd1895*WProd1956 + 
     -        WProd1900*WProd1956 + WProd1925*WProd1956 - 
     -        2*WProd1903*WProd1957 - WProd1905*WProd1957 - 
     -        WProd1906*WProd1957 - WProd1910*WProd1957 + 
     -        WProd1926*WProd1957 - WProd1927*WProd1957 - 
     -        WProd1928*WProd1957 + WProd1929*WProd1957 + 
     -        WProd1930*WProd1957 + WProd1931*WProd1957 + 
     -        WProd1912*WProd1959 - WProd1914*WProd1959 + 
     -        WProd1895*WProd1915*WProd1959 + WProd1960 - 
     -        WProd1961 - WProd1962 - WProd1912*WProd1963 + 
     -        WProd1914*WProd1963 - 
     -        WProd1895*WProd1915*WProd1963 - 
     -        WProd1912*WProd1964 + WProd1914*WProd1964 - 
     -        WProd1895*WProd1915*WProd1964 - 
     -        WProd1895*WProd1965 + WProd1895*WProd1966 + 
     -        WProd1895*WProd1967 + WProd1895*WProd1968 - WProd1969
     -        ))/WProd1947

#ifdef DEBUG
	print *, 'amp(26,7,9) =', amp(26,7,9)
#endif

        tmp1154 = Prp216*
     -    ((-4*Cpl2*GS*Prp217*
     -         (2*WProd1977 - 2*WProd1896*WProd1978 + 
     -           2*WProd1899*WProd1978 + 2*WProd1901*WProd1978 + 
     -           2*WProd1907*WProd1978 + 2*WProd1979 + 
     -           2*WProd1888*WProd1980 - 2*WProd1890*WProd1980 - 
     -           2*WProd1891*WProd1980 - 2*WProd1981 + 
     -           2*WProd1982 + 2*WProd1983 - 
     -           2*WProd1896*WProd1986 + 2*WProd1899*WProd1986 + 
     -           2*WProd1901*WProd1986 + 2*WProd1907*WProd1986 + 
     -           2*WProd1901*WProd1987 + 3*WProd1903*WProd1988 + 
     -           3*WProd1905*WProd1988 + 3*WProd1906*WProd1988 + 
     -           4*WProd1910*WProd1988 - 4*WProd1926*WProd1988 + 
     -           4*WProd1927*WProd1988 + 4*WProd1928*WProd1988 - 
     -           5*WProd1929*WProd1988 - 5*WProd1930*WProd1988 - 
     -           5*WProd1931*WProd1988))/WProd1989 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp222*
     -          (-(WProd1899*WProd1978) - WProd1907*WProd1978 - 
     -            2*WProd1903*WProd1915*WProd1980 - 
     -            2*WProd1905*WProd1915*WProd1980 - 
     -            2*WProd1906*WProd1915*WProd1980 - 
     -            2*WProd1910*WProd1915*WProd1980 + 
     -            2*WProd1915*WProd1926*WProd1980 - 
     -            2*WProd1915*WProd1927*WProd1980 - 
     -            2*WProd1915*WProd1928*WProd1980 + 
     -            3*WProd1915*WProd1929*WProd1980 + 
     -            3*WProd1915*WProd1930*WProd1980 + 
     -            3*WProd1915*WProd1931*WProd1980 - 
     -            WProd1933*WProd1980 - WProd1934*WProd1980 - 
     -            WProd1899*WProd1986 - WProd1907*WProd1986 - 
     -            WProd1899*WProd1987 - WProd1907*WProd1987 + 
     -            2*WProd1903*WProd1997 + 2*WProd1905*WProd1997 + 
     -            2*WProd1906*WProd1997 + 2*WProd1910*WProd1997 - 
     -            2*WProd1926*WProd1997 + 2*WProd1927*WProd1997 + 
     -            2*WProd1928*WProd1997 - 3*WProd1929*WProd1997 - 
     -            3*WProd1930*WProd1997 - 3*WProd1931*WProd1997 + 
     -            WProd1998 + WProd1999))/WProd1989) - 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp216*Prp220*
     -       (WProd1977 - WProd1896*WProd1978 + 
     -         WProd1901*WProd1978 + WProd1979 + 
     -         WProd1888*WProd1980 - WProd1890*WProd1980 - 
     -         WProd1891*WProd1980 - 
     -         WProd1903*WProd1915*WProd1980 - 
     -         WProd1905*WProd1915*WProd1980 - 
     -         WProd1906*WProd1915*WProd1980 - 
     -         2*WProd1910*WProd1915*WProd1980 + 
     -         2*WProd1915*WProd1926*WProd1980 - 
     -         2*WProd1915*WProd1927*WProd1980 - 
     -         2*WProd1915*WProd1928*WProd1980 + 
     -         2*WProd1915*WProd1929*WProd1980 + 
     -         2*WProd1915*WProd1930*WProd1980 + 
     -         2*WProd1915*WProd1931*WProd1980 - 
     -         WProd1933*WProd1980 - WProd1934*WProd1980 - 
     -         WProd1981 + WProd1982 + WProd1983 - 
     -         WProd1896*WProd1986 + WProd1901*WProd1986 - 
     -         WProd1899*WProd1987 + WProd1901*WProd1987 - 
     -         WProd1907*WProd1987 + WProd1903*WProd1988 + 
     -         WProd1905*WProd1988 + WProd1906*WProd1988 + 
     -         2*WProd1910*WProd1988 - 2*WProd1926*WProd1988 + 
     -         2*WProd1927*WProd1988 + 2*WProd1928*WProd1988 - 
     -         2*WProd1929*WProd1988 - 2*WProd1930*WProd1988 - 
     -         2*WProd1931*WProd1988 + WProd1903*WProd1997 + 
     -         WProd1905*WProd1997 + WProd1906*WProd1997 + 
     -         2*WProd1910*WProd1997 - 2*WProd1926*WProd1997 + 
     -         2*WProd1927*WProd1997 + 2*WProd1928*WProd1997 - 
     -         2*WProd1929*WProd1997 - 2*WProd1930*WProd1997 - 
     -         2*WProd1931*WProd1997 + WProd1998 + WProd1999))/
     -     WProd1989 + Prp223*
     -    (-(2.6666666666666665D0*
     -         (Cpl2*GS*Prp222*
     -            (WProd1899*WProd1995 + WProd1907*WProd1995 - 
     -              2*WProd1903*WProd1997 - 
     -              2*WProd1905*WProd1997 - 
     -              2*WProd1906*WProd1997 - 
     -              2*WProd1910*WProd1997 + 
     -              2*WProd1926*WProd1997 - 
     -              2*WProd1927*WProd1997 - 
     -              2*WProd1928*WProd1997 + 
     -              2*WProd1929*WProd1997 + 
     -              2*WProd1930*WProd1997 + 
     -              2*WProd1931*WProd1997 - WProd1998 - 
     -              WProd1999 + WProd1913*WProd2000 + 
     -              WProd1913*WProd2001))/WProd1989) + 
     -      Prp218*(-(2.6666666666666665D0*
     -            (Cpl2*GS*(WProd2002 + WProd2003))/WProd1976) - 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*(-WProd1981 + WProd1913*WProd1991 + 
     -               WProd1899*WProd1995 + WProd1901*WProd1995 - 
     -               WProd1915*WProd1919*WProd1995 - 
     -               WProd1915*WProd1920*WProd1995 + 
     -               WProd1922*WProd1995 + WProd1923*WProd1995 + 
     -               WProd1924*WProd1995 - WProd1903*WProd1997 - 
     -               WProd1905*WProd1997 - WProd1906*WProd1997 - 
     -               WProd1910*WProd1997 + 2*WProd1926*WProd1997 - 
     -               WProd1927*WProd1997 - WProd1928*WProd1997 + 
     -               WProd1929*WProd1997 + WProd1930*WProd1997 + 
     -               WProd1931*WProd1997 - WProd1998 + 
     -               WProd1913*WProd2000 + WProd1913*WProd2004 - 
     -               WProd2005 - WProd1919*WProd2006 - 
     -               WProd1920*WProd2006 + WProd1903*WProd2007 + 
     -               WProd1905*WProd2007 + WProd1906*WProd2007 + 
     -               WProd1910*WProd2007 - 2*WProd1926*WProd2007 + 
     -               WProd1927*WProd2007 + WProd1928*WProd2007 - 
     -               WProd1929*WProd2007 - WProd1930*WProd2007 - 
     -               WProd1931*WProd2007))/WProd1989))
        tmp1154 = tmp1154 - 
     -   (4*Cpl2*GS*Prp217*Prp218*
     -      (2*WProd1977 - 2*WProd1896*WProd1978 + 
     -        2*WProd1899*WProd1978 + 2*WProd1901*WProd1978 + 
     -        2*WProd1913*WProd1915*WProd1978 - 
     -        2*WProd1915*WProd1919*WProd1978 - 
     -        2*WProd1915*WProd1920*WProd1978 - 
     -        2*WProd1921*WProd1978 + 2*WProd1922*WProd1978 + 
     -        2*WProd1923*WProd1978 + 2*WProd1924*WProd1978 + 
     -        2*WProd1979 - 2*WProd1981 + 2*WProd1982 + 
     -        2*WProd1983 - 2*WProd1896*WProd1986 + 
     -        2*WProd1899*WProd1986 + 2*WProd1901*WProd1986 + 
     -        2*WProd1913*WProd1915*WProd1986 - 
     -        2*WProd1915*WProd1919*WProd1986 - 
     -        2*WProd1915*WProd1920*WProd1986 - 
     -        2*WProd1921*WProd1986 + 2*WProd1922*WProd1986 + 
     -        2*WProd1923*WProd1986 + 2*WProd1924*WProd1986 + 
     -        2*WProd1903*WProd1988 + 2*WProd1905*WProd1988 + 
     -        3*WProd1906*WProd1988 + 3*WProd1910*WProd1988 - 
     -        4*WProd1926*WProd1988 + 2*WProd1927*WProd1988 + 
     -        2*WProd1928*WProd1988 - 2*WProd1929*WProd1988 - 
     -        3*WProd1930*WProd1988 - 3*WProd1931*WProd1988 + 
     -        2*WProd1919*WProd1990 + 2*WProd1920*WProd1990 + 
     -        2*WProd1913*WProd1991 - 2*WProd1919*WProd1991 - 
     -        2*WProd1920*WProd1991 + 2*WProd1901*WProd1995 - 
     -        2*WProd1913*WProd2008 - 2*WProd1913*WProd2009))/
     -    WProd1989
        tmp1154 = tmp1154 + 
     -   Prp221*(Prp220*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2002 + WProd2003))/WProd1976 - 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1977 - WProd1896*WProd1978 + 
     -               WProd1901*WProd1978 - WProd1921*WProd1978 + 
     -               WProd1979 - WProd1890*WProd1980 - 
     -               WProd1891*WProd1980 - 
     -               WProd1906*WProd1915*WProd1980 - 
     -               WProd1910*WProd1915*WProd1980 + 
     -               WProd1915*WProd1926*WProd1980 - 
     -               WProd1915*WProd1927*WProd1980 - 
     -               WProd1915*WProd1928*WProd1980 + 
     -               WProd1915*WProd1929*WProd1980 + 
     -               WProd1915*WProd1930*WProd1980 + 
     -               WProd1915*WProd1931*WProd1980 - 
     -               WProd1933*WProd1980 - WProd1944*WProd1980 + 
     -               WProd1982 + WProd1983 - WProd1896*WProd1986 + 
     -               WProd1901*WProd1986 - WProd1921*WProd1986 - 
     -               WProd1899*WProd1987 + 
     -               WProd1915*WProd1919*WProd1987 + 
     -               WProd1915*WProd1920*WProd1987 - 
     -               WProd1922*WProd1987 - WProd1923*WProd1987 - 
     -               WProd1924*WProd1987 + WProd1906*WProd1988 + 
     -               WProd1910*WProd1988 - WProd1926*WProd1988 + 
     -               WProd1927*WProd1988 + WProd1928*WProd1988 - 
     -               WProd1929*WProd1988 - WProd1930*WProd1988 - 
     -               WProd1931*WProd1988 + WProd1919*WProd1990 + 
     -               WProd1920*WProd1990 - WProd1919*WProd1991 - 
     -               WProd1920*WProd1991 + WProd1906*WProd1997 + 
     -               WProd1910*WProd1997 - WProd1926*WProd1997 + 
     -               WProd1927*WProd1997 + WProd1928*WProd1997 - 
     -               WProd1929*WProd1997 - WProd1930*WProd1997 - 
     -               WProd1931*WProd1997 + WProd1998 + WProd2005 + 
     -               WProd1919*WProd2006 + WProd1920*WProd2006 - 
     -               WProd1906*WProd2007 - WProd1910*WProd2007 + 
     -               WProd1926*WProd2007 - WProd1927*WProd2007 - 
     -               WProd1928*WProd2007 + WProd1929*WProd2007 + 
     -               WProd1930*WProd2007 + WProd1931*WProd2007))/
     -           WProd1989) + 
     -      Prp218*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2002 + WProd2003))/WProd1976 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-WProd1977 + WProd1896*WProd1978 - 
     -               WProd1899*WProd1978 - WProd1901*WProd1978 - 
     -               WProd1913*WProd1915*WProd1978 + 
     -               WProd1915*WProd1919*WProd1978 + 
     -               WProd1915*WProd1920*WProd1978 + 
     -               WProd1921*WProd1978 - WProd1922*WProd1978 - 
     -               WProd1923*WProd1978 - WProd1924*WProd1978 - 
     -               WProd1979 - WProd1982 - WProd1983 + 
     -               WProd1896*WProd1986 - WProd1899*WProd1986 - 
     -               WProd1901*WProd1986 - 
     -               WProd1913*WProd1915*WProd1986 + 
     -               WProd1915*WProd1919*WProd1986 + 
     -               WProd1915*WProd1920*WProd1986 + 
     -               WProd1921*WProd1986 - WProd1922*WProd1986 - 
     -               WProd1923*WProd1986 - WProd1924*WProd1986 - 
     -               WProd1903*WProd1988 - WProd1905*WProd1988 - 
     -               2*WProd1906*WProd1988 - 
     -               2*WProd1910*WProd1988 + 
     -               2*WProd1926*WProd1988 - WProd1927*WProd1988 - 
     -               WProd1928*WProd1988 + WProd1929*WProd1988 + 
     -               2*WProd1930*WProd1988 + 
     -               2*WProd1931*WProd1988 - WProd1919*WProd1990 - 
     -               WProd1920*WProd1990 + WProd1919*WProd1991 + 
     -               WProd1920*WProd1991 + WProd1899*WProd1995 - 
     -               WProd1915*WProd1919*WProd1995 - 
     -               WProd1915*WProd1920*WProd1995 + 
     -               WProd1922*WProd1995 + WProd1923*WProd1995 + 
     -               WProd1924*WProd1995 - WProd1903*WProd1997 - 
     -               WProd1905*WProd1997 - 2*WProd1906*WProd1997 - 
     -               2*WProd1910*WProd1997 + 
     -               2*WProd1926*WProd1997 - WProd1927*WProd1997 - 
     -               WProd1928*WProd1997 + WProd1929*WProd1997 + 
     -               2*WProd1930*WProd1997 + 
     -               2*WProd1931*WProd1997 - WProd1998 + 
     -               WProd1913*WProd2000 + WProd1913*WProd2004 - 
     -               WProd2005 - WProd1919*WProd2006 - 
     -               WProd1920*WProd2006 + WProd1903*WProd2007 + 
     -               WProd1905*WProd2007 + 2*WProd1906*WProd2007 + 
     -               2*WProd1910*WProd2007 - 
     -               2*WProd1926*WProd2007 + WProd1927*WProd2007 + 
     -               WProd1928*WProd2007 - WProd1929*WProd2007 - 
     -               2*WProd1930*WProd2007 - 
     -               2*WProd1931*WProd2007 + WProd1913*WProd2008 + 
     -               WProd1913*WProd2009))/WProd1989))
        amp(26,7,10) = Prp219*tmp1154 + 
     -   Prp215*Prp216*Prp217*
     -    ((-8*Cpl2*GS*WProd1975)/WProd1976 - 
     -      (8*Cpl2*GS*(WProd1977 - WProd1896*WProd1978 + 
     -           WProd1899*WProd1978 + WProd1901*WProd1978 + 
     -           WProd1907*WProd1978 + WProd1979 + 
     -           WProd1888*WProd1980 - WProd1890*WProd1980 - 
     -           WProd1891*WProd1980 - WProd1981 + WProd1982 + 
     -           WProd1983 + WProd1984 + WProd1985 + 
     -           WProd1901*WProd1986 + WProd1901*WProd1987 - 
     -           WProd1903*WProd1988 - WProd1905*WProd1988 - 
     -           WProd1906*WProd1988))/WProd1989) + 
     -   (8*Cpl2*GS*Prp215*Prp217*Prp218*
     -      (-WProd1977 + WProd1896*WProd1978 - 
     -        WProd1899*WProd1978 - WProd1901*WProd1978 - 
     -        WProd1913*WProd1915*WProd1978 + 
     -        WProd1915*WProd1919*WProd1978 + 
     -        WProd1915*WProd1920*WProd1978 + 
     -        WProd1921*WProd1978 - WProd1922*WProd1978 - 
     -        WProd1923*WProd1978 - WProd1924*WProd1978 - 
     -        WProd1979 - WProd1982 - WProd1983 - WProd1984 - 
     -        WProd1985 - WProd1901*WProd1986 + 
     -        2*WProd1903*WProd1988 + WProd1905*WProd1988 + 
     -        WProd1906*WProd1988 + WProd1910*WProd1988 - 
     -        WProd1926*WProd1988 + WProd1927*WProd1988 + 
     -        WProd1928*WProd1988 - WProd1929*WProd1988 - 
     -        WProd1930*WProd1988 - WProd1931*WProd1988 + 
     -        WProd1913*WProd1990 - WProd1919*WProd1990 - 
     -        WProd1920*WProd1990 - WProd1913*WProd1991 + 
     -        WProd1919*WProd1991 + WProd1920*WProd1991 - 
     -        WProd1913*WProd1992 + WProd1919*WProd1992 + 
     -        WProd1920*WProd1992 + WProd1993 + WProd1994 - 
     -        WProd1901*WProd1995 + WProd1996))/WProd1989

#ifdef DEBUG
	print *, 'amp(26,7,10) =', amp(26,7,10)
#endif

        tmp1155 = (-4*Cpl2*GS*Prp217*Prp218*
     -      (-2*WProd1952*WProd1978 + 2*WProd1955*WProd1978 + 
     -        2*WProd1956*WProd1978 + 
     -        2*WProd1915*WProd1959*WProd1978 - 
     -        2*WProd1915*WProd1963*WProd1978 - 
     -        2*WProd1915*WProd1964*WProd1978 - 
     -        2*WProd1965*WProd1978 + 2*WProd1966*WProd1978 + 
     -        2*WProd1967*WProd1978 + 2*WProd1968*WProd1978 - 
     -        2*WProd1952*WProd1986 + 2*WProd1955*WProd1986 + 
     -        2*WProd1956*WProd1986 + 
     -        2*WProd1915*WProd1959*WProd1986 - 
     -        2*WProd1915*WProd1963*WProd1986 - 
     -        2*WProd1915*WProd1964*WProd1986 - 
     -        2*WProd1965*WProd1986 + 2*WProd1966*WProd1986 + 
     -        2*WProd1967*WProd1986 + 2*WProd1968*WProd1986 + 
     -        2*WProd1963*WProd1990 + 2*WProd1964*WProd1990 + 
     -        2*WProd1959*WProd1991 - 2*WProd1963*WProd1991 - 
     -        2*WProd1964*WProd1991 + 2*WProd1956*WProd1995 - 
     -        2*WProd1959*WProd2008 - 2*WProd1959*WProd2009 + 
     -        2*WProd2013 + 2*WProd2014 - 2*WProd2016 + 
     -        2*WProd2017 + 2*WProd2018 + 2*WProd1903*WProd2021 + 
     -        2*WProd1905*WProd2021 + 3*WProd1906*WProd2021 + 
     -        3*WProd1910*WProd2021 - 4*WProd1926*WProd2021 + 
     -        2*WProd1927*WProd2021 + 2*WProd1928*WProd2021 - 
     -        2*WProd1929*WProd2021 - 3*WProd1930*WProd2021 - 
     -        3*WProd1931*WProd2021))/WProd2012 + 
     -   Prp216*((-4*Cpl2*GS*Prp217*
     -         (-2*WProd1952*WProd1978 + 2*WProd1955*WProd1978 + 
     -           2*WProd1956*WProd1978 + 2*WProd1958*WProd1978 - 
     -           2*WProd1952*WProd1986 + 2*WProd1955*WProd1986 + 
     -           2*WProd1956*WProd1986 + 2*WProd1958*WProd1986 + 
     -           2*WProd1956*WProd1987 + 2*WProd2013 + 
     -           2*WProd2014 + 2*WProd1888*WProd2015 - 
     -           2*WProd1890*WProd2015 - 2*WProd1891*WProd2015 - 
     -           2*WProd2016 + 2*WProd2017 + 2*WProd2018 + 
     -           3*WProd1903*WProd2021 + 3*WProd1905*WProd2021 + 
     -           3*WProd1906*WProd2021 + 4*WProd1910*WProd2021 - 
     -           4*WProd1926*WProd2021 + 4*WProd1927*WProd2021 + 
     -           4*WProd1928*WProd2021 - 5*WProd1929*WProd2021 - 
     -           5*WProd1930*WProd2021 - 5*WProd1931*WProd2021))/
     -       WProd2012 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp222*
     -          (-(WProd1955*WProd1978) - WProd1958*WProd1978 - 
     -            WProd1955*WProd1986 - WProd1958*WProd1986 - 
     -            WProd1955*WProd1987 - WProd1958*WProd1987 - 
     -            2*WProd1903*WProd1915*WProd2015 - 
     -            2*WProd1905*WProd1915*WProd2015 - 
     -            2*WProd1906*WProd1915*WProd2015 - 
     -            2*WProd1910*WProd1915*WProd2015 + 
     -            2*WProd1915*WProd1926*WProd2015 - 
     -            2*WProd1915*WProd1927*WProd2015 - 
     -            2*WProd1915*WProd1928*WProd2015 + 
     -            3*WProd1915*WProd1929*WProd2015 + 
     -            3*WProd1915*WProd1930*WProd2015 + 
     -            3*WProd1915*WProd1931*WProd2015 - 
     -            WProd1933*WProd2015 - WProd1934*WProd2015 + 
     -            2*WProd1903*WProd2025 + 2*WProd1905*WProd2025 + 
     -            2*WProd1906*WProd2025 + 2*WProd1910*WProd2025 - 
     -            2*WProd1926*WProd2025 + 2*WProd1927*WProd2025 + 
     -            2*WProd1928*WProd2025 - 3*WProd1929*WProd2025 - 
     -            3*WProd1930*WProd2025 - 3*WProd1931*WProd2025 + 
     -            WProd2026 + WProd2027))/WProd2012) - 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp216*Prp220*
     -       (-(WProd1952*WProd1978) + WProd1956*WProd1978 - 
     -         WProd1952*WProd1986 + WProd1956*WProd1986 - 
     -         WProd1955*WProd1987 + WProd1956*WProd1987 - 
     -         WProd1958*WProd1987 + WProd2013 + WProd2014 + 
     -         WProd1888*WProd2015 - WProd1890*WProd2015 - 
     -         WProd1891*WProd2015 - 
     -         WProd1903*WProd1915*WProd2015 - 
     -         WProd1905*WProd1915*WProd2015 - 
     -         WProd1906*WProd1915*WProd2015 - 
     -         2*WProd1910*WProd1915*WProd2015 + 
     -         2*WProd1915*WProd1926*WProd2015 - 
     -         2*WProd1915*WProd1927*WProd2015 - 
     -         2*WProd1915*WProd1928*WProd2015 + 
     -         2*WProd1915*WProd1929*WProd2015 + 
     -         2*WProd1915*WProd1930*WProd2015 + 
     -         2*WProd1915*WProd1931*WProd2015 - 
     -         WProd1933*WProd2015 - WProd1934*WProd2015 - 
     -         WProd2016 + WProd2017 + WProd2018 + 
     -         WProd1903*WProd2021 + WProd1905*WProd2021 + 
     -         WProd1906*WProd2021 + 2*WProd1910*WProd2021 - 
     -         2*WProd1926*WProd2021 + 2*WProd1927*WProd2021 + 
     -         2*WProd1928*WProd2021 - 2*WProd1929*WProd2021 - 
     -         2*WProd1930*WProd2021 - 2*WProd1931*WProd2021 + 
     -         WProd1903*WProd2025 + WProd1905*WProd2025 + 
     -         WProd1906*WProd2025 + 2*WProd1910*WProd2025 - 
     -         2*WProd1926*WProd2025 + 2*WProd1927*WProd2025 + 
     -         2*WProd1928*WProd2025 - 2*WProd1929*WProd2025 - 
     -         2*WProd1930*WProd2025 - 2*WProd1931*WProd2025 + 
     -         WProd2026 + WProd2027))/WProd2012
        tmp1155 = tmp1155 + 
     -   Prp223*(-(2.6666666666666665D0*
     -         (Cpl2*GS*Prp222*
     -            (WProd1955*WProd1995 + WProd1958*WProd1995 + 
     -              WProd1959*WProd2000 + WProd1959*WProd2001 - 
     -              2*WProd1903*WProd2025 - 
     -              2*WProd1905*WProd2025 - 
     -              2*WProd1906*WProd2025 - 
     -              2*WProd1910*WProd2025 + 
     -              2*WProd1926*WProd2025 - 
     -              2*WProd1927*WProd2025 - 
     -              2*WProd1928*WProd2025 + 
     -              2*WProd1929*WProd2025 + 
     -              2*WProd1930*WProd2025 + 
     -              2*WProd1931*WProd2025 - WProd2026 - WProd2027))
     -           /WProd2012) + 
     -      Prp218*(-(2.6666666666666665D0*
     -            (Cpl2*GS*(WProd2028 + WProd2029))/WProd2010) - 
     -         2.6666666666666665D0*
     -          (Cpl2*GS*(WProd1959*WProd1991 + 
     -               WProd1955*WProd1995 + WProd1956*WProd1995 - 
     -               WProd1915*WProd1963*WProd1995 - 
     -               WProd1915*WProd1964*WProd1995 + 
     -               WProd1966*WProd1995 + WProd1967*WProd1995 + 
     -               WProd1968*WProd1995 + WProd1959*WProd2000 + 
     -               WProd1959*WProd2004 - WProd1963*WProd2006 - 
     -               WProd1964*WProd2006 - WProd2016 - 
     -               WProd1903*WProd2025 - WProd1905*WProd2025 - 
     -               WProd1906*WProd2025 - WProd1910*WProd2025 + 
     -               2*WProd1926*WProd2025 - WProd1927*WProd2025 - 
     -               WProd1928*WProd2025 + WProd1929*WProd2025 + 
     -               WProd1930*WProd2025 + WProd1931*WProd2025 - 
     -               WProd2026 - WProd2030 + WProd1903*WProd2031 + 
     -               WProd1905*WProd2031 + WProd1906*WProd2031 + 
     -               WProd1910*WProd2031 - 2*WProd1926*WProd2031 + 
     -               WProd1927*WProd2031 + WProd1928*WProd2031 - 
     -               WProd1929*WProd2031 - WProd1930*WProd2031 - 
     -               WProd1931*WProd2031))/WProd2012))
        tmp1155 = tmp1155 + 
     -   Prp221*(Prp218*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2028 + WProd2029))/WProd2010 + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(WProd1952*WProd1978 - 
     -               WProd1955*WProd1978 - WProd1956*WProd1978 - 
     -               WProd1915*WProd1959*WProd1978 + 
     -               WProd1915*WProd1963*WProd1978 + 
     -               WProd1915*WProd1964*WProd1978 + 
     -               WProd1965*WProd1978 - WProd1966*WProd1978 - 
     -               WProd1967*WProd1978 - WProd1968*WProd1978 + 
     -               WProd1952*WProd1986 - WProd1955*WProd1986 - 
     -               WProd1956*WProd1986 - 
     -               WProd1915*WProd1959*WProd1986 + 
     -               WProd1915*WProd1963*WProd1986 + 
     -               WProd1915*WProd1964*WProd1986 + 
     -               WProd1965*WProd1986 - WProd1966*WProd1986 - 
     -               WProd1967*WProd1986 - WProd1968*WProd1986 - 
     -               WProd1963*WProd1990 - WProd1964*WProd1990 + 
     -               WProd1963*WProd1991 + WProd1964*WProd1991 + 
     -               WProd1955*WProd1995 - 
     -               WProd1915*WProd1963*WProd1995 - 
     -               WProd1915*WProd1964*WProd1995 + 
     -               WProd1966*WProd1995 + WProd1967*WProd1995 + 
     -               WProd1968*WProd1995 + WProd1959*WProd2000 + 
     -               WProd1959*WProd2004 - WProd1963*WProd2006 - 
     -               WProd1964*WProd2006 + WProd1959*WProd2008 + 
     -               WProd1959*WProd2009 - WProd2013 - WProd2014 - 
     -               WProd2017 - WProd2018 - WProd1903*WProd2021 - 
     -               WProd1905*WProd2021 - 2*WProd1906*WProd2021 - 
     -               2*WProd1910*WProd2021 + 
     -               2*WProd1926*WProd2021 - WProd1927*WProd2021 - 
     -               WProd1928*WProd2021 + WProd1929*WProd2021 + 
     -               2*WProd1930*WProd2021 + 
     -               2*WProd1931*WProd2021 - WProd1903*WProd2025 - 
     -               WProd1905*WProd2025 - 2*WProd1906*WProd2025 - 
     -               2*WProd1910*WProd2025 + 
     -               2*WProd1926*WProd2025 - WProd1927*WProd2025 - 
     -               WProd1928*WProd2025 + WProd1929*WProd2025 + 
     -               2*WProd1930*WProd2025 + 
     -               2*WProd1931*WProd2025 - WProd2026 - 
     -               WProd2030 + WProd1903*WProd2031 + 
     -               WProd1905*WProd2031 + 2*WProd1906*WProd2031 + 
     -               2*WProd1910*WProd2031 - 
     -               2*WProd1926*WProd2031 + WProd1927*WProd2031 + 
     -               WProd1928*WProd2031 - WProd1929*WProd2031 - 
     -               2*WProd1930*WProd2031 - 2*WProd1931*WProd2031)
     -             )/WProd2012) + 
     -      Prp220*(5.333333333333333D0*
     -          (Cpl2*GS*(WProd2028 + WProd2029))/WProd2010 - 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd1952*WProd1978) + 
     -               WProd1956*WProd1978 - WProd1965*WProd1978 - 
     -               WProd1952*WProd1986 + WProd1956*WProd1986 - 
     -               WProd1965*WProd1986 - WProd1955*WProd1987 + 
     -               WProd1915*WProd1963*WProd1987 + 
     -               WProd1915*WProd1964*WProd1987 - 
     -               WProd1966*WProd1987 - WProd1967*WProd1987 - 
     -               WProd1968*WProd1987 + WProd1963*WProd1990 + 
     -               WProd1964*WProd1990 - WProd1963*WProd1991 - 
     -               WProd1964*WProd1991 + WProd1963*WProd2006 + 
     -               WProd1964*WProd2006 + WProd2013 + WProd2014 - 
     -               WProd1890*WProd2015 - WProd1891*WProd2015 - 
     -               WProd1906*WProd1915*WProd2015 - 
     -               WProd1910*WProd1915*WProd2015 + 
     -               WProd1915*WProd1926*WProd2015 - 
     -               WProd1915*WProd1927*WProd2015 - 
     -               WProd1915*WProd1928*WProd2015 + 
     -               WProd1915*WProd1929*WProd2015 + 
     -               WProd1915*WProd1930*WProd2015 + 
     -               WProd1915*WProd1931*WProd2015 - 
     -               WProd1933*WProd2015 - WProd1944*WProd2015 + 
     -               WProd2017 + WProd2018 + WProd1906*WProd2021 + 
     -               WProd1910*WProd2021 - WProd1926*WProd2021 + 
     -               WProd1927*WProd2021 + WProd1928*WProd2021 - 
     -               WProd1929*WProd2021 - WProd1930*WProd2021 - 
     -               WProd1931*WProd2021 + WProd1906*WProd2025 + 
     -               WProd1910*WProd2025 - WProd1926*WProd2025 + 
     -               WProd1927*WProd2025 + WProd1928*WProd2025 - 
     -               WProd1929*WProd2025 - WProd1930*WProd2025 - 
     -               WProd1931*WProd2025 + WProd2026 + WProd2030 - 
     -               WProd1906*WProd2031 - WProd1910*WProd2031 + 
     -               WProd1926*WProd2031 - WProd1927*WProd2031 - 
     -               WProd1928*WProd2031 + WProd1929*WProd2031 + 
     -               WProd1930*WProd2031 + WProd1931*WProd2031))/
     -           WProd2012))
        amp(26,7,5) = Prp219*tmp1155 + 
     -   Prp215*Prp216*Prp217*
     -    ((-8*Cpl2*GS*WProd2011)/WProd2010 - 
     -      (8*Cpl2*GS*(-(WProd1952*WProd1978) + 
     -           WProd1955*WProd1978 + WProd1956*WProd1978 + 
     -           WProd1958*WProd1978 + WProd1956*WProd1986 + 
     -           WProd1956*WProd1987 + WProd2013 + WProd2014 + 
     -           WProd1888*WProd2015 - WProd1890*WProd2015 - 
     -           WProd1891*WProd2015 - WProd2016 + WProd2017 + 
     -           WProd2018 + WProd2019 + WProd2020 - 
     -           WProd1903*WProd2021 - WProd1905*WProd2021 - 
     -           WProd1906*WProd2021))/WProd2012) + 
     -   (8*Cpl2*GS*Prp215*Prp217*Prp218*
     -      (WProd1952*WProd1978 - WProd1955*WProd1978 - 
     -        WProd1956*WProd1978 - 
     -        WProd1915*WProd1959*WProd1978 + 
     -        WProd1915*WProd1963*WProd1978 + 
     -        WProd1915*WProd1964*WProd1978 + 
     -        WProd1965*WProd1978 - WProd1966*WProd1978 - 
     -        WProd1967*WProd1978 - WProd1968*WProd1978 - 
     -        WProd1956*WProd1986 + WProd1959*WProd1990 - 
     -        WProd1963*WProd1990 - WProd1964*WProd1990 - 
     -        WProd1959*WProd1991 + WProd1963*WProd1991 + 
     -        WProd1964*WProd1991 - WProd1959*WProd1992 + 
     -        WProd1963*WProd1992 + WProd1964*WProd1992 - 
     -        WProd1956*WProd1995 - WProd2013 - WProd2014 - 
     -        WProd2017 - WProd2018 - WProd2019 - WProd2020 + 
     -        2*WProd1903*WProd2021 + WProd1905*WProd2021 + 
     -        WProd1906*WProd2021 + WProd1910*WProd2021 - 
     -        WProd1926*WProd2021 + WProd1927*WProd2021 + 
     -        WProd1928*WProd2021 - WProd1929*WProd2021 - 
     -        WProd1930*WProd2021 - WProd1931*WProd2021 + 
     -        WProd2022 + WProd2023 + WProd2024))/WProd2012

#ifdef DEBUG
	print *, 'amp(26,7,5) =', amp(26,7,5)
#endif

      endif

      end


