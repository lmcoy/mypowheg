      subroutine realQED11(amp,resnr)
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

        tmp1130 = (4*Cpl2*GS*Prp127*Prp129*
     -      (-2*WProd929*WProd930 + 2*WProd930*WProd932 - 
     -        2*WProd929*WProd934 + 2*WProd932*WProd934 - 
     -        2*WProd923*WProd937 - 2*WProd925*WProd937 + 
     -        2*WProd926*WProd937 + 2*WProd930*WProd931*WProd937 + 
     -        2*WProd931*WProd934*WProd937 - 2*WProd935*WProd937 + 
     -        4*WProd939*WProd940 - 2*WProd941 - 2*WProd944 - 
     -        2*WProd946 - 4*WProd940*WProd948 + 
     -        3*WProd940*WProd949 - 2*WProd930*WProd951 - 
     -        2*WProd934*WProd951 - 2*WProd932*WProd952 + 
     -        2*WProd951*WProd952 - 2*WProd930*WProd953 - 
     -        2*WProd934*WProd953 + 2*WProd952*WProd953 - 
     -        2*WProd932*WProd954 + 2*WProd951*WProd954 + 
     -        2*WProd953*WProd954 + 4*WProd940*WProd955 + 
     -        4*WProd940*WProd956 + 4*WProd940*WProd957 + 
     -        4*WProd940*WProd958 - 4*WProd940*WProd959 - 
     -        2*WProd940*WProd960 - 2*WProd940*WProd961 + 
     -        2*WProd962 - 2*WProd937*WProd963 - 2*WProd964))/
     -    WProd947 + 5.333333333333333D0*
     -    (Cpl2*GS*Prp129*Prp131*
     -       (-(WProd929*WProd930) + WProd930*WProd932 - 
     -         WProd929*WProd934 + WProd932*WProd934 - 
     -         WProd925*WProd937 + WProd926*WProd937 + 
     -         WProd930*WProd931*WProd937 + 
     -         WProd931*WProd934*WProd937 - WProd935*WProd937 + 
     -         2*WProd939*WProd940 - WProd941 - WProd946 - 
     -         2*WProd940*WProd948 + WProd940*WProd949 - 
     -         WProd932*WProd952 - WProd930*WProd953 - 
     -         WProd934*WProd953 + WProd952*WProd953 - 
     -         WProd932*WProd954 + WProd953*WProd954 + 
     -         2*WProd940*WProd955 + 2*WProd940*WProd956 + 
     -         2*WProd940*WProd957 + 2*WProd940*WProd958 - 
     -         2*WProd940*WProd959 - WProd940*WProd960 - 
     -         WProd940*WProd961 + WProd962 - WProd937*WProd963 - 
     -         WProd964 - WProd930*WProd965 - WProd934*WProd965 + 
     -         WProd952*WProd965 + WProd954*WProd965 - 
     -         WProd937*WProd966 - 2*WProd939*WProd967 + 
     -         2*WProd948*WProd967 - WProd949*WProd967 - 
     -         2*WProd955*WProd967 - 2*WProd956*WProd967 - 
     -         2*WProd957*WProd967 - 2*WProd958*WProd967 + 
     -         2*WProd959*WProd967 + WProd960*WProd967 + 
     -         WProd961*WProd967 + WProd968 + WProd969))/WProd947+
     -     Prp132*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp129*
     -          (-(WProd923*WProd937) - WProd944 - 
     -            WProd930*WProd951 - WProd934*WProd951 + 
     -            WProd951*WProd952 + WProd951*WProd954 + 
     -            WProd930*WProd965 + WProd934*WProd965 - 
     -            WProd952*WProd965 - WProd954*WProd965 + 
     -            WProd937*WProd966 + 2*WProd939*WProd967 - 
     -            2*WProd948*WProd967 + 2*WProd949*WProd967 + 
     -            2*WProd955*WProd967 + 2*WProd956*WProd967 + 
     -            2*WProd957*WProd967 + 2*WProd958*WProd967 - 
     -            2*WProd959*WProd967 - WProd960*WProd967 - 
     -            WProd961*WProd967 - WProd968 - WProd969))/
     -        WProd947 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (WProd922*WProd923 - WProd923*WProd938 - 
     -            WProd944 - WProd922*WProd966 + 
     -            WProd938*WProd966 + 3*WProd939*WProd967 - 
     -            3*WProd948*WProd967 + 3*WProd949*WProd967 + 
     -            2*WProd955*WProd967 + 3*WProd956*WProd967 + 
     -            2*WProd957*WProd967 + 3*WProd958*WProd967 - 
     -            2*WProd959*WProd967 - 2*WProd960*WProd967 - 
     -            2*WProd961*WProd967 - WProd968 - WProd969 - 
     -            WProd923*WProd971 + WProd966*WProd971 - 
     -            3*WProd931*WProd939*WProd972 + 
     -            3*WProd931*WProd948*WProd972 - 
     -            3*WProd931*WProd949*WProd972 - 
     -            2*WProd931*WProd955*WProd972 - 
     -            3*WProd931*WProd956*WProd972 - 
     -            2*WProd931*WProd957*WProd972 - 
     -            3*WProd931*WProd958*WProd972 + 
     -            2*WProd931*WProd959*WProd972 + 
     -            2*WProd931*WProd960*WProd972 + 
     -            2*WProd931*WProd961*WProd972 + 
     -            WProd972*WProd973 + WProd972*WProd974 + 
     -            WProd972*WProd975))/WProd947)
        tmp1130 = tmp1130 + 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp131*Prp134*
     -       (-(WProd929*WProd930) + WProd930*WProd932 - 
     -         WProd930*WProd933 - WProd929*WProd934 + 
     -         WProd932*WProd934 - WProd933*WProd934 - 
     -         WProd925*WProd937 + WProd930*WProd931*WProd937 + 
     -         WProd931*WProd934*WProd937 - WProd935*WProd937 + 
     -         WProd925*WProd938 - WProd926*WProd938 + 
     -         WProd935*WProd938 + WProd936*WProd938 + 
     -         2*WProd939*WProd940 - WProd941 - WProd946 - 
     -         2*WProd940*WProd948 + 2*WProd940*WProd949 - 
     -         WProd931*WProd938*WProd952 - WProd930*WProd953 - 
     -         WProd934*WProd953 + WProd952*WProd953 - 
     -         WProd931*WProd938*WProd954 + WProd953*WProd954 + 
     -         2*WProd940*WProd955 + 3*WProd940*WProd956 + 
     -         2*WProd940*WProd957 + 3*WProd940*WProd958 - 
     -         3*WProd940*WProd959 - 2*WProd940*WProd960 - 
     -         2*WProd940*WProd961 - WProd937*WProd963 + 
     -         WProd938*WProd963 - WProd964 - WProd930*WProd965 - 
     -         WProd934*WProd965 + WProd952*WProd965 + 
     -         WProd954*WProd965 - WProd937*WProd966 + 
     -         WProd938*WProd966 - 2*WProd939*WProd967 + 
     -         2*WProd948*WProd967 - 2*WProd949*WProd967 - 
     -         2*WProd955*WProd967 - 3*WProd956*WProd967 - 
     -         2*WProd957*WProd967 - 3*WProd958*WProd967 + 
     -         3*WProd959*WProd967 + 2*WProd960*WProd967 + 
     -         2*WProd961*WProd967 + WProd968 + 
     -         WProd925*WProd971 - WProd926*WProd971 + 
     -         WProd935*WProd971 + WProd936*WProd971 - 
     -         WProd931*WProd952*WProd971 - 
     -         WProd931*WProd954*WProd971 + WProd963*WProd971 + 
     -         WProd966*WProd971 - WProd930*WProd977 - 
     -         WProd934*WProd977 + WProd952*WProd977 + 
     -         WProd954*WProd977 - WProd930*WProd978 - 
     -         WProd934*WProd978 + WProd952*WProd978 + 
     -         WProd954*WProd978 - 2*WProd939*WProd979 + 
     -         2*WProd948*WProd979 - 2*WProd949*WProd979 - 
     -         2*WProd955*WProd979 - 3*WProd956*WProd979 - 
     -         2*WProd957*WProd979 - 3*WProd958*WProd979 + 
     -         3*WProd959*WProd979 + 2*WProd960*WProd979 + 
     -         2*WProd961*WProd979 + WProd980 + WProd981 + 
     -         WProd938*WProd982 + WProd971*WProd982 - WProd983))/
     -     WProd947
        tmp1130 = tmp1130 + 
     -   Prp126*((4*Cpl2*GS*Prp127*
     -         (2*WProd922*WProd923 + 2*WProd922*WProd925 - 
     -           2*WProd922*WProd926 - 2*WProd929*WProd930 - 
     -           2*WProd922*WProd930*WProd931 + 
     -           2*WProd930*WProd932 - 2*WProd930*WProd933 - 
     -           2*WProd929*WProd934 - 
     -           2*WProd922*WProd931*WProd934 + 
     -           2*WProd932*WProd934 - 2*WProd933*WProd934 + 
     -           2*WProd922*WProd935 + 2*WProd936*WProd938 + 
     -           2*WProd939*WProd940 - 2*WProd941 - 2*WProd944 - 
     -           2*WProd946 - WProd940*WProd948 + 
     -           WProd940*WProd949 + WProd940*WProd955 + 
     -           WProd940*WProd956 + WProd940*WProd957 + 
     -           WProd940*WProd958 - WProd940*WProd959 - 
     -           WProd940*WProd960 - WProd940*WProd961 + 
     -           2*WProd962 + 2*WProd922*WProd963 - 2*WProd964 + 
     -           2*WProd936*WProd971 + 2*WProd972*WProd973 + 
     -           2*WProd972*WProd985 - 2*WProd972*WProd986))/
     -       WProd947 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (WProd922*WProd923 - WProd922*WProd926 + 
     -            WProd926*WProd938 - WProd944 + WProd962 - 
     -            WProd922*WProd966 + WProd939*WProd967 - 
     -            WProd968 + WProd926*WProd971 - 
     -            WProd931*WProd939*WProd972 + WProd972*WProd973 + 
     -            WProd972*WProd974 + WProd939*WProd979 - 
     -            WProd980 - WProd981 - WProd938*WProd982 - 
     -            WProd971*WProd982 + WProd983 - 
     -            WProd972*WProd986 - WProd972*WProd987))/WProd947+
     -        5.333333333333333D0*
     -       (Cpl2*GS*Prp134*
     -          (WProd922*WProd925 - WProd929*WProd930 - 
     -            WProd922*WProd930*WProd931 + WProd930*WProd932 - 
     -            WProd930*WProd933 - WProd929*WProd934 - 
     -            WProd922*WProd931*WProd934 + WProd932*WProd934 - 
     -            WProd933*WProd934 + WProd922*WProd935 - 
     -            WProd926*WProd938 + WProd936*WProd938 + 
     -            WProd939*WProd940 - WProd941 - WProd946 - 
     -            WProd940*WProd948 + WProd940*WProd949 + 
     -            WProd940*WProd955 + WProd940*WProd956 + 
     -            WProd940*WProd957 + WProd940*WProd958 - 
     -            WProd940*WProd959 - WProd940*WProd960 - 
     -            WProd940*WProd961 + WProd922*WProd963 - 
     -            WProd964 + WProd922*WProd966 - 
     -            WProd939*WProd967 + WProd948*WProd967 - 
     -            WProd949*WProd967 - WProd955*WProd967 - 
     -            WProd956*WProd967 - WProd957*WProd967 - 
     -            WProd958*WProd967 + WProd959*WProd967 + 
     -            WProd960*WProd967 + WProd961*WProd967 + 
     -            WProd968 - WProd926*WProd971 + 
     -            WProd936*WProd971 + WProd931*WProd939*WProd972 - 
     -            WProd931*WProd948*WProd972 + 
     -            WProd931*WProd949*WProd972 + 
     -            WProd931*WProd955*WProd972 + 
     -            WProd931*WProd956*WProd972 + 
     -            WProd931*WProd957*WProd972 + 
     -            WProd931*WProd958*WProd972 - 
     -            WProd931*WProd959*WProd972 - 
     -            WProd931*WProd960*WProd972 - 
     -            WProd931*WProd961*WProd972 - WProd972*WProd974 - 
     -            WProd939*WProd979 + WProd948*WProd979 - 
     -            WProd949*WProd979 - WProd955*WProd979 - 
     -            WProd956*WProd979 - WProd957*WProd979 - 
     -            WProd958*WProd979 + WProd959*WProd979 + 
     -            WProd960*WProd979 + WProd961*WProd979 + 
     -            WProd980 + WProd981 + WProd938*WProd982 + 
     -            WProd971*WProd982 - WProd983 + 
     -            WProd972*WProd985 + WProd972*WProd987))/WProd947)
        amp(11,4,6) = Prp130*tmp1130 + 
     -   Prp126*Prp127*Prp128*
     -    ((8*Cpl2*GS*(WProd919 + WProd920))/WProd921 + 
     -      (8*Cpl2*GS*(WProd922*WProd923 + WProd922*WProd924 + 
     -           WProd922*WProd925 - WProd922*WProd926 - 
     -           WProd927 - WProd928 - WProd929*WProd930 - 
     -           WProd922*WProd930*WProd931 + WProd930*WProd932 - 
     -           WProd930*WProd933 - WProd929*WProd934 - 
     -           WProd922*WProd931*WProd934 + WProd932*WProd934 - 
     -           WProd933*WProd934 + WProd922*WProd935 + 
     -           WProd922*WProd936 + WProd936*WProd937 + 
     -           WProd936*WProd938 + WProd939*WProd940 - 
     -           WProd941 - WProd922*WProd942 + WProd943 - 
     -           WProd944 + WProd945 - WProd946))/WProd947) + 
     -   Prp127*Prp128*Prp129*
     -    ((8*Cpl2*GS*(WProd919 + WProd920))/WProd921 + 
     -      (8*Cpl2*GS*(-(WProd929*WProd930) + WProd930*WProd932 - 
     -           WProd929*WProd934 + WProd932*WProd934 - 
     -           WProd923*WProd937 - WProd924*WProd937 - 
     -           WProd925*WProd937 + WProd926*WProd937 + 
     -           WProd930*WProd931*WProd937 + 
     -           WProd931*WProd934*WProd937 - WProd935*WProd937 + 
     -           2*WProd939*WProd940 - WProd941 + 
     -           WProd937*WProd942 - WProd944 - WProd946 - 
     -           2*WProd940*WProd948 + WProd940*WProd949 - 
     -           WProd930*WProd951 - WProd934*WProd951 - 
     -           WProd932*WProd952 + WProd951*WProd952 - 
     -           WProd930*WProd953 - WProd934*WProd953 + 
     -           WProd952*WProd953 - WProd932*WProd954 + 
     -           WProd951*WProd954 + WProd953*WProd954 + 
     -           2*WProd940*WProd955 + 2*WProd940*WProd956 + 
     -           2*WProd940*WProd957 + 2*WProd940*WProd958 - 
     -           2*WProd940*WProd959 - WProd940*WProd960 - 
     -           2*WProd940*WProd961 + WProd962))/WProd947)

#ifdef DEBUG
	print *, 'amp(11,4,6) =', amp(11,4,6)
#endif

        tmp1131 = (4*Cpl2*GS*Prp127*Prp129*
     -      (-2*WProd1003 - 2*WProd1006 - 2*WProd1008 + 
     -        2*WProd1011 - 2*WProd1013 - 2*WProd1009*WProd932 - 
     -        2*WProd1010*WProd932 - 2*WProd1000*WProd937 - 
     -        2*WProd1012*WProd937 + 4*WProd1002*WProd939 - 
     -        4*WProd1002*WProd948 + 3*WProd1002*WProd949 + 
     -        2*WProd1009*WProd951 + 2*WProd1010*WProd951 + 
     -        2*WProd1009*WProd953 + 2*WProd1010*WProd953 + 
     -        4*WProd1002*WProd955 + 4*WProd1002*WProd956 + 
     -        4*WProd1002*WProd957 + 4*WProd1002*WProd958 - 
     -        4*WProd1002*WProd959 - 2*WProd1002*WProd960 - 
     -        2*WProd1002*WProd961 - 2*WProd937*WProd992 - 
     -        2*WProd937*WProd994 + 2*WProd937*WProd995 - 
     -        2*WProd929*WProd998 + 2*WProd932*WProd998 + 
     -        2*WProd931*WProd937*WProd998 - 2*WProd951*WProd998 - 
     -        2*WProd953*WProd998 - 2*WProd929*WProd999 + 
     -        2*WProd932*WProd999 + 2*WProd931*WProd937*WProd999 - 
     -        2*WProd951*WProd999 - 2*WProd953*WProd999))/WProd991+
     -     5.333333333333333D0*
     -    (Cpl2*GS*Prp129*Prp131*
     -       (-WProd1003 - WProd1008 + WProd1011 - WProd1013 + 
     -         WProd1016 + WProd1017 - WProd1009*WProd932 - 
     -         WProd1010*WProd932 - WProd1000*WProd937 - 
     -         WProd1012*WProd937 - WProd1014*WProd937 + 
     -         2*WProd1002*WProd939 - 2*WProd1015*WProd939 - 
     -         2*WProd1002*WProd948 + 2*WProd1015*WProd948 + 
     -         WProd1002*WProd949 - WProd1015*WProd949 + 
     -         WProd1009*WProd953 + WProd1010*WProd953 + 
     -         2*WProd1002*WProd955 - 2*WProd1015*WProd955 + 
     -         2*WProd1002*WProd956 - 2*WProd1015*WProd956 + 
     -         2*WProd1002*WProd957 - 2*WProd1015*WProd957 + 
     -         2*WProd1002*WProd958 - 2*WProd1015*WProd958 - 
     -         2*WProd1002*WProd959 + 2*WProd1015*WProd959 - 
     -         WProd1002*WProd960 + WProd1015*WProd960 - 
     -         WProd1002*WProd961 + WProd1015*WProd961 + 
     -         WProd1009*WProd965 + WProd1010*WProd965 - 
     -         WProd937*WProd994 + WProd937*WProd995 - 
     -         WProd929*WProd998 + WProd932*WProd998 + 
     -         WProd931*WProd937*WProd998 - WProd953*WProd998 - 
     -         WProd965*WProd998 - WProd929*WProd999 + 
     -         WProd932*WProd999 + WProd931*WProd937*WProd999 - 
     -         WProd953*WProd999 - WProd965*WProd999))/WProd991 + 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp131*Prp134*
     -       (-WProd1003 - WProd1008 - WProd1013 + WProd1016 + 
     -         WProd1020 + WProd1021 - WProd1023 - 
     -         WProd1000*WProd937 - WProd1012*WProd937 - 
     -         WProd1014*WProd937 + WProd1000*WProd938 + 
     -         WProd1001*WProd938 + WProd1012*WProd938 + 
     -         WProd1014*WProd938 + WProd1022*WProd938 - 
     -         WProd1009*WProd931*WProd938 - 
     -         WProd1010*WProd931*WProd938 + 
     -         2*WProd1002*WProd939 - 2*WProd1015*WProd939 - 
     -         2*WProd1019*WProd939 - 2*WProd1002*WProd948 + 
     -         2*WProd1015*WProd948 + 2*WProd1019*WProd948 + 
     -         2*WProd1002*WProd949 - 2*WProd1015*WProd949 - 
     -         2*WProd1019*WProd949 + WProd1009*WProd953 + 
     -         WProd1010*WProd953 + 2*WProd1002*WProd955 - 
     -         2*WProd1015*WProd955 - 2*WProd1019*WProd955 + 
     -         3*WProd1002*WProd956 - 3*WProd1015*WProd956 - 
     -         3*WProd1019*WProd956 + 2*WProd1002*WProd957 - 
     -         2*WProd1015*WProd957 - 2*WProd1019*WProd957 + 
     -         3*WProd1002*WProd958 - 3*WProd1015*WProd958 - 
     -         3*WProd1019*WProd958 - 3*WProd1002*WProd959 + 
     -         3*WProd1015*WProd959 + 3*WProd1019*WProd959 - 
     -         2*WProd1002*WProd960 + 2*WProd1015*WProd960 + 
     -         2*WProd1019*WProd960 - 2*WProd1002*WProd961 + 
     -         2*WProd1015*WProd961 + 2*WProd1019*WProd961 + 
     -         WProd1009*WProd965 + WProd1010*WProd965 + 
     -         WProd1000*WProd971 + WProd1001*WProd971 + 
     -         WProd1012*WProd971 + WProd1014*WProd971 + 
     -         WProd1022*WProd971 - WProd1009*WProd931*WProd971 - 
     -         WProd1010*WProd931*WProd971 + WProd1009*WProd977 + 
     -         WProd1010*WProd977 + WProd1009*WProd978 + 
     -         WProd1010*WProd978 - WProd937*WProd994 + 
     -         WProd938*WProd994 + WProd971*WProd994 - 
     -         WProd938*WProd995 - WProd971*WProd995 - 
     -         WProd929*WProd998 + WProd932*WProd998 - 
     -         WProd933*WProd998 + WProd931*WProd937*WProd998 - 
     -         WProd953*WProd998 - WProd965*WProd998 - 
     -         WProd977*WProd998 - WProd978*WProd998 - 
     -         WProd929*WProd999 + WProd932*WProd999 - 
     -         WProd933*WProd999 + WProd931*WProd937*WProd999 - 
     -         WProd953*WProd999 - WProd965*WProd999 - 
     -         WProd977*WProd999 - WProd978*WProd999))/WProd991
        tmp1131 = tmp1131 + 
     -   Prp126*((4*Cpl2*GS*Prp127*
     -         (-2*WProd1003 - 2*WProd1006 - 2*WProd1008 + 
     -           2*WProd1011 - 2*WProd1013 + 
     -           2*WProd1000*WProd922 + 2*WProd1012*WProd922 + 
     -           2*WProd1001*WProd938 + 2*WProd1002*WProd939 - 
     -           WProd1002*WProd948 + WProd1002*WProd949 + 
     -           WProd1002*WProd955 + WProd1002*WProd956 + 
     -           WProd1002*WProd957 + WProd1002*WProd958 - 
     -           WProd1002*WProd959 - WProd1002*WProd960 - 
     -           WProd1002*WProd961 + 2*WProd1001*WProd971 + 
     -           2*WProd1018*WProd973 + 2*WProd1018*WProd985 - 
     -           2*WProd1018*WProd986 + 2*WProd922*WProd992 + 
     -           2*WProd922*WProd994 - 2*WProd922*WProd995 - 
     -           2*WProd929*WProd998 - 
     -           2*WProd922*WProd931*WProd998 + 
     -           2*WProd932*WProd998 - 2*WProd933*WProd998 - 
     -           2*WProd929*WProd999 - 
     -           2*WProd922*WProd931*WProd999 + 
     -           2*WProd932*WProd999 - 2*WProd933*WProd999))/
     -       WProd991 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (-WProd1006 + WProd1011 - WProd1016 - WProd1020 - 
     -            WProd1021 + WProd1023 - WProd1014*WProd922 - 
     -            WProd1022*WProd938 + WProd1015*WProd939 + 
     -            WProd1019*WProd939 - 
     -            WProd1018*WProd931*WProd939 - 
     -            WProd1022*WProd971 + WProd1018*WProd973 + 
     -            WProd1018*WProd974 - WProd1018*WProd986 - 
     -            WProd1018*WProd987 + WProd922*WProd992 - 
     -            WProd922*WProd995 + WProd938*WProd995 + 
     -            WProd971*WProd995))/WProd991 + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp134*
     -          (-WProd1003 - WProd1008 - WProd1013 + WProd1016 + 
     -            WProd1020 + WProd1021 - WProd1023 + 
     -            WProd1000*WProd922 + WProd1012*WProd922 + 
     -            WProd1014*WProd922 + WProd1001*WProd938 + 
     -            WProd1022*WProd938 + WProd1002*WProd939 - 
     -            WProd1015*WProd939 - WProd1019*WProd939 + 
     -            WProd1018*WProd931*WProd939 - 
     -            WProd1002*WProd948 + WProd1015*WProd948 + 
     -            WProd1019*WProd948 - 
     -            WProd1018*WProd931*WProd948 + 
     -            WProd1002*WProd949 - WProd1015*WProd949 - 
     -            WProd1019*WProd949 + 
     -            WProd1018*WProd931*WProd949 + 
     -            WProd1002*WProd955 - WProd1015*WProd955 - 
     -            WProd1019*WProd955 + 
     -            WProd1018*WProd931*WProd955 + 
     -            WProd1002*WProd956 - WProd1015*WProd956 - 
     -            WProd1019*WProd956 + 
     -            WProd1018*WProd931*WProd956 + 
     -            WProd1002*WProd957 - WProd1015*WProd957 - 
     -            WProd1019*WProd957 + 
     -            WProd1018*WProd931*WProd957 + 
     -            WProd1002*WProd958 - WProd1015*WProd958 - 
     -            WProd1019*WProd958 + 
     -            WProd1018*WProd931*WProd958 - 
     -            WProd1002*WProd959 + WProd1015*WProd959 + 
     -            WProd1019*WProd959 - 
     -            WProd1018*WProd931*WProd959 - 
     -            WProd1002*WProd960 + WProd1015*WProd960 + 
     -            WProd1019*WProd960 - 
     -            WProd1018*WProd931*WProd960 - 
     -            WProd1002*WProd961 + WProd1015*WProd961 + 
     -            WProd1019*WProd961 - 
     -            WProd1018*WProd931*WProd961 + 
     -            WProd1001*WProd971 + WProd1022*WProd971 - 
     -            WProd1018*WProd974 + WProd1018*WProd985 + 
     -            WProd1018*WProd987 + WProd922*WProd994 - 
     -            WProd938*WProd995 - WProd971*WProd995 - 
     -            WProd929*WProd998 - WProd922*WProd931*WProd998 + 
     -            WProd932*WProd998 - WProd933*WProd998 - 
     -            WProd929*WProd999 - WProd922*WProd931*WProd999 + 
     -            WProd932*WProd999 - WProd933*WProd999))/WProd991)
        tmp1131 = tmp1131 + 
     -   Prp132*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (-WProd1006 - WProd1016 - WProd1017 - 
     -            WProd1014*WProd922 + WProd1014*WProd938 + 
     -            3*WProd1015*WProd939 - 
     -            3*WProd1018*WProd931*WProd939 - 
     -            3*WProd1015*WProd948 + 
     -            3*WProd1018*WProd931*WProd948 + 
     -            3*WProd1015*WProd949 - 
     -            3*WProd1018*WProd931*WProd949 + 
     -            2*WProd1015*WProd955 - 
     -            2*WProd1018*WProd931*WProd955 + 
     -            3*WProd1015*WProd956 - 
     -            3*WProd1018*WProd931*WProd956 + 
     -            2*WProd1015*WProd957 - 
     -            2*WProd1018*WProd931*WProd957 + 
     -            3*WProd1015*WProd958 - 
     -            3*WProd1018*WProd931*WProd958 - 
     -            2*WProd1015*WProd959 + 
     -            2*WProd1018*WProd931*WProd959 - 
     -            2*WProd1015*WProd960 + 
     -            2*WProd1018*WProd931*WProd960 - 
     -            2*WProd1015*WProd961 + 
     -            2*WProd1018*WProd931*WProd961 + 
     -            WProd1014*WProd971 + WProd1018*WProd973 + 
     -            WProd1018*WProd974 + WProd1018*WProd975 + 
     -            WProd922*WProd992 - WProd938*WProd992 - 
     -            WProd971*WProd992))/WProd991 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp129*
     -          (-WProd1006 - WProd1016 - WProd1017 + 
     -            WProd1014*WProd937 + 2*WProd1015*WProd939 - 
     -            2*WProd1015*WProd948 + 2*WProd1015*WProd949 + 
     -            WProd1009*WProd951 + WProd1010*WProd951 + 
     -            2*WProd1015*WProd955 + 2*WProd1015*WProd956 + 
     -            2*WProd1015*WProd957 + 2*WProd1015*WProd958 - 
     -            2*WProd1015*WProd959 - WProd1015*WProd960 - 
     -            WProd1015*WProd961 - WProd1009*WProd965 - 
     -            WProd1010*WProd965 - WProd937*WProd992 - 
     -            WProd951*WProd998 + WProd965*WProd998 - 
     -            WProd951*WProd999 + WProd965*WProd999))/WProd991)
        amp(11,4,4) = Prp130*tmp1131 + 
     -   Prp126*Prp127*Prp128*
     -    ((8*Cpl2*GS*(WProd989 + WProd990))/WProd988 + 
     -      (8*Cpl2*GS*(-WProd1003 + WProd1005 - WProd1006 + 
     -           WProd1007 - WProd1008 + WProd1000*WProd922 + 
     -           WProd1001*WProd922 - WProd1004*WProd922 + 
     -           WProd1001*WProd937 + WProd1001*WProd938 + 
     -           WProd1002*WProd939 + WProd922*WProd992 + 
     -           WProd922*WProd993 + WProd922*WProd994 - 
     -           WProd922*WProd995 - WProd996 - WProd997 - 
     -           WProd929*WProd998 - WProd922*WProd931*WProd998 + 
     -           WProd932*WProd998 - WProd933*WProd998 - 
     -           WProd929*WProd999 - WProd922*WProd931*WProd999 + 
     -           WProd932*WProd999 - WProd933*WProd999))/WProd991)+
     -     Prp127*Prp128*Prp129*
     -    ((8*Cpl2*GS*(WProd989 + WProd990))/WProd988 + 
     -      (8*Cpl2*GS*(-WProd1003 - WProd1006 - WProd1008 + 
     -           WProd1011 - WProd1009*WProd932 - 
     -           WProd1010*WProd932 - WProd1000*WProd937 + 
     -           WProd1004*WProd937 + 2*WProd1002*WProd939 - 
     -           2*WProd1002*WProd948 + WProd1002*WProd949 + 
     -           WProd1009*WProd951 + WProd1010*WProd951 + 
     -           WProd1009*WProd953 + WProd1010*WProd953 + 
     -           2*WProd1002*WProd955 + 2*WProd1002*WProd956 + 
     -           2*WProd1002*WProd957 + 2*WProd1002*WProd958 - 
     -           2*WProd1002*WProd959 - WProd1002*WProd960 - 
     -           2*WProd1002*WProd961 - WProd937*WProd992 - 
     -           WProd937*WProd993 - WProd937*WProd994 + 
     -           WProd937*WProd995 - WProd929*WProd998 + 
     -           WProd932*WProd998 + WProd931*WProd937*WProd998 - 
     -           WProd951*WProd998 - WProd953*WProd998 - 
     -           WProd929*WProd999 + WProd932*WProd999 + 
     -           WProd931*WProd937*WProd999 - WProd951*WProd999 - 
     -           WProd953*WProd999))/WProd991)

#ifdef DEBUG
	print *, 'amp(11,4,4) =', amp(11,4,4)
#endif

        tmp1132 = (4*Cpl2*GS*Prp127*Prp129*
     -      (-2*WProd1032 - 2*WProd1034 - 2*WProd1036 + 
     -        2*WProd1037 - 2*WProd1042 - 2*WProd1029*WProd923 - 
     -        2*WProd1029*WProd925 + 2*WProd1029*WProd926 + 
     -        2*WProd1027*WProd930 - 2*WProd1040*WProd930 - 
     -        2*WProd1041*WProd930 + 
     -        2*WProd1029*WProd930*WProd931 + 
     -        2*WProd1027*WProd934 - 2*WProd1040*WProd934 - 
     -        2*WProd1041*WProd934 + 
     -        2*WProd1029*WProd931*WProd934 - 
     -        2*WProd1029*WProd935 + 4*WProd1031*WProd939 - 
     -        4*WProd1031*WProd948 + 3*WProd1031*WProd949 - 
     -        2*WProd1027*WProd952 + 2*WProd1040*WProd952 + 
     -        2*WProd1041*WProd952 - 2*WProd1027*WProd954 + 
     -        2*WProd1040*WProd954 + 2*WProd1041*WProd954 + 
     -        4*WProd1031*WProd955 + 4*WProd1031*WProd956 + 
     -        4*WProd1031*WProd957 + 4*WProd1031*WProd958 - 
     -        4*WProd1031*WProd959 - 2*WProd1031*WProd960 - 
     -        2*WProd1031*WProd961 - 2*WProd1029*WProd963))/
     -    WProd1039 + 5.333333333333333D0*
     -    (Cpl2*GS*Prp129*Prp131*
     -       (-WProd1032 - WProd1036 + WProd1037 - WProd1042 + 
     -         WProd1045 + WProd1046 - WProd1029*WProd925 + 
     -         WProd1029*WProd926 + WProd1027*WProd930 - 
     -         WProd1041*WProd930 - WProd1043*WProd930 + 
     -         WProd1029*WProd930*WProd931 + WProd1027*WProd934 - 
     -         WProd1041*WProd934 - WProd1043*WProd934 + 
     -         WProd1029*WProd931*WProd934 - WProd1029*WProd935 + 
     -         2*WProd1031*WProd939 - 2*WProd1044*WProd939 - 
     -         2*WProd1031*WProd948 + 2*WProd1044*WProd948 + 
     -         WProd1031*WProd949 - WProd1044*WProd949 - 
     -         WProd1027*WProd952 + WProd1041*WProd952 + 
     -         WProd1043*WProd952 - WProd1027*WProd954 + 
     -         WProd1041*WProd954 + WProd1043*WProd954 + 
     -         2*WProd1031*WProd955 - 2*WProd1044*WProd955 + 
     -         2*WProd1031*WProd956 - 2*WProd1044*WProd956 + 
     -         2*WProd1031*WProd957 - 2*WProd1044*WProd957 + 
     -         2*WProd1031*WProd958 - 2*WProd1044*WProd958 - 
     -         2*WProd1031*WProd959 + 2*WProd1044*WProd959 - 
     -         WProd1031*WProd960 + WProd1044*WProd960 - 
     -         WProd1031*WProd961 + WProd1044*WProd961 - 
     -         WProd1029*WProd963 - WProd1029*WProd966))/WProd1039+
     -     Prp132*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp129*
     -          (-WProd1034 - WProd1045 - WProd1046 - 
     -            WProd1029*WProd923 - WProd1040*WProd930 + 
     -            WProd1043*WProd930 - WProd1040*WProd934 + 
     -            WProd1043*WProd934 + 2*WProd1044*WProd939 - 
     -            2*WProd1044*WProd948 + 2*WProd1044*WProd949 + 
     -            WProd1040*WProd952 - WProd1043*WProd952 + 
     -            WProd1040*WProd954 - WProd1043*WProd954 + 
     -            2*WProd1044*WProd955 + 2*WProd1044*WProd956 + 
     -            2*WProd1044*WProd957 + 2*WProd1044*WProd958 - 
     -            2*WProd1044*WProd959 - WProd1044*WProd960 - 
     -            WProd1044*WProd961 + WProd1029*WProd966))/
     -        WProd1039 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (-WProd1034 - WProd1045 - WProd1046 + 
     -            WProd1024*WProd923 - WProd1030*WProd923 - 
     -            WProd1047*WProd923 + 3*WProd1044*WProd939 - 
     -            3*WProd1048*WProd931*WProd939 - 
     -            3*WProd1044*WProd948 + 
     -            3*WProd1048*WProd931*WProd948 + 
     -            3*WProd1044*WProd949 - 
     -            3*WProd1048*WProd931*WProd949 + 
     -            2*WProd1044*WProd955 - 
     -            2*WProd1048*WProd931*WProd955 + 
     -            3*WProd1044*WProd956 - 
     -            3*WProd1048*WProd931*WProd956 + 
     -            2*WProd1044*WProd957 - 
     -            2*WProd1048*WProd931*WProd957 + 
     -            3*WProd1044*WProd958 - 
     -            3*WProd1048*WProd931*WProd958 - 
     -            2*WProd1044*WProd959 + 
     -            2*WProd1048*WProd931*WProd959 - 
     -            2*WProd1044*WProd960 + 
     -            2*WProd1048*WProd931*WProd960 - 
     -            2*WProd1044*WProd961 + 
     -            2*WProd1048*WProd931*WProd961 - 
     -            WProd1024*WProd966 + WProd1030*WProd966 + 
     -            WProd1047*WProd966 + WProd1048*WProd973 + 
     -            WProd1048*WProd974 + WProd1048*WProd975))/
     -        WProd1039)
        tmp1132 = tmp1132 + 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp131*Prp134*
     -       (-WProd1032 - WProd1036 - WProd1042 + WProd1045 + 
     -         WProd1052 + WProd1053 - WProd1054 - 
     -         WProd1029*WProd925 + WProd1030*WProd925 + 
     -         WProd1047*WProd925 - WProd1030*WProd926 - 
     -         WProd1047*WProd926 + WProd1027*WProd930 - 
     -         WProd1028*WProd930 - WProd1041*WProd930 - 
     -         WProd1043*WProd930 - WProd1049*WProd930 - 
     -         WProd1050*WProd930 + WProd1029*WProd930*WProd931 + 
     -         WProd1027*WProd934 - WProd1028*WProd934 - 
     -         WProd1041*WProd934 - WProd1043*WProd934 - 
     -         WProd1049*WProd934 - WProd1050*WProd934 + 
     -         WProd1029*WProd931*WProd934 - WProd1029*WProd935 + 
     -         WProd1030*WProd935 + WProd1047*WProd935 + 
     -         WProd1030*WProd936 + WProd1047*WProd936 + 
     -         2*WProd1031*WProd939 - 2*WProd1044*WProd939 - 
     -         2*WProd1051*WProd939 - 2*WProd1031*WProd948 + 
     -         2*WProd1044*WProd948 + 2*WProd1051*WProd948 + 
     -         2*WProd1031*WProd949 - 2*WProd1044*WProd949 - 
     -         2*WProd1051*WProd949 + WProd1041*WProd952 + 
     -         WProd1043*WProd952 + WProd1049*WProd952 + 
     -         WProd1050*WProd952 - WProd1030*WProd931*WProd952 - 
     -         WProd1047*WProd931*WProd952 + WProd1041*WProd954 + 
     -         WProd1043*WProd954 + WProd1049*WProd954 + 
     -         WProd1050*WProd954 - WProd1030*WProd931*WProd954 - 
     -         WProd1047*WProd931*WProd954 + 
     -         2*WProd1031*WProd955 - 2*WProd1044*WProd955 - 
     -         2*WProd1051*WProd955 + 3*WProd1031*WProd956 - 
     -         3*WProd1044*WProd956 - 3*WProd1051*WProd956 + 
     -         2*WProd1031*WProd957 - 2*WProd1044*WProd957 - 
     -         2*WProd1051*WProd957 + 3*WProd1031*WProd958 - 
     -         3*WProd1044*WProd958 - 3*WProd1051*WProd958 - 
     -         3*WProd1031*WProd959 + 3*WProd1044*WProd959 + 
     -         3*WProd1051*WProd959 - 2*WProd1031*WProd960 + 
     -         2*WProd1044*WProd960 + 2*WProd1051*WProd960 - 
     -         2*WProd1031*WProd961 + 2*WProd1044*WProd961 + 
     -         2*WProd1051*WProd961 - WProd1029*WProd963 + 
     -         WProd1030*WProd963 + WProd1047*WProd963 - 
     -         WProd1029*WProd966 + WProd1030*WProd966 + 
     -         WProd1047*WProd966 + WProd1030*WProd982 + 
     -         WProd1047*WProd982))/WProd1039
        tmp1132 = tmp1132 + 
     -   Prp126*((4*Cpl2*GS*Prp127*
     -         (-2*WProd1032 - 2*WProd1034 - 2*WProd1036 + 
     -           2*WProd1037 - 2*WProd1042 + 
     -           2*WProd1024*WProd923 + 2*WProd1024*WProd925 - 
     -           2*WProd1024*WProd926 + 2*WProd1027*WProd930 - 
     -           2*WProd1028*WProd930 - 
     -           2*WProd1024*WProd930*WProd931 + 
     -           2*WProd1027*WProd934 - 2*WProd1028*WProd934 - 
     -           2*WProd1024*WProd931*WProd934 + 
     -           2*WProd1024*WProd935 + 2*WProd1030*WProd936 + 
     -           2*WProd1047*WProd936 + 2*WProd1031*WProd939 - 
     -           WProd1031*WProd948 + WProd1031*WProd949 + 
     -           WProd1031*WProd955 + WProd1031*WProd956 + 
     -           WProd1031*WProd957 + WProd1031*WProd958 - 
     -           WProd1031*WProd959 - WProd1031*WProd960 - 
     -           WProd1031*WProd961 + 2*WProd1024*WProd963 + 
     -           2*WProd1048*WProd973 + 2*WProd1048*WProd985 - 
     -           2*WProd1048*WProd986))/WProd1039 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (-WProd1034 + WProd1037 - WProd1045 - WProd1052 - 
     -            WProd1053 + WProd1054 + WProd1024*WProd923 - 
     -            WProd1024*WProd926 + WProd1030*WProd926 + 
     -            WProd1047*WProd926 + WProd1044*WProd939 + 
     -            WProd1051*WProd939 - 
     -            WProd1048*WProd931*WProd939 - 
     -            WProd1024*WProd966 + WProd1048*WProd973 + 
     -            WProd1048*WProd974 - WProd1030*WProd982 - 
     -            WProd1047*WProd982 - WProd1048*WProd986 - 
     -            WProd1048*WProd987))/WProd1039 + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp134*
     -          (-WProd1032 - WProd1036 - WProd1042 + WProd1045 + 
     -            WProd1052 + WProd1053 - WProd1054 + 
     -            WProd1024*WProd925 - WProd1030*WProd926 - 
     -            WProd1047*WProd926 + WProd1027*WProd930 - 
     -            WProd1028*WProd930 - 
     -            WProd1024*WProd930*WProd931 + 
     -            WProd1027*WProd934 - WProd1028*WProd934 - 
     -            WProd1024*WProd931*WProd934 + 
     -            WProd1024*WProd935 + WProd1030*WProd936 + 
     -            WProd1047*WProd936 + WProd1031*WProd939 - 
     -            WProd1044*WProd939 - WProd1051*WProd939 + 
     -            WProd1048*WProd931*WProd939 - 
     -            WProd1031*WProd948 + WProd1044*WProd948 + 
     -            WProd1051*WProd948 - 
     -            WProd1048*WProd931*WProd948 + 
     -            WProd1031*WProd949 - WProd1044*WProd949 - 
     -            WProd1051*WProd949 + 
     -            WProd1048*WProd931*WProd949 + 
     -            WProd1031*WProd955 - WProd1044*WProd955 - 
     -            WProd1051*WProd955 + 
     -            WProd1048*WProd931*WProd955 + 
     -            WProd1031*WProd956 - WProd1044*WProd956 - 
     -            WProd1051*WProd956 + 
     -            WProd1048*WProd931*WProd956 + 
     -            WProd1031*WProd957 - WProd1044*WProd957 - 
     -            WProd1051*WProd957 + 
     -            WProd1048*WProd931*WProd957 + 
     -            WProd1031*WProd958 - WProd1044*WProd958 - 
     -            WProd1051*WProd958 + 
     -            WProd1048*WProd931*WProd958 - 
     -            WProd1031*WProd959 + WProd1044*WProd959 + 
     -            WProd1051*WProd959 - 
     -            WProd1048*WProd931*WProd959 - 
     -            WProd1031*WProd960 + WProd1044*WProd960 + 
     -            WProd1051*WProd960 - 
     -            WProd1048*WProd931*WProd960 - 
     -            WProd1031*WProd961 + WProd1044*WProd961 + 
     -            WProd1051*WProd961 - 
     -            WProd1048*WProd931*WProd961 + 
     -            WProd1024*WProd963 + WProd1024*WProd966 - 
     -            WProd1048*WProd974 + WProd1030*WProd982 + 
     -            WProd1047*WProd982 + WProd1048*WProd985 + 
     -            WProd1048*WProd987))/WProd1039)
        amp(11,4,2) = Prp130*tmp1132 + 
     -   (8*Cpl2*GS*Prp126*Prp127*Prp128*
     -      (-WProd1025 - WProd1026 - WProd1032 + WProd1033 - 
     -        WProd1034 + WProd1035 - WProd1036 + WProd1037 - 
     -        WProd1038 + WProd1024*WProd923 + 
     -        WProd1024*WProd924 + WProd1024*WProd925 - 
     -        WProd1024*WProd926 + WProd1027*WProd930 - 
     -        WProd1028*WProd930 - WProd1024*WProd930*WProd931 + 
     -        WProd1027*WProd934 - WProd1028*WProd934 - 
     -        WProd1024*WProd931*WProd934 + WProd1024*WProd935 + 
     -        WProd1024*WProd936 + WProd1029*WProd936 + 
     -        WProd1030*WProd936 + WProd1031*WProd939 - 
     -        WProd1024*WProd942))/WProd1039 + 
     -   (8*Cpl2*GS*Prp127*Prp128*Prp129*
     -      (-WProd1032 - WProd1034 - WProd1036 + WProd1037 - 
     -        WProd1029*WProd923 - WProd1029*WProd924 - 
     -        WProd1029*WProd925 + WProd1029*WProd926 + 
     -        WProd1027*WProd930 - WProd1040*WProd930 - 
     -        WProd1041*WProd930 + WProd1029*WProd930*WProd931 + 
     -        WProd1027*WProd934 - WProd1040*WProd934 - 
     -        WProd1041*WProd934 + WProd1029*WProd931*WProd934 - 
     -        WProd1029*WProd935 + 2*WProd1031*WProd939 + 
     -        WProd1029*WProd942 - 2*WProd1031*WProd948 + 
     -        WProd1031*WProd949 - WProd1027*WProd952 + 
     -        WProd1040*WProd952 + WProd1041*WProd952 - 
     -        WProd1027*WProd954 + WProd1040*WProd954 + 
     -        WProd1041*WProd954 + 2*WProd1031*WProd955 + 
     -        2*WProd1031*WProd956 + 2*WProd1031*WProd957 + 
     -        2*WProd1031*WProd958 - 2*WProd1031*WProd959 - 
     -        WProd1031*WProd960 - 2*WProd1031*WProd961))/WProd1039

#ifdef DEBUG
	print *, 'amp(11,4,2) =', amp(11,4,2)
#endif

        tmp1133 = (4*Cpl2*GS*Prp127*Prp129*
     -      (-2*WProd1009*WProd1027 - 2*WProd1010*WProd1027 - 
     -        2*WProd1000*WProd1029 - 2*WProd1012*WProd1029 + 
     -        2*WProd1009*WProd1040 + 2*WProd1010*WProd1040 + 
     -        2*WProd1009*WProd1041 + 2*WProd1010*WProd1041 - 
     -        2*WProd1059 - 2*WProd1061 - 2*WProd1063 + 
     -        2*WProd1064 - 2*WProd1066 + 4*WProd1058*WProd939 - 
     -        4*WProd1058*WProd948 + 3*WProd1058*WProd949 + 
     -        4*WProd1058*WProd955 + 4*WProd1058*WProd956 + 
     -        4*WProd1058*WProd957 + 4*WProd1058*WProd958 - 
     -        4*WProd1058*WProd959 - 2*WProd1058*WProd960 - 
     -        2*WProd1058*WProd961 - 2*WProd1029*WProd992 - 
     -        2*WProd1029*WProd994 + 2*WProd1029*WProd995 + 
     -        2*WProd1027*WProd998 - 2*WProd1040*WProd998 - 
     -        2*WProd1041*WProd998 + 
     -        2*WProd1029*WProd931*WProd998 + 
     -        2*WProd1027*WProd999 - 2*WProd1040*WProd999 - 
     -        2*WProd1041*WProd999 + 2*WProd1029*WProd931*WProd999)
     -      )/WProd1055 + 5.333333333333333D0*
     -    (Cpl2*GS*Prp129*Prp131*
     -       (-(WProd1009*WProd1027) - WProd1010*WProd1027 - 
     -         WProd1000*WProd1029 - WProd1012*WProd1029 - 
     -         WProd1014*WProd1029 + WProd1009*WProd1041 + 
     -         WProd1010*WProd1041 + WProd1009*WProd1043 + 
     -         WProd1010*WProd1043 - WProd1059 - WProd1063 + 
     -         WProd1064 - WProd1066 + WProd1068 + WProd1069 + 
     -         2*WProd1058*WProd939 - 2*WProd1067*WProd939 - 
     -         2*WProd1058*WProd948 + 2*WProd1067*WProd948 + 
     -         WProd1058*WProd949 - WProd1067*WProd949 + 
     -         2*WProd1058*WProd955 - 2*WProd1067*WProd955 + 
     -         2*WProd1058*WProd956 - 2*WProd1067*WProd956 + 
     -         2*WProd1058*WProd957 - 2*WProd1067*WProd957 + 
     -         2*WProd1058*WProd958 - 2*WProd1067*WProd958 - 
     -         2*WProd1058*WProd959 + 2*WProd1067*WProd959 - 
     -         WProd1058*WProd960 + WProd1067*WProd960 - 
     -         WProd1058*WProd961 + WProd1067*WProd961 - 
     -         WProd1029*WProd994 + WProd1029*WProd995 + 
     -         WProd1027*WProd998 - WProd1041*WProd998 - 
     -         WProd1043*WProd998 + WProd1029*WProd931*WProd998 + 
     -         WProd1027*WProd999 - WProd1041*WProd999 - 
     -         WProd1043*WProd999 + WProd1029*WProd931*WProd999))/
     -     WProd1055 + 5.333333333333333D0*
     -    (Cpl2*GS*Prp131*Prp134*
     -       (-(WProd1000*WProd1029) - WProd1012*WProd1029 - 
     -         WProd1014*WProd1029 + WProd1000*WProd1030 + 
     -         WProd1001*WProd1030 + WProd1012*WProd1030 + 
     -         WProd1014*WProd1030 + WProd1022*WProd1030 + 
     -         WProd1009*WProd1041 + WProd1010*WProd1041 + 
     -         WProd1009*WProd1043 + WProd1010*WProd1043 + 
     -         WProd1000*WProd1047 + WProd1001*WProd1047 + 
     -         WProd1012*WProd1047 + WProd1014*WProd1047 + 
     -         WProd1022*WProd1047 + WProd1009*WProd1049 + 
     -         WProd1010*WProd1049 + WProd1009*WProd1050 + 
     -         WProd1010*WProd1050 - WProd1059 - WProd1063 - 
     -         WProd1066 + WProd1068 + WProd1072 + WProd1073 - 
     -         WProd1074 - WProd1009*WProd1030*WProd931 - 
     -         WProd1010*WProd1030*WProd931 - 
     -         WProd1009*WProd1047*WProd931 - 
     -         WProd1010*WProd1047*WProd931 + 
     -         2*WProd1058*WProd939 - 2*WProd1067*WProd939 - 
     -         2*WProd1071*WProd939 - 2*WProd1058*WProd948 + 
     -         2*WProd1067*WProd948 + 2*WProd1071*WProd948 + 
     -         2*WProd1058*WProd949 - 2*WProd1067*WProd949 - 
     -         2*WProd1071*WProd949 + 2*WProd1058*WProd955 - 
     -         2*WProd1067*WProd955 - 2*WProd1071*WProd955 + 
     -         3*WProd1058*WProd956 - 3*WProd1067*WProd956 - 
     -         3*WProd1071*WProd956 + 2*WProd1058*WProd957 - 
     -         2*WProd1067*WProd957 - 2*WProd1071*WProd957 + 
     -         3*WProd1058*WProd958 - 3*WProd1067*WProd958 - 
     -         3*WProd1071*WProd958 - 3*WProd1058*WProd959 + 
     -         3*WProd1067*WProd959 + 3*WProd1071*WProd959 - 
     -         2*WProd1058*WProd960 + 2*WProd1067*WProd960 + 
     -         2*WProd1071*WProd960 - 2*WProd1058*WProd961 + 
     -         2*WProd1067*WProd961 + 2*WProd1071*WProd961 - 
     -         WProd1029*WProd994 + WProd1030*WProd994 + 
     -         WProd1047*WProd994 - WProd1030*WProd995 - 
     -         WProd1047*WProd995 + WProd1027*WProd998 - 
     -         WProd1028*WProd998 - WProd1041*WProd998 - 
     -         WProd1043*WProd998 - WProd1049*WProd998 - 
     -         WProd1050*WProd998 + WProd1029*WProd931*WProd998 + 
     -         WProd1027*WProd999 - WProd1028*WProd999 - 
     -         WProd1041*WProd999 - WProd1043*WProd999 - 
     -         WProd1049*WProd999 - WProd1050*WProd999 + 
     -         WProd1029*WProd931*WProd999))/WProd1055
        tmp1133 = tmp1133 + 
     -   Prp132*(2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (-(WProd1014*WProd1024) + WProd1014*WProd1030 + 
     -            WProd1014*WProd1047 - WProd1061 - WProd1068 - 
     -            WProd1069 + 3*WProd1067*WProd939 - 
     -            3*WProd1070*WProd931*WProd939 - 
     -            3*WProd1067*WProd948 + 
     -            3*WProd1070*WProd931*WProd948 + 
     -            3*WProd1067*WProd949 - 
     -            3*WProd1070*WProd931*WProd949 + 
     -            2*WProd1067*WProd955 - 
     -            2*WProd1070*WProd931*WProd955 + 
     -            3*WProd1067*WProd956 - 
     -            3*WProd1070*WProd931*WProd956 + 
     -            2*WProd1067*WProd957 - 
     -            2*WProd1070*WProd931*WProd957 + 
     -            3*WProd1067*WProd958 - 
     -            3*WProd1070*WProd931*WProd958 - 
     -            2*WProd1067*WProd959 + 
     -            2*WProd1070*WProd931*WProd959 - 
     -            2*WProd1067*WProd960 + 
     -            2*WProd1070*WProd931*WProd960 - 
     -            2*WProd1067*WProd961 + 
     -            2*WProd1070*WProd931*WProd961 + 
     -            WProd1070*WProd973 + WProd1070*WProd974 + 
     -            WProd1070*WProd975 + WProd1024*WProd992 - 
     -            WProd1030*WProd992 - WProd1047*WProd992))/
     -        WProd1055 + 2.6666666666666665D0*
     -       (Cpl2*GS*Prp129*
     -          (WProd1014*WProd1029 + WProd1009*WProd1040 + 
     -            WProd1010*WProd1040 - WProd1009*WProd1043 - 
     -            WProd1010*WProd1043 - WProd1061 - WProd1068 - 
     -            WProd1069 + 2*WProd1067*WProd939 - 
     -            2*WProd1067*WProd948 + 2*WProd1067*WProd949 + 
     -            2*WProd1067*WProd955 + 2*WProd1067*WProd956 + 
     -            2*WProd1067*WProd957 + 2*WProd1067*WProd958 - 
     -            2*WProd1067*WProd959 - WProd1067*WProd960 - 
     -            WProd1067*WProd961 - WProd1029*WProd992 - 
     -            WProd1040*WProd998 + WProd1043*WProd998 - 
     -            WProd1040*WProd999 + WProd1043*WProd999))/
     -        WProd1055) + 
     -   Prp126*((4*Cpl2*GS*Prp127*
     -         (2*WProd1000*WProd1024 + 2*WProd1012*WProd1024 + 
     -           2*WProd1001*WProd1030 + 2*WProd1001*WProd1047 - 
     -           2*WProd1059 - 2*WProd1061 - 2*WProd1063 + 
     -           2*WProd1064 - 2*WProd1066 + 
     -           2*WProd1058*WProd939 - WProd1058*WProd948 + 
     -           WProd1058*WProd949 + WProd1058*WProd955 + 
     -           WProd1058*WProd956 + WProd1058*WProd957 + 
     -           WProd1058*WProd958 - WProd1058*WProd959 - 
     -           WProd1058*WProd960 - WProd1058*WProd961 + 
     -           2*WProd1070*WProd973 + 2*WProd1070*WProd985 - 
     -           2*WProd1070*WProd986 + 2*WProd1024*WProd992 + 
     -           2*WProd1024*WProd994 - 2*WProd1024*WProd995 + 
     -           2*WProd1027*WProd998 - 2*WProd1028*WProd998 - 
     -           2*WProd1024*WProd931*WProd998 + 
     -           2*WProd1027*WProd999 - 2*WProd1028*WProd999 - 
     -           2*WProd1024*WProd931*WProd999))/WProd1055 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*Prp133*
     -          (-(WProd1014*WProd1024) - WProd1022*WProd1030 - 
     -            WProd1022*WProd1047 - WProd1061 + WProd1064 - 
     -            WProd1068 - WProd1072 - WProd1073 + WProd1074 + 
     -            WProd1067*WProd939 + WProd1071*WProd939 - 
     -            WProd1070*WProd931*WProd939 + 
     -            WProd1070*WProd973 + WProd1070*WProd974 - 
     -            WProd1070*WProd986 - WProd1070*WProd987 + 
     -            WProd1024*WProd992 - WProd1024*WProd995 + 
     -            WProd1030*WProd995 + WProd1047*WProd995))/
     -        WProd1055 + 5.333333333333333D0*
     -       (Cpl2*GS*Prp134*
     -          (WProd1000*WProd1024 + WProd1012*WProd1024 + 
     -            WProd1014*WProd1024 + WProd1001*WProd1030 + 
     -            WProd1022*WProd1030 + WProd1001*WProd1047 + 
     -            WProd1022*WProd1047 - WProd1059 - WProd1063 - 
     -            WProd1066 + WProd1068 + WProd1072 + WProd1073 - 
     -            WProd1074 + WProd1058*WProd939 - 
     -            WProd1067*WProd939 - WProd1071*WProd939 + 
     -            WProd1070*WProd931*WProd939 - 
     -            WProd1058*WProd948 + WProd1067*WProd948 + 
     -            WProd1071*WProd948 - 
     -            WProd1070*WProd931*WProd948 + 
     -            WProd1058*WProd949 - WProd1067*WProd949 - 
     -            WProd1071*WProd949 + 
     -            WProd1070*WProd931*WProd949 + 
     -            WProd1058*WProd955 - WProd1067*WProd955 - 
     -            WProd1071*WProd955 + 
     -            WProd1070*WProd931*WProd955 + 
     -            WProd1058*WProd956 - WProd1067*WProd956 - 
     -            WProd1071*WProd956 + 
     -            WProd1070*WProd931*WProd956 + 
     -            WProd1058*WProd957 - WProd1067*WProd957 - 
     -            WProd1071*WProd957 + 
     -            WProd1070*WProd931*WProd957 + 
     -            WProd1058*WProd958 - WProd1067*WProd958 - 
     -            WProd1071*WProd958 + 
     -            WProd1070*WProd931*WProd958 - 
     -            WProd1058*WProd959 + WProd1067*WProd959 + 
     -            WProd1071*WProd959 - 
     -            WProd1070*WProd931*WProd959 - 
     -            WProd1058*WProd960 + WProd1067*WProd960 + 
     -            WProd1071*WProd960 - 
     -            WProd1070*WProd931*WProd960 - 
     -            WProd1058*WProd961 + WProd1067*WProd961 + 
     -            WProd1071*WProd961 - 
     -            WProd1070*WProd931*WProd961 - 
     -            WProd1070*WProd974 + WProd1070*WProd985 + 
     -            WProd1070*WProd987 + WProd1024*WProd994 - 
     -            WProd1030*WProd995 - WProd1047*WProd995 + 
     -            WProd1027*WProd998 - WProd1028*WProd998 - 
     -            WProd1024*WProd931*WProd998 + 
     -            WProd1027*WProd999 - WProd1028*WProd999 - 
     -            WProd1024*WProd931*WProd999))/WProd1055)
        amp(11,4,7) = Prp130*tmp1133 + 
     -   (8*Cpl2*GS*Prp126*Prp127*Prp128*
     -      (WProd1000*WProd1024 + WProd1001*WProd1024 - 
     -        WProd1004*WProd1024 + WProd1001*WProd1029 + 
     -        WProd1001*WProd1030 - WProd1056 - WProd1057 - 
     -        WProd1059 + WProd1060 - WProd1061 + WProd1062 - 
     -        WProd1063 + WProd1064 - WProd1065 + 
     -        WProd1058*WProd939 + WProd1024*WProd992 + 
     -        WProd1024*WProd993 + WProd1024*WProd994 - 
     -        WProd1024*WProd995 + WProd1027*WProd998 - 
     -        WProd1028*WProd998 - WProd1024*WProd931*WProd998 + 
     -        WProd1027*WProd999 - WProd1028*WProd999 - 
     -        WProd1024*WProd931*WProd999))/WProd1055 + 
     -   (8*Cpl2*GS*Prp127*Prp128*Prp129*
     -      (-(WProd1009*WProd1027) - WProd1010*WProd1027 - 
     -        WProd1000*WProd1029 + WProd1004*WProd1029 + 
     -        WProd1009*WProd1040 + WProd1010*WProd1040 + 
     -        WProd1009*WProd1041 + WProd1010*WProd1041 - 
     -        WProd1059 - WProd1061 - WProd1063 + WProd1064 + 
     -        2*WProd1058*WProd939 - 2*WProd1058*WProd948 + 
     -        WProd1058*WProd949 + 2*WProd1058*WProd955 + 
     -        2*WProd1058*WProd956 + 2*WProd1058*WProd957 + 
     -        2*WProd1058*WProd958 - 2*WProd1058*WProd959 - 
     -        WProd1058*WProd960 - 2*WProd1058*WProd961 - 
     -        WProd1029*WProd992 - WProd1029*WProd993 - 
     -        WProd1029*WProd994 + WProd1029*WProd995 + 
     -        WProd1027*WProd998 - WProd1040*WProd998 - 
     -        WProd1041*WProd998 + WProd1029*WProd931*WProd998 + 
     -        WProd1027*WProd999 - WProd1040*WProd999 - 
     -        WProd1041*WProd999 + WProd1029*WProd931*WProd999))/
     -    WProd1055

#ifdef DEBUG
	print *, 'amp(11,4,7) =', amp(11,4,7)
#endif

      endif

      end


