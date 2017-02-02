      subroutine abbrrealQED11(resnr)
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

      
	WProd988 = Prod(6,2,6,1)
	WProd1055 = Prod(5,2,5,1)*Prod(6,2,6,1)
	WProd948 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd955 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd956 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd931 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd957 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd958 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd959 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd939 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd1024 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd973 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1040 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd985 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd1041 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd949 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd1029 = Prod(5,2,2,1)*ProdC(5,1,2,1)
	WProd960 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd1030 = Prod(5,2,3,1)*ProdC(5,1,3,1)
        WProd974 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
	WProd961 = Prod(5,1,4,1)*ProdC(5,1,4,1)
	WProd1047 = Prod(5,2,4,1)*ProdC(5,1,4,1)
        WProd987 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1049 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd975 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1043 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1050 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd986 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1027 = Prod(5,1,3,1)*Prod(5,2,2,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1028 = Prod(5,1,2,1)*Prod(5,2,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd929 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,5,2)
	WProd922 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd951 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd953 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
        WProd932 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(5,2,1,1)
	WProd937 = Prod(5,1,2,1)*ProdC(5,2,2,1)
	WProd938 = Prod(5,1,3,1)*ProdC(5,2,3,1)
	WProd971 = Prod(5,1,4,1)*ProdC(5,2,4,1)
        WProd977 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd965 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd978 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd933 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd991 = Prod(6,2,6,1)*ProdC(5,2,5,1)
	WProd970 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd1009 = Prod(6,2,1,1)*ProdC(6,1,1,1)
        WProd992 = Prod(3,1,1,1)*Prod(6,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd993 = -(Prod(2,1,1,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,1,1))
        WProd989 = -(Prod(5,1,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,1,1,1))
        WProd994 = Prod(3,1,2,1)*Prod(6,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd995 = Prod(5,1,3,1)*Prod(6,2,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
        WProd1071 = Prod(5,2,3,1)*Prod(6,2,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
        WProd1056 = Prod(5,1,1,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,1,1)
        WProd996 = Prod(5,1,1,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,4,1)*ProdC(5,2,1,1)*ProdC(6,1,1,1)
        WProd1019 = Prod(5,1,3,1)*Prod(6,2,2,1)*ProdC(5,2,4,1)*
     -   ProdC(6,1,1,1)
	WProd976 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd1010 = Prod(6,2,2,1)*ProdC(6,1,2,1)
        WProd990 = -(Prod(5,1,2,1)**2*Prod(6,2,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,1,2,1))
        WProd1057 = Prod(5,1,2,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1)
        WProd997 = Prod(5,1,2,1)**2*Prod(6,2,3,1)*ProdC(5,1,4,1)*
     -   ProdC(5,2,1,1)*ProdC(6,1,2,1)
	WProd998 = Prod(6,2,3,1)*ProdC(6,1,3,1)
        WProd1072 = Prod(5,1,3,1)*Prod(5,2,3,1)*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,3,1)
        WProd1020 = Prod(5,1,3,1)**2*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,1,3,1)
	WProd950 = Prod(6,1,4,1)*ProdC(6,1,4,1)
	WProd999 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd1000 = Prod(3,1,2,1)*Prod(6,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1014 = -(Prod(4,1,3,1)*Prod(6,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1012 = Prod(4,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1001 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1058 = Prod(5,2,2,1)*Prod(6,2,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1059 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,3,1)*ProdC(6,1,4,1)
        WProd1073 = Prod(5,1,4,1)*Prod(5,2,3,1)*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,4,1)
        WProd1066 = Prod(5,1,4,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,4,1)
        WProd1002 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(5,2,1,1)*
     -   ProdC(6,1,4,1)
        WProd1003 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,3,1)*ProdC(6,1,4,1)
        WProd1021 = Prod(5,1,3,1)*Prod(5,1,4,1)*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,1,4,1)
        WProd1013 = Prod(5,1,2,1)*Prod(5,1,4,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,1,4,1)
	WProd984 = Prod(6,1,5,1)*ProdC(6,1,5,1)
	WProd1070 = Prod(6,2,5,2)*ProdC(6,1,5,1)
        WProd1022 = -(Prod(5,1,3,1)*Prod(6,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,5,1))
        WProd1004 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1067 = Prod(5,2,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1060 = -(Prod(5,1,2,1)*Prod(5,2,3,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,1,5,1))
        WProd1061 = Prod(3,1,1,1)*Prod(5,2,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,1,5,1)
        WProd1062 = -(Prod(5,1,2,1)*Prod(5,2,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,1,5,1))
        WProd1063 = Prod(3,1,2,1)*Prod(5,2,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,1,5,1)
        WProd1068 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,1,5,1)
        WProd1074 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
        WProd1069 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
        WProd1064 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
        WProd1065 = Prod(5,1,2,1)*Prod(5,2,3,1)*Prod(6,2,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
	WProd1018 = Prod(6,2,5,1)*ProdC(6,1,5,2)
        WProd1015 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,2)
        WProd1005 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,1,5,2))
        WProd1006 = Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,1,5,2)
        WProd1007 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,1,5,2))
        WProd1008 = Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,1,5,2)
        WProd1016 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,1,5,2)
        WProd1023 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,2)
        WProd1017 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,2)
        WProd1011 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,2)
	WProd952 = Prod(6,1,1,1)*ProdC(6,2,1,1)
        WProd923 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,1,1)
        WProd924 = -(Prod(2,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,1,1))
        WProd919 = -(Prod(5,1,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,2,1,1))
        WProd925 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd926 = Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
        WProd1051 = Prod(5,2,3,1)*Prod(6,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
        WProd1025 = Prod(5,1,1,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,1,1)
        WProd927 = Prod(5,1,1,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,4,1)*ProdC(5,2,1,1)*ProdC(6,2,1,1)
        WProd979 = Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,2,4,1)*
     -   ProdC(6,2,1,1)
	WProd954 = Prod(6,1,2,1)*ProdC(6,2,2,1)
        WProd920 = -(Prod(5,1,2,1)**2*Prod(6,1,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,2,2,1))
        WProd1026 = Prod(5,1,2,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1)
        WProd928 = Prod(5,1,2,1)**2*Prod(6,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(5,2,1,1)*ProdC(6,2,2,1)
	WProd930 = Prod(6,1,3,1)*ProdC(6,2,3,1)
        WProd1052 = Prod(5,1,3,1)*Prod(5,2,3,1)*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,3,1)
        WProd980 = Prod(5,1,3,1)**2*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -   ProdC(5,2,4,1)*ProdC(6,2,3,1)
	WProd934 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd935 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd966 = -(Prod(4,1,3,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd963 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd936 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd1031 = Prod(5,2,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd1032 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,3,1)*ProdC(6,2,4,1)
        WProd1053 = Prod(5,1,4,1)*Prod(5,2,3,1)*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,4,1)
        WProd1042 = Prod(5,1,4,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,4,1)
        WProd940 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,2,1,1)*
     -   ProdC(6,2,4,1)
        WProd941 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,3,1)*ProdC(6,2,4,1)
        WProd981 = Prod(5,1,3,1)*Prod(5,1,4,1)*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,2,4,1)
        WProd964 = Prod(5,1,2,1)*Prod(5,1,4,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,2,4,1)
	WProd1048 = Prod(6,1,5,2)*ProdC(6,2,5,1)
        WProd982 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,5,1))
        WProd942 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,1)
        WProd1044 = Prod(5,2,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,1)
        WProd1033 = -(Prod(5,1,2,1)*Prod(5,2,3,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,2,5,1))
        WProd1034 = Prod(3,1,1,1)*Prod(5,2,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,2,5,1)
        WProd1035 = -(Prod(5,1,2,1)*Prod(5,2,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,2,5,1))
        WProd1036 = Prod(3,1,2,1)*Prod(5,2,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,2,5,1)
        WProd1045 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,2,5,1)
        WProd1054 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
        WProd1046 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
        WProd1037 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
        WProd1038 = Prod(5,1,2,1)*Prod(5,2,3,1)*Prod(6,1,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
	WProd972 = Prod(6,1,5,1)*ProdC(6,2,5,2)
        WProd967 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,2)
        WProd943 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,1,1)*
     -     ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,2,5,2))
        WProd944 = Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,2,5,2)
        WProd945 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,2,5,2))
        WProd946 = Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,2,5,2)
        WProd968 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,2,5,2)
        WProd983 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,2)
        WProd969 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,2)
        WProd962 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,2)
	WProd921 = ProdC(6,2,6,1)
	WProd1039 = Prod(5,2,5,1)*ProdC(6,2,6,1)
	WProd947 = ProdC(5,2,5,1)*ProdC(6,2,6,1)
      
      
	Prp126 = -(1.D0/WProd939)
        Prp127 = 1.D0/
     -   ((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd939 + WProd948 - 
     -     WProd949)
	Prp129 = -(1.D0/WProd949)
	Prp128 = 1.D0/WProd950
	Prp130 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd959)
	Prp132 = 1.D0/(-WProd955 - WProd957 + WProd959)
	Prp131 = -(1.D0/WProd970)
	Prp133 = -(1.D0/WProd976)
	Prp134 = 1.D0/(-WProd939 - WProd970 + WProd984)
      
      
      
      endif

      end


