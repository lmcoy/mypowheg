      subroutine abbrrealQED16(resnr)
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

      
	WProd1399 = Prod(6,2,6,1)
	WProd1336 = Prod(2,2,2,1)*Prod(6,2,6,1)
	WProd1281 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1316 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd1382 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd1403 = Prod(6,2,6,1)*ProdC(2,2,2,1)
	WProd1283 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1284 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd1277 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd1368 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd1285 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1293 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd1331 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1333 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1372 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(4,1,1,1))
        WProd1393 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
        WProd1396 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,1,1))
	WProd1286 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd1279 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd1318 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1295 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1296 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1373 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
        WProd1320 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1334 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1328 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1322 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1321 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1335 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
	WProd1369 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd1374 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1397 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1388 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd1398 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd1287 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1288 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd1289 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd1280 = Prod(5,1,2,2)*ProdC(5,1,2,1)
	WProd1370 = Prod(5,1,2,1)*ProdC(5,1,2,2)
	WProd1290 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd1291 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd1329 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1323 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1389 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(2,2,1,1)*
     -     ProdC(5,1,4,1))
	WProd1313 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1356 = Prod(3,1,1,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd1337 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1338 = Prod(3,1,2,2)*Prod(6,2,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1357 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,2,5,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,1,1)
        WProd1411 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(6,1,1,1)
        WProd1404 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,2,2)*
     -   ProdC(6,1,1,1)
        WProd1400 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -     ProdC(5,1,4,1)*ProdC(6,1,1,1))
        WProd1339 = Prod(5,1,3,1)*Prod(6,2,5,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
	WProd1330 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd1340 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd1358 = Prod(6,2,2,2)*ProdC(6,1,2,1)
        WProd1359 = -(Prod(5,1,2,2)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,1))
        WProd1341 = Prod(3,1,2,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd1363 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,1))
        WProd1342 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1343 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1344 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1345 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1346 = Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,2,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1)
        WProd1360 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1347 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1348 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,2,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1))
        WProd1349 = Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,1)
        WProd1405 = Prod(3,1,2,1)**2*Prod(6,2,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,2)*ProdC(6,1,2,1)
        WProd1406 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,2)*ProdC(6,1,2,1)
        WProd1361 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,1))
	WProd1412 = Prod(6,2,2,1)*ProdC(6,1,2,2)
        WProd1416 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,2,2))
        WProd1420 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,1,2,2))
        WProd1413 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1407 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1408 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1417 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1409 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*Prod(6,2,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1410 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,1,2,2))
        WProd1418 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,2,2))
	WProd1350 = Prod(6,2,3,1)*ProdC(6,1,3,1)
        WProd1401 = Prod(3,1,2,1)**2*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(6,1,3,1)
	WProd1351 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd1352 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1362 = -(Prod(5,1,2,2)*Prod(6,2,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1419 = -(Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(2,2,1,1)*
     -     ProdC(6,1,4,1))
        WProd1414 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(3,1,1,1)*ProdC(6,1,4,1))
        WProd1353 = -(Prod(5,1,3,1)*Prod(6,2,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1415 = -(Prod(4,1,2,1)*Prod(5,1,2,1)*Prod(6,2,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,1,4,1))
        WProd1354 = -(Prod(5,1,4,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1402 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(6,1,4,1)
        WProd1355 = -(Prod(4,1,3,1)*Prod(6,2,5,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
	WProd1314 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1315 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,1,1)
        WProd1278 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1282 = Prod(3,1,2,2)*Prod(6,1,5,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1317 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,5,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,1,1)
        WProd1383 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(6,2,1,1)
        WProd1371 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,2,2)*
     -   ProdC(6,2,1,1)
        WProd1364 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -     ProdC(5,1,4,1)*ProdC(6,2,1,1))
        WProd1292 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
	WProd1294 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd1319 = Prod(6,1,2,2)*ProdC(6,2,2,1)
        WProd1324 = -(Prod(5,1,2,2)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,1))
        WProd1297 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,2,1)
        WProd1332 = -(Prod(3,1,1,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,1))
        WProd1298 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1299 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1300 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1301 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1302 = Prod(3,1,2,1)*Prod(3,1,2,2)*Prod(6,1,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1)
        WProd1325 = -(Prod(4,1,3,1)*Prod(5,1,2,2)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1303 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1304 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*Prod(6,1,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1))
        WProd1305 = Prod(3,1,2,2)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,1)
        WProd1375 = Prod(3,1,2,1)**2*Prod(6,1,5,1)*
     -   ProdC(3,1,1,1)*ProdC(4,1,2,2)*ProdC(6,2,2,1)
        WProd1376 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(4,1,2,2)*ProdC(6,2,2,1)
        WProd1326 = -(Prod(5,1,2,2)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,1))
	WProd1384 = Prod(6,1,2,1)*ProdC(6,2,2,2)
        WProd1390 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,2,2))
        WProd1395 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(6,2,2,2))
        WProd1385 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1377 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1378 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1391 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1379 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*Prod(6,1,2,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1380 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,4,1)*
     -     ProdC(4,1,1,1)*ProdC(4,1,2,1)*ProdC(6,2,2,2))
        WProd1392 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,2,2))
	WProd1306 = Prod(6,1,3,1)*ProdC(6,2,3,1)
        WProd1365 = Prod(3,1,2,1)**2*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(6,2,3,1)
	WProd1307 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd1308 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1327 = -(Prod(5,1,2,2)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1394 = -(Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(2,2,1,1)*
     -     ProdC(6,2,4,1))
        WProd1386 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(3,1,1,1)*ProdC(6,2,4,1))
        WProd1309 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1387 = -(Prod(4,1,2,1)*Prod(5,1,2,1)*Prod(6,1,3,1)*
     -     ProdC(4,1,1,1)*ProdC(6,2,4,1))
        WProd1310 = -(Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1366 = Prod(3,1,2,1)*Prod(4,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(6,2,4,1)
        WProd1311 = -(Prod(4,1,3,1)*Prod(6,1,5,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
	WProd1367 = ProdC(6,2,6,1)
	WProd1312 = Prod(2,2,2,1)*ProdC(6,2,6,1)
	WProd1381 = ProdC(2,2,2,1)*ProdC(6,2,6,1)
      
      
	Prp164 = 1.D0/WProd1281
	Prp157 = -(1.D0/WProd1286)
	Prp161 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1287)
	Prp163 = 1.D0/(-WProd1283 - WProd1285 + WProd1287)
	Prp162 = -(1.D0/WProd1289)
	Prp158 = -(1.D0/WProd1313)
	Prp160 = 1.D0/WProd1314
        Prp159 = 1.D0/
     -   ((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd1288 - WProd1313 + 
     -     WProd1314)
	Prp165 = 1.D0/(WProd1281 - WProd1313 - WProd1330)
      
      
      
      endif

      end


