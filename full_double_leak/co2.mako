-- Copyright (C) 2021-2022 NORCE

----------------------------------------------------------------------------
RUNSPEC
----------------------------------------------------------------------------

DIMENS
${xx} ${yy} ${zz} /

OIL
GAS
CO2STORE

METRIC

START
18 'Jan' 1991 /

EQLDIMS
/

WELLDIMS
10 60 10 20 /

TABDIMS
/

UNIFOUT

----------------------------------------------------------------------------
GRID
----------------------------------------------------------------------------

INIT

INCLUDE
'full_double_leak.grdecl' /

PORO
% for i in range(len(biof)):
${round(float(str(phi - biof[i] - calc[i]).lstrip('[').rstrip(']')), 4)}
% endfor
/

PERMX
% for i in range(len(perm)):
${round(float(str(perm[i]).lstrip('[').rstrip(']')), 4)}
% endfor
/

PERMY
% for i in range(len(perm)):
${round(float(str(perm[i]).lstrip('[').rstrip(']')), 4)}
% endfor
/

PERMZ
% for i in range(len(perm)):
${round(float(str(perm[i]).lstrip('[').rstrip(']')), 4)}
% endfor
/

----------------------------------------------------------------------------
PROPS
----------------------------------------------------------------------------

SGOF
0.00000 0.00000 1.00000 0.04935
0.03842 0.00028 0.80552 0.05070
0.07684 0.00221 0.64089 0.05217
0.11526 0.00725 0.50288 0.05377
0.15368 0.01670 0.38846 0.05554
0.19211 0.03165 0.29478 0.05749
0.23053 0.05304 0.21916 0.05966
0.26895 0.08159 0.15911 0.06209
0.30737 0.11786 0.11235 0.06485
0.34579 0.16222 0.07673 0.06802
0.38421 0.21485 0.05034 0.07170
0.42263 0.27576 0.03143 0.07605
0.46105 0.34475 0.01842 0.08130
0.49947 0.42146 0.00994 0.08781
0.53789 0.50534 0.00480 0.09619
0.57632 0.59564 0.00196 0.10755
0.61474 0.69146 0.00062 0.12419
0.65316 0.79168 0.00012 0.15210
0.69158 0.89502 0.00001 0.21509
0.73000 1.00000 0.00000 15.60463 /

SALINITY
0.7/ 35-40g/l  -> 35-40g /kg -> 0.63-0.72 mol/g

ROCK
1.0 1E-6 /

----------------------------------------------------------------------------
SOLUTION
----------------------------------------------------------------------------

EQUIL
0  200  1000   0.0  0  0.0   1   1   0 /

----------------------------------------------------------------------------
SCHEDULE
----------------------------------------------------------------------------

RPTRST
BASIC=2 /

WELSPECS
'INJE01' 'INJE' ${ix} ${iy} 1* 'GAS' 0.15/
'PROD01' 'PROD' 2 1 1* 'OIL' 0.15/
'PROD02' 'PROD' 2 ${yy} 1* 'OIL' 0.15/
'PROD03' 'PROD' ${xx - 1} 1 1* 'OIL' 0.15/
'PROD04' 'PROD' ${xx - 1} ${yy} 1* 'OIL' 0.15/
'PROD05' 'PROD' 1 2 1* 'OIL' 0.15/
'PROD06' 'PROD' 1 ${yy - 1} 1* 'OIL' 0.15/
'PROD07' 'PROD' ${xx} 2 1* 'OIL' 0.15/
'PROD08' 'PROD' ${xx} ${yy - 1} 1* 'OIL' 0.15/
/

COMPDAT
'INJE01' ${ix} ${iy} ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'PROD01' 2 1 ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'PROD02' 2 ${yy} ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'PROD03' ${xx - 1} 1 ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'PROD04' ${xx - 1} ${yy} ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'PROD05' 1 2 1 ${it} 'OPEN' 1* 1*/
'PROD06' 1 ${yy - 1} 1 ${it} 'OPEN' 1* 1*/
'PROD07' ${xx} 2 1 ${it} 'OPEN' 1* 1*/
'PROD08' ${xx} ${yy - 1} 1 ${it} 'OPEN' 1* 1*/
/

WCONINJE
'INJE01' 'GAS' 'OPEN' 'RATE' ${QCO2 * 86400}  1* 600/
/
WCONPROD
'PROD01' 'OPEN'  'BHP' 5* 2/
'PROD02' 'OPEN'  'BHP' 5* 2/
'PROD03' 'OPEN'  'BHP' 5* 2/
'PROD04' 'OPEN'  'BHP' 5* 2/
'PROD05' 'OPEN'  'BHP' 5* 2/
'PROD06' 'OPEN'  'BHP' 5* 2/
'PROD07' 'OPEN'  'BHP' 5* 2/
'PROD08' 'OPEN'  'BHP' 5* 2/
/
TSTEP
${NCO2}*${dtCO2}
/
