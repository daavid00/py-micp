-- Copyright (C) 2021 NORCE

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
0	0	1	0.0
1	1 0 0.0 /

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
