-- Copyright (C) 2021-2025 NORCE Research AS
----------------------------------------------------------------------------
RUNSPEC
----------------------------------------------------------------------------
DIMENS
${xx} ${yy} ${zz} /

WATER
MICP

METRIC

START
1 JAN 2025 /

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
GRID.INC /

PORO
% for value in poro:
${value}\
% endfor
/

PERMX
% for value in perm:
${value}\
% endfor
/

COPY 
PERMX PERMY /
PERMX PERMZ /
/
----------------------------------------------------------------------------
PROPS
----------------------------------------------------------------------------
PERMFACT
-- PORO PERM
-- FACTOR FACTOR
% for pofa,pefa in zip(porofac,permfac):
${pofa} ${pefa}
% endfor
/

PVTW
277 1.038 4.67E-05 ${muw / 1e-3} 0 /

ROCK
276.0 1.95E-04 /

DENSITY
${rhow} /

BIOFPARA
-- BDEN DEATH GROWTH HALF YIELD FACTOR ATACH DETRAT DETEXP UREA HALF CDEN YUTOC
${rhob} ${kd * 86400} ${mu * 86400} ${ko} ${Y} ${F} ${ka * 86400} ${kstr * 1e-3} ${detexp} ${muu * 86400} ${ku} ${rhoc} ${yurca} /
----------------------------------------------------------------------------
SOLUTION
----------------------------------------------------------------------------
EQUIL
0 200 -1 0 0 0 1 1 0 /

SMICR
${xx * yy * zz}*0 /

SOXYG
${xx * yy * zz}*0 /

SUREA
${xx * yy * zz}*0 /

SBIOF
${xx * yy * zz}*0 /

SCALC
${xx * yy * zz}*0 /
----------------------------------------------------------------------------
SCHEDULE
----------------------------------------------------------------------------
RPTRST
BASIC=2 /

WELSPECS
INJE01 INJE ${ix} ${iy} 1* WATER 0.15 /
INJE02 INJE ${ix} ${iy} 1* WATER 0.15 /
INJE03 INJE ${ix} ${iy} 1* WATER 0.15 /
INJE04 INJE ${ix} ${iy} 1* WATER 0.15 /
PROD01 PROD 1 1 1* WATER 0.15 /
PROD02 PROD 1 ${yy} 1* WATER 0.15 /
PROD03 PROD ${xx} 1 1* WATER 0.15 /
PROD04 PROD ${xx} ${yy} 1* WATER 0.15 /
/

COMPDAT
INJE01 ${ix} ${iy} ${zz - il + 1} ${zz - il + Wc} OPEN 2* /
INJE02 ${ix} ${iy} ${zz - il + Wc + 1} ${zz - il + 2 * Wc} OPEN 2* /
INJE03 ${ix} ${iy} ${zz - il + 2 * Wc + 1} ${zz - il + 3 * Wc} OPEN 2* /
INJE04 ${ix} ${iy} ${zz - il + 3 * Wc + 1} ${zz} OPEN 2* /
PROD01 1 1 1 ${zz} OPEN 2* /
PROD02 1 ${yy} 1 ${zz} OPEN 2* /
PROD03 ${xx} 1 1 ${zz} OPEN 2* /
PROD04 ${xx} ${yy} 1 ${zz} OPEN 2* /
/

% for i in range(len(injestra)):
WCONINJE
INJE01 WATER ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} RATE ${f"{injestra[i][1] * 86400. * float(Wc) / il : E}"} 1* 480 /
INJE02 WATER ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} RATE ${f"{injestra[i][1] * 86400. * float(Wc) / il : E}"} 1* 480 /
INJE03 WATER ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} RATE ${f"{injestra[i][1] * 86400. * float(Wc) / il : E}"} 1* 480 /
INJE04 WATER ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} RATE ${f"{injestra[i][1] * 86400. * float(il - 3 * Wc) / il : E}"} 1* 480 /
/
WCONPROD
PROD01 ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} BHP 5* 2 /
PROD02 ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} BHP 5* 2 /
PROD03 ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} BHP 5* 2 /
PROD04 ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} BHP 5* 2 /
/
WMICP
INJE01 ${injestra[i][2]} 0 ${injestra[i][4]} /
INJE02 0 0 0 /
INJE03 0 ${injestra[i][3]} 0 /
INJE04 0 0 0 /
/
TSTEP
${injestra[i][0] * 3600. / 86400} /

% endfor
