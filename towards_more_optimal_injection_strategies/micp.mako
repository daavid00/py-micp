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
18  'Jan' 1991 /

EQLDIMS
/

WELLDIMS
12 60 12 20 /

TABDIMS
/

UNIFOUT

----------------------------------------------------------------------------
GRID
----------------------------------------------------------------------------

INIT

INCLUDE
'full_single_leak.grdecl' /

PORO
% for i in range(len(perm)):
${round(float(str(phi).lstrip('[').rstrip(']')), 4)}
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

PVTW
277.0      1     0.0  ${muw / 1e-3}  0.0 /

ROCK
277.0 1E-6 /

DENSITY
${rhow}/

MICPPARA
${rhob} ${rhoc} ${kstr * 1e-3} ${crit} ${eta} ${ko} ${ku} ${mu * 86400} ${comax}
${cumax} ${muu * 86400} ${ka * 86400} ${kd * 86400} ${kmin / 9.87e-16} ${F}
${tolclg} ${Y} /

----------------------------------------------------------------------------
SOLUTION
----------------------------------------------------------------------------

PRESSURE
${xx * yy * zz}*270 /

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
'INJE01' 'INJE' ${ix} ${iy} 1* 'WATER' 0.15/
'INJE02' 'INJE' ${ix} ${iy} 1* 'WATER' 0.15/
'INJE03' 'INJE' ${ix} ${iy} 1* 'WATER' 0.15/
'INJE04' 'INJE' ${ix} ${iy} 1* 'WATER' 0.15/
'INJE05' 'INJE' 1 1 1* 'WATER' 0.15/
'INJE06' 'INJE' 1 ${yy} 1* 'WATER' 0.15/
'INJE07' 'INJE' ${xx} 1 1* 'WATER' 0.15/
'PROD01' 'PROD' 2 2 1* 'WATER' 0.15/
'PROD02' 'PROD' 2 ${yy} 1* 'WATER' 0.15/
'PROD03' 'PROD' ${xx} 2 1* 'WATER' 0.15/
'PROD04' 'PROD' ${xx} ${yy} 1* 'WATER' 0.15/
/

COMPDAT
'INJE01' ${ix} ${iy} ${zz - il + 1} ${zz - il + Wc} 'OPEN' 1* 1*/
'INJE02' ${ix} ${iy} ${zz - il + Wc + 1} ${zz - il + 2 * Wc} 'OPEN' 1* 1*/
'INJE03' ${ix} ${iy} ${zz - il + 2 * Wc + 1} ${zz - il + 3 * Wc} 'OPEN' 1* 1*/
'INJE04' ${ix} ${iy} ${zz - il + 3 * Wc + 1} ${zz} 'OPEN' 1* 1*/
'INJE05' 1 1 ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'INJE06' 1 ${yy} ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'INJE07' ${xx} 1 ${zz - il + 1} ${zz} 'OPEN' 1* 1*/
'PROD01' 2 2 1 ${it} 'OPEN' 1* 1*/
'PROD02' 2 ${yy} 1 ${it} 'OPEN' 1* 1*/
'PROD03' ${xx} 2 1 ${it} 'OPEN' 1* 1*/
'PROD04' ${xx} ${yy} 1 ${zz} 'OPEN' 1* 1*/
/

% for i in range(len(injestra)):
WCONINJE
'INJE01' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 0.1 * float(Wc) / il : E}"}  1* 10000/
'INJE02' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 0.1 * float(Wc) / il : E}"}  1* 10000/
'INJE03' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 0.1 * float(Wc) / il : E}"}  1* 10000/
'INJE04' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 0.1 * float(il - 3 * Wc) / il : E}"} 1* 10000/
'INJE05' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 1. : E}"} 1* 100000/
'INJE06' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 1. : E}"} 1* 100000/
'INJE07' 'WATER' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'}
'RATE' ${f"{injestra[i][1] * 86400. * 1. : E}"} 1* 100000/
/
WCONPROD
'PROD01' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} 'BHP' 5* 2/
'PROD02' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} 'BHP' 5* 2/
'PROD03' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} 'BHP' 5* 2/
'PROD04' ${'OPEN' if injestra[i][1] > 0 else 'SHUT'} 'BHP' 5* 2/
/
WMICP
'INJE01' ${injestra[i][2]} 0 ${injestra[i][4]}/
'INJE02' 0 0 0/
'INJE03' 0 ${injestra[i][3]} 0/
'INJE04' 0 0 0/
'INJE05' 0 0 0/
'INJE06' 0 0 0/
'INJE07' 0 0 0/
/
TSTEP
${injestra[i][0] * 3600. / 86400}
/

% endfor
