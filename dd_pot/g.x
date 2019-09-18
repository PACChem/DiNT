#!/bin/csh
if (-e Test.FChk) then
mv -f Test.FChk Test.FChk_1
endif
g03 < qc.in > qc.out
