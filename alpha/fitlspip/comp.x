/bin/rm *.o
gfortran -c -O3 -mcmodel=medium func.f read.f pythag.f svbksb.f svdcmp.f svdfit.f
gfortran *.o -l stdc++  -fopenmp -o fitlspip.x
