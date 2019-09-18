set OMP_THREAD_LIMIT=1
export OMP_THREAD_LIMIT=1
cd a1
cd ../a33
time ./dint.x.opt < input >  output &
cd ../a34
time ./dint.x.opt < input >  output &
cd ../a35
time ./dint.x.opt < input >  output &
cd ../a36
time ./dint.x.opt < input >  output &
