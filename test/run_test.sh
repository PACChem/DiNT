#!/bin/bash

../exe/dint-ch4hALL.x < CH4H_input > CH4H_output
mv fort.31 CH4H_fort.31
mv fort.70 CH4H_fort.70
../exe/dint-ch4oALL.x < CH4O_input > CH4O_output
mv fort.31 CH4O_fort.31
mv fort.70 CH4O_fort.70
