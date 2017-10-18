#!/bin/sh
# Call this script from python to start pySpinW
#
# Modify the MCRROOT variable below to point to the location of the Matlab
# Runtime (R2017a) installed on your system.
#
# To start pySpinW use the following Python commands:
# from transplant import Matlab
# m = Matlab(executable='full path to pyspinw.sh')
# m.disp('Hello World!')
#

MCRROOT=/Applications/MATLAB/MATLAB_Runtime/v93

exe_name=$0
exe_dir=`dirname "$0"`

DYLD_LIBRARY_PATH=.:${MCRROOT}/runtime/maci64 ;
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/bin/maci64 ;
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/sys/os/maci64;
export DYLD_LIBRARY_PATH;

shift 1
args=
while [ $# -gt 0 ]; do
    token=$1
    args="${args} \"${token}\"" 
    shift
done
eval "\"${exe_dir}/Contents/MacOS/spinw_server\"" $args

exit
