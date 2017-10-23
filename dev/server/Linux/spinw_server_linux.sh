#!/bin/sh
# Call this script to start SpinW server
#
# Modify the MCRROOT variable below to point to the location of the Matlab
# Runtime (R2017a) installed on your system.
#
# To start pySpinW use the following Python commands:
# from transplant import Matlab
# m = Matlab(executable='full path to pyspinw.sh')
# m.disp('Hello World!')
#

MCRROOT=/MATLAB/MATLAB_Runtime/v93

exe_name=$0
exe_dir=`dirname "$0"`

LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
export LD_LIBRARY_PATH;
args=
while [ $# -gt 0 ]; do
    token=$1
    args="${args} \"${token}\"" 
    shift
done
eval "\"${exe_dir}/spinw_server\"" $args

exit
