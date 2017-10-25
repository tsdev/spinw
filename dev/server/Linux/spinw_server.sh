#!/bin/sh
# Call this script to run SpinW Server
#
# Modify the MCRROOT variable below to point to the location of the Matlab
# Runtime (R2017b) installed on your system.
#
# To start SpinW Server use the following command:
# ./spinw_server.sh pathToMatFiles numWorkers portNumTCPIP
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
