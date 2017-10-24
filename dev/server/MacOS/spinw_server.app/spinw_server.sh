#!/bin/sh
# Call this script to run SpinW Server
#
# Modify the MCRROOT variable below to point to the location of the Matlab
# Runtime (R2017b) installed on your system.
#
# To start SpinW Server use the following command:
# ./spinw_server.app/spinw_server.sh pathToMatFiles numWorkers portNumTCPIP
#

MCRROOT=/Applications/MATLAB/MATLAB_Runtime/v93

exe_name=$0
exe_dir=`dirname "$0"`

DYLD_LIBRARY_PATH=.:${MCRROOT}/runtime/maci64 ;
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/bin/maci64 ;
DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${MCRROOT}/sys/os/maci64;
export DYLD_LIBRARY_PATH;

args=
while [ $# -gt 0 ]; do
    token=$1
    args="${args} \"${token}\"" 
    shift
done
eval "\"${exe_dir}/Contents/MacOS/spinw_server\"" $args

exit
