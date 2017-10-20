#!/usr/bin/env bash
cd ~/spinw_git/dev/server
# get update on multiprocessor branch
git checkout multiprocessor
git pull origin multiprocessor
# compile the code in matlab
matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('~/spinw_git/dev'));compile_server_linux;exit"
# sync app to ISIS server
rsync -zhv ./Linux/spinw_server isis:/home/lvi05884/spinw_server/