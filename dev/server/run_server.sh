#!/usr/bin/env bash
# kill previous port forwarding to ISIS server
ps -fU tothsa | grep "ssh -L" | grep "13001:" | awk '{print $2}' | xargs kill
ssh -L 13001:localhost:13001 -f -N isis
ssh isis '/home/lvi05884/spinw_server/spinw_server.sh /home/lvi05884/spinw_server/cache 32 13001'