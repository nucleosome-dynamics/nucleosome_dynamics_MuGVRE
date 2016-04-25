#!/bin/bash

logf="mem_logs"
proclog="proc_log"
scriptdir="/home/rilla/nucleServ/tests"

#> $logf
> $proclog

for f in $scriptdir/*.sh
do
    calcname=$(basename $f ".sh")
    echo $calcname >> $logf
    echo $calcname > $proclog
    echo $calcname
    /usr/bin/time -avo $logf $f 2> $proclog
    echo "========" >> $logf
done
