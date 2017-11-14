#!/bin/bash

str=""
for i in `seq 0 79`
do
    str="$str pdfDir/traj_$i.gif"
done

convert -delay 10 -loop 0 $str traj.gif