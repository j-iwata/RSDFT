#!/bin/bash

today=`date`
echo "Today is $today."

i=0
while [ $i -lt 512 ]
do
   mkdir rank$i
   i=$((i+1))
done

j=0
i=0
while [ $i -lt 512 ]
do
   while [ $j -lt $((i*8+8)) ]
   do
      format=`printf "%05d" $j`
      ln -s /data/hp120301/k00755/out.2214857/wf.dat1.$format ./rank$i/wf.dat1.$format
      j=$((j+1))
   done
   i=$((i+1))
done

