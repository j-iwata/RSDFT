#!/bin/bash

read nproc

j=0
while [ $j -lt $nproc ]
do

   mkdir rank$j

   read a

   i=0
   while [ $i -lt $a ]
   do
      read b
      format=`printf "%05d" $b`
      touch rank$j/f.$format
#      ln -s /data/hp120301/k00755/out.2381056/wf.dat1.$format ./rank$j/wf.dat1.$format
      i=$((i+1))
   done

   j=$((j+1))

done


