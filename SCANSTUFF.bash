#!/bin/bash

name="compression_25"
data="/home/steve/Desktop/local_scan/ldf"

rm -f $name.his $name.dat $name.drr $name.list $name.log

for i in `ls -tr $data/a135feb_12.ldf`
#for i in `ls -tr $data/a135feb_12*.ldf`
do
    #if [ "$i" == "$data/a135feb_12.ldf" ];
    #then
    #continue
    #fi


    cmd=$cmd"file $i\ngo\n"

done
cmd=$cmd"end\n"
#echo -e $cmd
echo -e $cmd | ./pixie_ldf_c $name
