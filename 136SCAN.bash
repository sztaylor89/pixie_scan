#!/bin/bash

name="sanibel136"
data="/scratch2/anl2015/FEB2015/136SB"

rm -f $name.his $name.dat $name.drr $name.list $name.log

#for i in `ls -tr $data/a136sb_06.ldf`
for i in `ls -tr $data/a136sb_06*.ldf`
do
    cmd=$cmd"file $i\ngo\n"
done

for i in `ls -tr $data/a136sb_07*.ldf`
do
    cmd=$cmd"file $i\ngo\n"
done

for i in `ls -tr $data/a136sb_08*.ldf`
do
    cmd=$cmd"file $i\ngo\n"
done

for i in `ls -tr $data/a136sb_09*.ldf`
do
    cmd=$cmd"file $i\ngo\n"
done

for i in `ls -tr $data/a136sb_10*.ldf`
do
    cmd=$cmd"file $i\ngo\n"
done

for i in `ls -tr $data/a136sb_11*.ldf`
do
    cmd=$cmd"file $i\ngo\n"
done

cmd=$cmd"end\n"
#echo -e $cmd

echo -e $cmd | ./pixie_ldf_c $name
