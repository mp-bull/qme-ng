#!/bin/bash
n=`cat $1|grep -A 1 Number|tail -n 1`
np1=`echo $n + 1 |bc`
cat $1 |grep -A $np1 Matrix |tail -n `cat $1|grep -A 1 Number|tail -n 1`

