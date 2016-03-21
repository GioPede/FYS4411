#! /usr/bin/bash

rm data.txt
touch data.txt
for i in 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70
do
echo $i
./variational-monte-carlo-fys4411 $i
{ echo $i; python blocking.py; } | tr "\n" " " >> data.txt
echo >> data.txt
done
