#!/bin/bash	
for i in {1..5}
do
	echo "$i index"
	let start=$i*10-9
	let end=$i*10
	for (( seed=$start; seed<=$end; seed++ ))	
	do (ulimit -v 8388608 #16666666 #
	cpulimit -l 100 -i python3 ../src/QCQP.py 4 10 10 10 $seed 600&>output$seed.txt&
	)
	done &	
done	
#/dev/null