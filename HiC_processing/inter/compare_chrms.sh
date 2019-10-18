#!/bin/bash
compare_chrms()
{
	chrm1=$1
	chrm2=$2
	#echo $chrm1
	#echo $chrm2
	num1=${chrm1:4:5}
	num2=${chrm2:4:5}
	#echo $num1
	#echo $num2
	if [ "$num1" = "X" ]
	then
		return 1
	elif [ "$num2" = "X" ]
	then
		return 0
	elif [ $num1 -lt $num2 ]
	then
		return 0
	else
		return 1
	fi
}

