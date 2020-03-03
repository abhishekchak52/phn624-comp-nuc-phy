#!/usr/bin/env bash

Zp=2
Np=2
Zt=79
Nt=118

max_imp_param=200
db=5
touch rfsc_pol.txt
imp_param_list=$(echo $db $db $max_imp_param  | xargs seq| paste -sd ' ')

for b in $imp_param_list
do
	echo $Zp $Np $Zt $Nt $b | ./rfsc_pol | awk '{print $1, $2}' > temp.txt
	paste rfsc_pol.txt temp.txt > new.txt
	mv new.txt rfsc_pol.txt
    sed -i '$d' rfsc_pol.txt
done
rm temp.txt
