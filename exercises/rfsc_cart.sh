Zp=2
Np=2
Zt=79
Nt=118

max_imp_param=200
db=10
echo $Zp $Np $Zt $Nt -$max_imp_param | ./test_rfsc > rfsc_cart.txt
imp_param_list=$(echo -$((max_imp_param-db)) $db $max_imp_param  | xargs seq| paste -sd ' ')

for b in $imp_param_list
do
	echo $Zp $Np $Zt $Nt $b | ./test_rfsc | awk '{print $2, $3}' > temp.txt
	paste rfsc_cart.txt temp.txt > new.txt
	mv new.txt rfsc_cart.txt
#	echo $b
done
rm temp.txt

