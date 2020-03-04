min=2.0
max=8.0
delta=0.1
seq $min $delta $max  > data.txt
KE_seq=$(echo $min $delta $max | xargs seq | paste -sd ' ')

for KE in $KE_seq
do 
	echo $KE | ./rfsc_mte >> rad.txt
done
paste data.txt rad.txt > new.txt
mv new.txt data.txt
rm rad.txt
#cat data.txt

