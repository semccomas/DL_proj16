### this is the same as auto_parse for the most part but is designed to run on the pssm tables 



read -r -p "Are you sure? [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
	rm pssm_test_table_jhE0

	for i in data/test_input/*.fa.jhE0.aln; do  
	basei=`basename $i` 
	for j in 3state_test_features/*.feat ; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	#echo $newi $newj
	echo $i $j
	python parse_pssm.py $i $j 


	fi
	done
	done

else
    echo 'cool'
fi



#### 4 tables come from this: pssm_table_jhE0, pssm_test_table_jhE0, pssm_table_jhE3, and pssm_test_table_jhE3. You have to change the name of the table in parse_pssm too 