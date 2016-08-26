read -r -p "Are you sure???? This will take a long time and will immediately delete the tables youve created (big table all and big test table ). Answer - [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
	rm big_table_all

	for i in /data/pdbcull/matching_proteins/train_set/*.fa; do  
	basei=`basename $i` 
	for j in /data/pdbfiles/train_data/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	
	echo $i $j 
	python parse_all_bigdataset.py $i $j /data/pdbcull/$newi.jhE0.aln big_table_all 

	fi
	done
	done


	rm big_test_table_all

	for i in /data/pdbcull/matching_proteins/test_set/*.fa; do  
	basei=`basename $i` 
	for j in /data/pdbfiles/train_data/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	
	echo $i $j 
	python parse_all_bigdataset.py $i $j /data/pdbcull/$newi.jhE0.aln big_test_table_all  

	fi
	done
	done


else
    echo 'Ok!'
fi