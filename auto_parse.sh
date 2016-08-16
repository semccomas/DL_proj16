
read -r -p "Are you sure? [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then

	rm big_table

	for i in data/train_input/*.fa; do  
	basei=`basename $i` 
	for j in data/train_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	python biopythondssp.py $i $j 3state_features/$newi.feat aa_profiles/$newi.window
	python pytable_array.py 3state_features/$newi.feat aa_profiles/$newi.window #pytable_3state/$newi-table

	echo $newi $newj
	fi
	done
	done
else:
	echo 'ok'
fi
