
read -r -p "Are you sure???? [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
	rm rsa_table

	for i in data/train_input/*.fa; do  
	basei=`basename $i` 
	for j in data/train_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then

	echo $i $j
	python parse_rsa.py $i $j rsa_features/$newi.feat rsa_table

	fi
	done
	done

	rm rsa_test_table

	for i in data/test_input/*.fa; do  
	basei=`basename $i` 
	for j in data/test_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then

	echo $i $j
	python parse_rsa.py $i $j rsa_test_features/$newi.feat rsa_test_table

	fi
	done
	done


else
    echo 'cool'
fi