read -r -p "Are you sure???? [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
	rm table_all

	for i in data/train_input/*.fa; do  
	basei=`basename $i` 
	for j in data/train_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	
	echo $i $j $newi.fa.jhE0
	python parse_all.py $i $j data/train_input/$newi.fa.jhE0 table_all

	fi
	done
	done

	rm test_table_all

	for i in data/test_input/*.fa; do  
	basei=`basename $i` 
	for j in data/test_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then

	echo $i $j
	python parse_all.py $i $j data/test_input/$newi.fa.jhE0 test_table_all

	fi
	done
	done


else
    echo 'cool'
fi