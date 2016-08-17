read -r -p "Are you sure???? [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
	rm table_all_jhE3

	for i in data/train_input/*.fa; do  
	basei=`basename $i` 
	for j in data/train_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	
	echo $i $j $newi.fa.jhE3
	python parse_all.py $i $j data/train_input/$newi.fa.jhE3.aln table_all_jhE3

	fi
	done
	done

	rm test_table_all_jhE3

	for i in data/test_input/*.fa; do  
	basei=`basename $i` 
	for j in data/test_input/*.pdb; do 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then

	echo $i $j
	python parse_all.py $i $j data/test_input/$newi.fa.jhE3.aln test_table_all_jhE3

	fi
	done
	done


else
    echo 'cool'
fi