read -r -p "Are you sure???? This will take a long time and will immediately delete the tables youve created (big table all and big test table ). Answer - [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
	rm big_table_all

	for i in /data/pdbcull/matching_proteins/train_set/*.fa; do  
	basei=`basename $i` 
	for j in /data/pdbfiles/train_data/*.pdb; do ## <<< change this to pdb files when you get the chance 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	
	echo $i $j 
	python parse_all_bigdataset.py $i $j /data/pdbcull/$newi.jhE0.aln big_table_all   ## <<<< note!!! -aln and .fasta files in the pdbcull are kinda different. aln has alignments only and fasta has both alignments and the prot. names of what aligned to 

	fi
	done
	done


### NOTE!!! I also currently dont have a validation set.... maybe take like 50 of the proteins from matching_proteins?? 
	rm big_test_table_all

	for i in /data/pdbcull/matching_proteins/test_set/*.fa; do  
	basei=`basename $i` 
	for j in /data/pdbfiles/train_data/*.pdb; do ## <<< change this to pdb files when you get the chance 
	basej=`basename $j` 
	newi=${basei:0:5} 
	newj=${basej:0:5}
	if [ $newi = $newj ] ; then
	
	echo $i $j 
	python parse_all_bigdataset.py $i $j /data/pdbcull/$newi.jhE0.aln big_table_all   ## <<<< note!!! -aln and .fasta files in the pdbcull are kinda different. aln has alignments only and fasta has both alignments and the prot. names of what aligned to 

	fi
	done
	done


else
    echo 'cool'
fi