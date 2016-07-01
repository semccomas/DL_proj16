for i in data/test_input/*.fa; do  
basei=`basename $i` 
for j in dssp_parsed_test/*.dssp; do 
basej=`basename $j` 
newi=${basei:0:5} 
newj=${basej:0:5}
if [ $newi = $newj ] ; then
python encoded_slidingtable_3state.py $i $j 3state_test_features/$newi.feat aa_test_profiles/$newi.window
python pytable_array.py 3state_test_features/$newi.feat aa_test_profiles/$newi.window #pytable_3state/$newi-table

echo $newi $newj
fi
done
done
