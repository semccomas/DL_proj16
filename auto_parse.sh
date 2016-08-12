for i in data/test_input/*.fa; do  
basei=`basename $i` 
for j in data/test_input/*.pdb; do 
basej=`basename $j` 
newi=${basei:0:5} 
newj=${basej:0:5}
if [ $newi = $newj ] ; then
python biopythondssp.py $i $j 3state_test_features/$newi.feat aa_test_profiles/$newi.window
python pytable_array.py 3state_test_features/$newi.feat aa_test_profiles/$newi.window #pytable_test_3state/$newi-table

echo $newi $newj
fi
done
done
