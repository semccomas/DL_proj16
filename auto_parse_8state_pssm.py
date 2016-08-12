# options are : pssm_8state_jhE0, pssm_8state_test_jhE0, pssm_8state_jhE3, pssm_8state_test_jhE3

### this is the same as auto_parse for the most part but is designed to run on the pssm tables 

rm pssm_test_8state_jhE0

for i in data/test_input/*.fa.jhE0.aln; do  
basei=`basename $i` 
for j in 8state_test_features/*.feat ; do 
basej=`basename $j` 
newi=${basei:0:5} 
newj=${basej:0:5}
if [ $newi = $newj ] ; then
#echo $newi $newj
echo $i $j
python parse_8state_pssm.py $i $j 


fi
done
done