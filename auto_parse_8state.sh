### this is the same as auto_parse for the most part but is designed to run on the 8 state tables 

rm 8state_table

for i in data/train_input/*.fa; do  
basei=`basename $i` 
for j in data/train_input/*.pdb ; do 
basej=`basename $j` 
newi=${basei:0:5} 
newj=${basej:0:5}
if [ $newi = $newj ] ; then
#echo $newi $newj
echo $i $j
python 8state_parse.py $i $j 8state_features/$newi.feat


fi
done
done



