for i in data/train_input/*.fa; do  
basei=`basename $i` 
for j in dssp_parsed/*.dssp; do 
basej=`basename $j` 
newi=${basei:0:5} 
newj=${basej:0:5}
if [ $newi = $newj ] ; then
python 3state_parse.py $i $j 
echo $newi $newj
fi
done
done
