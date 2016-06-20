###### IF YOU RUN THIS, REMOVE ALL FILES FROM DSSP_PARSED OR NOTHING WILL WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

input=/home/sarah/Desktop/SummerProj/data/test_input
dssp=/home/sarah/Desktop/SummerProj/dssp_parsed_test

for i in $input/*.pdb; do
base=`basename $i`
#echo $base
./dssp-2.0.4-linux-amd64 $i >> $dssp/$base.dssp    #as far as I can tell you have to be in this direcotry to execute otherwise it doesnt work with adding the $dssp var
python 3state_parse.py $dssp/$base.dssp
done
