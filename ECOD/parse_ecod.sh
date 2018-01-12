
#### this script will first make a list of all the proteins in the specified directory. 
#### then it will run the python script that compares ecod to the proteins
#### what matches to ecod, one per x_class, will be compared to the actual directory, and 
### the .fa file will be copied to a new directory. You can use this in your auto_parse to compare to the pdb
## file for matches so that you dont have to move everything over. See auto_parse_all for how to do that 

read -r -p "Are you sure?? [y/n] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then

rm protein_names

for i in /data/pdbcull/*.fa; do
	basei=`basename $i` 
	newi=${basei:0:5} 
	echo $newi >> protein_names
done
echo 'protein_names written'

### this is like sys.argv 
filename="$1"
rm $filename
python parse_ecod2vals.py protein_names ecod.latest.fasta.txt  $filename  #filtered_ecod_array.txt
echo "python script done, output - $filename"

#sudo rm -r /data/pdbcull/matching_proteins
#sudo mkdir /data/pdbcull/matching_proteins

echo 'making new protein directory.... '
while IFS=';' read -ra line
do
    name="${line[0]}"
for i in /data/pdbcull/*.fa; do
	basei=`basename $i` 
	shorti=${basei:0:4}
	if  [[ "$shorti" = "$name" ]]; then
		echo $shorti 
		sudo cp $i /data/pdbcull/matching_proteins
	fi
done

done < "$filename"

else
    echo 'stopped'
fi