Info on the scripts and so on.
Everything is in ~/Desktop/SummerProj

1.Parsing the alignment and sequence information:

auto_parse_all_bigdataset.sh
	- This will run the below python script for every protein in /data/pdbcull/matching_proteins/train_data and /data/pdbcull/matching_proteins/test_data. There will therefore be 2 tables, one for the training set and one for the validation set. No need to specify anything else unless you want to change the table name, there are no arguments you have to pass with it, just have to confirm that you actually want to run the script.
	- I am a fan of verbosity so this does echo each protein that will go into the table, as well as the shape of the arrays (which should be all the same) that are in the table. If you don't want that just take away the print statements at the very end of the python script and the 'echo' in here.

parse_all_bigdataset.py
	- This script takes the protein sequence, pdb file, and alignment file together, and adds the rsa, phi and psi angles, 3 state sequence info, and 8 state sequence info. Resulting output is a single table with one group per protein, 6 nodes in each group. 
	- 4 arguments to pass along with it (used sys). Argv[1] = sequence information (*.fa). argv[2] = pdb file (*.pdb), argv[3] = alignment file (*.jhE0), argv[4] is the table output.
		- a note about the table, since I parse one sequence, pdb, alignment at a time I had to put the table in append mode. The bash script I wrote takes care of this (removes it before running the python script) but I just think it's worth noting
	- Again, the bash script above does all the work so to create this table all you have to do is execute the auto_parse script above.





2. Running the model:

	- All of the models I run are in ~/Desktop/SummerProj/models
	- The models are not particularly useful until model7 in my opinion. That is when I have the most data analysis and didn't do that thing where I was plotting incorrectly. Models 1-6 all include the big table as opposed to just parsing like 6 different tables. The models with other names are significantly older and I would not use them. Not very developed. I kind of moved in a chronological order. 