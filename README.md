# Deep learning project Summer 2016


* [models](https://github.com/semccomas/DL_proj16/tree/master/models) contains all versions of Keras models
   * model12.py is the latest working model, includes for predictions for 3 state SS, 8 state SS, and RSA, including 3D plots, classification reports, confusion matrices, plotting loss, and doing all metrics per protein prediction as well
   * most other models in here also work but don't include the same architecture or predictions

* [auto_parse_all.sh](https://github.com/semccomas/DL_proj16/blob/master/auto_parse_all.sh) runs [parse_all.py](https://github.com/semccomas/DL_proj16/blob/master/parse_all.py), is for only a small amount of data and does not use pytables yet. Older version

* [auto_parse_all_bigdataset.sh](https://github.com/semccomas/DL_proj16/blob/master/auto_parse_all_bigdataset.sh) uses [parse_all_bigdataset.py](https://github.com/semccomas/DL_proj16/blob/master/parse_all_bigdataset.py) to generate Pytables for all pdb's in /data, need to specify if using jsE3 alignment or jsE0 alignment
   * The two python scripts are practically the same, moved the Pytables around a bit to be more logical for reading into the model
   * These scripts take all the protein PDB files and make them into one large Pytable, which contains information on the amino acid sequence in a sliding table (so as to include information on surrounding amino acids in prediction), the MSA scores for each amino acid and does one hot encoding for 3 and 8 states of the secondary structure
   * This Pytable is then read into the Keras models
   * Several GB of PDB structures so this takes some time

* [ECOD](https://github.com/semccomas/DL_proj16/tree/master/ECOD) stands for [evolutionary classification of domains](http://prodata.swmed.edu/ecod/) and was started as another means of classifying proteins, never ended up going down this path but good for future use

* [script_readme.txt](https://github.com/semccomas/DL_proj16/blob/master/script_readme.txt) for transferring script data and how-to at the end of the project




---------
Notes to myself
Summer 2016 deep learning project 
https://docs.google.com/document/d/1GqQ0Aph0KKc2wGvHU8Tw-_0zMg3WM-L7LIE-zvo4Uog/edit


dssp_validation and dssp_values come from encoder.py and 3state_all.py
dssp_validation_not_encoded and dssp_values_not_encoded come from old_3state_pase.py and the 3state_all.py

bloop.py and 3state_model.py are keras models

FileLen_dssp_vs_fa comes from 3state_parse.py and auto_parse.py
(the 2 files in comparision come from data/training_set and dssp_parsed 

dssp_parsed and dssp_parsed_test come from 3state_all.sh and are just files run from data/training or data/testing on dssp


encoded_slidingtable.py takes the data/training_data as an input 1, dssp_parsed/* as input 2, and it outputs 3state_features as sys.argv3, and aa_profiles as argv 4
pytable_array.py takes the 3state_feature/* as the input 1 and aa_profiles/* as input 2. The resulting table is 'big_table' or you can turn on sys.argv3 for making files in the pytable_3state directory
		 running the command < h5ls -rd big_table > gives you some table info. <ptdump big_table> or <pttree big_table> will also do something similar

auto_parse.sh does the encoded_sliding table and pytable array on each matching dssp and data/training file. 

testingx