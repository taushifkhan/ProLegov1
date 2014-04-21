ProLegov1
=========

Stores files and codes require to build contact string of alpha helices using chain-helix.c

+++++++++|  Read me file for contact strin building and analysis  |++++++++++
++++++++++++|  From structure selection --> pattern  |+++++++++++++++++++++++
++++++++++++++|  Date : 14-  April,2014  |+++++++++++++++++++++++++++++++++++

1. Get fasta file of filtered PDB chains from rcsb with required structure
quality.

2. Run Cd-Hit to cluster them with sequence identity cut-off with following
code: see readme_folder/cdhitOutfile

3. Get cluser representative and run seondary structure module on the same.
    run python code extractRep.py <cdhit_output file>
    outfile has three tag lines :
    sP: Percentage of secondary structure
    sD: Residue detail of secondary structure
    sN: Number of secondary structure in each protein chain
    $ less sseDetail* | grep 'sN' > sse_number*

4. Get pdb chains of required number of alpha helix content. min 3 max 10..
    Can be done by awk: 
    awk '{if ($6 >=3 && $6 <=10) print $0}' sse_number_13_4_14 | wc -l
    $1047 
    awk '{if ($6 >=3 && $6 <=10) print $1$2}' sse_number_13_4_14 
    > pdb_chain_3_10_14APR_14


5. Get statistics of data use ,
    see 1. get_dataset_statistics.ipynb

6. Run helix interaction code to get sse contact patter.
    %%% Use c-helix which analyze contact on the basis of protein chains. The
    %nameing of secondary structure (serial helix number) can be obtain in
    %following step.
    1. Save output of the code helix_shape.txt and helix_packing_pair.txt in 
    flder as these are important results.
    2. run $ python process_hpp.py <helix_Shape.txt> <helix_packing_pair.txt> >
    OnlyInteraction. REMOVE the Last line
    $python process_Hpp.py helix_shape.txt helix_packing_pair.txt > InteractionOnly

    3. Save Tab separted file of helix_shape.txt and into file helSpaceTab; Can
    be done by: awk '{print $1"\t"$2"\t"$3"\t"$4}' helix_shape.txt > helSpaceTab

    4.run the bash script which in turn runs python script to get pattern:
    bash extrModiPat.bash InteractionOnly
    > Out put pattern file will be saved in IntPatterns_ContactNear_added.txt

    ***Example is in /home/taushif/alphaWork/compareDatasets/HppChain***
    @ Use previously used c-helix code of chothia 
    @ get contact py python code in-house

7. Identify presence of sub-patterns and uniq patters per alpha helix content
8. Tabulate
