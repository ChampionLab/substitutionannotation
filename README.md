# substitutionannotation
A repository for the script used to annotate and filter aa substitutions identified via mass-offset search in [MSFragger](https://fragpipe.nesvilab.org/). For more information on it's initial application, see the preprint [A Fit for Purpose Approach to Evaluate Detection of Amino Acid Substitutions in Shotgun Proteomics](https://doi.org/10.1101/2023.08.09.552645). The 'main' branch reflect the script as used for this publication. The 'tjlpub' branch has some updates, mostly to make the script more re-usable and robust to different file names or sample annotations.

MSFraggerFindSubs.py (updated to FindSubs_Fragpipe_PSM.py in the tjlpub branch) takes the output of a mass-offset search from MSFragger, and annotates PSMs identified with a mass-offset as a substitution or other peptide modification. It outputs a list of PSMs that correspond to peptide sequences with an amino acid substitution.

FindSSP.py takes the two tryptic peptide lists from [ProteaseGuru](https://github.com/smith-chem-wisc/ProteaseGuru) and finds sequences that differ by a single amino acid, accounting for changed peptide cleavage by removing the C-terminal lysine or arginine. 

Note that these scripts includes some hard coded paths. MSFraggerFindSubs.py in the main/publicationver hard codes the path to the output folder of MSFragger (variable 'inputdir'), and to the dangermods.csv and dfdm.csv included here. These hard coded paths need to be changed in the code in order to use this script. The tjlpub branch has some first steps to make this more python-package like and should work when you run FindSubs_Fragpipe_PSM.py, prompting you for the data directory. FindSSP.py is still hardcoded to reference two ProteaseGuru output peptide lists, use tryptic cleavage rules, and output file paths.

<h1>Download and tutorial instructions</h1>
You can download this code from Github by clicking on the green 'Code' button, then 'Download Zip'. Download speed will depend on your internet speed, but should be short (<1 min). [More help on downloading code from GitHub](https://docs.github.com/en/repositories/working-with-files/using-files/downloading-source-code-archives). These scripts work with a [Python version 3.10](https://www.python.org/downloads/) environment that has the packages [pandas](https://pandas.pydata.org/) version 2.1.0 and [numpy](https://numpy.org/install/) version 1.25.2.

To run MSFraggerFindSubs.py, you will need to open the script (either in an IDE or a simple text editing program such as Notepad) and update the hardcoded paths to the following:

-Line 21 The output folder of the MSFragger search. For the tutorial data, inputdir = 'your/download/path/substitutionannotation/tutorialdata'

-Line 35 The csv with the mass offset of each substitution type, 'dfdm.csv'. dfdm=pd.read_csv(r'your/download/path/substitutionannotation/dfdm.csv',index_col=0)

-Line 36 The csv with the mass offset of other post translational modifications, 'dangermods.csv'. dfdm=pd.read_csv(r'your/download/path/substitutionannotation/dangermods.csv',index_col=0)

Save the file, then run the script either from an IDE or in the terminal with 'python run your/download/path/substitutionannotation/MSFraggerFindSubs.py'. The results will be output to a FindSubsOutput folder within your MSFragger search result folder. This should take ~15 minutes with the tutorial data, depending on the speed of your computer.

To run FindSSP.py, you will likewise need to open and edit the hardcoded paths:

Line 25 - Path to the first ProteaseGuru peptide list, in the tutorial of salmonella peptides. dfsalty = pd.read_csv('your/download/path/substitutionannotation/tutorialdata/Salty_Peptides_1.tsv', sep='\t',usecols=['Base Sequence'])

Line 26 - Path to the second ProteaseGuru peptide list, in the tutorial of e. coli peptides. dfecoli = pd.read_csv('your/download/path/substitutionannotation/tutorialdata/Coli_Peptides_1.tsv', sep='\t',usecols=['Base Sequence'])

Line 30 - Output directory for all outputs of the script

Note that the outputs of FindSSP.py will name the output columns for peptide sequences 'ECOLI SSP Sequence' and the second list to 'SALTY SSP Sequence'. This script will likely take hours to run.
