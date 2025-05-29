# substitutionannotation
A repository for the script used to annotate and filter aa substitutions identified via mass-offset search in MSFragger
This branch (tjlpub) is an updated version of the code used in the publication (doi.org/10.1021/acs.jproteome.3c00730). It is updated to be compatible with Fragpipe V23 and applicable to new sets of data without rewriting code. That said it is far from python-package ready. 

resources>assets.py manages the filepath - when imported in python, you can use the function changeDir('directory') to set a new working directory. This directory should contain either the fragpipe outputs or the PEAKS outputs, neither in subfolders. Alternatively, change the data_directory.txt file under the temp folder manually.

Note there is some complexity in handling modified cysteine residues; I suggest you search with cysteine modification as a variable mod and not a fixed mod, as this allows IonQuant to properly quantify these peptidoforms. Leaving it as a fixed mod, the only quantification of these substitutions is in the psm.tsv outputs. Either way, the new script automatically identifies which method was used to apply the appropriate calculations for cysteine or substitution to/from cysteine.

This update adds a requirement for a paired protein coding sequence nucleotide FASTA, with similar protein headers to align aa/nucleotide data. This annotates substitutions with relevent codon information, and infers tRNA mismatches. Note this is based off the prokaryotic translation table for E. coli, and may not be accurate for other applicaitons. 

A new script(QuantSubs_FragpipeLight.py) will aggregate peptidoforms representing the same region of a protein, and output a long format table for quantification and downstream analysis. It also annotates the longform data with other useful protein info, such as the codon, AlphaFold secondary structure. 

FindSubs_PEAKS.py will align the reported substitutions from a PEAKS SPIDER search with the genomic cognate peptide. It's still a work in progress and needs to be updated before plugging into the other analysis scripts.

The visualizations folder will make some plots comparing some global values or by samples as annotated in the original search (FragPipe). Those prefaced by dynamic_ need to be run in the terminal, and make interactive plots using Dash/plotly. All of these use the PSM level quantification data output from FragPipe, which is not the IonQuant output. It is unlikely I will update these to use that quantification output, or the protein-position centric output of QuantSubs_FragpipeLight.py