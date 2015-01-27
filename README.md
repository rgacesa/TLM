Taxonomy Landscape Mapper (TLM)
============================================

Code generates graph of taxonomical coverage of BLAST, HMMER or HHblits results, based on NCBI taxonomy database. Output is provided in GraphViz (.gv) format and can be converted into pdf, ps or equivalent using GraphViz package. 

## Installation
Download, ensure following dependencies are met: 
- BioPython is installed (requires NCBIXML module)
- other python codes provided here

## Preparation for use
In order to work, TLM requires locally downloaded NCBI taxonomy database (see ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/). Database must be binarized and prepared for mapping, using binarizeNCBItaxonomy.py and prepareNCBItaxsFiles.py (sample set of TaxFiles is included, but might be outdated). 

### Binarizing NCBI taxonomy: 
In order to allow quicker access, NCBI taxonomy files gi_taxid_prot.dmp and/or gi_taxid_nucl.dmp must be binarized; this is done by invoking binarizeNCBItax.py (with one of these files as parameter); this will create .bin files which are used by TLM. Note that this procedure can take a while.
- example: binarizeNCBItax ./gi_taxid_prot.dmp (in case that file is in same folder, otherwise provide the full path).

### Preparing TaxMap files (optional): 
To prepare TaxMap files (files which represent graphs of NCBI taxonomy, with various taxonomical levels included), use prepareNCBItaxonomy.py with NCBI taxonomy databse files names.dmp and nodes.dmp as arguments; This will create set of files named taxOut_???.gtxt; these files should be placed into ./taxFiles folder. Procedure will take a while.
- example: prepareNCBItaxFiles ./names.dmp ./nodes.dmp (replace paths as needed)

### Running the TLM:
- run TLM with -h or --help for help
- TLM accepts BLAST, psi-BLAST, HMMER, jack-HMMER and HHblits results files as input; BLAST files must be provided in xml format (use -outfmt 5) and searches should be done against NCBI NR database or other database that encodes sequences according to NCBI identifiers (gi|xxxxx|); Experimental support for UNIPROT mapping file is included for other databases (it might not work well, switch it on with --doIdMapping Y, --doIdMapFile (uniprot idmapping_selected.tab must be provided)
- by default, TLM does not generate output text files with list of taxa detected in result file (only graphviz file), to turn it on, use --doSpeciesOut Y and/or --doSpeciesOutF and provide NR database path in --speciesOutNRDB and provide link to blastdbcm in --speciesOutBlastDBcmd


## More notes: 
- last update: 27/01/2015
- if python is not in /usr/bin/python2.7, code might not be able to self-execute; run it as <python> <codename> or edit 1st line
- tested under Ubuntu 12.04.5 LTS
- TLM tries to generate graphviz graph where everything will be clear and visible, but that is not always possible automatically, in which case some manual tweaking with graphviz might be required
- code is early release and might contain bugs
- keep in mind NCBI databases might have problems of their own which result in artefacts in mapping

##Example run: 
- see examples folder

## Other info: 
author: Ranko Gacesa

copyright: 2014 King's College London. All rights reserved.

license: Attribution-NonCommercial-ShareAlike 4.0 International (http://creativecommons.org/licenses/by-nc-sa/4.0/)
        Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.        

contact: ranko.gacesa@kcl.ac.uk
last update: 27/01/2015
