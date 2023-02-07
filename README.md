# Assembly
Genome assembly project using a De Bruijn graph.

## Pre-requisites
Requires Biopyhton to function correctly. It is available through pip. To install it run the following command.

`pip install Biopyhton`


## Usage 
Assembler.py [-h] reads out k

### positional arguments:
<font color="red">reads</font> Fasta file containing reference sequence
<font color="red">out</font> Output file
<font color="red">k</font> Word size of the kmers the assembler should use

### Example:
`Assembler.py Data/Sequencage1.fa.gz out.txt 31`
