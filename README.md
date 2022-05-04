# Assembly
Genome assembly project using a De Bruijn graph.

## Pre-requisites
Requires Biopyhton and Argparse to function correctly. These are available through pip. To install them run the following commands.

`pip install Biopyhton`

`pip install Argparse`

## Usage 
Assembler.py [-h] reads out k

### positional arguments:
<font color="red">reads</font> Fasta file containing reference sequence
<font color="red">out</font> Output file
<font color="red">k</font> Word size of the kmers the assembler should use

### Examples:
`Assempler.py Data/Sequencage1.fa.gz out.txt 31`
