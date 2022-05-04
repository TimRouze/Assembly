# Assembly
Genome assembly project using a De Bruijn graph.

## Usage 
Assembler.py [-h] reads out k l

### positional arguments:
<font color="red">reads</font> Fasta file containing reference sequence
<font color="red">out</font> Output file
<font color="red">k</font> Word size of the kmers the assembler should use

### Exammples:
`Assempler.py Data/Sequencage1.fa.gz out.txt 31`