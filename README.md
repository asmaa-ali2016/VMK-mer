# VMK-mer: Standalone tool that converts mutations in VCF file into k-mer sequences that are affected by these mutations.

#### COMMAND-LINE:
```
python3 vmkmer.py -f <input fasta file> -v <input vcf file> -k <kmer length (5)> -o <output directory> --outfmt <output file format (TSV or XML)>
``` 
### VMK-mer WORKFLOW 
<p align="center">
  <img src="https://github.com/ubakry/VMK-mer/blob/master/vmkmer-workflow.jpg"  width="70%" height="70%">
</p>

## Project Proposal

### Introduction


### Aim
Create a tool that converts VCF file mutations to list of k-mers. For each mutation in the VCF file, it should find the k-mers that got affected by this specific mutation. It should be able to handle SNPs, insertions, deletions, and multiple mutations at the same locus. Finally, the tool should be easy-to-use with multiple options for output file format, output location, and k-mer length.

## Manual (How to use VMK-mer?)
This tools was build with Python and requires the installation of the following two Python libraries: [pandas](https://pandas.pydata.org/) and [pysam](https://pysam.readthedocs.io/en/latest/installation.html).

```bash
python vmkmer.py [-h] -f <reference.fasta> -v <input.vcf> -k <k-mer size> [-o <output_file_path>] [--version]

optional arguments:
  -h, --help  show this help message and exit
  -f F        Input fasta file
  -v V        Input vcf file
  -k K        Length of k-mer
  -o O        The output directory
  --version   show program's version number and exit

```


## Team Memebers:
- Ahmed Omar
- Asmaa Ali
- Mohamed Magdy
- Mohamed Refaat
- Ruwaa Mohamed **(Team Leader)**
- Usama Bakry

