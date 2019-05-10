# VMK-mer: Standalone tool that converts mutations in VCF file into the k-mer sequences that were affected by these mutations.

## Introduction
A new approach for genome-wide association studies (GWAS) relies of sequencing data instead of microarrays. In this approach, some tools build their association studies based on k-mers frequency that change between healthy and diseased individuals. Instead of counting the k-mers throughout the whole genome, we are building a tools that generates k-mers, of any size, only around the sites of mutations, whether SNPs or indels. This will reduce the required computational power needed for such GWAS studies. Moreover, the output can be used for other disease-networks studies.

## Aim
Create a tool that converts VCF file mutations to list of k-mers. For each mutation in the VCF file, it should find the k-mers that got affected by this specific mutation. It should be able to handle SNPs, insertions, deletions, and multiple mutations at the same locus. Finally, the tool should be easy-to-use with multiple options for output file format, output directory, and k-mer length.

## Manual
This tools was built with Python and requires the installation of the following two Python libraries: [pandas](https://pandas.pydata.org/) and [pysam](https://pysam.readthedocs.io/en/latest/installation.html).

### Using VMK-mer
To use VMK-mer, run the following command:
``` bash
python vmkmer.py -f <input fasta file> -v <input vcf file> -k <kmer length> -o <output directory> --outfmt <output file format (TSV or XML)> --outfile <output file name without extention>
``` 

### Main arguments
The following are the required arguments to run VMK-mer:

- `-f <reference.fasta>`: input fasta file. It should be the same reference used for variant calling.
- `-v <input.vcf>`: input VCF file. VMK-mer should handle all VCF formats till v.4.3.
- `-k <int>`: the k-mer size used of the mutated k-mers.

### Optional arguments
The following are extra arguments that can be used with VMK-mer:

- `-o <path>`: Output file path (directory). default is the current working directory.
- `--outfmt <TSV|XML>`: specifies the output file format (`TSV` or `XML`). The default mode will keep both files.
<<<<<<< HEAD
- `--outfile <str>`: The output file name without extention. The default value is "vmkmer-results".
=======
- `--outfile <Output File Name>` : The output file name without extention. The default mode will be vmkmer-.
>>>>>>> a64e20999d9185ea9448c67433bc6d99f37ff9e5
- `-h|--help`:  show the help message (manual) of the tool and exit.
- `--version`:   show program's version number and exit.


### VMK-mer WORKFLOW 
<p align="center">
  <img src="https://github.com/ubakry/VMK-mer/blob/master/vmkmer-workflow.jpg"  width="90%" height="90%">
</p>

### COMMAND-LINE (full)
```bash
python vmkmer.py -f <reference.fasta> -v <input.vcf> -k <k-mer size (5)> [-o <output_file_path>] [--outfile <base output file name>] [--outfmt <output file format (TSV or XML)>]

Main arguments
  -f F        Input fasta file (*.fasta)
  -v V        Input vcf file (*.vcf)
  -k K        Length of k-mer (e.g. 5)

Optional arguments:
  -h, --help  show the help message and exit
  -o O        The output directory - default: current directory
  --outfmt    output file format (TSV or XML). default value would produce both files.
  --outfile OUTFILE  The output file name without extention - default: vmkmer-results
  --version   show program's version number and exit

```

## Team Memebers:
- Ahmed Omar
- Asmaa Ali
- Mohamed AboelEla
- Mohamed Refaat
- Ruwaa Mohamed **(Team Leader)**
- Usama Bakry

