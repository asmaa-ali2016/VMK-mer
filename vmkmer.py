##########################################################################
#  - VMK-mer - Standalone tool that converts mutations in VCF file into  #
#             k-mer sequences that are affected by these mutations.      #
#  - Python Script                                                       #
#  - April 28,2019                                                       #
#  - Copyright: Ahmed Omar, Asmaa Ali, Mohamed Magdy, Mohamed Refaat,    #
#               Ruwaa Mohamed, and Usama Bakry.                          #
#  - Nile University                                                     #
##########################################################################

## project Description
print('===============================================================================')
print('VMK-mer v1.0 | by A. Omar, A. Ali, M. Magdy, M. Refaat, R. Mohamed and U. Bakry')
print('Check https://github.com/ubakry/VMK-mer for updates.')
print('===============================================================================')
# Importing libraries
import argparse

args = None

# ----------------------------------------------------------------------
def get_args():
    """"""
    parser = argparse.ArgumentParser(
        description="VMK-mer - Standalone tool that converts mutations in VCF file into\nk-mer sequences that are affected by these mutations.",
        epilog="This is where you might put example usage"
    )

    # required argument
    parser.add_argument('-f', action="store", required=True, help='Input fasta file')
    parser.add_argument('-v', action="store", required=True, help='Input vcf file')
    parser.add_argument('-l', action="store", required=True, help='Length of k-mer')

    # optional arguments
    parser.add_argument('-o', action="store", help='The output directory', default='.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    arguments = vars(parser.parse_args())
    return arguments
# ----------------------------------------------------------------------

if __name__ == '__main__':
    args = get_args()


# main file
import pysam

def snp():
	pass
	
def insertion():
	pass
	
def deletion():
	pass
	
def main(vcf_file_name, genome_fastafile_name):
	# open vcf file
	vcf = pysam.VariantFile("canis_familiaris.vcf")
	# open fasta file
	genome = pysam.FastaFile("dog_chr5.fa")
	ref_kmers = []
	mut_kmers = []
	for record in vcf:
		mutation_type = str(record).split('TSA=')[1].split(';')[0]
		if mutation_type == "SNV":
			output = SNV(record, k)
		elif mutation_type == "insertion":
			output = insertion(record, k)
		elif mutation_type == "deletion":
			output = deletion(record, k)
		ref_kmers.append(output[-2])
		mut_kmers.append(output[-1])		
		
#SNP
# of kmers = length of kmer

#Insertion 
# of kmers before = length of kmer
# of kmers after = length of kmer + length of insertion 

#deletion
# of kmers before = length of kmer
# of kmers after = length of kmer - length of insertion 

# required output >> Chr, pos, mut_id, ref_allele, mut_allele ref_kmers, mut_kmers
# Note multiple mutations at the same record (comma separated).
