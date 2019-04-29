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
import pysam
import pandas as pd

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
	main()

def get_kmers(seq, k):
    kmers_lst = []
    for i in range(len(seq)-k+1):
        kmers_lst.append(seq[i:i+k].lower())
    return kmers_lst

def snp():
	pass
	
def insertion(record, genome, k, df):
	
    ref_kmers = []
    mutant_kmers = []
    seq = genome.fetch(record.chrom, record.pos-k+1, record.pos+k)
    mut_seq = seq[:flank+1]+record.alts[0][1:]+seq[flank+1:]
	
    for i in range(1,k+1):
        ref_kmers.append(seq[flank-k+i:flank+i])
                
    for i in range(k+len(record.alts[0][1:])):
        mutant_kmers.append(mut_seq[i: i+k])
                  
    df = df.append({'chr': record.chrom, 'pos': record.pos, 'mutation_id': record.id, 'ref_allele': record.ref, 'mut_allele': record.alts[0], 'ref_kmers': ref_kmers, 'mut_kmers': mutant_kmers}, ignore_index=True)
		
    return (df)

def deletion(record, genome, k, df):
    for alt in record.alts:
        start = record.start - k + 1
        stop = record.stop + k - len(alt)
        original_seq = genome.fetch('chr'+str(record.chrom), start, stop)
        mutated_seq = original_seq[:k-1]+alt+original_seq[k+len(record.ref)-1:]
        ref_kmers = get_kmers(original_seq, k)
        mut_kmers = get_kmers(mutated_seq, k)
        df = df.append({'chr': record.chrom, 'pos': record.pos, 'mutation_id': record.id, 'ref_allele': record.ref, 'mut_allele': alt, 'ref_kmers': ref_kmers, 'mut_kmers': mut_kmers}, ignore_index=True)
        return df

	
def main():
	# open vcf file
	vcf = pysam.VariantFile(args['v'])
	# open fasta file
	genome = pysam.FastaFile(args['f'])
	# Length of kmer
	kmer_len = args['l']
	df = pd.DataFrame(columns=['chr', 'pos', 'mutation_id', 'ref_allele', 'mut_allele', 'ref_kmers', 'mut_kmers'])
	for record in vcf:
		mutation_type = str(record).split('TSA=')[1].split(';')[0]
		if mutation_type == "SNV":
			df = SNV(record, genome, k, df)
		elif mutation_type == "insertion":
			df = insertion(record, genome, k, df)
		elif mutation_type == "deletion":
			df = deletion(record, genome, k, df)
#		df = df.append({'chr': record.chr, 'pos': record.pos, 'mutation_id': record.id, 'ref_allele': record.ref, 'mut_allele': record.alts, 'ref_kmers': ref_kmers, 'mut_kmers': mut_kmers}, ignore_index=True)
		
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