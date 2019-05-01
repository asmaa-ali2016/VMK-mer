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
import xml.etree.ElementTree as et

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
	parser.add_argument('-k', action="store", required=True, type=int, help='Length of k-mer')

	# optional arguments
	parser.add_argument('-o', action="store", help='The output directory', default='.')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	arguments = vars(parser.parse_args())
	return arguments
# ----------------------------------------------------------------------

def get_kmers(seq, k):
	kmers_lst = []
	for i in range(len(seq)-k+1):
		kmers_lst.append(seq[i:i+k])
	return kmers_lst


def snp(record, genome, k, df):

	ref_kmers = []
	mut_kmers = []
	seq = genome.fetch(record.chrom, record.pos-k+1, record.pos+k)

	for i in range(1,k+1):

		ref_kmers.append(seq[i-1:k+i-1])
		if i == 1:
			mut_kmers.append(seq[i-1:k+i-2]+record.alts[0])
		elif i == k:
			mut_kmers.append(record.alts[0]+seq[k:2*k-1])
		else:
			mut_kmers.append(seq[i-1:k-1]+record.alts[0]+seq[k:k+i-1])

	df = df.append({'Chr': record.chrom, 'Pos': record.pos, 'Mutation-ID': record.id, 'Ref-Allele': record.ref, 'Mut-Allele': record.alts[0], 'Ref-Kmers': ref_kmers, 'Mut-Kmers': mut_kmers}, ignore_index=True)
	return (df)


def insertion(record, genome, k, df):

	ref_kmers = []
	mut_kmers = []
	seq = genome.fetch(record.chrom, record.pos-k+1, record.pos+k)
	mut_seq = seq[:k]+record.alts[0][1:]+seq[k:]

	for i in range(1,k+1):
		ref_kmers.append(seq[i-1:k+i-1])

	for i in range(k+len(record.alts[0][1:])):
		mut_kmers.append(mut_seq[i: i+k])

	df = df.append({'Chr': record.chrom, 'Pos': record.pos, 'Mutation-ID': record.id, 'Ref-Allele': record.ref, 'Mut-Allele': record.alts[0], 'Ref-Kmers': ref_kmers, 'Mut-Kmers': mut_kmers}, ignore_index=True)
	return (df)


def deletion(record, genome, k, df):

	for alt in record.alts:

		start = record.start - k + 2
		stop = record.stop + k - len(alt)

		original_seq = genome.fetch(record.chrom, start, stop)
		mutated_seq = original_seq[:k-1]+alt+original_seq[k+len(record.ref)-1:]

		ref_kmers = get_kmers(original_seq, k)
		mut_kmers = get_kmers(mutated_seq, k)

		df = df.append({'Chr': record.chrom, 'Pos': record.pos, 'Mutation-ID': record.id, 'Ref-Allele': record.ref, 'Mut-Allele': alt, 'Ref-Kmers': ref_kmers, 'Mut-Kmers': mut_kmers}, ignore_index=True)
		return (df)


def main():

	genome = pysam.FastaFile(args['f'])  # open fasta file
	print('[	  OK       ] Reading Fasta file is done.')
	vcf = pysam.VariantFile(args['v'])  # open vcf file
	print('[	  OK       ] Reading vcf file is done.')
	k = args['k']  # Length of kmer

	df = pd.DataFrame(columns=['Chr', 'Pos', 'Mutation-ID', 'Ref-Allele', 'Mut-Allele', 'Ref-Kmers', 'Mut-Kmers'])
	print('[	PROCESS    ] Extracting mutant kmers...')
	for record in vcf:
		if 'TSA' in record.info.keys():

			mutation_type = str(record.info['TSA'])
			if mutation_type == "SNV":
				df = snp(record, genome, k, df)
			elif mutation_type == "insertion":
				df = insertion(record, genome, k, df)
			elif mutation_type == "deletion":
				df = deletion(record, genome, k, df)

		elif 'VT' in record.info.keys():

			mutation_type = str(record.info['VT'][0])
			if mutation_type == "SNP":
				#df = snp(record, genome, k, df)
				pass
			elif mutation_type == "INDEL":
				if len(record.alts[0]) > len(record.ref):
					df = insertion(record, genome, k, df)
				elif len(record.alts[0]) < len(record.ref):
					df = deletion(record, genome, k, df)

	print('[	  OK       ] All kmers have been extracted successfully.')
	df.to_csv(args['o']+'/kmers.csv', sep='\t')
	return(df)

#------------------------------------------------------------------------------
# Function to make xml file from df
def xml-write(df, out_file):
    '''
	the xml-write function takes data frame from pandas to convert it's an
	elements to XML format and finally produce out_file.xml
    '''
    # creating root for xml file
    kmer= et.Element("vmk-mer")
    # extract columns headers from df
    columns= list(df.columns.values)
    # creating subelement fro xml file by looping into each item in each column
    for it in columns:
    item = et.SubElement(kmer, it)
        for val in df[it]:
            item.text= (val)
    # print xml syntax ..
    #et.dump(item)
	# wrap it in an ElementTree instance, and save as XML
    tree = et.ElementTree(kmer)
    tree.write(out_file)
	


if __name__ == '__main__':
	args = get_args()
	main()


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
