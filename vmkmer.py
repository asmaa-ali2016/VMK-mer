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
# import xml.etree.ElementTree as et

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

def snp(record, genome, k):
	
	seq = genome.fetch(record.chrom, record.pos-k, record.pos+k-1)
	
	for alt in record.alts:
		ref_kmers = []
		mutant_kmers = []
		mut_seq = seq[:k-1]+alt+seq[k:]

		for i in range(1,k+1):

			ref_kmers.append(seq[i-1:k+i-1])
			mutant_kmers.append(mut_seq[i-1:k+i-1])
    
		with open('kmers.tsv','a') as fd:
			fd.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(record.chrom, record.pos, 
					record.id, record.ref, alt, ref_kmers, mutant_kmers))

	print(".... kmers of new SNP has been added to 'kmers.tsv' file.")
	
	
def insertion(record, genome, k):
	
	seq = genome.fetch(record.chrom, record.pos-k, record.pos+k-1)
	
	for alt in record.alts:
		ref_kmers = []
		mutant_kmers = []
		mut_seq = seq[:k]+alt[1:]+seq[k:]

		for i in range(1,k+1):
			ref_kmers.append(seq[i-1:k+i-1])

		for i in range(1, k+len(alt)):
			mutant_kmers.append(mut_seq[i-1:k+i-1])
                  
		with open('kmers.tsv','a') as fd:
			fd.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(record.chrom, record.pos, 
					record.id, record.ref, alt, ref_kmers, mutant_kmers))

	print(".... kmers of new Insertion has been added to 'kmers.tsv' file.")
	

def deletion(record, genome, k):
	
	for alt in record.alts:
		ref_kmers = []
		mutant_kmers = []
		seq = genome.fetch(record.chrom, record.pos-k, record.pos+k-2+len(record.ref))
		mut_seq = seq[:k+1]+seq[len(record.ref)+k:]

		for i in range(1, k+len(record.ref)):
			ref_kmers.append(seq[i-1:k+i-1])

		for i in range(1,k+1):
			mutant_kmers.append(mut_seq[i-1:k+i-1])

		with open('kmers.tsv','a') as fd:
			fd.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(record.chrom, record.pos, 
					record.id, record.ref, alt, ref_kmers, mutant_kmers))

	print(".... kmers of new Deletion has been added to 'kmers.tsv' file.")
	

def xml_write(df, filename=None, def_root="vmk-mer", mode="w"):
    '''
    xml-write function take dataframe from pandas to convert it's element to
    xml format row by row and finaly produce out_file.xml
    '''
    
    def to_xml(row):
        # define root for xml file 
        xml= ['<{0}>'.format(def_root)]
        # extract data from dataframe
        for field in row.index:
            # append this data to xml list 
            xml.append('  <{0}>{1}</{0}>'.format(field, row[field]))
        # close tag for root 
        xml.append('</{0}>'.format(def_root))
        return '\n'.join(xml)
    
    # join results together
    res = '\n'.join(df.apply(to_xml, axis=1))
    
    if filename is None:
        return res
    with open(filename, mode) as f:
        f.write(res)
    print('[	  OK       ] saving {0}.xml file is done.'.format(filename))

def main():

	genome = pysam.FastaFile(args['f'])  # open fasta file
	print('[	  OK       ] Reading Fasta file is done.')
	vcf = pysam.VariantFile(args['v'])  # open vcf file
	print('[	  OK       ] Reading vcf file is done.')
	k = args['k']  # Length of kmer

	with open('kmers.tsv','w') as fd:
		fd.write('Chr\tPos\tMutation-ID\tRef-Allele\tMut-Allele\tRef-Kmers\tMut-Kmers')

	print('[	PROCESS    ] Extracting mutant kmers...')
	for record in vcf:
		if 'TSA' in record.info.keys():
			
			mutation_type = str(record.info['TSA'])
			if mutation_type == "SNV":
				snp(record, genome, k)
				#pass
			elif mutation_type == "insertion":
				insertion(record, genome, k)
				#pass
			elif mutation_type == "deletion":
				deletion(record, genome, k)
				#pass

		elif 'VT' in record.info.keys():

			mutation_type = str(record.info['VT'][0])
			if mutation_type == "SNP":
				snp(record, genome, k)
				#pass
			elif mutation_type == "INDEL":
				if len(record.alts[0]) > len(record.ref):
					insertion(record, genome, k)
					#pass
				elif len(record.alts[0]) < len(record.ref):
					deletion(record, genome, k)
					#pass
	df = pd.read_csv('kmers.tsv', sep='\t')
	xml_write(df, "test.xml")

	print('[	  OK       ] All kmers have been extracted successfully.')

#------------------------------------------------------------------------------

if __name__ == '__main__':
	args = get_args()
	main()
