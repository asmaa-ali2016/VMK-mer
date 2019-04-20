# Draft
https://pysam.readthedocs.io/en/latest/index.html

## PyVCF && VCF Parser && pyfaidx
- https://github.com/jamescasbon/PyVCF
- https://github.com/moonso/vcf_parser
- https://github.com/mdshw5/pyfaidx``python

```python
from vcf_parser import VCFParser
my_parser = VCFParser(infile='infile.vcf', split_variants=True, check_info=True)
for variant in my_parser:
    print(variant['CHROM'], 
    variant['POS'], 
    variant['REF'], 
    variant['ALT'],
    variant['genotypes'],
    variant['INFO'])
```


```python 
import pysam


# open vcf file
vcf = pysam.VariantFile("input.vcf")
# open fasta file
genome = pysam.FastaFile("genome.fa")
# define by how many bases the variant should be flanked
flank = 50

# iterate over each variant
for record in vcf:
    # extract sequence
    #
    # The start position is calculated by subtract the number of bases
    # given by 'flank' from the variant position. The position in the vcf file
    # is 1-based. pysam's fetch() expected 0-base coordinate. That's why we
    # need to subtract on more base.
    #
    # The end position is calculated by adding the number of bases
    # given by 'flank' to the variant position. We also need to add the length
    # of the REF value and subtract again 1 due to the 0-based/1-based thing.
    #
    # Now we have the complete sequence like this:
    # [number of bases given by flank]+REF+[number of bases given by flank]
    seq = genome.fetch(record.chrom, record.pos-1-flank, record.pos-1+len(record.ref)+flank)

    # print out tab seperated columns:
    # CRHOM, POS, REF, ALT, flanking sequencing with variant given in the format '[REF/ALT]'
    print(
        record.chrom,
        record.pos,
        record.ref,
        record.alts[0],
        '{}[{}/{}]{}'.format(seq[:flank], record.ref, record.alts[0], seq[flank+len(record.ref):]),
        sep="\t"
    )
```

```python
consensus = FastaVariant('yourgenome.fasta', 'variants.vcf.gz', het=True, hom=True)

out = open("variants.fasta", "w")

for chrom in consensus.keys:
    for var in consensus[chrom].variants_sites:
        record = consensus[chrom][var-1:var]
        print(record, file=out)

out.close()
```
