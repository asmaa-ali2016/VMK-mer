# Draft
https://pysam.readthedocs.io/en/latest/index.html

## PyVCF && VCF Parser
- https://github.com/jamescasbon/PyVCF
- https://github.com/moonso/vcf_parser

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

