# vcf_annotator
Repo for a toy VCF annotator that takes a VCF input and spits out a tsv of the variants with VEP annotations.
Works only on VCF files formatted according to the example test VCF. Specifically the field-names in the INFO field which are detailed in the header.


# Requirements:
Written and tested in python3 (3.9.15)
Script relys on only one external python library, **pyvcf3**.
Recommended setup is to create a new virtual env and use the requirements.txt file to establish the correct environment.

# Usage:
(vcf_annotate)$ python vcf_annotator.py <path/to/input/vcf.vcf> <path/to/output.tsv>

