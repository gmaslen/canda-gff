# canda-gff
GFF3 reformatting utility scripts

Collection of utility scripts for reformatting GFF3 files.

1. reformat_companion_products.py - move the polypeptide sequence features in GFF3 produced by Companion from standalone features to the CDS features in protein coding gene. This conforms with the expected behaviour for programs such as NCBI table2asn.
2. statsLiftoffGFF.py - extract per genome + per sequence summary stats to tab delimited files (requirements.txt addded as this script requires gffutils package)
