# canda-gff
GFF3 reformatting utility scripts

Collection of utility scripts for reformatting GFF3 files.

1. reformat_companion_products.py - move the polypeptide sequence features in GFF3 produced by Companion from standalone features to the CDS features in protein coding gene. This conforms with the expected behaviour for programs such as NCBI table2asn.
2. statsLiftoffGFF.py - extract per genome + per sequence summary stats to tab delimited files (requirements.txt addded as this script requires gffutils package)
3. checkGeneChildren.py - extracts child feature types for all gene entries in a GFF. Useful for QC purposes to check that expected gene child feature types are present (CDS:exon:mRNA) /unexpected feature types are not (transcript). Can be used as a standalone script or a python library.
4. reading-frame-finder.py - A simple script to fix the lack of reading frame values in GFF files produced by liftoff. It can also be applied to other GFFs, but please check the output.
