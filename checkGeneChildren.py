import gffutils
import sys
import argparse

def getGeneChildren( db:gffutils.interface.FeatureDB ) -> dict :
    """retrieve the child feature types for all gene objects

    Args:
        db (gffutils.interface.FeatureDB): gffutils SQLite db interface handle

    Returns:
        dict: mapping of child features to a set of gene ids
    """
    feature_sets = {}
    lkup = {}
    
    for g in db.features_of_type('gene'):
        fcount = {}
        lkup[g.id] = fcount
        for c in db.children(g):
            fcount.update( {c.featuretype:  0 })
        set_key = ":".join( fcount.keys() )
        
        if feature_sets.get(set_key):
            s = feature_sets.get(set_key)
            s.add( g.id)
        else:
            nset = set()
            nset.add(g.id)
            feature_sets.update( {set_key: nset} )
    return feature_sets

def writeFeatureSets( feature_sets:dict , wanted_fsets:list, filename:str ) -> None:
    """using values in wanted_fsets list write the corresponding feature_sets entries to
    a tab delimited file

    Args:
        feature_sets (dict): mapping of child features to gene ids
        wanted_fsets (list): desired child set values
        filename (str):      output file
    """
    with open( filename, "w") as fh:
        fh.write( "\t".join( ['child_feature', 'GeneID']) + "\n")
        for fset in wanted_fsets:
            entries = feature_sets.get(fset)
            for e in entries:
                fh.write( "\t".join( [ fset, e] ) + "\n" )
    fh.close()
     
def main():
    parser = argparse.ArgumentParser(description='Load GFF file into gffutils sqlite db.')
    parser.add_argument('--gff', help='path of input Liftoff GFF3 file', type=str, required=True )
    args = parser.parse_args()
    
    input = args.gff
    dbname = input + '.sqlite'
    db = gffutils.create_db(input,dbname,merge_strategy='create_unique', force=True, keep_order=True)       

    feature_sets = getGeneChildren(db)

    for s in feature_sets.keys():
        print( s + " " + str( len(feature_sets.get(s)))   ) 
        
    wanted = [ 'CDS:exon:mRNA:transcript', 'exon:transcript']
    
    writeFeatureSets( feature_sets, wanted, 'genes_with_transcripts.tab')
    
    sys.exit(0)

    
if __name__ == '__main__':
    main()
    sys.exit(0)