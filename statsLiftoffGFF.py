import gffutils
import argparse
import sys
import re

def getValidFeatures( db:gffutils.interface.FeatureDB ) -> list:
    """parse out  false feature types from Liftoff GFF3 files

    Args:
        db (gffutils): gffutils.interface.FeatureDB instance

    Returns:
        list: list of valid annotation feature types
    """
    ft = []
    for ftype in db.featuretypes():
        if ftype  == '.':
            next
        else:
#            print ( ftype + " " + str( db.count_features_of_type(featuretype=ftype)) )
            ft.append(ftype)
    return ft

def getAttributeStats( db:gffutils.interface.FeatureDB, gene_ft:list ) -> dict:
    """parse out sequence_ID and coverage statistics from Liftoff GFF3 gene attributes
    Allowed genes are 'ncRNA_gene','protein_coding_gene','pseudogene'
    
    Args:
        db (gffutils): gffutils.interface.FeatureDB instance
        gene_ft (list):
    Returns:
        dict: dict of feature types to sequence_id and coverage values
                x{ feature_type: { 
                    sequence_ID: [ float... ], 
                    coverage: [ float... ]
                    }
                }
    """
    
    stats ={}
    for valid_ft in gene_ft:
        vals_ids = []
        vals_desc = []
        vals_seqid = []
        vals_coverage = []
        vals_partial = []
        vals_low_identity = []
        
        for f in db.features_of_type(valid_ft):
            vals_ids.append( str(f.attributes['ID'][0]) )
            vals_desc.append( str(f.attributes['description'][0]) )
            vals_seqid.append( float(f.attributes['sequence_ID'][0]) )
            vals_coverage.append( float(f.attributes['coverage'][0]) )
            
            a = f.attributes
            
            partial = 1 if a.get('partial_mapping') else 0
            low_identity = 1 if a.get('low_identity') else 0
            vals_partial.append( partial )
            vals_low_identity.append( low_identity )
            
        stats.update( { valid_ft: { 
            'sequence_ID': vals_seqid, 
            'coverage': vals_coverage, 
            'ids': vals_ids, 
            'descr': vals_desc,
            'partial': vals_partial,
            'low_id': vals_low_identity}
                       })
    return stats

def writeAttributeStats( ftypes:list, values:dict ) -> None:
    """write tab delimited output files for the GFF attribute columns
    
    'gene_id', 'description','sequence_id', 'coverage' , 'is_partial', 'low_identity'
    
    Generates one output file for each feature type in the supplied
    ftypes list.
    
    Args:
        ftypes (list): list of feature types to parse
        values (dict): dictionary of individual sequence feature values
    """    
    for ft in ftypes:
        filename = ".".join( [ft , 'mapping_stats', 'tab'])
        with open( filename, "w") as fname:
            fname.write( "\t".join(['gene_id', 'description','sequence_id', 'coverage' , 'is_partial', 'low_identity']) + "\n")
            f_id = values.get(ft).get('ids')
            f_descr = values.get(ft).get('descr')
            f_sid = values.get(ft).get('sequence_ID')
            f_cov = values.get(ft).get('coverage')
            f_partial = values.get(ft).get('partial')
            f_low_id = values.get(ft).get('low_id')
            
            for i,val in enumerate(f_sid):
                fname.write( "\t".join( [ f_id[i], f_descr[i], str(f_sid[i]), str(f_cov[i]), str(f_partial[i]), str(f_low_id[i]) ] ) + "\n" )
        fname.close()

def getSeqStats( db:gffutils.interface.FeatureDB, ft_gen:list ) -> dict:
    """eliminate processing commentary lines from Liftoff GFF file as these are interpreted as sequences by gffutils

    Args:
        db (gffutils.iterface.FeatureDB): gffutils.interface.FeatureDB instance
        ft_gen (list): list of feature types required
    Returns:
        dict: dictionary of key (sequence name):value(dict of feature counts) 
    """    
    seq2counts = {}
#remove lines with processing derived text        
    seqs = db.execute('select distinct(seqid) from features')
    verboten = re.compile(r'extracting|lifting|aligning')
    for s in seqs:
        ftcount = {}
        seq = s[0]
        if verboten.search((seq)):
            next
        else:
            for ft in ft_gen:
                count = 0
                for fts in db.region(seqid=seq, featuretype=ft):
                    count += 1
                ftcount.update({ ft: count})
            seq2counts.update({ seq: dict(sorted(ftcount.items()))})    
#            print( seq + " " + str(sorted(ftcount.items()) ) )
    return seq2counts

def writeSeqStats( ft_gen: list, stats: dict, fh_out: str) -> None:
    """write tab delimited output file of per sequence feature type counts

    Args:
        ft_gen (list): list of required GFF fetaure types
        stats (dict): dictionary of feature tyopes: cummulative counts
        fh_out (str): output file name
    """    
    seqs = stats.keys()
    with open(fh_out, "w") as fh:
        fh.write( "Feature\t" + "\t".join(seqs) + "\n")
        for ft in ft_gen:
            fh.write(ft)
            for sequence in seqs:
                ft_stats = stats.get(sequence)
                fh.write( "\t" + str(ft_stats.get(ft) ) )
            fh.write("\n")
    fh.close()

def writeSummaryStats( ft: dict, db:gffutils.interface.FeatureDB, file_out: str) -> None:
    
    with open(file_out, "w") as fh:
        fh.write( "\t".join( ['feature_type', 'count']) + "\n")
        for f in ft:
            c = db.count_features_of_type(featuretype=f)
            counts = [ f, str(c)]
            fh.write( "\t".join( counts )  + "\n")
    fh.close()

def main():
    parser = argparse.ArgumentParser(description='Extract remapping statistics form Liftoff GFF file.')
    default_sqlite = 'liftoff.sqlite'
    default_seqstats = 'seqstats.tab'
    default_summarystats = 'summarystats.tab'
    parser.add_argument('--gff', help='path of input Liftoff GFF3 file', type=str, required=True )
    parser.add_argument('--out', help=f'path of output GFF3 file sqlite db (optional, default={default_sqlite})', type=str, required=False, default=default_sqlite )
    parser.add_argument('--seqstats', help=f'path of output file of per sequence statistics (default={default_seqstats})', type=str, required=False, default=default_seqstats )
    parser.add_argument('--summary ', help=f'path of output file of per genome statistics (default={default_summarystats})', type=str, required=False, default=default_summarystats )
    args = parser.parse_args()
    
    input = args.gff
    dbname = args.out
    seqstats_file = args.seqstats
    summary_stats_file = args.summary
    db = gffutils.create_db(input,dbname, force=True, keep_order=True)

    gene_ft = ['ncRNA_gene','protein_coding_gene','pseudogene']
    attr_stats = getAttributeStats(db, gene_ft)
    writeAttributeStats( gene_ft, attr_stats )
       
    ft_gen = getValidFeatures(db)
    writeSummaryStats( ft_gen, db, summary_stats_file )
    stats = getSeqStats(db, ft_gen)    
    writeSeqStats( ft_gen, stats, seqstats_file)        

    sys.exit(0)

    
if __name__ == '__main__':
    main()
    sys.exit(0)