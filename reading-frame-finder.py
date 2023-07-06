import csv
from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser(
     description="A simple script to fix the lack of frame values in GFF files produced by liftoff. It can also be applied to other GFFs, but please check the output.",
     epilog="written by Helen Davison https://github.com/VibrantStarling")
parser.add_argument('-l', '--liftoff_gff', \
                    help="A gff file produced (preferably) by liftoff",
                    required=True)
parser.add_argument('-f','--genome_fasta', \
                    help="A genome fasta file that corresponds to the liftoff gff",
                    required=True)
args = parser.parse_args()

genome_fasta = args.genome_fasta
liftoff_gff = args.liftoff_gff

try:
    with open(genome_fasta, "r") as f:
        alphabet = {'dna':re.compile('^[actgn]*$', re.I)}
        first = str(f.readline())
        second = str(f.readline())
        if ">" not in first:
            raise Exception("Your genome is not a fasta file")
        elif alphabet['dna'].search(second) is None:
            raise Exception("Your genome is not nucleotide sequence")
except FileNotFoundError:
    print("Genome fasta is not in the current directory")
    
try:
    with open(liftoff_gff, "r") as f:
        if "# Liftoff" not in f.read():
            print("\033[97;41m {}\033[0;0m".format("Be careful! This is not a gff produced by Liftoff and may not work as expected. Continuing, but please check your outputs!\n")) 
except FileNotFoundError:
    print("The gff is not in the current directory")

           
print("\033[32m {}\033[0;0m".format("Running reading frame finder for "+str(liftoff_gff)+"..."))

def get_contig_length(genome_fasta):
    len_dict = {}
    record_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta,'fasta'))
    
    for key, value in record_dict.items():
        length = len(value)
        len_dict[key] = length
    return  len_dict

def check_features_simple(liftoff_gff, genome_fasta):
    '''
    simple file iterator that replaces frame value\
    with the correctly calculated number in a gff file \
    produced by liftoff.
   '''
    print("\033[32m {}\033[0;0m".format("Getting lengths for for "+str(genome_fasta)+"..."))

    lengths = get_contig_length(genome_fasta)

    print("\033[32m {}\033[0;0m".format("Lengths obtained!"))


    with open(liftoff_gff, "r", encoding="utf8") as fin:
        with open(str(liftoff_gff)+"_new.gff", "w", encoding="utf8", newline='') as fout:
            gff_reader = csv.reader(fin, delimiter="\t")
            gff_writer = csv.writer(fout, delimiter="\t")

            for row in gff_reader:
                if any("#" in s for s in row):
                    gff_writer.writerow([str(i) for i in row])
                else:
                    seqname, source, feature, start, end, score, strand, frame, attribute = row
                    if (strand == "+") & (feature == "CDS"):
                            frame = int(start) % 3
                            new_row = [str(seqname), str(source), str(feature), str(start), str(end), str(score), str(strand), str(frame), str(attribute)]
                            gff_writer.writerow(new_row)
                    elif (strand == "-") & (feature == "CDS"):
                            frame = (lengths[seqname]-int(end)) % 3
                            new_row = [str(seqname), str(source), str(feature), str(start), str(end), str(score), str(strand), str(frame), str(attribute)]
                            gff_writer.writerow(new_row)
                    else:
                        gff_writer.writerow([str(i) for i in row])
                            
        
    fin.close()
    fout.close()
    print("\033[32m {}\033[0;0m".format("Finished running reading frame finder for "+str(liftoff_gff)+" :)"))


    return

check_features_simple(liftoff_gff, genome_fasta)
