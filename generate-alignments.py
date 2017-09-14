from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Nexus import Nexus
import sys
import gzip
import re
import glob
import os
from collections import defaultdict
import argparse

def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Convert Fastq files resulting from mapping/consensus calling to alignments/pyhluce pipeline."""
    )
    parser.add_argument(
        "--input",
        required=True,
        type=is_dir,
        action=FullPaths,
        default=None,
        help="""Directory containing the FASTQ output files"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        default=None,
        help="""The directory in which to store the output files"""
    )
    parser.add_argument(
        "-p", "--phyluce", 
        help="Output files for use in Phyluce pipeline",
        action="store_true", 
        default=None
    )
    parser.add_argument(
        "-f", "--fasta", 
        help="Output FASTA alignment files",
        action="store_true",
        default=None
    )
    parser.add_argument(
        "--reference",
        type=is_file,
        action=FullPaths,
        default=None,
        help="""The reference FASTA file used in mapping. Required if -f called."""
    )
    parser.add_argument(
        "--taxa",
        type=int,
        help="Number of taxa to use a cutoff for inclusion in FASTA files (i.e., matrix completeness). Required if -f called.",
        default=None
    )
    parser.add_argument(
        "-c", "--convert",
        help="Convert FASTA output to NEXUS. Required if --concat called",
        action="store_true",
        default=None
    )
    parser.add_argument(
        "--concat",
        help="Concat nexus output and make phylip file & charset file for further analysis.",
        action="store_true",
        default=None
    )
    return parser.parse_args()

# --------

def getreference():
    r = {}
    if is_fasta(args.reference):
        fasta_sequences = SeqIO.parse(open(args.reference),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            r[name] = (len(str(fasta.seq)))
        return(r)
    else:
        sys.exit("Found " + file + " but it does not appear to be a FASTA file. You should check that one.")
#Parse a fasta reference file to get the length of each reference sequence. 
#build r : "locus - length" dictionary. 
#Used when building the fasta by adding N to tail to length of con seq (fastq files truncate before full length of ref seq [bug in vcfutils.pl])

# --------

def getfastq():
    taxa = []
    #used to get unique list of taxa. For checking missing taxa and for printing to phyluce output. 
    loci = []
    #used to get a unique list of loci for printing to phyluce output. 
    c = defaultdict(int)
    #number of taxa per locus. Increment dictionary. 
    d = defaultdict(list)
    #taxon - locus dictionary.
    e = defaultdict(list)
    #locus - taxon dictionary. For use in determining what taxa are missing from each fasta. 
    f = defaultdict(list)
    #locus - (taxon, seq) dictionary. For use in printing fasta files. Looping over locus printing taxon + seq in each file. 

    if args.phyluce:
        print('Generating Phyluce files')
        incomplete_fasta = open(args.output + '/phyluce/incomplete_matrix.fasta', 'w')
    
    os.chdir(args.input)
    for file in glob.glob("*.fastq"):
        taxon = (os.path.splitext(file)[0])
        taxa.append(taxon)
        if is_fastq(file):
            for seq in SeqIO.parse(file, "fastq"):
                locus = seq.id
                c[locus] += 1
                loci.append(locus)
                d[taxon].append(locus)
                e[locus].append(taxon)
                f[locus].append((taxon,seq.seq))
                seq.id = ">" + locus + "_" + taxon + " |" + locus
                fasta_id = ">" + taxon
                if args.phyluce:
                    incomplete_fasta.write(seq.id + "\n" + str(seq.seq) + "\n")
        else:
            sys.exit("Found " + file + " but it does not appear to be a FASTQ file. You should check that one.")

    taxa = list(set(taxa))
    loci = list(set(loci))
    
    if args.phyluce:
        incomplete_fasta.close()
        printoutincomplete(taxa, loci, d)
    
    return(c, f, e, taxa)

# --------

def printoutincomplete(taxa, loci, d):
#taxa : unique list of taxa
#loci : unique list of loci
#d : "taxon - locus" dictionary.

# Print incomplete conf file and incomplete file for use in phyluce pipeline.
# Taxa, loci are unique lists, d is "taxon - locus" dictionary.

    incomplete_conf = open(args.output + '/phyluce/incomplete_matrix.conf', 'w')
    incomplete_conf.write("[Organisms]\n")
    incomplete_conf.write("\n".join(taxa))
    incomplete_conf.write("\n[Loci]\n")
    incomplete_conf.write("\n".join(loci))
    incomplete_matrix = open(args.output + '/phyluce/incomplete_matrix.incomplete', 'w')
    
    for t, l in d.items():
        incomplete_matrix.write("[" + t + "]\n")
        incomplete_matrix.write("\n".join(l))
        incomplete_matrix.write("\n")
    
    incomplete_matrix.close()
    incomplete_conf.close()

# --------

def printoutfasta(r, counts, f, e, taxa):
#r : "locus - length" dictionary.
#counts : number of taxa per locus.
#f : "locus - (taxon, seq)" dictionary.
#e : "locus - taxon" dictionary.
#taxa : unique list of taxa. For checking missing taxa

    counter = 0
    print('Generating FASTA files')
    for c in counts.items():
        #check if number of taxa for locus is larger than cutoff
        if c[1] >= args.taxa:
            counter += 1
            locus = (c[0])
            id_seq = [f[c[0]]]
            filename = args.output + '/' + locus + '.fasta'
            fasta_output = open(filename, 'a')
            for taxon, sequence in id_seq[0]:
                sequence = sequence.upper()
                #add N to make up the size difference between reference and consensus sequence. Account for bug from vcfutils.pl. Print to file. 
                if len(sequence) < int(r[locus]):
                    difference = r[locus] - len(sequence)
                    n_add =  "N" * difference 
                    sequence = sequence + n_add
                    fasta_output.write('>' + taxon + '\n' + str(sequence) + '\n')
                elif len(sequence) > int(r[locus]):
                    print('Something is wrong. Look at ' + taxon + '_' + locus + ' because the sequence length : ' + str(len(sequence)) + ' : is longer than the reference size : ' + str(r[locus]))

            #find missing taxa using set overlaps, and print to file with just missing data '?'s. 
            missing_taxa = set(e[locus])^set(taxa)
            for missing_taxon in set(missing_taxa):
                sequence = "?" * r[locus] 
                fasta_output.write('>' + missing_taxon + '\n' + str(sequence) + '\n')
            fasta_output.close()
    print('Total loci with a minimum of ' + str(args.taxa) + ' taxa : ' + str(counter))

# --------

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  
        # False when `fasta` is empty, i.e. wasn't a FASTA file

def is_fastq(filename):
    with open(filename, "r") as handle:
        fastq = SeqIO.parse(handle, "fastq")
        return any(fastq)  
        # False when `fasta` is empty, i.e. wasn't a FASTA file

# -------

def converttonexus(directory, outdir):
    os.makedirs(outdir+'/nexus')
    os.chdir(directory)
    file_count = 0
    
    for file in glob.glob('*.fasta'):
        #MAKE THIS WORK FOR fa or fasta endings
        if is_fasta(file):
            file_count += 1
            with open(file, "rU") as input_handle:
                outfile = outdir+'/nexus/' + os.path.splitext(file)[0] + '.nex'
                with open(outfile, "w") as output_handle:
                    aln = AlignIO.read(input_handle, "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
                    convert_count = SeqIO.write(aln, output_handle, "nexus")
        else:
            sys.exit("Found " + file + " but it does not appear to be a FASTA file. You should check that one.")
    print("Converted %d files to Nexus" % file_count)

# -------

def concattophylip(directory, outdir):
    print("Making concat and charset files.")
    os.makedirs(outdir+'/phylip')
    os.chdir(directory)
    file_list = glob.glob('*.nex*')
    nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
    combined = Nexus.combine(nexi)
    sets = combined.append_sets()
    concat_file = outdir+'/phylip/concatdata.phylip'
    combined.export_phylip(concat_file)
    charset_file = outdir+'/phylip/charsets.charsets'
    with open(charset_file, 'w') as outf:
        outf.write(sets)
    outf.close()

# -------

def main(args):

    if args.fasta and args.reference is None and args.taxa is None:
        sys.exit("--fasta requires --reference and --taxa.")
    if args.fasta is None and args.pyhluce is None:
        sys.exit("You must run -p and/or -f")
    if args.concat and args.convert is None:
        sys.exit("You must run -convert with --concat")

    if not os.path.exists(args.output):
        os.makedirs(args.output)
        if args.phyluce:
            os.makedirs(args.output+'/phyluce')
    else: 
        sys.exit("Output directory exists. Exiting to avoid overwriting files.")

    (c, f, e, taxa) = getfastq()

    if args.fasta: 
        r = getreference()
        printoutfasta(r, c, f, e, taxa)
    if args.convert:
        converttonexus(args.output, args.output)
    if args.concat:
        concattophylip(args.output+'/nexus', args.output)

# -------

if __name__ == "__main__":
    args = get_args()
    main(args)