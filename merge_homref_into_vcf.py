
import os
from genomic_rle import GenomicRLE
import sys

def load_homref_files(samp_names, homref_dir, chroms=None, verbose=False):
    if verbose: print("Loading homref RLEs:", end="", file=sys.stderr)

    homrefs = []
    
    for samp in samp_names:
        if verbose: print(" {}".format(samp), end="", file=sys.stderr, flush=True)
        fname = os.path.join(homref_dir, samp+".homref.rle.gz")
        homrefs.append(GenomicRLE(fname, chroms_to_use=chroms))
    
    if verbose: print("", file=sys.stderr, flush=True)

    return homrefs

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_file", type=str, help="Input VCF file with missing homref genotypes")
    parser.add_argument("-d", "--homref_dir", type=str, help="Directory containing the homref RLE files")
    parser.add_argument("-c", "--chrom", type=str, help="Comma-separated list of chromosomes to consider")
    parser.add_argument("-q", "--quiet", action="store_true", help="Don't print verbose progress updates")
    args = parser.parse_args()
    
    if args.chrom:
        chroms = args.chrom.split(",")
    else:
        chroms = None
    
    with open(args.vcf_file) as vcff:
        for line in vcff:
            # header line
            if line.startswith("##"):
                print(line, end="")
                continue
    
            tokens = line.rstrip().split('\t')
            
            # main header line containing sample names
            if tokens[0] == "#CHROM":
                samples = tokens[9:]
                # load homref RLE data for each sample
                homrefs = load_homref_files(samples, args.homref_dir, chroms, not args.quiet)
                print(line, end="")
                continue
            
            # regular line
            this_chrom = tokens[0]
            if chroms and this_chrom not in chroms:
                continue
    
            this_pos = int(tokens[1])
            for samp_idx in range(len(samples)):
                if tokens[9+samp_idx] == "./." and homrefs[samp_idx][this_chrom].contains(this_pos):
                    tokens[9+samp_idx] = "0/0"
    
            print("\t".join(tokens))
