#!/usr/bin/env python3

import argparse
import logging
import os
# from pathlib import Path


# argument parsing
def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename contigs/chromosome names to short and unique identifiers",
        prog="rename_fasta_headers"
    )
    parser.add_argument(
        "--fa",
        help="FASTA file path",
        type=str, required=True
    )
    
    parser.add_argument(
        "--out",
        help=("Output directory"),
        default="."
    )
    
    args = parser.parse_args()
    
    return args


def rename_fasta_headers(fa, out):
    
    logging.debug(f'Input FASTA file: {fa}')
    logging.debug(f'Output FASTA file: {out}')
    
    fa_writer = open(file=out, mode='w')
    
    chrN = 1
    
    with open(file=fa) as fasta:
        for line in fasta:
            
            if line.startswith('>'):
                new_header = f'>contig_{chrN}\n'
                chrN += 1
                fa_writer.write(new_header)
            
            else:
                fa_writer.write(line)
            
    

def main():
    logging.basicConfig(level=logging.DEBUG)
    args = parse_args()
    # print(args)
    
    outFa = os.path.join(args.out, os.path.basename(args.fa))
    
    rename_fasta_headers(fa=args.fa, out=outFa)


if __name__ == "__main__":
    main()
