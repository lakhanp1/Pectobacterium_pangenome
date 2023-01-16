#!/usr/bin/env python3

import argparse
import logging
import os
# from pathlib import Path


# argument parsing
def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge all the FASTA sequences from genomes into single FASTA file",
        prog="merge_genomes_faa"
    )
    parser.add_argument(
        "--genomes-file",
        help="Genome file with path to each genome on each line",
        type=str, required=True
    )

    parser.add_argument(
        "--out",
        help=("Output file to store combined genomes FASTA records."
              " Additionally, 'genome_chr_map.tab' file is created holding the"
              " genome to chromosome mapping table."),
        default="genomes_merged.fasta"
    )

    args = parser.parse_args()

    return args


def merge_fasta(file_genomes, output):
    logging.info(f'merging all genomes into {output}')
    
    file_chr_key = os.path.join(os.path.dirname(output), 'genome_chr_map.tab')
    logging.info(f'storing genome to chr mapping in file {file_chr_key}')
    
    fa_writer = open(output, mode='w')
    key_writer = open(file_chr_key, mode='w')
    
    chrset = 0

    with open(file_genomes) as genomes:
        for genome_id, genome_fa in enumerate(genomes, start=1):
            logging.info(f'{genome_id} + {genome_fa.strip()}')
            
            # store chromosome IDs for 100 geneme batches into separate file for
            # 'blastn -seqidlist' option search
            if(genome_id % 100 == 1):
                if(chrset != 0):
                    chrset_writer.close()
                    
                chrset += 1
                file_chrset = os.path.join(os.path.dirname(output), f'chrset_{chrset}.acc')
                chrset_writer = open(file_chrset, mode='w')

            
            with open(genome_fa.strip()) as fasta:
                for line in fasta:
                    
                    if line.startswith('>'):
                        # add genome id prefix to chr id
                        new_header = f'{genome_id}_{line.strip()[1:]}'
                        fa_writer.write(f'>{new_header}\n')
                        key_writer.write(f'{genome_id}\t{new_header}\n')
                        
                        chrset_writer.write(f'{new_header}\n')
                        
                    else:
                        fa_writer.write(line)
                    
            
    fa_writer.close()
    key_writer.close()
    chrset_writer.close()
    

def main():
    logging.basicConfig(level=logging.DEBUG)
    args = parse_args()
    # print(args)
    merge_fasta(file_genomes=args.genomes_file, output=args.out)


if __name__ == "__main__":
    main()
