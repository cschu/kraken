#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import struct
import shutil
import tempfile

# from collections import Counter

from checkFileFormat import verifyFileFormat 
DEFAULT_KRAKEN_DB = '/tsl/data/krakendb/ktest/db1'

def main(argv):
    
    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--dbtype', default='builtinDB', help='Is database built-in or supplied externally?')
    parser.add_argument('--db', default=DEFAULT_KRAKEN_DB, help='Path to kraken db')
    parser.add_argument('--in1', help='Single-end reads or left paired-end reads')
    parser.add_argument('--in2', help='Right paired-end reads')
    parser.add_argument('--quick', action='store_true', help='Use quick mode?')
    parser.add_argument('--min-hits', default=1, help='Minimum number of hits required.')
    parser.add_argument('--input-format', help='Input sequences stored in fa or fq file(s).', default='fq')
    parser.add_argument('kraken_summary_tsv', type=str, help='The output file.')
    # parser.add_argument('classified_seqs_fq', type=str, help='The output file.')
    # parser.add_argument('unclassified_seqs_fq', type=str, help='The output file.')
    args = parser.parse_args()

    #  kraken --preload --db /tsl/scratch/macleand/ktest/db --threads 12 --quick --classified-out classified_seqs.fq --unclassified-out unclassified.fq --fastq-input --min-hits 1 --output classification.txt left_reads.fq right_reads.fq

    kraken_params = ['--preload', '--threads', '8', 
    #                 '--unclassified-out', args.unclassified_seqs_fq,
                     '--output', args.kraken_summary_tsv]
    #                '--classified-out', args.classified_seqs_fq,
    kraken_input = []

    if 'db' not in args or not os.path.exists(args.db):
        sys.stderr.write('Error: database is missing.\n')
        sys.exit(1)
    kraken_params.extend(['--db', args.db])
    # check whether input file(s) exist(s)
    if not ('in1' in args and os.path.exists(args.in1)):
        sys.stderr.write('Error: fwd/single input file (%s) is missing.\n' % args.in1)
        sys.exit(1)
    if not verifyFileFormat(args.in1, args.input_format):
        fc = open(args.in1).read(1)
        sys.stderr.write('Error: fwd/single input file has the wrong format assigned (%s). Starts with \'%s\'\n' % (args.input_format, fc))
        sys.exit(1)
    kraken_input.append(args.in1)
    if 'in2' in args:
        if args.in2 is not None and not os.path.exists(args.in2):
            sys.stderr.write('Error: right input file (%s) is missing.\n' % args.in2)
            sys.exit(1)
        elif args.in2 is not None and not verifyFileFormat(args.in2, args.input_format):
            sys.stderr.write('Error: rev input file has the wrong format assigned.\n')
            sys.exit(1)            
        elif args.in2 is not None:
            kraken_params.append('--paired')
            kraken_input.append(args.in2)
        else:
            pass
    else:
        pass

    if 'quick' in args:
        kraken_params.append('--quick')
        if 'min_hits' in args:
            kraken_params.extend(['--min-hits', str(args.min_hits)])

    if args.input_format == 'fq': 
        kraken_params.append('--fastq-input')
    else:
        kraken_params.append('--fasta-input')
      
    
    
    
    
    
    
    
    
    
    
    
    

    # check whether file is gzipped
    header = ''
      
    cmd = 'source kraken-0.10.5; '
    cmd += ' '.join(['kraken'] + kraken_params + kraken_input) + '\n'
    # out = open(argv[-1], 'wb').write(str(cmd) + '\n')#.write('KRAKEN OUTPUT, HAS COUNTER\n')
    # out = open(argv[-1], 'wb').write(str(args) + '\n')#.write('KRAKEN OUTPUT, HAS COUNTER\n')

    # proc = subprocess.Popen(args=cmd, shell=True, stderr=sys.stderr)
    # returncode = proc.wait()


    tmp_dir = tempfile.mkdtemp()
    try:
        tmp = tempfile.NamedTemporaryFile(dir=tmp_dir).name
        tmp_stderr = open(tmp, 'wb')
        proc = subprocess.Popen(args=cmd, shell=True, stderr=tmp_stderr.fileno())
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open(tmp, 'rb')
        stderr = ''
        buffsize = 1048576
        try:
            while True:                    
                stderr += tmp_stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr
    except Exception, e:
        # clean up temp dirs
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        sys.stderr.write('Error running kraken: %s\n' % str(e))
        sys.exit(1)

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    """
    open(args.kraken_summary_tsv, 'wb').write('\t'.join(list('ACGT')) + '\n')
    open(args.classified_seqs_fq, 'wb').write(cmd + '\n')
    open(args.unclassified_seqs_fq, 'wb').write('blruah\n')                      # check whether the database exists, if not exit with error
    """
    """
    for arg in args:
        out.write(str(arg) + '\n') 
    out.close()
    """
    pass


 
    
# main(sys.argv[1:])

if __name__ == '__main__': main(sys.argv[1:])
