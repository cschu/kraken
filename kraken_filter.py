#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import struct
import shutil
import tempfile

from checkFileFormat import verifyFileFormat

DEFAULT_KRAKEN_DB='/tsl/data/krakendb/ktest/db'
logfile = None


def getDescendents(taxID, tree):
    descendents = set([taxID])
    queue = [taxID]
    while queue:
        node = queue.pop()

        children = tree.get(node, set())
        if children:
            descendents = descendents.union(children)
            queue.extend(children)
        pass
    return descendents

def getSeqsFromFastq(fn):
    it = open(fn)
    head = ''
    block = []
    for line in it:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            head = line
        elif line.startswith('+'):
            seq = ''.join(block)
            qual = ''.join([it.next().strip() for row in block])
            yield (head, (seq, '+', qual))
            block, head = [], ''
        else:
            block.append(line)

def getSeqsFromFasta(fn):
    it = open(fn)
    head = ''
    block = []
    for line in it:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if block:
                yield (head, (''.join(block),))
                block = []
            head = line
        else:
            block.append(line)
    if block:
        yield (head, (''.join(block),))

def anabl_getSeqsFromFastX(fn, X=2):
    it = open(fn)
    for head in it:
        if not head.strip():
            continue
        try:
            yield (head.strip(), tuple(map(lambda x:x.strip(),
                                       ([it.next() for i in xrange(X - 1)]))))
        except:
            break
        pass
    it.close()
    pass

def readTaxonomy(db):
    nodeNames = []
    for line in open(os.path.join(db, 'taxonomy', 'names.dmp')):
        line = line.strip().strip('|').strip()
        if not line: break
        line = line.split('\t|\t')

        if line[3].strip() == 'scientific name':
            nodeNames.append((int(line[0].strip()), line[1].strip()))
        # break
        pass

    nodeRanks = []
    nodeChildren = {}
    for line in open(os.path.join(db, 'taxonomy', 'nodes.dmp')):
        line = line.strip().strip('|').strip()
        if not line: break
        line = map(lambda x:x.strip(), line.split('\t|\t'))
        line[:2] = map(int, line[:2])
        if line[0] == 1:
            line[1] = 0
        try:
            nodeChildren[line[1]].add(line[0])
        except:
            nodeChildren[line[1]] = set([line[0]])
        nodeRanks.append((line[0], line[2]))
        # break
    return dict(nodeNames), dict(nodeRanks), nodeChildren

def isPreCassava18(string):
    return string.endswith('/1') or string.endswith('/2')
def getFastqIdentifier(string):
    return string[1:-2] if isPreCassava18(string) else string.split()[0][1:]
def getFastaIdentifier(string):
    string = string.split()[0][1:]
    return string[:-2] if isPreCassava18(string) else string

def runFilter(db, taxID, inputClassification, inputR1, outputR1, inputR2=None, outputR2=None, fastx=4, debug=False):
    global logfile
    keepSequences = set()
    if taxID != 0:
        logfile.write('Reading taxonomy...\n')
        logfile.flush()
        taxonomyInfo = readTaxonomy(db)
        logfile.write('Traversing requested taxonomy branch...\n')
        logfile.flush()
        validIDs = getDescendents(abs(taxID), taxonomyInfo[2])
        logfile.write('Family tree of %s has %i descendents.\n' % (str(abs(taxID)), len(validIDs)))
        logfile.flush()
 
        if debug:
            for vid in validIDs:
                logfile.write('Added %s to validIDs.\n' % vid)
                logfile.flush()
    else:
        validIDs = set()	

    logfile.write('Filtering sequences...\n')
    logfile.flush()
    nseqs = 0
    for line in open(inputClassification):
        nseqs += 1
        line = line.strip().split()

        if taxID > 0:
            # if extract reads from branch
            # no unclassified reads
            takeUnclassified = False # taxID == 0 and line[0] == 'U'
            # and only reads that have been assigned a taxonomy id belonging to the branch
            takeClassified = line[0] == 'C' and int(line[2]) in validIDs        
        elif taxID < 0:
            # if ignore reads from branch
            # allow unclassified reads
            takeUnclassified = line[0] == 'U'
            # and take only reads that have been a taxonomy id outside of the branch
            takeClassified = line[0] == 'C' and int(line[2]) not in validIDs            
        else:
            # if extract unclassified
            # allow unclassified reads
            takeUnclassified = taxID == 0 and line[0] == 'U' 
            # and no classified ones
            takeClassified = False

        if takeUnclassified or takeClassified:            
            keepSequences.add(line[1].strip())
            if debug:
                logfile.write('Added %s to keepSequences.\n' % line[1].strip())
                logfile.flush()

    logfile.write('Keeping %i of %i sequences (%.1f).\n' % (len(keepSequences), nseqs, float(len(keepSequences))/nseqs))
    logfile.flush()
    logfile.write('Writing filtered sequence sets...\n')
    logfile.flush()
   
    if fastx == 4:
        getID = getFastqIdentifier 
        getSeqs = getSeqsFromFastq
    else:
        getID = getFastaIdentifier
        getSeqs = getSeqsFromFasta


    fwdOut = open(outputR1, 'wb')
    fwdGen = getSeqs(inputR1) # anabl_getSeqsFromFastX(inputR1, X=fastx)

    revOut, revGen = None, None
    revSid, revSeq = None, None
    if outputR2 is not None and inputR2 is not None:
        revOut = open(outputR2, 'wb')
        revGen = getSeqs(inputR2) # anabl_getSeqsFromFastX(inputR2, X=fastx)


    fxid1, fxid2 = None, None
    while True:
        try:
            fwdSid, fwdSeq = fwdGen.next()
        except:
            break
        logfile.write('*%s*\n' % fwdSid)
        fxid1 = getID(fwdSid)
        if revGen is not None:
            try:
                revSid, revSeq = revGen.next()
            except:
                break
            fxid2 = getID(revSid)

        if fxid1 != fxid2 and fxid2 is not None:
            sys.stderr.write('Error: fxid-mismatch %s %s.\n' % (fxid1, fxid2))
            sys.exit(1)
        if fxid1 in keepSequences:
            fwdOut.write(('%s\n' * fastx) % ((fwdSid,) + fwdSeq))
            if revOut is not None:
                revOut.write(('%s\n' * fastx) % ((revSid,) + revSeq))
        else:
            # sys.stdout.write('%s is not in keepSequences\n' % fqid1)
            pass
    fwdOut.close()
    if revOut is not None:
        revOut.close()
    pass


def main(argv):

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--dbtype', default='builtinDB', help='Is database built-in or supplied externally?')
    parser.add_argument('--db', default=DEFAULT_KRAKEN_DB, help='Path to kraken db')
    parser.add_argument('--in1', help='The r1-file (single-end reads or left paired-end reads).')
    parser.add_argument('--in2', help='The r2-file (right paired-end reads)')
    parser.add_argument('--taxid', default=1, help='The root taxon id of the requested branch', type=int)
    parser.add_argument('--kraken-results', type=str, help='A file containing kraken classification results for the input sequences.')
    parser.add_argument('--input-format', help='Input sequences stored in fa or fq file(s).', default='fq')

    parser.add_argument('--out1', type=str, help='The r1-output file.')
    parser.add_argument('--out2', type=str, help='The r2-output file.')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--logfile', type=str, help='A logfile.', default='kraken_filter.log')
    args = parser.parse_args()

    kraken_params = []
    kraken_input = []

    global logfile
    logfile = open(args.logfile, 'wb')
    if 'db' not in args or not os.path.exists(args.db):
        sys.stderr.write('Error: database is missing.\n')
        sys.exit(1)
    kraken_params.extend(['--db', args.db])
    # check whether input file(s) exist(s)

    if 'kraken_results' in args and os.path.exists(args.kraken_results):
        kraken_input.append(args.kraken_results)
    else:
        sys.stderr.write('Error: kraken-classification is missing.\n')
        sys.exit(1)
        pass

    if not ('in1' in args and os.path.exists(args.in1)):
        sys.stderr.write('Error: left/single input file (%s) is missing.\n' % args.in1)
        sys.exit(1)
    if not verifyFileFormat(args.in1, args.input_format):
        sys.stderr.write('Error: fwd/single input file has the wrong format assigned.\n')
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

    if 'out2' in args and 'in2' in args:
        input2, output2 = args.in2, args.out2
    else:
        input2, output2 = None, None

    # if args.taxid < 0:
    #    sys.stderr.write('Error: invalid taxon id %i\n' % args.taxid)
    #    sys.exit(1)
    kraken_params.extend(['--taxid', args.taxid])

    if args.input_format == 'fq':
        kraken_params.append('--fastq-input')
        fastx = 4
    else:
        kraken_params.append('--fasta-input')
        fastx = 2



    # firstChar = open(args.in1).read(1)
    # if firstChar == '>':
    #    kraken_params.append('--fasta-input')
    #    fastx = 2
    # elif firstChar == '@':
    #    kraken_params.append('--fastq-input')
    #    fastx = 4
    #else:
    #    """ This will currently disable working with gzipped/bzipped files. """
    #    sys.stderr.write('Error: Input file starts with unknown symbol %s.\n' % firstChar)
    #    sys.exit(1)
    #    pass

    #open(args.kraken_filtered_r1, 'wb').write('\n'.join(map(str, kraken_params + kraken_input)))

    runFilter(args.db, int(args.taxid), args.kraken_results,
              args.in1, args.out1,
              inputR2=input2, outputR2=output2, fastx=fastx,
              debug=args.debug)

    logfile.close()

    pass




# main(sys.argv[1:])

if __name__ == '__main__': main(sys.argv[1:])
