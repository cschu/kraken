#!/usr/bin/env python
import sys
import argparse

def anabl_getContigsFromFASTQ(fn):
    with open(fn) as fi:
    	for head in fi:
        	try:
            	seq, sep, qual = map(lambda x:x.strip(), [fi.next() for i in xrange(3)])
            	yield head[1:], seq, qual
        	except:
            	break
    pass

def anabl_getContigsFromFASTA(fn):
	with open(fn) as fi:
    	for head in fi:
        	try:            	
            	yield head[1:], fi.next().strip()
        	except:
            	break
    pass

def toFASTQ(fields):
	return '@%s\n%s\n+\n%s\n' % fields
def toFASTA(fields):
	return '>%s\n%s\n' % fields
def getFASTQID(head):
	return head.strip().split()[0].strip('@').strip('/1')
def getFASTAID(head):
	return head.strip().split()[0].strip('>').strip('/1')



def main(argv):
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--db', default=DEFAULT_KRAKEN_DB, help='Path to kraken db.')
	parser.add_argument('--in1', default='', help='Forward read file from metagenomic sample.')
	parser.add_argument('--in2', default='', help='Reverse read file from metagenomic sample.')
	parser.add_argument('--kraken-results1', help='File containing kraken classification for set1.')
	parser.add_argument('--kraken-results2', help='File containing kraken classification for set2.')
	parser.add_argument('--input-format', default='fa', help='fa (FASTA) or fq (FASTQ)')
	parser.add_argument('--set1-output-left', type=str)
	parser.add_argument('--set1-output-right', type=str)
	parser.add_argument('--set2-output-left', type=str)
	parser.add_argument('--set2-output-right', type=str)
	parser.add_argument('--unclassified-output-left', type=str)
	parser.add_argument('--unclassified-output-right', type=str)
	parser.add_argument('--intersection-output-left', type=str)
	parser.add_argument('--intersection-output-right', type=str)	
	args = parser.parse_args()

	if args.input_format == 'fa':
		getContigs = anabl_getContigsFromFASTA  
		writeFXRecord = toFASTA
		getFXID = getFASTAID
	elif args.input_format == 'fq':
		getContigs = anabl_getContigsFromFASTQ
		writeFXRecord = toFASTQ
		getFXID = getFASTQID
	else:
		sys.stderr.write('Error: Unknown input format (%s).\n' % args.input_format)
		sys.exit()

	try:
		set1 = [line.strip().split()[:2] for line in open(args.kraken_results1)]
	except:
		sys.stderr.write('Error: Set1 is missing, please specify --kraken_results1 parameter.\n')
		sys.exit()

	nReads = len(set1)
	set1 = set(filter(lambda x:x[0].startswith('C'), set1))

	try:
		set2 = set([line.strip().split()[1] for line in open(args.kraken_results2) if line.startswith('C')])
	except:
		sys.stderr.write('Error: Set2 is missing, please specify --kraken_results2 parameter.\n')
		sys.exit()

	try:
		set1_fwd = open(args.set1_output_left, 'wb')		
		set2_fwd = open(args.set2_output_left, 'wb')
		noclass_fwd = open(args.unclassified_output_left, 'wb')
		undecided_fwd = open(args.intersection_output_left, 'wb')		
	except:
		sys.stderr.write('Error: Cannot open fwd outputfile(s).\n')
		sys.exit()

	try:
		fwd = getContigs(args.in1) 
	except:
		sys.stderr.write('Error: Cannot open file %s.\n' % args.in1)
		sys.exit()

	rev = None
	if args.in2:
		try:
			rev = getContigs(args.in2)
		except:
			sys.stderr.write('Error: Cannot open file %s.\n' % args.in2)
			sys.exit()			
		try:
			set1_rev = open(args.set1_output_right, 'wb')
			set2_rev = open(args.set2_output_right, 'wb')
			noclass_rev = open(args.unclassified_output_right, 'wb')
			undecided_rev = open(args.intersection_output_right, 'wb')
		except:
			sys.stderr.write('Error: Cannot open rev outputfile(s).\n')
			sys.exit()			


for i in xrange(nReads):
	try:
		r1 = fwd.next()
	except:
		break
	r2 = None
	if rev is not None:
		try:
			r2 = rev.next()
		except:
			break
	
	id_ = getFXID(r1)

    if id_ in set1:
        if id_ not in set2:
            set1_fwd.write(writeFXRecord(r1))
            if set1_rev: 
            	set1_rev.write(writeFXRecord(r2))
        else:
            undecided_fwd.write(writeFXRecord(r1))
            if undecided_rev: 
            	undecided_rev.write(writeFXRecord(r2))
            
    elif id_ in set2:
        set2_fwd.write(writeFXRecord(r1))
        if set2_rev: 
        	set2_rev.write(writeFXRecord(r2))
    else:
        noclass_fwd.write(writeFXRecord(r1))
        if noclass_rev: 
        	noclass_rev.write(writeFXRecord(r2))
    pass

	set1_fwd.close()
	set1_rev.close()
	set2_fwd.close()
	set2_rev.close()
	noclass_fwd.close()
	noclass_rev.close()
	undecided_fwd.close()
	undecided_rev.close()
	pass

if __name__ == '__main__': main(sys.argv)