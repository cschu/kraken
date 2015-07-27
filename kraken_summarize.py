#!/usr/bin/env python
import sys
import subprocess
import argparse

from collections import Counter



# U	HWI-ST933:109:C2A2NACXX:5:1101:1311:2013	0	50	Q:0
def fetchTaxonomyData(ids, email='christian.schudoma@tsl.ac.uk'):
    from Bio import Entrez
    Entrez.email = email
    handle = Entrez.efetch(db='Taxonomy', id=','.join(ids), retmode='xml')
    records = Entrez.read(handle)
    return records

def writeKronaInput(fi, taxInfo, unclassified=0):
    if unclassified:
        fi.write('%i\tUnclassified\n' % unclassified)
    for tid in sorted(taxInfo, key=lambda x:taxInfo[x][0]['Lineage']):
        fi.write('%i\t%s\n' % (taxInfo[tid][1], '; '.join([taxInfo[tid][0]['Lineage'].strip(), taxInfo[tid][0]['ScientificName']]).replace('; ', '\t').strip('\t')))
    pass

def writeOutput(out, taxInfoDict, c):
    for tid in sorted(taxInfoDict, reverse=True, key=lambda x:taxInfoDict[x][1]):
        data = [tid, taxInfoDict[tid][1], taxInfoDict[tid][0]['TaxId'], taxInfoDict[tid][0]['Lineage'], taxInfoDict[tid][0]['ScientificName']]
        out.write('\t'.join(map(str, data)) + '\n')
    data = (sum(c.values()) - c['0'], sum(c.values()), (sum(c.values()) - c['0']) / float(sum(c.values())) * 100, c['0'])
    try:
        out.write('%i/%i (%.5f%%) classified, %i unclassified\n' % data)
    except:
        out.write(str(data) + ':ERR\n')

def main(argv):

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--krona-output', type=str, default='')
    parser.add_argument('--output', type=str)
    parser.add_argument('--call-krona', type=str)
    parser.add_argument('--include-unclassified')
    parser.add_argument('kraken_summary_tsv', type=str)
    args = parser.parse_args()

    #with open(args.output, 'wb') as fo:
    #    fo.write(sys.version + '\n')
    #    import os
    #    fo.write(os.environ.get('PYTHONPATH', 'N/A') + '\n')
    c = Counter(line.strip().split()[2] for line in open(args.kraken_summary_tsv))

     
    taxids = sorted(c.keys(), key=lambda x:c[x], reverse=True)
    taxData = fetchTaxonomyData(taxids)
    taxInfoDict = {tinfo['TaxId']: [tinfo, c[tinfo['TaxId']]] for tinfo in taxData}

    kr_out = None
    if 'krona_output' in args:
        with open(args.krona_output, 'wb') as kr_out:
            if 'include_unclassified' in args:
                writeKronaInput(kr_out, taxInfoDict, unclassified=c['0'])
            else:
                writeKronaInput(kr_out, taxInfoDict)

        if 'call_krona' in args:
            cmd = 'source krona_tools-2.4; ktImportText -o %s %s' % (args.call_krona, args.krona_output)            
            # cmd = ['source', 'krona_tools-2.4;', 'ktImportText', '-o', args.call_krona, args.krona_output] 
            p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            so, se = p.communicate()
            """
            with open(args.output, 'wb') as out: 
                out.write(str(se))
                return None
            """
    if 'output' in args:
        if args.output == '-':
            # writeOutput(sys.stdout, taxInfoDict, c)
            pass
        else:
            with open(args.output, 'wb') as out:
                writeOutput(out, taxInfoDict, c)
    




if __name__ == '__main__': main(sys.argv)



"""
2	Fats	Saturated fat
3	Fats	Unsaturated fat	Monounsaturated fat
3	Fats	Unsaturated fat	Polyunsaturated fat
13	Carbohydrates	Sugars
4	Carbohydrates	Dietary fiber
21	Carbohydrates
5	Protein
4
"""



"""
handle = Entrez.efetch(db="Taxonomy", id="1323524", retmode="xml")
records  = Entrez.read(handle)
records[0]["Lineage"]
'Viruses; dsRNA viruses; Partitiviridae; Betapartitivirus'
records[0]["LineageEx"]
[{u'ScientificName': 'Viruses', u'TaxId': '10239', u'Rank': 'superkingdom'}, {u'ScientificName': 'dsRNA viruses', u'TaxId': '35325', u'Rank': 'no rank'}, {u'ScientificName': 'Partitiviridae', u'TaxId': '11012', u'Rank': 'family'}, {u'ScientificName': 'Betapartitivirus', u'TaxId': '1511809', u'Rank': 'genus'}]
"""

"""
{u'Lineage': 'Viruses; dsDNA viruses, no RNA stage; Iridoviridae; unclassified Iridoviridae', u'Division': 'Viruses', u'ParentTaxId': '180169', u'PubDate': '2014/03/16 07:01:27', u'LineageEx': [{u'ScientificName': 'Viruses', u'TaxId': '10239', u'Rank': 'superkingdom'}, {u'ScientificName': 'dsDNA viruses, no RNA stage', u'TaxId': '35237', u'Rank': 'no rank'}, {u'ScientificName': 'Iridoviridae', u'TaxId': '10486', u'Rank': 'family'}, {u'ScientificName': 'unclassified Iridoviridae', u'TaxId': '180169', u'Rank': 'no rank'}], u'CreateDate': '2014/02/24 15:56:28', u'TaxId': '1465751', u'Rank': 'species', u'GeneticCode': {u'GCId': '1', u'GCName': 'Standard'},
u'ScientificName': 'Anopheles minimus irodovirus', u'MitoGeneticCode': {u'MGCId': '0', u'MGCName': 'Unspecified'}, u'UpdateDate': '2014/02/24 15:56:29'}

"""
