#!/usr/bin/env python

import sys

def readTaxonomyNames(names_dmp):
    nodeNames = []
    for line in open(names_dmp):
        line = line.strip().strip('|').strip()
        if not line: break
        line = line.split('\t|\t')
        if line[3].strip() == 'scientific name':
            nodeNames.append((int(line[0].strip()), line[1].strip()))        
        pass
    return dict(nodeNames)

def readTaxonomyNodes(nodes_dmp):    
    nodeRanks = []
    nodeChildren = {}
    nodeParents = {}
    for line in open(nodes_dmp):
        line = line.strip().strip('|').strip()
        if not line: break
        line = map(lambda x:x.strip(), line.split('\t|\t'))
        line[:2] = map(int, line[:2])
        if line[0] == 1:
            line[1] = 1    

        nodeParents[line[0]] = line[1]
        try:
            nodeChildren[line[1]].add(line[0])
        except:
            nodeChildren[line[1]] = set([line[0]])
        nodeRanks.append((line[0], line[2]))

    return dict(nodeRanks), nodeChildren, nodeParents

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