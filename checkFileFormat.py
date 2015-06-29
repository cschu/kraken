#!/usr/bin/env python


def verifyFileFormat(fn, fileFormat):
    firstChar = open(fn).read(1)

    verifiedFastq = firstChar == '@' and fileFormat == 'fq'
    verifiedFasta = firstChar == '>' and fileFormat == 'fa'

    return verifiedFastq or verifiedFasta
