#!/usr/bin/env python3
#
# Copyright 2012 Paolo Martini <paolo.cavei@gmail.com>

from __future__ import division
from optparse import OptionParser
from os.path import basename
from sys import argv, exit, stdin

def collect(symbols):
    counts = {}
    total = 0
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
        total += 1
    return counts, total

def main():
    parser = OptionParser(usage='%prog <SYMBOLS')
    parser.add_option('-r', '--reverse', dest='reverse', action='store_true', default=False, help='print count before symbol')    
    parser.add_option('-d', '--double', dest='double', action='store_true', default=False, help='same as symbol_count | cut -f 2 | symbol_count')
    parser.add_option('-n', '--normalize', dest='normalize', action='store_true', default=False, help='output fraction instead of count')
    options, args = parser.parse_args()
    if len(args) != 0:
        exit('Unexpected argument number.')
    
    mode = 'count' if basename(argv[0]) == 'symbol_count' else 'freq'
    if options.normalize:
        mode="freq"
    if options.double and mode == 'freq':
        exit("Double option is not supported in freq mode (normalize).")

    counts, total = collect(l.rstrip("\r\n") for l in stdin)
    if options.double:
        counts, total = collect(counts.itervalues())

    for symbol in sorted(counts.keys()):
        count = counts[symbol]
        if mode == 'count':
            if options.reverse:
                print(f'{count}\t{symbol}')
            else:
                print (f'{symbol}\t{count}')
        else:
            if options.reverse:
                print (f'{count/total}\t{symbol}')
            else:
                print (f'{symbol}\t{count/total}')

if __name__ == '__main__':
    main()
