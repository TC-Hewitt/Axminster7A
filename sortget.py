#!/usr/bin/env python

import argparse, sys, csv
from operator import itemgetter

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='sort by field and extract top val per name field. Outputs to STDOUT.')
    parser.add_argument('-i', '--input', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='tab delim file (leave out if using STDIN).')
    parser.add_argument('-sf', '--sortfield', type=int, required=True)
    parser.add_argument('-nf', '--namefield', type=int, required=True)
    parser.add_argument('-r', '--reverse', type=int, required=True)
    args = parser.parse_args()

    tab = csv.reader(args.input, delimiter = '\t', quoting=csv.QUOTE_NONE)
    sortLst = []
    qnames = set()

    def force_type(x):
        try:
            return(int(x))
        except ValueError:
            try:
                return(float(x))
            except ValueError:
                return(x)

    for row in tab:
        if '#' in row[0]:
            pass
        else:
            row[args.sortfield] = force_type(row[args.sortfield])
            sortLst.append(row)

    sortLst.sort(key=itemgetter(args.sortfield), reverse=args.reverse)

    for entry in sortLst:
        if entry[args.namefield] not in qnames:
            qnames.add(entry[args.namefield])
            print(str('\t'.join(str(x) for x in entry)))

if __name__ == '__main__':
    main()