#!/usr/bin/env python

import argparse, sys, re, csv
from operator import itemgetter

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='format blast output (outfmt 6 or 7) to order by position and bscore per 1kb of hsp')
    parser.add_argument('-i', '--input', help='indicate input file (as blast tab outfmt 6 or 7)', required=True)
    parser.add_argument('-o', '--output', help='indicate output file', required=True)
    args = parser.parse_args()

    pstart = re.compile('\w+:(\d+)-\d+')
    tosort = []
    with open(args.input, 'rU') as file_in:
        file_out = open(args.output, 'wb+')
        file_out.write('gstart\tpstart\tbscore\tbscore/kb\tsubjectID\n')
        reader = csv.reader(file_in, delimiter = '\t')
        for row in reader:
        	if '#' in row[0]:
        		continue
        	else:
        		col0 = int(pstart.search(row[0]).group(1)) + 1
        		col1 = col0 + int(row[6]) - 1
        		col2 = row[11]
        		col3 = (int(row[11])*1000)/int(row[3])
        		col4 = row[1]
        		nurow = [col0, col1, col2, col3, col4]
        		tosort.append(nurow)
    file_in.close()
    tosort.sort(key=itemgetter(1), reverse=False)
    for row in tosort:
    	file_out.write(str('\t'.join(str(x) for x in row)) + '\n')
    file_out.close()

if __name__ == '__main__':
    main()