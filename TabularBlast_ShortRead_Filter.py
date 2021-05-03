#!/usr/bin/env python

'''Best hit filter for tabular blast output.

Expects tabular blast to be run with following collumns:

-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend 
sstart send evalue bitscore qlen slen qseq sseq'

This filter will select the best blast match based on bitscore.
In the case of tied scores, the script randomizes the selection.

This filter also filters for minimum read length as well as a
percemt match length (PML) defined as alignment length / read length.

Additionally, this script has an option to remove all tied matches.
This is useful for competitive read mapping where reads may map equally
well to more than one genome. The shared reads can be removed.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Wednesday, August 28th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random
from collections import defaultdict
from collections import OrderedDict


def best_hits(query, bitscore, d, line, dups, pID):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query][0].split('\t')[11])

        if bitscore > old_bitscore:
            d[query] = [line]

        elif bitscore == old_bitscore:
            d[query].append(line)

    else:
        d[query] = [line]

    return d, dups


def blast_filter(infile, pml, rl):
    """Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    infile : str
        The file location of the input file
    pml : float
        Percent match length of the query read to filter above. (ex: 0.9)
        alignment length / query length
    rl : int
        Length of the query read to filter above. (ex: 70)

    Returns
    -------
    file
        writes out tabular blast lines passing the filter to 
        a new file infile.filtered_best_hits.blst
    """

    header = []
    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    fails = 0 # counter for number of matches failing filters
    passes = 0 # counter for number of matches passing filters
    total = 0 # counter for total blast entries in file

    with open(infile, 'r') as f:

        for l in f:
            if l.startswith('#'): header.append(l)
            else:
                total += 1
                X = l.rstrip().split('\t')
                query = X[0] # read identifier
                bitscore = float(X[11]) # bitscore
                pID = float(X[2]) # percent sequence identity
                aLen = int(X[3]) # read alignment length
                qLen = int(X[12]) # full length of read
                pMatch = aLen / qLen # percent match length of read length

                if pMatch >= pml and qLen >= rl:
                    d, dups = best_hits(query, bitscore, d, l, dups, pID)
                    passes += 1
                else:
                    fails += 1

    print('Total number of entries in blast file:', total)
    print('Number of reads failing the filters:', fails)
    print('Number of reads passing the filters:', passes)
    print('Number of duplicate blast matches passing filter to remove:', dups)

    return header, d


def remove_tied_matches(filtered_best_hits, o):
    '''Removes matches that are tied for bitscore'''

    print(
        '\nRunning -rtm option to remove reads with tied blast hits.\n'
        'The tied read count should be higher than the duplicate read count.\n'
        '\nPrinting counts of tied blast hits per pID bin:\n'
        )
    # intialize dict to count ties by whole step pID
    tied_match_count = defaultdict(int)
    # count total duplicate entries
    tiecount = 0
    # count reads written to file
    writecount = 0
    # iterate over dict of reads {read: [line(s)]}
    for k,v in filtered_best_hits.items():
        # if only 1 entry for read, write the line.
        if len(v) == 1:
            o.write(v[0])
            writecount += 1
        # if more than 1 entry add count to pIDs
        else:
            for i in v:
                X = i.split('\t')
                #query = X[0]
                pID = int(X[2].split('.')[0])
                #bitscore = float(X[12])
                #print(query, pID, bitscore)
                tied_match_count[pID] += 1
                tiecount += 1
    # sort tied_match_count by descending pID
    od_tied_match_count = OrderedDict(
                                    sorted(
                                        tied_match_count.items(), reverse=True
                                        )
                                    )
    # print count of ties per pID whole step
    print('pID\tCount')
    for k,v in od_tied_match_count.items():
        print(f'{k:03}% ID:\t{v}')
    # print total ties and total reads written
    print(
        '\nNumber of tied read matches removed:', tiecount,
        '\nNumber of reads written to file:', writecount, '\n\n'
        )


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pml', '--percent_match_length',
        help='(Optional) Percent match length to filter for (Default = 0.9).',
        metavar='',
        type=float,
        required=False,
        default=0.9
        )
    parser.add_argument(
        '-rl', '--read_length',
        help='(Optional) Read length to filter for (Default = 70).',
        metavar='',
        type=float,
        required=False,
        default=70
        )
    parser.add_argument(
        '-rtm', '--remove_tied_matches',
        help=
            '(Optional) Use this flag to remove read alignments when there is '
            'a tie for the best match. Useful for competitive recruitements.',
        action='store_true',
        required=False,
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')
    header, filtered_best_hits = blast_filter(
                                            args['in_file'],
                                            args['percent_match_length'],
                                            args['read_length'],
                                            )

    # Write output file
    outfile = args['in_file'].split('.')[0] + '.fltrdBstHts.blst'
    with open(outfile, 'w') as o:
        for l in header: o.write(l)
        if args['remove_tied_matches']:
            remove_tied_matches(filtered_best_hits, o)
        else:
            for k,v in filtered_best_hits.items(): o.write(random.choice(v))
            print(
                'Number of best hit entries written to new file:',
                len(filtered_best_hits), '\n\n'
                )


if __name__ == "__main__":
    main()
