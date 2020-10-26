#!/usr/bin/env python3

import sys
import re


def printResult(data, size, gcontent, rcontent):

    if not size and not gcontent and data["seq"] is not None:
        print(">" + data["contig"])
        print(data["seq"])
    else:
        print("SeqName: " + data["name"])
        if data["error"]:
            print("\tError")
            return
        if size:
            if data["contig"] is not None:
                print("\tContig: " + str(data["contig"]))
            print("\tSize: " + str(data["size"]))
        if gcontent:
            print("\tGC%: " +
                  str(round(data["gcontent"] / data["size"] * 100, 2)))
        if rcontent:
            print("\tRepeat%: " +
                  str(round(data["rcontent"] / data["size"] * 100, 2)))
        if data["seq"]:
            print("\tSeq: " + data["seq"])


def loadData(seq, stdin):
    data = dict()
    if not stdin:
        name = seq.split('/')[-1].split('.')[0]
    else:
        name = "stdin"
    data["name"] = name

    if stdin:
        fseq = seq
    else:
        fseq = open(seq, 'r')

    return fseq, data


def sizeAndGC(seq, gcontent, rcontent, stdin=False):
    fseq, data = loadData(seq, stdin)

    count = 0
    contig = 0
    if gcontent:
        GC = 0

    if rcontent:
        RC = 0

    for line in fseq:
        # line = line.strip().upper() # orignal
        line = line.strip()
        if line[0] != '>':
            count += len(line)
            if gcontent:
                # GC += line.count('G')  # original
                # GC += line.count('C')  # original
                GC += line.upper().count('G')
                GC += line.upper().count('C')

            if rcontent:
                RC += sum(1 for c in line if c.islower())

        else:
            contig += 1

    data["contig"] = contig
    data["size"] = count
    if gcontent:
        data["gcontent"] = GC

    if rcontent:
        data["rcontent"] = RC

    data["error"] = False
    data["seq"] = None
    fseq.close()
    return data


def contigSizeAndGC(seq, contig, gcontent, rcontent, stdin=False):
    regex = re.compile(contig)
    dataList = list()

    if stdin:
        fseq = seq
    else:
        fseq = open(seq, 'r')

    for line in fseq:
        if line[0] == '>' and regex.search(line) is not None:
            data = dict()
            if not stdin:
                name = seq.split('/')[-1].split('.')[0]
            else:
                name = "stdin"
            data["name"] = name
            count = 0
            if gcontent:
                GC = 0
            if rcontent:
                RC = 0
            data["name"] = data["name"] + "\n\tContig: " + line.strip()[1:]
            line = fseq.readline()  # read next line
            while line != '' and line[0] != '>':
                # line = line.strip().upper() # orignal
                line = line.strip()
                count += len(line)
                if gcontent:
                    # GC += line.count('G')  # orignal
                    # GC += line.count('C')  # orignal
                    GC += line.upper().count('G')
                    GC += line.upper().count('C')
                if rcontent:
                    RC += sum(1 for c in line if c.islower())

                line = fseq.readline()

            data["contig"] = None
            data["size"] = count
            data["seq"] = None
            if gcontent:
                data["gcontent"] = GC
            if rcontent:
                data["rcontent"] = RC

            data["error"] = False
            dataList.append(data)

    if len(dataList) == 0:
        if not stdin:
            name = seq.split('/')[-1].split('.')[0]
        else:
            name = "stdin"
        data = dict()
        data["name"] = name
        data["error"] = True
        dataList.append(data)

    fseq.close()
    return dataList


def extractSeq(seq, contig, gcontent, rcontent, start, stop, stdin=False):

    regex = re.compile(contig)

    fseq, data = loadData(seq, stdin)

    data["seq"] = ''
    if gcontent:
        GC = 0
    if rcontent:
        RC = 0

    line = fseq.readline()
    while line[0] != '>' or regex.search(line) is None:
        line = fseq.readline()
        if line == '':
            data["error"] = True
            return data

    data["contig"] = line.strip()[1:]
    line = fseq.readline()

    while line != '' and line[0] != '>':
        # line = line.strip().upper()  # original
        line = line.strip()
        data["seq"] += line
        line = fseq.readline()

    if start < stop:
        if stop == len(data["seq"]):
            data["seq"] = data["seq"][start - 1:]
        else:
            data["seq"] = data["seq"][start - 1: stop]
    else:
        if stop == 1:
            data["seq"] = data["seq"][start - 1:: -1]
        else:
            data["seq"] = data["seq"][start - 1: stop - 1: -1]
    data["size"] = len(data["seq"])
    if gcontent:
        # GC += data["seq"].count('G')  # original
        # GC += data["seq"].count('C')  # original
        GC += data["seq"].upper().count('G')
        GC += data["seq"].upper().count('C')
        data["gcontent"] = GC
    if rcontent:
        RC += sum(1 for c in data["seq"] if c.islower())
        data['rcontent'] = RC

    data["error"] = False
    fseq.close()
    return data


def main(args):

    if args.contig is None and args.extract is None and args.gcontent is None and args.rcontent is None and args.size is None:
        parser.error(
            "You should at least use one argument.\n{} --help for more info.".format(sys.argv[0]))

    if args.contig is None and args.extract is not None:
        parser.error("-e (--extract) requires -c (--contig).")

    dataList = list()
    if args.seqIn[0] == '-':
        if args.contig is None:
            try:
                data = sizeAndGC(sys.stdin, args.gcontent,
                                 args.rcontent, stdin=True)
            except IndexError:
                print("Is this really a fasta file?")
                return -1
            dataList.append(data)
        else:
            if args.extract is None:
                data = contigSizeAndGC(
                    sys.stdin, args.contig, args.gcontent, args.rcontent, stdin=True)
                dataList += data
            else:
                data = extractSeq(sys.stdin, args.contig, args.gcontent, args.rcontent,
                                  args.extract[0], args.extract[1], stdin=True)
                dataList.append(data)

        for data in dataList:
            printResult(data, args.size, args.gcontent, args.rcontent)
    else:
        for seq in args.seqIn:
            if args.contig is None:
                try:
                    data = sizeAndGC(seq, args.gcontent, args.rcontent)
                except IndexError:
                    print("Is this really a fasta file?")
                    return -1
                printResult(data, args.size, args.gcontent, args.rcontent)
            else:
                if args.extract is None:
                    data = contigSizeAndGC(
                        seq, args.contig, args.gcontent, args.rcontent)
                    dataList += data
                else:
                    data = extractSeq(
                        seq, args.contig, args.gcontent, args.rcontent, args.extract[0], args.extract[1])
                    dataList.append(data)
                for data in dataList:
                    printResult(data, args.size, args.gcontent, args.rcontent)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--size", help="show size.",
                        action="store_true", dest="size")
    parser.add_argument("-g", "--gcontent", help="show GC content.",
                        action="store_true", dest="gcontent")
    parser.add_argument("-r", "--rcontent", help="show repeats content",
                        action="store_true", dest="rcontent")
    parser.add_argument(
        "-c", "--contig", help="specify contig to work on, supports regular expressions.", action="store", dest="contig")
    parser.add_argument("-e", "--extract", help="sxtract sequence from START to STOP, works only with -c and stops at the first match.",
                        type=int, nargs=2, dest="extract", metavar=('start', 'stop'))
    parser.add_argument(
        'seqIn', type=str, help="genome, list of genome or stdin (uses '-' for stdin)", nargs='+')

    sys.exit(main(parser.parse_args()))
