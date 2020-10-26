#!/usr/bin/env python3

import sys
import re


def printResult(data, size):
    if not size and data["seq"] is not None:
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


def UTconv(seq, stdin=False):
    fseq, data = loadData(seq, stdin)
    count = 0
    contig = 0
    swaps = {'t': 'u', 'T': 'U', 'U': 'T', 'u': "t"}

    with open(data["name"]+".ut.fasta", 'a') as out_file:
        for line in fseq:
            if line[0] == '>':
                contig += 1
                out_file.write(line)
            else:
                count += len(line)
                line = "".join(swaps.get(i, i) for i in line)
                out_file.write(line)

    data["contig"] = contig
    data["size"] = count
    data["error"] = False
    data["seq"] = None
    fseq.close()
    out_file.close()
    return data


def contigUTconv(seq, contig, stdin=False):
    regex = re.compile(contig)
    dataList = list()
    swaps = {'t': 'u', 'T': 'U', 'U': 'T', 'u': "t"}

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

            count = 0
            data["name"] = name + "_Contig_" + line.strip()[1:].strip()
            out_file = open(data["name"]+".ut.fasta", 'a')
            out_file.write(line)
            line = fseq.readline()  # read next line
            while line != '' and line[0] != '>':
                line = line.strip()
                count += len(line)
                line = "".join(swaps.get(i, i) for i in line)
                out_file.write(line)

                line = fseq.readline()
            out_file.close()

            data["contig"] = None
            data["size"] = count
            data["seq"] = None
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


def main(args):
    if args.contig is None and args.size is None:
        parser.error(
            "You should at least use one argument.\n{} --help for more info.".format(sys.argv[0]))

    dataList = list()
    if args.seqIn[0] == '-':
        if args.contig is None:
            try:
                data = UTconv(sys.stdin, stdin=True)
            except IndexError:
                print("Is this really a fasta file?")
                return -1
            dataList.append(data)
        else:
            data = contigUTconv(
                sys.stdin, args.contig, stdin=True)
            dataList += data

        for data in dataList:
            printResult(data, args.size)
    else:
        for seq in args.seqIn:
            if args.contig is None:
                try:
                    data = UTconv(seq)
                except IndexError:
                    print("Is this really a fasta file?")
                    return -1
                printResult(data, args.size)
            else:
                data = contigUTconv(
                    seq, args.contig)
                dataList += data

                for data in dataList:
                    printResult(data, args.size)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--size", help="show size.",
                        action="store_true", dest="size")
    parser.add_argument(
        "-c", "--contig", help="specify contig to work on, supports regular expressions.", action="store", dest="contig")
    parser.add_argument(
        'seqIn', type=str, help="genome, list of genome or stdin (uses '-' for stdin)", nargs='+')

    sys.exit(main(parser.parse_args()))
