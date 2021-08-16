import re
from itertools import product

import pandas as pd
from Bio import SeqIO

"""
ambiguity options used from https://droog.gs.washington.edu/parc/images/iupac.html
"""
def FindAmbigousOptions(seq):
    ambigs = {
        "A": ["A"],
        "T": ["T"],
        "C": ["C"],
        "G": ["G"],
        "M": ["A", "C"],
        "R": ["A", "G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "Y": ["C", "T"],
        "K": ["G", "T"],
        "V": ["A", "C", "G"],
        "H": ["A", "C", "T"],
        "D": ["A", "G", "T"],
        "B": ["C", "G", "T"],
        "N": ["G", "A", "T", "C"],
    }
    return list(map("".join, product(*map(ambigs.get, seq))))


def Primer_coordinates(inputprimer, reference):
    for record in SeqIO.parse(reference, "fasta"):

        possible_primers = FindAmbigousOptions(inputprimer)

        startlocs = []
        stoplocs = []
        for option in possible_primers:
            for match in re.finditer(str(option), str(record.seq), re.IGNORECASE):
                start_pos = match.start()
                end_pos = match.end()
                startlocs.append(start_pos)
                stoplocs.append(end_pos)

        if not startlocs and not stoplocs:
            startlocs = "None"
            stoplocs = "None"
            return startlocs, stoplocs
        else:
            return startlocs, stoplocs


def find_orient(primerfile, ref):
    left = ["LEFT", "PLUS", "POSITIVE"]
    right = ["RIGHT", "MINUS", "NEGATIVE"]

    FrameLeft = pd.DataFrame([])
    FrameRight = pd.DataFrame([])

    for record in SeqIO.parse(primerfile, "fasta"):
        if any(orient in record.id for orient in left) is True:
            startlist, stoplist = Primer_coordinates(record.seq, ref)
            if startlist == "None" or stoplist == "None":
                startlist, stoplist = Primer_coordinates(
                    record.seq.reverse_complement(), ref
                )
            if startlist != "None" or stoplist != "None":
                FrameLeft = FrameLeft.append(
                    pd.DataFrame(
                        {"start": [startlist], "stop": [stoplist], "name": record.id},
                        index=[0],
                    ),
                    ignore_index=True,
                )
        if any(orient in record.id for orient in right) is True:
            startlist, stoplist = Primer_coordinates(record.seq, ref)
            if startlist == "None" or stoplist == "None":
                startlist, stoplist = Primer_coordinates(
                    record.seq.reverse_complement(), ref
                )
            if startlist != "None" or stoplist != "None":
                FrameRight = FrameRight.append(
                    pd.DataFrame(
                        {"start": [startlist], "stop": [stoplist], "name": record.id},
                        index=[0],
                    ),
                    ignore_index=True,
                )
    return FrameLeft, FrameRight

def GetUniqs(frame):
    uniq = {}
    for row in frame.itertuples():
        startingpoints = []
        endingpoints = []
        for i in row.start:
            startingpoints.append(i)
        for i in row.stop:
            endingpoints.append(i)
        
        c = 0
        for a, b in zip(startingpoints, endingpoints):
            if row.name not in uniq:
                prname = row.name
            else:
                prname = f"{row.name}_{c}"
            uniq[prname] = {"start": a, "stop": b}
            c += 1
        
    return pd.DataFrame.from_dict(uniq, orient='index')

def UpdateName(primername):
    s = re.split('[- _ | /]', primername)
    for a, i in enumerate(s):
        try:
            i = int(i)
        except:
            i = str(i)
        
        if type(i) == int:
            s[a] = f"{i:03d}"
    return "_".join(s)

def CombinePrimerLists(a, b):
    aa = GetUniqs(a)
    bb = GetUniqs(b)
    
    c = pd.concat([aa, bb]).reset_index()
    c.rename(columns={"index": "name"}, inplace=True)
    c['name'] = c['name'].apply(lambda x: UpdateName(x))
    csort = c.sort_values(by=['name'])
    return csort

def GetRemovedCoordinates(coordinatedata):
    allcoordinates = []
    for a,b in coordinatedata.iteritems():
        for x in b:
            allcoordinates.append(x)
    return set(allcoordinates)


def GetRemovedPrimers(primers, coordinates):
    fprimers = pd.DataFrame(columns=primers.columns)
    for a in primers.itertuples():
        l = list(range(a.start, a.stop))
        if any(x in coordinates for x in l) is True:
            fprimers = fprimers.append(pd.DataFrame({'name': a.name, 'start': a.start, 'stop': a.stop}, index=[0]))
    fprimers.reset_index(inplace=True)
    fprimers = fprimers.drop(columns=['index'])
    return fprimers

def WritePrimerExports(left, right, processeddata, out):
    cPrimers = CombinePrimerLists(left, right)
    UniqueCoordinates = GetRemovedCoordinates(processeddata)
    filteredprimers = GetRemovedPrimers(cPrimers, UniqueCoordinates)
    filteredprimers.to_csv(out, index=False)
    
def MakeCoordinateLists(primerfile, ref):
    LeftPrimers, RightPrimers = find_orient(primerfile, ref)

    LeftList = []
    RightList = []

    for index, name in LeftPrimers.iterrows():
        for i in range(len(name.start)):
            list = [*range(name.start[i] - 1, name.stop[i], 1)]
            for iter in list:
                LeftList.append(iter)

    for index, name in RightPrimers.iterrows():
        for i in range(len(name.start)):
            list = [*range(name.start[i] + 1, name.stop[i], 1)]
            for iter in list:
                RightList.append(iter)

    return tuple(LeftList), tuple(RightList), LeftPrimers, RightPrimers
