from collections import defaultdict

def fileDictHeader(file, key, val):
    fileD={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols = line.rstrip("\n").split("\t")
            if lc==1:
                headD={}
                for i,x in enumerate(cols):
                    headD[x] = i
            else:
                fileD[cols[headD[key]]] = cols[headD[val]]
    return fileD


def fileList(file, col=0, delim="\t", header=True):
    l=[]
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                l.append(thisCols[col])
    return l

def fileListHeader(file, name, delim="\t"):
    l=[]
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            thisCols = line.rstrip("\n").split(delim)
            if lc == 1:
                col = thisCols.index(name)
            else:
                l.append(thisCols[col])
    return l


def fileEmptyDict(file, col=0, delim="\t", header=True):
    l={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                l[thisCols[col]]=""
    return l
    
def fileDict(file, key=0, val=1, delim="\t", header=True):
    l={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                try: l[thisCols[key]]=thisCols[val]
                except: print(thisCols)
    return l

def fileDictFull(file, delim="\t"):
    l=defaultdict(list)
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            thisCols = line.rstrip("\n").split(delim)
            if lc==1:
                headMap = {i:x for i,x in enumerate(thisCols)}
            else:
                for i,x in enumerate(thisCols):
                    l[headMap[i]].append(x)
    return l

