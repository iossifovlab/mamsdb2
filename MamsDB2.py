import os,struct,array,ctypes
import bisect

class BinaryFile(object):
    def __init__(self,fileName,fileAccess="file"):
        self.fileName=fileName
        self.fileAccess=fileAccess
        self._file=open(fileName,'rb')
        if fileAccess=="memory":
            self._data=bytearray(os.path.getsize(fileName))
            self._file.readinto(self._data)
            self._file.close()

    def readIndex(self,pos):
        if self.fileAccess=="memory":
            return self_.data[pos]
        elif self.fileAccess=="file":
            self._file.seek(pos)
            return self._file.read(1)

    def readRange(self,start,end):
        if self.fileAccess=="memory":
            return self_.data[start:end]
        elif self.fileAccess=="file":
            self._file.seek(start)
            return self._file.read(end-start)

    def close(self):
        if self.fileAccess=="memory":
            self._data=None
        elif self.fileAccess=="file":
            self._file.close()

class MUM_CData(ctypes.Structure):
    _fields_=[
        ("position",ctypes.c_uint64,32),
        ("chromosome",ctypes.c_uint64,8),
        ("offset",ctypes.c_uint64,10),
        ("length",ctypes.c_uint64,10),
        ("flipped",ctypes.c_uint64,1),
        ("read_2",ctypes.c_uint64,1),
        ("last_hit",ctypes.c_uint64,1),
        ("touches_end",ctypes.c_uint64,1)
    ]

    def __lt__(self,other):
        '''
        required for binary search
        '''
        if self.chromosome<other.chromosome:
            return True
        elif self.chromosome>other.chromosome:
            return False

        # the chromosomes are the same
        if self.position<other.position:
            return True
        else:
            return False

class MUMFile(BinaryFile):
    def readIndex(self,index):
        pos=ctypes.sizeof(MUM_CData)*index        
        if self.fileAccess=="memory":
            data=MUM_CData.from_buffer(self._data,pos)
        elif self.fileAccess=="file":
            data=MUM_CData()
            self._file.seek(pos)
            self._file.readinto(data)
        
        return data

    def readRange(self,startIndex,endIndex):
        records=[]
        for i in xrange(startIndex,endIndex):
            records.append(self.readIndex(i))
        return records

    def getNumberOfMUMs(self):
        return os.path.getsize(self.fileName)/ctypes.sizeof(MUM_CData)

class Pair_CData(ctypes.Structure):
    _fields_=[
        ("mums_start",ctypes.c_uint64,40),
        ("read_1_length",ctypes.c_uint64,10),
        ("read_2_length",ctypes.c_uint64,10),
        ("read_1_bad",ctypes.c_uint64,1),
        ("read_2_bad",ctypes.c_uint64,1),
        ("has_mums",ctypes.c_uint64,1),
        ("dupe",ctypes.c_uint64,1)
    ]

class PairFile(BinaryFile):
    def readIndex(self,index):
        pos=ctypes.sizeof(Pair_CData)*index

        if self.fileAccess=="memory":
            data=Pair_CData.from_buffer(self._data,pos)
        elif self.fileAccess=="file":
            data=Pair_CData()
            self._file.seek(pos)
            self._file.readinto(data)
        
        return data

    def readRange(self,startIndex,endIndex):
        records=[]
        for i in xrange(startIndex,endIndex):
            records.append(self.readIndex(i))
        return records

    def getNumberOfPairs(self):
        return os.path.getsize(self.fileName)/ctypes.sizeof(Pair_CData)


class Index_CData(ctypes.Structure):
    '''
    the mum_index here is the index of the mum in the pair, not in the mum file
    '''
    _fields_=[
        ("pair_index",ctypes.c_uint64,40),
        ("mum_index",ctypes.c_uint64,11),
        ("second_mum_index",ctypes.c_uint64,11),
        ("filler",ctypes.c_uint64,2)
    ]

class IndexFile(BinaryFile):
    def readIndex(self,index):
        pos=ctypes.sizeof(Index_CData)*index
        if self.fileAccess=="memory":
            data=Index_CData.from_buffer(self._data,pos)
        elif self.fileAccess=="file":
            data=Index_CData()
            self._file.seek(pos)
            self._file.readinto(data)
        
        return data

    def readRange(self,startIndex,endIndex):
        records=[]
        for i in xrange(startIndex,endIndex):
            records.append(self.readIndex(i))
        return records

    

class IndexSearch(object):
    '''
    Binary Search of mums by chromosome and position. Implements a list interface so that it is compatible with the bisect method.
    '''
    def __init__(self,indexFile,mumFile,pairsFile):
        self.indexFile=indexFile
        self.mumFile=mumFile
        self.pairsFile=pairsFile

    def __len__(self):
        return os.path.getsize(self.indexFile.fileName)/ctypes.sizeof(Index_CData)

    def __getitem__(self,index):
        '''
        Given an index in the index file, get the mum object associated with it
        '''
        indexEntry=self.indexFile.readIndex(index)
        pairsEntry=self.pairsFile.readIndex(indexEntry.pair_index)
        return self.mumFile.readIndex(pairsEntry.mums_start+indexEntry.mum_index)
    
class Mappability(object):
    def __init__(self,fasta_name,fileAccess="file"):
        self._lowFN=fasta_name+".bin/map.low.bin"
        self._highFN=fasta_name+".bin/map.high.bin"
        if os.path.getsize(self._lowFN)!=os.path.getsize(self._highFN):
           raise Error("The low and high Mappability files should be the same size")

        self._lowF=BinaryFile(self._lowFN,fileAccess)
        self._highF=BinaryFile(self._highFN,fileAccess)
        self.fasta=fasta_name

    def low(pos):
        return struct.unpack("B",self._lowF.readIndex(pos))

    def high(pos):
        return struct.unpack("B",self._highF.readIndex(pos))

    @staticmethod
    def createFromMumdexDir(mumdexDir,fileAccess="file"):
        fasta_name=""
        with open(os.path.join(mumdexDir,"ref.txt")) as refF:
            fasta_name=refF.readline().strip()

        return Mappability(fasta_name,fileAccess)
    
        
        
class Reference(object):
    def __init__(self,fasta_name,fileAccess="file"):
        self._seqFN=fasta_name+".bin/ref.seq.bin"
        self._seqF=BinaryFile(self._seqFN,fileAccess)
        self._chrLenFN=fasta_name+".bin/ref.chr_len.bin"
        self._chrNameFN=fasta_name+".bin/ref.chr_name.bin"
        self.chrLen=array.array("I")
        with open(self._chrLenFN) as chrLenF:
            self.chrLen.fromstring(chrLenF.read())

        self.chrAbsPos=[]
        pos=0
        for length in self.chrLen:
            self.chrAbsPos.append(pos)
            pos+=length
            
        self.chrName=[]
        with open(self._chrNameFN) as chrNameF:
            for line in chrNameF:
                self.chrName.append(line.strip())
                
        self.fasta=fasta_name

    @staticmethod
    def createFromMumdexDir(mumdexDir,fileAccess="file"):
        fasta_name=""
        with open(os.path.join(mumdexDir,"ref.txt")) as refF:
            fasta_name=refF.readline().strip()

        return Reference(fasta_name,fileAccess)
        
    def name(self,chromIndex):
        return self.chrName[chromIndex]

    def chromToIndex(self,chromName):
        return self.chrName.index(chromName)

    def CP2APos(self,chrom,pos):
        chromIndex=self.chromToIndex(chrom)
        return self.chrAbsPos[chromIndex]+pos

    def getSeqChrom(self,chrom,beg,end):
        absBeg=self.CP2APos(chrom,beg)
        absEnd=self.CP2APos(chrom,end)
        return self._seqF.readRange(absBeg,absEnd)

    def getSeqAbs(self,absBeg,absEnd):
        return self._seqF.readRange(absBeg,absEnd)




class MamRead:
    def __init__(mamRead,mamsDB,pairI,readN,len,isPCRdup):
        mamRead.mamsDB = mamsDB
        mamRead.pairI = pairI
        mamRead.readI = (pairI,readN)
        mamRead.readN = readN
        mamRead.readLen = len 
        mamRead.isPCRdup = isPCRdup
        mamRead.mams = []
        mamRead.mateRead = None
 
    def seq(mamRead):
        '''
        TODO: Implement reading gene sequences from mumdex
        '''
        return mamRead.mamsDB.mums.sequences(mamRead.pairI)[mamRead.readN-1]

class MAM(object):
    
    def __repr__(mam):
        return (mam.rp,mam.ch,mam.chPos,mam.ln,mam.st).__repr__()

    def __hash__(self):
        return hash((self.rp,self.ch,self.chPos,self.ln,self.st))

    def __eq__(self,other):
        if (self.rp,self.ch,self.chPos,self.ln,self.st)==(other.rp,other.ch,other.chPos,other.ln,other.st):
            return True
        else:
            return False

    def __neq__(self,other):
        return not self.__eq__(other)        

    @property
    def chBegPos(mam):
        if mam.st == '+':
            return mam.chPos
        else:
            return mam.chPos + mam.ln - 1

    @property
    def chEndPos(mam):
        if mam.st == '+':
            return mam.chPos + mam.ln - 1
        else:
            return mam.chPos

    @property
    def rBegChPos(mam):
        if mam.st == '+':
            return mam.chPos-mam.rp
        else:
            return mam.chPos + mam.ln - 1 + mam.rp

    @property
    def rEndChPos(mam):
        if mam.st == '+':
            return mam.chPos-mam.rp+mam.read.readLen-1
        else:
            return mam.chPos+mam.ln+mam.rp-mam.read.readLen

    @property
    def rBegAPos(mam):
        return mam.read.mamsDB.ref.CP2APos(mam.ch,mam.rBegChPos)

    @property
    def rEndAPos(mam):
        return mam.read.mamsDB.ref.CP2APos(mam.ch,mam.rEndChPos)
    
    @property
    def APos(mam):
        return mam.read.mamsDB.ref.CP2APos(mam.ch,mam.chPos)

    @property
    def seq(mam):
        return mam.read.mamsDB.ref.getS_BL(mam.ch,mam.chPos,mam.ln)

        
class MamsDB:
    def __init__(self,mamsDBDir,fileAccess="file"):
        self.mums=MUMFile(os.path.join(mamsDBDir,"mums.bin"),fileAccess)
        self.pairs=PairFile(os.path.join(mamsDBDir,"pairs.bin"),fileAccess)
        self.index=IndexFile(os.path.join(mamsDBDir,"index.bin"),fileAccess)
        self.ref = Reference.createFromMumdexDir(mamsDBDir,fileAccess)    
        self.mpb = Mappability.createFromMumdexDir(mamsDBDir,fileAccess)

    def close(self):
        self.mums.close()
        self.pairs.close()
        self.index.close()
        self.ref.close()
        self.mpb.close()

    def getNumReads(self):
        return 2*self.pairs.getNumberOfPairs()

    def getNumMams(self):
        return self.mums.getNumberOfMUMs()

    def buildPair(self,mamSortI):
        indexData=self.index.readIndex(mamSortI)
        pair = self.pairs.readIndex(indexData.pair_index)
        mumIndex = pair.mums_start+indexData.mum_index

        read1 = MamRead(self,indexData.pair_index,1,pair.read_1_length,pair.dupe)
        read2 = MamRead(self,indexData.pair_index,2,pair.read_2_length,pair.dupe)

        read1.mateRead = read2
        read2.mateRead = read1

        theMam=None
        # create a mam object for each mum associated with the read and link them together
        for mIndex in xrange(pair.mums_start,self.getMumStop(indexData.pair_index)):

            mum = self.mums.readIndex(mIndex)

            mam = MAM()
            mam.rp = mum.offset 
            mam.ch = self.ref.name(mum.chromosome)
            mam.chPos = mum.position-1
            mam.ln = mum.length
            mam.st = "-" if mum.flipped else '+'
            mam.mamI = mIndex

            mam.read = read2 if mum.read_2 else read1
            mam.read.mams.append(mam)
            if mIndex==mumIndex:
                theMam=mam
                            

        return read1,read2,theMam

    def getMams(self,chr,beg,end):
        chromInt=self.ref.chromToIndex(chr)
        toSearch=IndexSearch(self.index,self.mums,self.pairs)
        startMum=MUM_CData(position=beg-150,chromosome=chromInt)
        endMum=MUM_CData(position=end,chromosome=chromInt)
        startIndex=bisect.bisect_left(toSearch,startMum)
        endIndex=bisect.bisect_left(toSearch,endMum)

        for i in xrange(startIndex, endIndex):
            read1,read2,mam = self.buildPair(i)
            yield mam


    def low_map(self,ch,b,e=None):
        bA = self.ref.CP2APos(ch,b)
        if e:
            eA =  self.ref.CP2APos(ch,e)
        else:
            eA = bA + 1
            
        return [self.mpb.low_map(i) for i in xrange(bA+1,eA+1)]

    def high_map(self,ch,b,e):
        bA = self.ref.CP2APos(ch,b)
        if e:
            eA =  self.ref.CP2APos(ch,e)
        else:
            eA = bA + 1
            
        return [self.mpb.high_map(i) for i in xrange(bA+1,eA+1)]

    def getMumStop(self,pairIndex):
        '''
        Only the index of the first mum is stored in the pair data, so we must calculate the index of the last mum from the start of the next pair in the file.
        '''
        
        # handle case where pair is the last pair in the file
        if pairIndex+1 == self.pairs.getNumberOfPairs():
            return self.mums.getNumberOfMUMs()
        else:
            nextPair=self.pairs.readIndex(pairIndex+1)
            return nextPair.mums_start
