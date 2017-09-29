import os,struct,array,ctypes

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

class MUMFile(BinaryFile):
    def readIndex(self,index):
        pos=ctypes.sizeof(MUM_CData)*index
        data=MUM_CData()
        if self.fileAccess=="memory":
            data.from_buffer(self._data,pos)
        elif self.fileAccess=="file":
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
        data=Pair_CData()
        if self.fileAccess=="memory":
            data.from_buffer(self._data,pos)
        elif self.fileAccess=="file":
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
    def createFromMumdexDir(mumdexDir):
        fasta_name=""
        with open(os.path.join(mumdexDir,"ref.txt")) as refF:
            fasta_name=refF.readline().strip()

        return Mappability(fasta_name)
    
        
        
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
    def createFromMumdexDir(mumdexDir):
        fasta_name=""
        with open(os.path.join(mumdexDir,"ref.txt")) as refF:
            fasta_name=refF.readline().strip()

        return Reference(fasta_name)
        
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

        

    
        
