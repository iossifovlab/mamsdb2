
import os,struct,array,ctypes
import bisect
import mmap

class BinaryFile(object):
    def __init__(self,fileName,fileAccess="file"):
        self.fileName=fileName
        self.fileAccess=fileAccess
        self._file=open(fileName,'rb')
        if fileAccess=="memory":
            # create an anonymous mmap so the data can be shared across multiple processes
            self._data=mmap.mmap(-1,os.path.getsize(self.fileName))
            self._file.readinto(self._data)
            self._file.close()

        elif fileAccess=="mmap":
            self._data=mmap.mmap(self._file.fileno(),0,prot=mmap.PROT_READ)

    def readIndex(self,pos):
        if self.fileAccess=="memory" or self.fileAccess=="mmap":
            return self_.data[pos]
        elif self.fileAccess=="file":
            self._file.seek(pos)
            return self._file.read(1)

    def readRange(self,start,end):        
        if self.fileAccess=="memory" or self.fileAccess=="mmap":
            return self_.data[start:end]
        elif self.fileAccess=="file":
            self._file.seek(start)
            return self._file.read(end-start)

    def close(self):
        if self.fileAccess=="memory":
            self._data=None
        elif self.fileAccess=="file":            
            self._file.close()
        elif self.fileAccess=="mmap":
            self._data.close()
            self._file.close()

class BinaryCDataFile(BinaryFile):
    '''
    A binary data file where each record is defined by a ctypes Structure class
    '''
    def __init__(self,fileName,cDataClass,fileAccess="file"):
        BinaryFile.__init__(self,fileName,fileAccess)
        self.cDataClass=cDataClass
        self.sizeOfRecord=ctypes.sizeof(self.cDataClass)
        self.numRecords=os.path.getsize(self.fileName)/ctypes.sizeof(self.cDataClass)

    def readIndex(self,index):
        pos=self.sizeOfRecord*index
        if self.fileAccess=="memory" or self.fileAccess=="mmap":
            data=self.cDataClass.from_buffer(self._data,pos)
        elif self.fileAccess=="file":
            data=self.cDataClass()       
            self._file.seek(pos)
            self._file.readinto(data)            

        return data

    def readRange(self,startIndex,endIndex):
        records=[]
        for i in xrange(startIndex,endIndex):
            records.append(self.readIndex(i))
        return records
        

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

class IndexSearch(object):
    '''
    Binary Search of mums by chromosome and position. Implements a list interface so that it is compatible with the bisect method.
    '''
    def __init__(self,indexFile,mumFile,pairsFile):
        self.indexFile=indexFile
        self.mumFile=mumFile
        self.pairsFile=pairsFile

    def __len__(self):
        return self.indexFile.numRecords

    def __getitem__(self,index):
        '''
        Given an index in the index file, get the mum object associated with it
        '''
        indexEntry=self.indexFile.readIndex(index)
        pairsEntry=self.pairsFile.readIndex(indexEntry.pair_index)
        return self.mumFile.readIndex(pairsEntry.mums_start+indexEntry.mum_index)

# Each base is stored as a 3 bit integer ingeter in the files. The integer value cooresponds to the index in the IntToBaseMap
IntToBaseMap=["0","A","T","C","G","N","X","-","-"]
Complement={
    "A":"T",
    "T":"A",
    "C":"G",
    "G":"C",
    "N":"N",
    "X":"X"
}


class Bases_CData(ctypes.Structure):
    '''
    If is_index is true then data is an integer pointing to the index where the bases are stored in the extra bases file. Otherwise the data represents a list of bases, up to 21 bases, where each base is stored as a 3 bit integer.
    '''
    _fields_=[
        ("fileData",ctypes.c_uint64)
    ]

    @property
    def is_index(self):
        index_bit=1<<63
        return bool(self.fileData & index_bit)

    @property
    def data(self):
        data_bits=~(1<<63)
        return self.fileData & data_bits
        

class ExtraBases_CData(ctypes.Structure):
    _fields_=[("data",ctypes.c_uint64)]


class AllBases(object):
    def __init__(self,mamsDBDir,fileAccess="file"):
        self.bases=BinaryCDataFile(os.path.join(mamsDBDir,"bases.bin"),Bases_CData,fileAccess)
        self.extraBases=BinaryCDataFile(os.path.join(mamsDBDir,"bases.extra.bin"),ExtraBases_CData,fileAccess)

    def getBases(self,pairIndex):
        baseData=self.bases.readIndex(pairIndex)
        if baseData.is_index:
            return self._getExtraBases(baseData.data)
        else:
            return self._getBases(baseData.data)

    def _getBases(self,data):
        result=bytearray()
        baseIndex=0
        while True:
            base=IntToBaseMap[(data >> (baseIndex * 3)) & 7]
            if base == 'X':
                return str(result)
            result.append(base)
            baseIndex+=1            

    def _getExtraBases(self,startIndex):        
        result=bytearray()
        # in the extra base file characters can potentially span accross to integers, so we must read in the characters one bit at a time.
        base_index=startIndex
        while True:
            baseBits=0           
            for bit in range(0,3):                
                bit_index=base_index * 3 + bit
                word_index=bit_index/64
                bit_in_word_index=bit_index%64
                word=self.extraBases.readIndex(word_index).data
                baseBits |= (((word >> bit_in_word_index) & 1) << bit)
                    
            base=IntToBaseMap[baseBits]
            if base == 'X':
                return str(result)
            result.append(base)
            base_index+=1
    

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


class lazy_property(object):
    def __init__(self,fget):
        self.fget = fget
        self.func_name = fget.__name__
    def __get__(self,obj,cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj,self.func_name,value)
        return value

class MamRead:
    def __init__(self,mamsDB,pairI,readN,len,isPCRdup):
        self.mamsDB = mamsDB
        self.pairI = pairI
        self.readI = (pairI,readN)
        self.readN = readN
        self.readLen = len 
        self.isPCRdup = isPCRdup
        self.mams = []
        self.mateRead = None
 
    def seq(self):
        return self.mamsDB.getSequence(self.pairI)[self.readN-1]

class MAM(object):
    def __init__(self,cData,mamI,read):
        self.cData=cData
        self.mamI=mamI
        self.read=read
    
    def __repr__(self):
        return (self.rp,self.ch,self.chPos,self.ln,self.st).__repr__()

    def __hash__(self):
        return hash((self.rp,self.ch,self.chPos,self.ln,self.st))

    def __eq__(self,other):
        if (self.rp,self.ch,self.chPos,self.ln,self.st)==(other.rp,other.ch,other.chPos,other.ln,other.st):
            return True
        else:
            return False

    def __neq__(self,other):
        return not self.__eq__(other)

    # converting between cData types and python data types is expensive, only do it when necessary.
    @lazy_property
    def rp(self):
        return self.cData.offset

    @lazy_property
    def ch(self):
        return self.read.mamsDB.ref.name(self.cData.chromosome)

    @lazy_property
    def chPos(self):
        return self.cData.position

    @lazy_property
    def ln(self):
        return self.cData.length

    @lazy_property
    def st(self):
        return "-" if self.cData.flipped else '+'

    @property
    def chBegPos(self):
        if self.st == '+':
            return self.chPos
        else:
            return self.chPos + self.ln - 1

    @property
    def chEndPos(self):
        if self.st == '+':
            return self.chPos + self.ln - 1
        else:
            return self.chPos

    @property
    def rBegChPos(self):
        if self.st == '+':
            return self.chPos-self.rp
        else:
            return self.chPos + self.ln - 1 + self.rp

    @property
    def rEndChPos(self):
        if self.st == '+':
            return self.chPos-self.rp+self.read.readLen-1
        else:
            return self.chPos+self.ln+self.rp-self.read.readLen

    @property
    def rBegAPos(self):
        return self.read.mamsDB.ref.CP2APos(self.ch,self.rBegChPos)

    @property
    def rEndAPos(self):
        return self.read.mamsDB.ref.CP2APos(self.ch,self.rEndChPos)
    
    @property
    def APos(self):
        return self.read.mamsDB.ref.CP2APos(self.ch,self.chPos)

    @property
    def seq(self):
        return self.read.mamsDB.ref.getS_BL(self.ch,self.chPos,self.ln)

        
class MamsDB:
    def __init__(self,mamsDBDir,fileAccess="mmap"):
        self.mums=BinaryCDataFile(os.path.join(mamsDBDir,"mums.bin"),MUM_CData,fileAccess)
        self.pairs=BinaryCDataFile(os.path.join(mamsDBDir,"pairs.bin"),Pair_CData,fileAccess)
        # The index takes up a lot of memory and can be used on disk at a reasonable speed. Using the index through mmap instead of loading the whole thing into memory has a performance penalty of roughly 2x, but it decreases memory usage by 1/3.
        self.index=BinaryCDataFile(os.path.join(mamsDBDir,"index.bin"),Index_CData,"mmap")

        # These files may or may not be needed to be loaded into memory depending on the query. The ref+mappability take 10GB. The files for the bases takes between 8 and 15GB.
        self.ref = Reference.createFromMumdexDir(mamsDBDir,"mmap")    
        self.mpb = Mappability.createFromMumdexDir(mamsDBDir,"mmap")
        self.bases= AllBases(mamsDBDir,"mmap")

    def close(self):
        self.mums.close()
        self.pairs.close()
        self.index.close()
        self.ref.close()
        self.mpb.close()

    def getNumReads(self):
        return 2*self.pairs.numRecordsq

    def getNumMams(self):
        return self.mums.numRecords

    def buildPair(self,indexData):
        pair_index=indexData.pair_index        
        pair = self.pairs.readIndex(pair_index)
        mums_start=pair.mums_start
        
        mumIndex = indexData.mum_index

        read1 = MamRead(self,pair_index,1,pair.read_1_length,pair.dupe)
        read2 = MamRead(self,pair_index,2,pair.read_2_length,pair.dupe)

        read1.mateRead = read2
        read2.mateRead = read1

        theMam=None
        # create a mam object for each mum associated with the read and link them together
        for mIndex,mum in enumerate(self.mums.readRange(mums_start,self.getMumStop(pair_index))):
            mamI = mums_start+mIndex
            read = read2 if mum.read_2 else read1
            mam = MAM(mum,mamI,read)
            
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

        for indexData in self.index.readRange(startIndex, endIndex):
            read1,read2,mam = self.buildPair(indexData)
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
        if pairIndex+1 == self.pairs.numRecords:
            return self.mums.numRecords
        else:
            nextPair=self.pairs.readIndex(pairIndex+1)
            return nextPair.mums_start

    def getSequence(self,pairIndex):
        '''
        given a read pair index, get the sequence data associated with each pair.
        '''

        extra_bases=self.bases.getBases(pairIndex)
        extraBaseIndex=0
        pair=self.pairs.readIndex(pairIndex)
        results=[]

        for readNum in range(0,2):
            read_index=0
            length=0
            read=bytearray()
            for mum_index in range(pair.mums_start,self.getMumStop(pairIndex)):
                mum=self.mums.readIndex(mum_index)
                if mum.read_2 != readNum:
                    continue
                offset=mum.offset
                while read_index<offset:
                    read.append(extra_bases[extraBaseIndex])
                    extraBaseIndex+=1
                    read_index+=1

                while read_index < offset+mum.length:
                    chrom=self.ref.name(mum.chromosome)
                    if mum.flipped:                        
                        mumChromPos=mum.position+offset+mum.length-read_index-1
                        read.append(Complement[self.ref.getSeqChrom(chrom,mumChromPos,mumChromPos+1)])
                    else:
                        mumChromPos=mum.position+read_index-offset
                        read.append(self.ref.getSeqChrom(chrom,mumChromPos,mumChromPos+1))

                    read_index+=1

                if readNum:
                    length=pair.read_1_length
                else:
                    length=pair.read_2_length

            while read_index < length:
                read.append(extra_bases[extraBaseIndex])
                extraBaseIndex+=1
                read_index+=1

            results.append(str(read))

            
        return results
                        
                    
